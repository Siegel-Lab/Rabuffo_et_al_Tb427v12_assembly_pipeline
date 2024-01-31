#!/usr/bin/env python

import argparse
from collections import defaultdict
import sys
import gffutils
import bokeh.plotting
import bokeh.models
import wiggelen
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches
import matplotlib.backends.backend_pdf
import matplotlib.gridspec
from matplotlib import markers
from matplotlib.path import Path
#Raul additins####
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
##added transparent=True to plt.savefig()

"""
##### Todo ####
- Add barplotting for html version
- set tick distance in html version
- remove redundancy in _coverage_bar_subplot_matplotlib
  and _coverage_line_subplot_matplotlib

- Split class into 
  - AnnotationSelector
  - VizHTML
  - VizPDF
  - helpers

bokeh.io.vplot was deprecated in Bokeh 0.12.0; please use bokeh.models.layouts.Column instead
"""


"""

Copyright (c) 2016, Konrad Foerstner <konrad@foerstner.org>

Permission to use, copy, modify, and/or distribute this software for
any purpose with or without fee is hereby granted, provided that the
above copyright notice and this permission notice appear in all
copies.

THE SOFTWARE IS PROVIDED 'AS IS' AND THE AUTHOR DISCLAIMS ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
         
"""
__description__ = ""
__author__ = "Konrad Foerstner <konrad@foerstner.org>"
__copyright__ = "2016 by Konrad Foerstner <konrad@foerstner.org>"
__license__ = "ISC license"
__email__ = "konrad@foerstner.org"
__version__ = "0.6dev"

def align_marker(marker, halign='center', valign='middle'):
    """
    create markers with specified alignment.

    Parameters
    ----------

    marker : a valid marker specification.
      See mpl.markers

    halign : string, float {'left', 'center', 'right'}
      Specifies the horizontal alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'center',
      -1 is 'right', 1 is 'left').

    valign : string, float {'top', 'middle', 'bottom'}
      Specifies the vertical alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'middle',
      -1 is 'top', 1 is 'bottom').

    Returns
    -------

    marker_array : numpy.ndarray
      A Nx2 array that specifies the marker path relative to the
      plot target point at (0, 0).

    Notes
    -----
    The mark_array can be passed directly to ax.plot and ax.scatter, e.g.::

        ax.plot(1, 1, marker=align_marker('>', 'left'))

    """

    if isinstance(halign, str):
        halign = {'right': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'left': 1.,
                  }[halign]

    if isinstance(valign, str):
        valign = {'top': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'bottom': 1.,
                  }[valign]

    # Define the base marker
    bm = markers.MarkerStyle(marker)

    # Get the marker path and apply the marker transform to get the
    # actual marker vertices (they should all be in a unit-square
    # centered at (0, 0))
    m_arr = bm.get_path().transformed(bm.get_transform()).vertices

    # Shift the marker vertices for the specified alignment.
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2

    return Path(m_arr, bm.get_path().codes)

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("--gff_file", required=True)
    parser.add_argument("--output_file", default="output.html")
    parser.add_argument("--replicons", help="Replicon IDs. If not given all "
                        "replicons are used.", nargs="+")
    parser.add_argument("--feature_types", help="Features. If not given all "
                        "feature types will be used.", nargs="+")
    parser.add_argument("--constant_width", default=False, action="store_true",
                        help="Plot with independt of replicon length")
    parser.add_argument("--wiggle_files", nargs="+", type=str)
    parser.add_argument("--wiggle_y_max", nargs="+", type=float,
                        help="Maximum y values for the wiggle plotting. For "
                        "each of the given wiggle file there must be value "
                        "given.")
    parser.add_argument("--wiggle_y_min", nargs="+", type=float,
                        help="Minimum y values for the wiggle plotting. For "
                        "each of the given wiggle file there must be value "
                        "given.")
    parser.add_argument("--default_color")
    parser.add_argument("--attribute_exclude",
                        help="Exclude entries if the given attribute contains "
                        "the given string. E.g. 'product=VSG domain;ID=Boing'")
    parser.add_argument("--attribute_include",
                        help="Require that the iven attribute contains "
                        "the given string. E.g. 'product=VSG domain;ID=Boing'")
    parser.add_argument("--feature_color_mapping", help="Allocate color to "
                        "features. E.g. 'rRNA=#FF0000;tRNA=#00FF00'")
    parser.add_argument("--attribute_color_mapping",
                        help="Allocate color to entries in the given "
                        "attribute contain the given string. E.g. "
                        "'product|VSG domain|#FF0000|||ID|Boing|Blue'")
    parser.add_argument("--report_feature_selection", default=False,
                        action="store_true", help="Report the feature used.")
    parser.add_argument("--pdf_multi_page", default=False,
                        action="store_true", help="Plot each graph on single "
                        "PDF page (neglected for HTML output.")
    parser.add_argument("--coverage_bar_plot", default=False,
                        action="store_true", help="Plot coverage as bar plot.")
    parser.add_argument("--bar_width", type=float,  default=1,
                        help="With of the bar in a bar plot.")
    parser.add_argument("--x_tick_distance", type=int,  default=None,
                        help="Distance of the ticks of thex x-axis")
    parser.add_argument("--alpha", type=float, help="Alpha value.")
    parser.add_argument("--font_size", type=int, default=10, help="Font size.")
    parser.add_argument("--version", "-v", action="version",
                        version=__version__, help="Print version.")
    args = parser.parse_args()

    anno_viz = AnnoViz(args.gff_file, args.output_file, args.replicons,
                       args.feature_types, args.feature_color_mapping,
                       args.attribute_color_mapping,
                       args.default_color, args.constant_width,
                       args.report_feature_selection,
                       args.attribute_include, args.attribute_exclude,
                       args.alpha, args.wiggle_files,
                       args.wiggle_y_min, args.wiggle_y_max,
                       args.pdf_multi_page, args.coverage_bar_plot,
                       args.bar_width, args.x_tick_distance,
                       args.font_size)
    anno_viz.get_replicon_lengths()
    anno_viz.read_annotations()
    anno_viz.read_wiggle_files()
    anno_viz.viz()


class AnnoViz(object):

    def __init__(self, gff_file, output_file, replicons, feature_types,
                 feature_color_mapping, attribute_color_mapping,
                 default_color, constant_width,
                 report_feature_selection, attribute_include,
                 attribute_exclude, alpha, wiggle_files,
                 wiggle_y_min, wiggle_y_max, pdf_multi_page,
                 coverage_bar_plot, bar_width, x_tick_distance, font_size):

        self._gff_file = gff_file
        self._output_file = output_file
        self._feature_types = feature_types
        self._replicons = replicons
        self._constant_width = constant_width
        self._max_plot_width = 1200
        self._plot_heigth = 200
        # The width should not become too small as otherwise due to
        # rounding the width ratios are affected.
        self._max_plot_width_matplotlib = 120
        self._plot_heigth_matplotlib_multipage = 30
        self._font_size_matplotlib = 10
        self._default_color = "#FF8080"
        if default_color is not None:
            self._default_color = default_color
        self.__create_feature_color_mapping(feature_color_mapping)
        self.__create_atttribute_color_mapping(attribute_color_mapping)
        self.__create_include_attributes(attribute_include)
        self.__create_exclude_attributes(attribute_exclude)
        self._report_feature_selection = report_feature_selection
        self._alpha = 0.5
        if alpha is not None:
            if not 0 < alpha < 1.0:
                sys.stderr.write("Alpha value must be between 0 and 1.\n")
                sys.exit()
            self._alpha = alpha
        self._wiggle_files = wiggle_files
        if wiggle_files is not None:
            if wiggle_y_min is not None:
                if len(wiggle_files) != len(wiggle_y_min):
                    sys.stderr.write("--wiggle_y_min must have the same "
                                     "number of values as --wiggle_files.")
                    sys.exit(1)
            
            if wiggle_y_max is not None:
                if len(wiggle_files) != len(wiggle_y_max):
                    sys.stderr.write("--wiggle_y_max must have the same "
                                     "number of values as --wiggle_files.")
                    sys.exit(1)
            self._wiggle_y_min = wiggle_y_min
            self._wiggle_y_max = wiggle_y_max
        self._pdf_multi_page = pdf_multi_page
        self._coverage_bar_plot = coverage_bar_plot
        self._bar_width = bar_width
        self._x_tick_distance = x_tick_distance
        self._font_size = font_size

    def __parse_mapping_string(self, mapping_string):
        dictionary = {}
        if mapping_string is None:
            return dictionary
        for key_and_value in mapping_string.split(";"):
            key, value = key_and_value.split("=")
            dictionary[key] = value
        return dictionary

    def __create_include_attributes(self, attribute_include):
        self._include_attributes = self.__parse_mapping_string(
            attribute_include)

    def __create_exclude_attributes(self, attribute_exclude):
        self._exclude_attributes = self.__parse_mapping_string(
            attribute_exclude)

    def __create_feature_color_mapping(self, feature_color_mapping_str):
        self._feature_color_mapping = self.__parse_mapping_string(
            feature_color_mapping_str)

    def __create_atttribute_color_mapping(self, attribute_color_mapping_str):
        self._atttribute_color_mapping = []
        if attribute_color_mapping_str is not None:
            for attr_string_color in attribute_color_mapping_str.split("|||"):
                attribute, string, color = attr_string_color.split("|")
                self._atttribute_color_mapping.append(
                    (attribute, string, color))

    def get_replicon_lengths(self):
        self._replicons_and_lengths = {}
        for line in open(self._gff_file):
            if line.startswith("##sequence-region"):
                self._replicons_and_lengths[line.split()[1]] = int(
                    line.split()[3])
        if self._replicons is None:
            self._replicons = self._replicons_and_lengths.keys()

    def read_annotations(self):
        self._annotation_glyphs_by_replicon = {}
        annotations = gffutils.create_db(self._gff_file, ':memory:')
        feature_counter = 0
        for replicon in self._replicons:
            glyph_data = defaultdict(list)
            for feature in annotations.region(replicon):
                if not self.__use_feature(feature):
                    continue
                feature_counter += 1
                assert feature.strand in ["+", "-", "."]
                bottom = 5
                top = 10
                if feature.strand == "-":
                    bottom = 0
                    top = 5
                if feature.strand == ".":
                    bottom = 0
                    top = 10
                glyph_data["bottom"].append(bottom)
                glyph_data["top"].append(top)
                glyph_data["left"].append(feature.start)
                glyph_data["right"].append(feature.end)
                glyph_data["start"].append(feature.start)
                glyph_data["end"].append(feature.end)
                glyph_data["zorder"].append(self._feature_types.index(feature.featuretype))
                glyph_data["ID"].append(
                    feature.attributes.get("ID", [""])[0])
                glyph_data["product"].append(
                    feature.attributes.get("product", [""])[0][:50])
                glyph_data["feature_type"].append(feature.featuretype)
                glyph_data["color"].append(self._color_of_feature(feature))
            # sort by zorder
            zipped_sorted = sorted(zip(glyph_data["zorder"], glyph_data["bottom"], glyph_data["top"], 
                glyph_data["left"], glyph_data["right"], glyph_data["start"], glyph_data["ID"], glyph_data["product"], glyph_data["feature_type"], glyph_data["color"]))
            self._annotation_glyphs_by_replicon[replicon] = {
                "zorder": [x[0] for x in zipped_sorted],
                "bottom": [x[1] for x in zipped_sorted],
                "top": [x[2] for x in zipped_sorted],
                "left": [x[3] for x in zipped_sorted],
                "right": [x[4] for x in zipped_sorted],
                "start": [x[5] for x in zipped_sorted],
                "ID": [x[6] for x in zipped_sorted],
                "product": [x[7] for x in zipped_sorted],
                "feature_type": [x[8] for x in zipped_sorted],
                "color": [x[9] for x in zipped_sorted],
            }
        if self._report_feature_selection:
            print("\nNumber of plotted features: {}".format(feature_counter))

    def _color_of_feature(self, feature):
        resulting_color = self._feature_color_mapping.get(
            feature.featuretype, self._default_color)
        if len(self._atttribute_color_mapping) > 0:
            for attribute, string, color in self._atttribute_color_mapping:
                if string in feature.attributes.get(attribute, [""])[0]:
                    resulting_color = color
        return resulting_color

    def __use_feature(self, feature):
        if self._feature_types is None:
            use_feature = True
        elif feature.featuretype not in self._feature_types:
            return False
        use_feature = False
        if len(self._include_attributes) > 0:
            for attribute, string in self._include_attributes.items():
                if string in feature.attributes.get(attribute, [""])[0]:
                    use_feature = True
                    break
        else:
            use_feature = True
        if len(self._exclude_attributes) > 0:
            for attribute, string in self._exclude_attributes.items():
                if string in feature.attributes.get(attribute, [""])[0]:
                    use_feature = False
        if self._report_feature_selection and use_feature:
            print(feature)
        return use_feature

    def viz(self):
        if self._output_file.endswith(".html"):
            self._viz_bokeh()
        elif (self._output_file.endswith(".pdf") or
              self._output_file.endswith(".svg")):
            #matplotlib.rcParams.update(
            #    {"font.size": self._font_size_matplotlib})
            if self._pdf_multi_page:
                self._viz_matplotlib_multipage()
            else:
                self._viz_matplotlib_singlepage()

    def _viz_matplotlib_singlepage(self):
        total_number_subplots = len(self._replicons)
        figure = plt.figure(figsize=(6.4, total_number_subplots))
        if self._wiggle_files is not None:
            total_number_subplots += len(
                self._replicons) + len(self._wiggle_files)
        subplot_index = 1
        max_width = max([self._plot_width_matplotlib(
            replicon) for replicon in self._replicons])
        grid_spec = matplotlib.gridspec.GridSpec(
            total_number_subplots, max_width, hspace=4)
        for replicon in self._replicons:
            width = self._plot_width_matplotlib(replicon)
            if self._wiggle_files is not None:
                for wiggle_file_index, wiggle_file in enumerate(
                        self._wiggle_files):
                    self._generate_coverage_subplot_matplotlib_singlepage(
                        figure, self._coverages[replicon][wiggle_file],
                        subplot_index, wiggle_file_index, replicon, width,
                        grid_spec)
                    subplot_index += 1
            self._generate_replicon_subplot_matplotlib_singlepage(
                figure, replicon, subplot_index, width, grid_spec)
            subplot_index += 1
        #plt.tight_layout()
        
        plt.savefig(self._output_file, transparent=True)
#        plt.savefig(self._output_file)

    def _generate_coverage_subplot_matplotlib_singlepage(
            self, figure, positions_and_coverages,
            subplot_index, wiggle_file_index, replicon, width, grid_spec):
        ax = plt.subplot(grid_spec[subplot_index-1, 0:width])
        if not self._coverage_bar_plot:
            self._coverage_line_subplot_matplotlib(
                wiggle_file_index, positions_and_coverages, replicon)
        else:
            self._coverage_bar_subplot_matplotlib(
                wiggle_file_index, positions_and_coverages, replicon)
        ax.set_title(replicon)
        self._set_ticks(ax)

    def _set_ticks(self, ax):
        if self._x_tick_distance == -1:
            ax.axes.get_xaxis().set_visible(False)
        if self._x_tick_distance is not None:
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, self._x_tick_distance))
            ax.tick_params(axis='x', which='both', labelbottom=False)
            ax.ticklabel_format(style='plain')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_position('center')
        ax.spines['bottom'].set_position('center')
        ax.xaxis.set_ticks_position('both')


    def _generate_replicon_subplot_matplotlib_singlepage(
            self, figure, replicon, subplot_index,
            width, grid_spec):
        ax = plt.subplot(grid_spec[subplot_index-1, 0:max(width, 1)])
        self._replicon_subplot_matplotlib(replicon, ax)
        ax.set_title(replicon, horizontalalignment="left", x=0)
        self._set_ticks(ax)

    def _generate_replicon_subplot_matplotlib_multipage(
            self, figure, replicon, total_number_subplots):
        ax = figure.add_subplot(
            total_number_subplots, 1, total_number_subplots)
        self._replicon_subplot_matplotlib(replicon, ax)
        plt.suptitle(replicon)
        self._set_ticks(ax)
        
    def _replicon_subplot_matplotlib(self, replicon, ax):
        if len(self._annotation_glyphs_by_replicon[replicon]["left"]) > 0:
            plt.xlim(min(self._annotation_glyphs_by_replicon[replicon]["left"]),
                    max(self._annotation_glyphs_by_replicon[replicon]["right"]))
            plt.ylim(min(self._annotation_glyphs_by_replicon[replicon]["bottom"]),
                    max(self._annotation_glyphs_by_replicon[replicon]["top"]))
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_yaxis().set_ticks([])
            for left, bottom, right, top, color in zip(
                    self._annotation_glyphs_by_replicon[replicon]["left"],
                    self._annotation_glyphs_by_replicon[replicon]["bottom"],
                    self._annotation_glyphs_by_replicon[replicon]["right"],
                    self._annotation_glyphs_by_replicon[replicon]["top"],
                    self._annotation_glyphs_by_replicon[replicon]["color"]):
                if color[1] == ":":
                    marker, color = color.split(":")
                else:
                    marker = None
                if marker == "-":
                    rect = matplotlib.patches.Rectangle(
                        (left, 4.975), right-left, 0.05, color=color,
                        alpha=self._alpha)
                    ax.add_patch(rect)
                else:
                    rect = matplotlib.patches.Rectangle(
                        (left, bottom), right-left, top-bottom, color=color,
                        alpha=self._alpha)
                    ax.add_patch(rect)
                    if marker is not None:
                        ax.plot((left + right) / 2, -1 if marker in "^d" else 11, 
                                marker=align_marker(marker.replace("p", "d"), 
                                valign='top' if marker in "^d" else "bottom"), 
                                color=color, clip_on=False, markersize=10)
        # Genome middle bar
        # bar = matplotlib.patches.Rectangle(
        #     (0, 4.975), self._replicons_and_lengths[replicon], 0.05,
        #     color="#000000")
        # ax.add_patch(bar)
        plt.xlim(0, self._replicons_and_lengths[replicon])

    def _viz_matplotlib_multipage(self):
        pp = matplotlib.backends.backend_pdf.PdfPages(self._output_file)
        for replicon in self._replicons:
            total_number_subplots = 1
            if self._wiggle_files is not None:
                total_number_subplots = len(self._wiggle_files) + 1
            figure = plt.figure(figsize=(
                self._plot_width_matplotlib(replicon),
                self._plot_heigth_matplotlib_multipage *
                total_number_subplots))
            if self._wiggle_files is not None:
                for wiggle_file_index, wiggle_file in enumerate(
                        self._wiggle_files):
                    self._generate_coverage_subplot_matplotlib_multipage(
                        figure, self._coverages[replicon][wiggle_file],
                        total_number_subplots, wiggle_file_index,
                        replicon)
            self._generate_replicon_subplot_matplotlib_multipage(
                figure, replicon, total_number_subplots)
            pp.savefig()
        pp.close()

    def _generate_coverage_subplot_matplotlib_multipage(
            self, figure, positions_and_coverages, total_number_subplots,
            wiggle_file_index, replicon):
        ax = figure.add_subplot(
            total_number_subplots, 1, wiggle_file_index+1)
        if not self._coverage_bar_plot:
            self._coverage_line_subplot_matplotlib(
                wiggle_file_index, positions_and_coverages, replicon)
        else:
            self._coverage_bar_subplot_matplotlib(
                wiggle_file_index, positions_and_coverages, replicon)
        self._set_ticks(ax)

    def _coverage_line_subplot_matplotlib(
            self, wiggle_file_index, positions_and_coverages, replicon):
        positions = [positions_and_coverage[0]
                     for positions_and_coverage in positions_and_coverages]
        coverages = [positions_and_coverage[1]
                     for positions_and_coverage in positions_and_coverages]
#        plt.plot(positions, coverages, color="#000000")
        plt.plot(positions, coverages, color="#0079bf")
        plt.xlim(0, max(positions))
        y_min = min(coverages)
        y_max = max(coverages)
        if self._wiggle_y_min is not None:
            y_min = self._wiggle_y_min[wiggle_file_index]
        if self._wiggle_y_max is not None:
            y_max = self._wiggle_y_max[wiggle_file_index]
        plt.ylim(y_min, y_max)
        plt.xlim(0, self._replicons_and_lengths[replicon])

    def _coverage_bar_subplot_matplotlib(
            self, wiggle_file_index, positions_and_coverages, replicon):
        positions = [positions_and_coverage[0]
                     for positions_and_coverage in positions_and_coverages]
        coverages = [positions_and_coverage[1]
                     for positions_and_coverage in positions_and_coverages]
        width = [self._bar_width] * len(positions)
        plt.bar(positions, coverages, width, color="#000000", linewidth=0)
        plt.xlim(0, max(positions))
        y_min = min(coverages)
        y_max = max(coverages)
        if self._wiggle_y_min is not None:
            y_min = self._wiggle_y_min[wiggle_file_index]
        if self._wiggle_y_max is not None:
            y_max = self._wiggle_y_max[wiggle_file_index]
        plt.ylim(y_min, y_max)
        plt.xlim(0, self._replicons_and_lengths[replicon])

    def _viz_bokeh(self):
        bokeh.plotting.output_file(self._output_file)
        sub_plots = []
        # if self._wiggle_files is not None:
        #     self._set_y_axis_ranges()
        for replicon in self._replicons:
            # Coverage plot
            if self._wiggle_files is not None:
                for wiggle_file_index, wiggle_file in enumerate(
                        self._wiggle_files):
                    sub_plots.append(self._generate_coverage_subplot(
                        self._coverages[replicon][wiggle_file], replicon,
                        wiggle_file_index))
            # Gene plot
            sub_plots.append(self._generate_replicon_subplot(replicon))
        figure_combo = bokeh.layouts.column(*sub_plots)
        bokeh.plotting.save(figure_combo)

    def _generate_coverage_subplot(self, positions_and_coverages, replicon,
                                   wiggle_file_index):
        figure = bokeh.plotting.figure(
            plot_width=self._plot_width(replicon),
            plot_height=self._plot_heigth,
            title=replicon, logo=None)
        positions = [positions_and_coverage[0]
                     for positions_and_coverage in positions_and_coverages]
        coverages = [positions_and_coverage[1]
                     for positions_and_coverage in positions_and_coverages]
        figure.line(x=positions, y=coverages)
        self._set_y_axis_ranges(figure, positions_and_coverages,
                                wiggle_file_index)
        return figure

    def _set_y_axis_ranges(self, figure, positions_and_coverages,
                           wiggle_file_index):
        # Both not specified => Nothing to due - just use the
        # automatically set range
        if self._wiggle_y_min is None and self._wiggle_y_max is None:
            return
        if self._wiggle_y_min is not None:
            y_min = self._wiggle_y_min[wiggle_file_index]
        else:
            y_min = min([cov for pos, cov in positions_and_coverages])
        if self._wiggle_y_max is not None:
            y_max = self._wiggle_y_max[wiggle_file_index]
        else:
            y_max = max([cov for pos, cov in positions_and_coverages])
        figure.y_range = bokeh.models.Range1d(y_min, y_max)

    def _plot_width_matplotlib(self, replicon):
        if self._constant_width:
            return self._max_plot_width_matplotlib
        replicon_length = [
            self._replicons_and_lengths[replicon]
            for replicon in self._replicons]
        width = int(self._replicons_and_lengths[replicon] /
                    max(replicon_length) *
                    self._max_plot_width_matplotlib)
        return width

    def _plot_width(self, replicon):
        if self._constant_width:
            return self._max_plot_width
        replicon_length = [
            self._replicons_and_lengths[replicon]
            for replicon in self._replicons]
        width = int(
            self._replicons_and_lengths[replicon] /
            max(replicon_length) *
            self._max_plot_width)
        return width

    def _generate_replicon_subplot(self, replicon):
        source = bokeh.plotting.ColumnDataSource(
            data=self._annotation_glyphs_by_replicon[replicon])
        hover = bokeh.models.HoverTool(
            tooltips=[("ID", "@ID"), ("Start", "@start"), ("End", "@end"),
                      ("Feature", "@feature_type"), ("Product", "@product")])
        figure = bokeh.plotting.figure(
            plot_width=self._plot_width(replicon),
            plot_height=self._plot_heigth,
            title=replicon, #logo=None,
            tools=[hover,
                   bokeh.models.BoxZoomTool(),
                   bokeh.models.ResetTool(),
                   bokeh.models.PanTool(),
                   #bokeh.models.ResizeTool(),
                   bokeh.models.WheelZoomTool()])

        #  bokeh.models.PreviewSaveTool() # TODO - seems to cause an error
        
        figure.quad("left", "right", "top", "bottom",
                    source=source, color="color",
                    fill_alpha=self._alpha, line_alpha=self._alpha)
        # Genome middle bar
        figure.quad(left=0, right=self._replicons_and_lengths[replicon],
                    top=5.1, bottom=4.9, color="#000000",
                    fill_alpha=self._alpha, line_alpha=self._alpha)
        figure.xgrid.grid_line_color = None
        figure.ygrid.grid_line_color = None
        figure.xaxis.axis_line_color = None
        figure.yaxis.minor_tick_line_color = None
        figure.xaxis.axis_label_text_color = None
        figure.yaxis.visible = False
        return figure

    def read_wiggle_files(self):
        self._coverages = defaultdict(lambda: defaultdict(list))
        if self._wiggle_files is not None:
            for wiggle_file in self._wiggle_files:
                self._read_wiggle_file(wiggle_file)

    def _read_wiggle_file(self, wiggle_file):
        for replicon, pos, coverage in wiggelen.walk(open(wiggle_file)):
            if replicon not in self._replicons:
                continue
            self._coverages[replicon][wiggle_file].append((pos, coverage))

if __name__ == "__main__":
    main()
