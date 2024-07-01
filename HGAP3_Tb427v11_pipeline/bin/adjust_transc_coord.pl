#!/bin/perl
# usage = perl adjust_transc_coord.pl gene_mRNA_ordered.txt > adjusted_mRNA_coord.gff3
# Input is
# Output is

use strict;
use warnings;

my $gene_mRNA_ordered_file = $ARGV[0];

$" = "\t"; # Setting the output field separator to the tab character.	

my $last;
my $out;
open(my $gene_in, '<', $gene_mRNA_ordered_file) or die "Could not open file '$gene_mRNA_ordered_file' $!";

while(my $gene_row = <$gene_in>){

my ($Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes, $ID) = split("\t",$gene_row);
	
	if($Feature =~ m/gene/){# If feature is gene or pseudogene print the lines and store the coordinates
	
	$last = [$Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes, $ID];
	$out = [$Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes];
	print "@$out\n";
	}
	elsif(($Feature =~ m/mRNA/ || $Feature =~ m/pseudogenic_transcript/) && $ID eq $last->[9] ){#adjust the coordinates of the mRNA and pseudogenic_transcript entries to much the one of their respectives gene and pseudogene entries
	$out = [$Chrom, $Source, $Feature, $last->[3], $last->[4], $Score, $Strand, $Frame, $Attributes];
	print "@$out\n";
	}
	else{
	print "Something wrong happened, i am not expecting this scenario\n"
	}

}

close($gene_in) or die "Could not close filehandle '$gene_in' $!";
exit;


