#!/bin/perl
# usage = perl adapt_gff.pl input.gff3 chrom_name cut_start cut_end gap_length> adapted.gff3
# Input file is the gff3 file from the genome that was modified, and the coordinates and chromosome name of the region that was removed
# The output is a gff3 file removing all annotation that were within the region removed and shifting the coordinates of the things downstream the region removed

use strict;
use warnings;

my $input_gff = $ARGV[0];
$" = "\t"; # Setting the output field separator to the tab character.

my $chrom_name = $ARGV[1];
my $cut_start = $ARGV[2];
my $cut_end = $ARGV[3];
my $gap = $ARGV[4];

my $shift = ($cut_end - $cut_start) - $gap;

open(my $gff_in, '<', $input_gff) or die "Could not open file '$input_gff' $!";

while(my $gff_row = <$gff_in>){
	# if it is a header or a comment line
	if ($gff_row =~ m/^#/ ){
		# if it is a header or a comment line of the length of a chromosome not involved in the cut, print the line and move to the next
		if($gff_row !~ m/$chrom_name/ ){
			print $gff_row;
			next;
		}
		# if it is the header of the length of the chromosome where the cut was done, adapt the end coordinate
		else{
			chomp $gff_row;
			my ($Prefix, $Chrom, $Start, $End) = split(/\s+/,$gff_row);
			my $New_end = $End - $shift -1;
			print "$Prefix   $Chrom $Start $New_end\n";
			next;
		}

	}

chomp $gff_row;
my ($Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes) = split("\t",$gff_row);

	
	if($Chrom !~ m/^$chrom_name/ ) {

		print "$Chrom\t$Source\t$Feature\t$Start\t$End\t$Score\t$Strand\t$Frame\t$Attributes\n";
	next;
	}
	# for contig entries, print them even if they partially overlap the cutting region, but shift the coordinates if they do
	elsif($Feature =~ m/^contig/ ){
		# if the feature overlaps with only the start position --> modify the end coordinate to one less than the cut_start
		if($Start < $cut_start and $End ~~ [$cut_start..$cut_end]){
			my $New_end = $cut_start - 1;
			print "$Chrom\t$Source\t$Feature\t$Start\t$New_end\t$Score\t$Strand\t$Frame\t$Attributes\n";
		next;
		}
		# if the feature overlaps only with the end position --> modify the start coordinate to one more than the cut_end
		elsif($Start ~~ [$cut_start..$cut_end] and $End > $cut_end ){
		my $New_start = $cut_start + $gap;
		my $New_end = $End - $shift - 1;
		print "$Chrom\t$Source\t$Feature\t$New_start\t$New_end\t$Score\t$Strand\t$Frame\t$Attributes\n";
		next;
		}
		# if the feature includes the cutting region, shorten the end by cut_end
		elsif($Start < $cut_start and $End > $cut_end ){
		my $New_end =$End - $shift;
		print "$Chrom\t$Source\t$Feature\t$Start\t$New_end\t$Score\t$Strand\t$Frame\t$Attributes\n";
                next;
		}
	}
	# do not print those entries in the cutting region (even those that start OR end within the cutting region, except for contigs entries that where already handled before)
	elsif(($Start ~~ [$cut_start..$cut_end]) or ($End ~~ [$cut_start..$cut_end]) ) {
	next;
	}
	# if both start and end are smaller than cut_start print them as they are
	elsif(($Start < $cut_start) and ($End < $cut_start)) {
		print "$Chrom\t$Source\t$Feature\t$Start\t$End\t$Score\t$Strand\t$Frame\t$Attributes\n";
	next;
	}
	# if the feature includes the cutting region, remove it
	elsif(($Start < $cut_start) and ($End > $cut_end)) {
	next;
	}
	# if the feature starts after the cutting region shift the start and end coordinates by cut_end - cut_start + 1000
	elsif($Start > $cut_end) {
		my $New_start = $Start - $shift - 1;
		my $New_end = $End - $shift - 1;
		print "$Chrom\t$Source\t$Feature\t$New_start\t$New_end\t$Score\t$Strand\t$Frame\t$Attributes\n";
	next;
	}
	else{
		print "I DID NOT CONSIDER THIS SCENARIO\n";
	next;
	}

}
close($gff_in) or die "Could not close filehandle '$gff_in' $!";
#print a line with the GAP
#format: BES17_Tb427v10	.	gap	65824	66823	.	.	.	estimated_length=1000;gap_type=within scaffold
my $gap_start = $cut_start;
my $gap_end = $cut_start + $gap -1;
print "$chrom_name\t.\tgap\t$gap_start\t$gap_end\t.\t.\t.\testimated_length=$gap;gap_type=within scaffold";

exit;

