#!/bin/perl
# usage = perl join_overlapping_entries_VSGs_gff3.pl VSG_gff3 > VSG_overlapping_joined_output.gff3
# Input file is the gff3 output generated with Blast_output_to_gff3_v2.pl from the BLAST output file from the Cross et al. (2014) VSG list to our genome. --> Sorted by chromosome and then by start position!
# Output is a gff3 file in which overlapping VSG entries were merged into one, and the only text in the attributes column is "description=Tb427VSG-XX"

use strict;
use warnings;

my $VSG_gff = $ARGV[0];

$" = "\t"; # Setting the output field separator to the tab character.	
my $last;
my $inc;
## Merge overlapping VSG features.. are there some? Yes, there are quite a bit. Like ~350.

open(my $VSG_in, '<', $VSG_gff) or die "Could not open file '$VSG_gff' $!";

while(my $VSG_row = <$VSG_in>){

chomp $VSG_row;
my ($VSG_Chrom, $VSG_Source, $VSG_Feature, $VSG_Start, $VSG_End, $VSG_Score, $VSG_Strand, $VSG_Frame, $VSG_Attributes) = split("\t",$VSG_row);
my $VSG_Name;
my $VSG_Description;
	if($VSG_Chrom =~ m/^#/ ) {
	next;
	}

	if($last) {# Test if last overlaps with the new entry and merge them in last.
		if ($last->[0] ne $VSG_Chrom || $last->[4] < $VSG_Start) {# If the chromosome of the current entry is not equal to the previous one OR the end of the previous entry is smaller than the start of the current one, print the previous entry.
		print "@$last\n";
		$inc = 0; #turn $inc to zero every time a merging ended
		}
		else {# If the ranges between the previous entry and the current one overlap, assign the start of the previous to $VSG_Start, If the current entry END is smaller than the previous, assign the previous one to $VSG_End,and add the annotation of the current one to $VSG_Attributes 
			## Modify this part to do two things: 
			#	(1): Put VSG_ID ("Tb427VSG-XX") alone in the attributes section, as description, in the following way "description=Tb427VSG-XX"
			#	(2): Assure that in case of merging entries, the coordinates and name of the one  with higher similarity is being used. --> i will just keep the one with the biggest alignment
			#Break the Attributes of the last entry ($last->[8]) storing only the first part into $VSG_Name
			if ($inc == 0 ) {# first two overlapping genes
			$last->[8] =~ m/(Name=(.*);)/;
			my $Last_Name = $1;
			$Last_Name =~ s/Name=//;
			$Last_Name =~ s/;$//;
				
			$VSG_Attributes =~ m/(Name=(.*);)/; # Match the "Name= until the ;" and store it in $1 variable
			$VSG_Name = $1 =~ s/Name=//r; # Assign to $VSG_Name the $1 variable but removing "Name="
			$VSG_Name =~ s/;$//;
				# if the current VSG alignment is bigger, keep its VSG_ID for the naming 
				if (($VSG_End - $VSG_Start) > ($last->[4] - $last->[3])) {
				$VSG_Description = $VSG_Name =~ s/Similar to //r;
				}
				# if the previous was bigger, keep the VSG_ID of the previous for the naming
				else {
				$VSG_Description = $Last_Name =~ s/Similar to //r;
				}
			}	
			else { # merging an entry on already merged genes
			$VSG_Attributes =~ m/(Name=(.*);)/; # Match the "Name= until the ;" and store it in $1 variable
			$VSG_Name = $1 =~ s/Name=//r; # Assign to $VSG_Name the $1 variable but removing "Name="
			$VSG_Name =~ s/;$//;	
				if (($VSG_End - $VSG_Start) > ($last->[4] - $last->[3])) {
				$VSG_Description = $VSG_Name =~ s/Similar to //r;
				}
				else {
				$VSG_Description = $last->[8] =~ s/description=//r;
				}



			}
			$VSG_Attributes = "description=$VSG_Description";

			# Assign the start of the previous to $VSG_Start and if the current entry END is smaller than the previous, assign the previous one to $VSG_End
			$VSG_Start = $last->[3];
			$VSG_End = $last->[4] if $VSG_End < $last->[4];
			$inc++;
		}

	}
	$last = [$VSG_Chrom, $VSG_Source, $VSG_Feature, $VSG_Start, $VSG_End, $VSG_Score, $VSG_Strand, $VSG_Frame, $VSG_Attributes];
}
print "@$last\n";
close($VSG_in) or die "Could not close filehandle '$VSG_in' $!";
exit;
