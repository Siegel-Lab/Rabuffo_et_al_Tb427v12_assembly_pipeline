#!/bin/perl
# usage = perl merge_VSG_Companion_gff3.pl joined_VSG_Companion.gff3 > merged_VSG_Companion.gff3
# Input file is the gff3 file joining VSG_gff3 file with Companion, sorted.
# The output is a gff3 file were the VSG entries and the gene/pseudogene entries from Companion that are overlapping are being merged.

use strict;
use warnings;

my $joined_gff = $ARGV[0];

$" = "\t"; # Setting the output field separator to the tab character.
my $last;

open(my $joined_in, '<', $joined_gff) or die "Could not open file '$joined_gff' $!";

while(my $joined_row = <$joined_in>){

chomp $joined_row;
my ($Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes) = split("\t",$joined_row);

	if($Chrom =~ m/^#/ ){
	next;
	}
	# Print out every line which does not have in $Feature gene or pseudogene, and go to the next entry.
	if($Feature !~ m/gene/ ) {
		print "$Chrom\t$Source\t$Feature\t$Start\t$End\t$Score\t$Strand\t$Frame\t$Attributes\n";
	next;
	}
	elsif($last){# Test if last overlaps with the new entry and merge them in last.      
		if ($last->[0] ne $Chrom || $last->[4] < $Start) {# If the chromosome of the current entry is not equal to the previous one OR the end of the previous entry is smaller than the start of the current one, print the previous entry.
		print "@$last\n";
		}
		elsif($last->[8] =~m/ID=/) {# If the previous entry has ID (Is Companion or Companion merged with VSG), 
			if($Attributes =~ m/ID=/) {# If current entry also has ID (is Companion), print the last one (I am not going to merge Companion entries
			print "@$last\n";
			}

			else{# Current entry does not have ID (it is not Companion-derived, then it is BLAST-derived), then join.
			# Here I would be joining a current VSG-BLAST result together with a previous Companion-only result or with a Companion already merged with one(or more?) VSG entries.
				if($last->[8] =~ m/Name=Similar to Tb427VSG/){# If the previous matches "Name=Similar to Tb427VSG", then, it has already been joined with a BLAST-based annotation, and special parsing should be done.
					
					$Attributes =~ m/(Name=(.*);)/;
					my $VSG_name = $1;
					$VSG_name =~ s/;$//;

					$Attributes =~ m/(Note=(.*))/;
					my $VSG_note = $1;
					$last->[8] =~ s/Name=/$VSG_name,/;
					$last->[8] =~ s/Note=/$VSG_note /;
					$Attributes = $last->[8];
				}

				else{# if it hadn't been joined with a BLAST-based annotation already, classical merging of annotation should be done. Situation: Previous = Companion-only-derived ; Actual = BLAST-only-derived
					if($last->[8] =~ m/Name=/){#If the previous Companion-only-derived entry already has a "Name=", I need to join them properly.

						$Attributes =~ m/(Name=(.*);)/;
						my $VSG_name = $1;
						$VSG_name =~ s/;$//;
						
						$Attributes =~ m/(Note=(.*))/;
						my $VSG_note = $1;
						$last->[8] =~ s/Name=/$VSG_name,/;
						$Attributes = "$last->[8];$VSG_note";
					}
					else{#the previous Companion-only-derived entry does not have a "Name=",I just need to join it with the current BLAST-only-derived

						$Attributes = "$last->[8];$Attributes";
					}
				}
			
			$Source = $last->[1];
			$Feature = $last->[2];
			$Start = $last->[3];
			$End = $last->[4] if $End < $last->[4]; # Assign $last->[4] to $End if $End is smaller than $last->[4]
			$Score = $last->[5];
			$Strand = $last->[6];
			$Frame = $last->[7];

			}
		}
		else{# If the previous entry does not have ID, join them.--> Here, the previous would be a BLAST-only-derived and the current a Companion-only-derived
			$End = $last->[4] if $End < $last->[4]; # Assign $last->[4] to $End if $End is smaller than $last->[4]
			$Start = $last->[3];
			#Need to handle the situation if the current Companion-derived entry already has a "Name="
			if($Attributes =~m/Name=/){#If the actual Companion-derived entry already has a "Name=" join them

				$last->[8] =~ m/(Name=(.*);)/;
				my $VSG_name = $1;
				$VSG_name =~ s/;$//;	
				$Attributes =~ s/Name=/$VSG_name,/; # substitute "Name=" by "Name=Similar to ... ,"

				$last->[8] =~ m/(Note=(.*))/;
				my $VSG_note = $1;
				
				$Attributes = "$Attributes;$VSG_note"; #Add ";Note=..." to the end of Attributes

			}
			else{# If the actual Companion-derived entry does not have a "Name=", just paste the attributes
				$Attributes = "$Attributes;$last->[8]";
			}
		}

	}

	$last = [$Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes];
}
print "@$last\n";
close($joined_in) or die "Could not close filehandle '$joined_in' $!";
exit;

