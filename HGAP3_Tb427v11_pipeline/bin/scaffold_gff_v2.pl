#!/bin/perl
# this script goes through a list with "headers" tab separated (those that should be joined) an input multifasta file and a gff file. It scaffolds the gffs with 1000 Ns in between those sequences from the headers tab separated.
# usage = perl input.fasta input.gff data.txt > output.gff

# I am supposing the presence of 22 alleles (11 chromosomes), and that the chromosomes in the guide file are ordered in decreasing number (Chr11, Chr10, Chr9, etc..) and that the "A" allele is first, and then the "B" allele

use strict;
use Bio::DB::Fasta;

my $fastafile = $ARGV[0];
my $gfffile = $ARGV[1];
my $datafile = $ARGV[2];

$" = "\t"; # Setting the output field separator to the tab character.

my $db = Bio::DB::Fasta->new($fastafile);

#store the $gfffile into an array of lines named @gff_lines
open(my $gff_in, '<', $gfffile) or die "Could not open file '$gfffile' $!";
my @gff_lines = <$gff_in>;
chomp (@gff_lines);
close($gff_in) or warn $!;

my $Allele = 22;
my $Char;
my $New_Chrom_name;

open(FILE,$datafile);

#print "##gff-version 3\n";

while(<FILE>){
chomp;
my @data = split("\t");
my @data2 = split("\t");
my $dataSize = @data;
my $dataLast = $dataSize - 1;
	if ($Allele > 0){
		# To get Chromosome number from allele, supposing the alleles are ordered in the guide file decreasingly from Chr11 to Chr1, and the alleles are together, which will lead to Chr11, Chr11 , Chr10, Chr10, Chr9, Chr9, etc.
		my $Chrom = int(($Allele / 2) + 0.5);
#		I add an if statement, if $Allele is even i print "_A" after the chromosome number, if $Allele is odd i print "_B" after the chromosome number.. this supposes that the alleles are ordered in the guide file, "A" allele first, "B" allele then,
		if ($Allele % 2 == 0){
			$Char = "A";
		}
		else {

			$Char = "B";
		}
				
#		$New_Chrom_name = "Chr".$Chrom."_".$Char."_Tb427v10";
		$New_Chrom_name = "Chr".$Chrom."_hap".$Char."_Tb427v11";
	}

	else { # after the scaffolded chromosomes, print the rest unscaffolded BESs
	$New_Chrom_name = $data2[0];
	}
# Until here is to generate to $New_Chrom_name variable

	if ($Allele > 0){
		$Allele -=1;
		my $addition = 0;
		my $i;
		for ( $i = 0; $i<= $dataLast; $i++){
		
		my $length = $db->length($data[$i]); # the length of the current contig
	
			if ($i == 0 && $data[$i] =~/^BES/){
			# reverse complement gff
				foreach my $gff_line (@gff_lines){
					my ($Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes) = split("\t",$gff_line);
					if($Chrom eq $data[$i]){
						my $New_start = $length - int($End) + 1;
						my $New_end = $length - int($Start) + 1;
						my $New_strand;
						if($Strand eq "+"){
							$New_strand = "-";
						}
						elsif($Strand eq "-"){
							$New_strand = "+";
						}
						else{
							$New_strand = $Strand;
						}

						print "$New_Chrom_name\t$Source\t$Feature\t$New_start\t$New_end\t$Score\t$New_strand\t$Frame\t$Attributes\n";
					}
				}
		

			}

			else{
			#print all the gff entries from that chromosome, modifiying the chrom_name Start and End
				foreach my $gff_line (@gff_lines) {
					my ($Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes) = split("\t",$gff_line);
					if($Chrom eq $data[$i]){
						my $New_start = int($Start) + $addition;
						my $New_end = int($End) + $addition;
						print "$New_Chrom_name\t$Source\t$Feature\t$New_start\t$New_end\t$Score\t$Strand\t$Frame\t$Attributes\n";
					}

				}
			}
			$addition = $addition + $length + 1000;
		}#close the for loop through the tab separated contigs to be joined

	} # close the if loop through the 22 chromosomes




#			my $subseqstr; # I don't need this
#				if ($i == 0 && $data[$i] =~/^BES/){
#					$subseqstr = $db->seq($data[$i],$length,1); # I don't need this
#					# I need to reverse complement the annotation
#				}# Reverse-complement the BES that go on the 5'
#				else {
#					$subseqstr = $db->seq($data[$i],1,$length); # I don't need this
#				}
#		print "$subseqstr";
#	
#			if ($i < $dataLast){
#				print 'N' x 1000;
#			}
#		}
#		print "\n";
#	}
#	else {
#		my $length = $db->length($data[0]);
#		my $seqstr = $db->seq($data[0],1,$length);
#		print "$seqstr\n";
#	}
}
#close($gff_in) or die "Could not close filehandle '$gff_in' $!";
close FILE;

exit;
