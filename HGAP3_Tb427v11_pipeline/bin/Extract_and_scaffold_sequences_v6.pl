#!/bin/perl
# this script goes through a list with "headers" tab separated (those that should be joined) and an input multifasta file. It scaffolds with 1000 Ns in between those sequences from the headers tab separated.
# usage = perl input.fasta data.txt > output.fasta

# I am supposing the presence of 22 alleles (11 chromosomes), and that the chromosomes in the guide file are ordered in decreasing number (Chr11, Chr10, Chr9, etc..) and that the "A" allele is first, and then the "B" allele

use strict;
use Bio::DB::Fasta;

my $fastafile = $ARGV[0];
my $datafile = $ARGV[1];

my $db = Bio::DB::Fasta->new($fastafile);

open(FILE,$datafile);

my $Allele = 22;

my $Char;

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
		print ">Chr".$Chrom."_hap".$Char;


	print "_Tb427v11\n";

	}

	else { # after the eleven scaffolded chromosomes, print the rest unscaffolded forks
	print ">$data2[0]\n";
	}

# Until here is to generate to combined header.

	if ($Allele > 0){
		$Allele -=1;

		for (my $i = 0; $i<= $dataLast; $i++){
			my $subseqstr;
			my $length = $db->length($data[$i]);
				if ($i == 0 && $data[$i] =~/^BES/){
					$subseqstr = $db->seq($data[$i],$length,1);
				}# Reverse-complement the BES that go on the 5'
				else {
					$subseqstr = $db->seq($data[$i],1,$length);
				}
		print "$subseqstr";
	
			if ($i < $dataLast){
				print 'N' x 1000;
			}
		}
		print "\n";
	}
	else {
		my $length = $db->length($data[0]);
		my $seqstr = $db->seq($data[0],1,$length);
		print "$seqstr\n";
	}
}

close FILE;

exit;
