#!/bin/perl
# usage = perl Add_IDs_to_VSGs.pl file.gff GeneID_prefix 1st_gene_ID_to_use > file_with_IDs_in_VSGs.gff
# Input file is the gff file with VSG annotated based on BLAST results without geneID, the geneID_prefix (e.g. "Tb427v9_") and the first number to use in the gene ID (e.g. "001844900")
# Output is the gff file with the VSGs annotated based on BLAST results, now having a gene ID.

use strict;
use warnings;

my $Input_gff = $ARGV[0];
my $Prefix = $ARGV[1];
my $ID = $ARGV[2];
my $Suffix =":pseudogene";

open(my $in, '<', $Input_gff) or die "Could not open file '$Input_gff' $!";

while(my $in_row = <$in>){
chomp $in_row;
	if($in_row =~m/^#/ ){
		print "$in_row\n";
	}
	else{
	my ($Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes) = split("\t",$in_row);
		if($Source eq "BLAST"){
			# Add the GeneID to the Attributes.
		$Attributes = "ID=".$Prefix.$ID.$Suffix.";$Attributes";
		print "$Chrom\t$Source\tpseudogene\t$Start\t$End\t$Score\t$Strand\t$Frame\t$Attributes\n";
		# ID=Tb427_000006000.1;Parent=Tb427_000006000:pseudogene
		my $trans_attributes = "ID=".$Prefix.$ID.".1".";Parent=$Prefix".$ID.$Suffix;
		print "$Chrom\t$Source\tpseudogenic_transcript\t$Start\t$End\t$Score\t$Strand\t$Frame\t$trans_attributes\n";

		$ID= $ID + 100;

		}
		else{
		print "$Chrom\t$Source\t$Feature\t$Start\t$End\t$Score\t$Strand\t$Frame\t$Attributes\n";
		}
	}
	
}
close($in) or die "Could not close filehandle '$in' $!";
exit;

