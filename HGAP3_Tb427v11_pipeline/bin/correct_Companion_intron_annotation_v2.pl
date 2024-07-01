#!/bin/perl
# usage = perl correct_Companion_intron_annotation.pl input.gff3 ID_list_file.txt(a text file with the IDs of the entries with introns, one ID per line) > output.gff3
#
# Input file is the gff3 file from Companion ( or Companion-derived) for which wrong assigned introns need to be corrected
# The output is a gff3 file were the intron containing entries are corrected and the gene converted in pseudogene
#
#
# Example on how the annotation looks and how it should look
################################################################################
##WRONG INTRON EXAMPLES
# v10
#Chr9_core_Tb427v10	.	polypeptide	867298	867852	.	+	.	ID=Tb427_090038100.1:pep;Derives_from=Tb427_090038100.1;orthologous_to=Tb927.9.6830;ortholog_cluster=ORTHOMCL5130;product=term%3Dhypothetical protein%2C conserved
#Chr9_core_Tb427v10	RATT	CDS	867298	867504	.	+	0	ID=Tb427_090038100.1:CDS:1;Parent=Tb427_090038100.1
#Chr9_core_Tb427v10	RATT	mRNA	867298	867852	.	+	.	ID=Tb427_090038100.1;Parent=Tb427_090038100
#Chr9_core_Tb427v10	RATT	CDS	867508	867852	.	+	0	ID=Tb427_090038100.1:CDS:2;Parent=Tb427_090038100.1
#Chr9_core_Tb427v10	RATT	gene	867298	867852	.	+	.	ID=Tb427_090038100;ratt_ortholog=Tb927.9.6830
#
#v9.2
#
#Chr10_core_Tb427v9	RATT	gene	67142	69209	.	+	.	ID=Tb427v9_000116500;ratt_ortholog=Tb927.10.530
#Chr10_core_Tb427v9	RATT	mRNA	67142	69209	.	+	.	ID=Tb427v9_000116500.1;Parent=Tb427v9_000116500
#Chr10_core_Tb427v9	RATT	CDS	67142	67534	.	+	0	ID=Tb427v9_000116500.1:CDa:1;Parent=Tb427v9_000116500.1
#Chr10_core_Tb427v9	RATT	CDS	67539	69209	.	+	0	ID=Tb427v9_000116500.1:CDS:2;Parent=Tb427v9_000116500.1
#Chr10_core_Tb427v9	.	polypeptide	67142	69209	.	+	.	ID=Tb427v9_000116500.1:pep;Parent=Tb427v9_000116500.1:CDS:2;Derives_from=Tb427v9_000116500.1;orthologous_to=Tb927.10.530;ortholog_cluster=ORTHOMCL8743;product=term%3DPlasma-membrane choline transporter%2C putative%3Bevidence%3DIEA%3Bwith%3DGeneDB:Tb927.10.530;Dbxref=UniProtKB:Q38CI7
#Chr10_core_Tb427v9	Pfam	protein_match	68030	69088	.	+	.	Name=PF04515;Parent=Tb427v9_000116500.1:pep;signature_desc=Plasma-membrane choline transporter
#
#
##PSEUDOGENE EXAMPLE
#
#v10
#
#Chr9_core_Tb427v10	.	polypeptide	946971	947985	.	-	.	ID=Tb427_090041900.1:pep;Derives_from=Tb427_090041900.1;orthologous_to=Tb427_090041500.1;ortholog_cluster=ORTHOMCL5100;product=term%3DESAG protein%2C putative%3Bevidence%3DIEA%3Bwith%3DPfam:PF03238
#Chr9_core_Tb427v10	.	pseudogenic_exon	946971	947222	.	-	.	ID=Tb427_090041900.1:pseudogenic_exon:1;Parent=Tb427_090041900.1
#Chr9_core_Tb427v10	.	pseudogenic_transcript	946971	947985	.	-	.	ID=Tb427_090041900.1;Parent=Tb427_090041900:pseudogene
#Chr9_core_Tb427v10	Pfam	protein_match	946975	947652	.	-	.	Parent=Tb427_090041900.1:pep;Name=PF03238;signature_desc=ESAG protein
#Chr9_core_Tb427v10	.	pseudogenic_exon	947224	947985	.	-	.	ID=Tb427_090041900.1:pseudogenic_exon:2;Parent=Tb427_090041900.1
#Chr9_core_Tb427v10	.	pseudogene	946971	947985	861	-	.	ID=Tb427_090041900:pseudogene;has_internal_stop=true;has_frameshift=true;Target=Tb927.1.5120:mRNA 1 325;original_prot_length=328
#
#v9.2
#
#Chr10_core_Tb427v9	.	pseudogene	113830	117082	6.04e+03	+	.	ID=Tb427v9_000118400:pseudogene;has_internal_stop=false;has_frameshift=true;has_start=true;has_stop=true;Target=Tb927.10.720:mRNA 1 1084;original_prot_length=1084
#Chr10_core_Tb427v9	.	pseudogenic_transcript	113830	117082	.	+	.	ID=Tb427v9_000118400.1;Parent=Tb427v9_000118400:pseudogene
#Chr10_core_Tb427v9	.	pseudogenic_exon	113830	115794	.	+	.	ID=Tb427v9_000118400.1:pseudogenic_exon:1;Parent=Tb427v9_000118400.1
#Chr10_core_Tb427v9	.	pseudogenic_exon	115796	117082	.	+	.	ID=Tb427v9_000118400.1:pseudogenic_exon:2;Parent=Tb427v9_000118400.1
#Chr10_core_Tb427v9	.	polypeptide	113830	117082	.	+	.	ID=Tb427v9_000118400.1:pep;Derives_from=Tb427v9_000118400.1;orthologous_to=Tb927.10.720;ortholog_cluster=ORTHOMCL8728;product=term%3DPB1 domain containing protein%2C putative%3Bevidence%3DIEA%3Bwith%3DGeneDB:Tb927.10.720;Dbxref=UniProtKB:Q38CG8;Ontology_term=GO:0005515
#Chr10_core_Tb427v9	Pfam	protein_match	115663	115923	.	+	.	Name=PF00564;Parent=Tb427v9_000118400.1:pep;Ontology_term=GO:0005515;signature_desc=PB1 domain
#
##CORRECTED EXAMPLES
#
#v10
#
#Chr9_core_Tb427v10	.	polypeptide	867298	867852	.	+	.	ID=Tb427_090038100.1:pep;Derives_from=Tb427_090038100.1;orthologous_to=Tb927.9.6830;ortholog_cluster=ORTHOMCL5130;product=term%3Dhypothetical protein%2C conserved
#Chr9_core_Tb427v10	RATT	pseudogenic_exon	867298	867504	.	+	0	ID=Tb427_090038100.1:pseudogenic_exon:1;Parent=Tb427_090038100.1
#Chr9_core_Tb427v10	RATT	pseudogenic_exon	867508	867852	.	+	0	ID=Tb427_090038100.1:pseudogenic_exon:2;Parent=Tb427_090038100.1
#Chr9_core_Tb427v10	RATT	pseudogenic_transcript	867298	867852	.	+	.	ID=Tb427_090038100.1;Parent=Tb427_090038100:pseudogene
#Chr9_core_Tb427v10	RATT	pseudogene	867298	867852	.	+	.	ID=Tb427_090038100:pseudogene;has_internal_stop=true;has_frameshift=true;ratt_ortholog=Tb927.9.6830
#
#
#v9.2
#
#Chr10_core_Tb427v9	.	polypeptide	67142	69209	.	+	.	ID=Tb427v9_000116500.1:pep;Derives_from=Tb427v9_000116500.1;orthologous_to=Tb927.10.530;ortholog_cluster=ORTHOMCL8743;product=term%3DPlasma-membrane choline transporter%2C putative%3Bevidence%3DIEA%3Bwith%3DGeneDB:Tb927.10.530;Dbxref=UniProtKB:Q38CI7
#Chr10_core_Tb427v9	RATT	pseudogenic_exon	67142	67534	.	+	0	ID=Tb427v9_000116500.1:pseudogenic_exon:1;Parent=Tb427v9_000116500.1
#Chr10_core_Tb427v9	RATT	pseudogenic_exon	67539	69209	.	+	0	ID=Tb427v9_000116500.1:pseudogenic_exon:2;Parent=Tb427v9_000116500.1
#Chr10_core_Tb427v9	RATT	pseudogenic_transcript	67142	69209	.	+	.	ID=Tb427v9_000116500.1;Parent=Tb427v9_000116500:pseudogene
#Chr10_core_Tb427v9	RATT	pseudogene	67142	69209	.	+	.	ID=Tb427v9_000116500:pseudogene;has_internal_stop=true;has_frameshift=true;ratt_ortholog=Tb927.10.530
#
#
##THINGS TO DO TO CORRECT
#
# polypeptide entry --> For v9.2 Remove "Parent=*CDS:[Num];" 
# # CDS entries --> CHANGE "\tRATT\t" to "." on Tool, CHANGE "\tCDS\t" to "\tpseudogenic_exon\t" on Type, change ":CDS:" to ":pseudogenic_exon:"
# # mRNA entry --> CHANGE "\tRATT\t" to "." on Tool, CHANGE "mRNA" to "pseudogenic_transcript" and add ":pseudogene" at the end of the Note column
# # gene entry --> CHANGE "\tRATT\t" to "." on Tool, CHANGE "\tgene\t" to "\tpseudogene\t" and change ";ratt_ortholog=" by ":pseudogene;has_internal_stop=true;has_frameshift=true;Target=" 
#
################################################################################

use strict;
use warnings;

my $input_gff_file = $ARGV[0];
my $ID_list_file = $ARGV[1];

$" = "\t"; # Setting the output field separator to the tab character.

#Storing the ID_list into an array
open(my $ID_list_in, '<', $ID_list_file) or die "Could not open file '$ID_list_file' $!";

chomp(my @IDs = <$ID_list_in>);
#print "@IDs\n";
close $ID_list_in;


open(my $input_gff_in, '<', $input_gff_file) or die "Could not open file '$input_gff_file' $!";

while(my $input_gff_row = <$input_gff_in>){

chomp $input_gff_row;

#Counter
my $i=0;

	# Go throught the array with IDs and check if any one is a substring of the line in the gff file
	foreach my $ID (@IDs) {
		if ( grep ( /$ID/, $input_gff_row ) ) {
		#add one to the match counter
		$i++;
	
		# Split the line 
		my ($Chrom, $Source, $Feature, $Start, $End, $Score, $Strand, $Frame, $Attributes) = split("\t",$input_gff_row);

		# Actions to perform if there is a line match
#		$Source =~ s/RATT/./;

			if($Feature eq "polypeptide"){
				$Attributes =~ s/Parent=.*:CDS:[0-9];//;
			
			}

			if($Feature eq "CDS"){
			$Feature =~ s/CDS/pseudogenic_exon/;
			$Attributes =~ s/:CDS:/:pseudogenic_exon:/;
			}
		
			if($Feature eq "mRNA"){
			$Feature =~ s/mRNA/pseudogenic_transcript/;
			$Attributes .= ":pseudogene";

			}
				

			if($Feature eq "gene"){
			$Feature =~ s/gene/pseudogene/;
			$Attributes =~ s/;ratt_ortholog=/:pseudogene;has_internal_stop=true;has_frameshift=true;ratt_ortholog=/;
		
			}

		print "$Chrom\t$Source\t$Feature\t$Start\t$End\t$Score\t$Strand\t$Frame\t$Attributes\n";

		}
	}
	# If there was no match, print the line as it was
	if ($i == 0) {
	#If there is no match, print the line as it is
	print "$input_gff_row\n";

	}
}

close($input_gff_in) or die "Could not close filehandle '$input_gff_in' $!";
exit;

