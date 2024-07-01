#!/bin/perl
# usage = perl Blast_output_to_gff3_v2.pl Blast_output_file VSG-multifasta-file.fasta > output.gff3
# Input files are the multifasta file containing all the CDS of the VSGs identified by Cross GA et al. (2014) in fasta format and the BLAST output file from these VSG list to our genome, generated with BLAST+ using blastn with "-outfmt 7" 
use strict;
use warnings;
use Bio::DB::Fasta;

my $datafile = $ARGV[0];
my $fastafile = $ARGV[1];
my $start;
my $end;
my @subject_ids;
my $j;
my $order;
my $db = Bio::DB::Fasta->new($fastafile);
my @secondary_hits;
open(FILE,$datafile);
print "##gff-version 3\n";
while(<FILE>){
chomp;

my ($queryId, $subjectId, $percIdentity, $alnLength, $mismatchCount, $gapOpenCount, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $eVal, $bitScore) = split("\t");

# skip the comment lines at the beginning of the BLAST output file
if($queryId =~ /^#/ ) {
     next;
     } 
my $seq_length = $db->length($queryId);
my $seq_header = $db->header($queryId);



# select those hits that have high percentage of identity (more than 95%) and query coverage (more than 90%)
if(($percIdentity > 95) && ($alnLength > ($seq_length)*90/100)){


# Trying to update the ordering of the different hits that go to the same VSG. So, I would like to order the hits to the same VSG that go above the threshold by their alignment score, because only the first is always the best one (maybe not the only best one) but for the next this no guaranteed because then the hits (that are above the threshold) that belong to the same chromosome as the best hit are shown.
# Print the first one. For the followings, if the query is equal to the ones before store the important values in an array of arrays (?). do not print. if the next entry has a different queryId, order the array of arrays by the e-value, add a variable to the array based in the position in the array (starting from 2) , and print them in order with the information of this extra variable. Then process the actual on

#---------------- consider the coding strand ------------------------

my $flow = "+";
$start = $subjectStart;
$end = $subjectEnd;
if($subjectStart > $subjectEnd){
$flow = "-";
$start = $subjectEnd;
$end = $subjectStart;
}
#--------------- consider the coding strand ---------------------

# A loop to see if this is the first positive match to a VSG from the list or if it is a secondary hit. If it is a secondary hit, store the important variables in an array. And store this array in another array containing all the secondary hits for this VSG
if(grep {$_ eq $queryId} @subject_ids){

my @secondary_hit = ($queryId, $subjectId, $percIdentity, $alnLength, $eVal, $bitScore, $seq_length, $seq_header, $start, $end, $flow);

push @secondary_hits, [ @secondary_hit ];

}
	else {# If it is the first positive hit to a VSG from the list, add this VSG to the list to test first hit. This also means that the hits to the previous VSG are already finished, so now we can sort them by Bitscore and print them, and also print the result line of the best hit for this new VSG. 
# add the queryId to the @subject_ids array.
unshift @subject_ids, $queryId; #add the VSG ($queryId) to the list of already seen VSGs (@subject_ids)
# order and print the previous secondary hits
my @ordered_secondary_hits = sort { $b->[5] <=> $a->[5] } @secondary_hits; # order the secondary hits of the previous VSG by decreasing bitScore.


for ($j=0; $j<=$#ordered_secondary_hits; $j++) {# for loop to go through the ordered list of secondary hits for the previous VSG and print them. If $ordered_secondary_hits is empty (like it will happen in the first positive hit of the first VSG, or after a VSG that has only one positive hit, it is not going inside this loop.

	if ($j == 0){

	$order = "the second most similar gene";
	}

	elsif ($j == 1){

	$order = "the third most similar gene";
	}

	elsif ($j == 2){

	$order = "the forth most similar gene";
	}

	else {

	$order = "not one of the most similar genes";
	}

print "$ordered_secondary_hits[$j][1]\tBLAST\tgene\t$ordered_secondary_hits[$j][8]\t$ordered_secondary_hits[$j][9]\t$ordered_secondary_hits[$j][4]\t$ordered_secondary_hits[$j][10]\t.\tName=Similar to $ordered_secondary_hits[$j][0];Note=From $ordered_secondary_hits[$j][8] to $ordered_secondary_hits[$j][9], $ordered_secondary_hits[$j][2] percent identity to VSG $ordered_secondary_hits[$j][0]. The alignment length is $ordered_secondary_hits[$j][3] nt, while the query VSG gene sequence is $ordered_secondary_hits[$j][6] nt long. This is $order to $ordered_secondary_hits[$j][0] ($ordered_secondary_hits[$j][7]) in the genome.\n"; 

} # order and print to secondary hits of the previous VSG


undef(@secondary_hits);# returns @secondary_hits to be empty.

# print the resulting line of the best hit for the present VSG.
print "$subjectId\tBLAST\tgene\t$start\t$end\t$eVal\t$flow\t.\tName=Similar to $queryId;Note=From $start to $end, $percIdentity percent identity to VSG $queryId. The alignment length is $alnLength nt, while the query VSG gene sequence is $seq_length nt long. This is the most similar gene to $queryId ($seq_header) in the genome.\n"; 

}

} # close the if related to what to do if it goes above the threshold to be considered a hit
} # close the while
close(FILE);
exit;
