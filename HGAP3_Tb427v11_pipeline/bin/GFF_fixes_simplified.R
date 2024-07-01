# To create a valid gff, gene and mRNA features need to be updated to include the UTRs

library(rtracklayer)

# Load UTRme output 'best score' and the original genome gff, which was used to run UTRme

#setwd("~/LabProjects/2022-03-01-Test-generation-of-UTR-containing-annotation-file/")
args <- commandArgs(trailingOnly = TRUE)
five_gff <- import(args[1],format="GFF3")
three_gff <- import(args[2],format="GFF3")
original_gff <- import(args[3],format="GFF3")

# five_gff <- import("UTRme-BSF_Run-5-best-score_subset.gff")
# three_gff <- import("UTRme-BSF_Run-3-best-score_subset.gff",format="GFF3")
# original_gff <- import("~/Genomes/Tbrucei/Lister427/HGAP3_Tb427v10/HGAP3_Tb427v10.gff3",format="GFF3")


#this gets the index numbers of the entries in the three prime gff that holds the three_prime_UTR annotation
three_primes <- which(three_gff$type == "three_prime_UTR")
id <- ""
print(paste("Starting the gff correction for", length(three_primes), "3'UTRs...please wait"))

# this loops through the three_prime_UTR indexes
for (utr in three_primes) {
  # get the base gene_id , store it in id, if it was not already stored there, if it was move to the next utr (?when that happen? if there are more than one UTR annotation for a gene, not it my case where i am using only the "best")

  id <- strsplit(strsplit(unlist(three_gff[utr]$Parent), ":")[[1]][1], ".", fixed = T)[[1]][1]
  # if it positive strand
  if (runValue(strand(three_gff[utr]) == "+")) {
    # get all the features related to the same gene (have the same value in ID as id) -->couldn't i get undesired hits with this grep?
    related_features <- grep(id, three_gff$ID)
    #stores the the position of the longest UTR predicted --> not necessary for the "best" because there is only one
    utr_end <- max(end(three_gff[related_features[which(three_gff[related_features]$type == "three_prime_UTR")]]))

    # get the indexes of the entries with the same id in the original gff
    original_features <- grep(id, original_gff$ID)
    # asigns to the end coordinate of the longest UTR to the entries sharing the same id that are feature "gene" or "mRNA" 
    end(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA" , "pseudogene", "pseudogenic_transcript"))]]) <- utr_end
  } else {
    related_features <- grep(id, three_gff$ID)
    utr_start <- min(start(three_gff[related_features[which(three_gff[related_features]$type == "three_prime_UTR")]]))

    original_features <- grep(id, original_gff$ID)

    start(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA" , "pseudogene", "pseudogenic_transcript"))]]) <- utr_start

  }
}
print("...done!")

five_primes <- which(five_gff$type == "five_prime_UTR")
id <- ""
print(paste("Starting the gff correction for", length(five_primes), "5'UTRs...please wait"))
for (utr in five_primes) {
  # get the base gene_id
  if (id == strsplit(strsplit(unlist(five_gff[utr]$Parent), ":")[[1]][1], ".", fixed = T)[[1]][1]) next
  id <- strsplit(strsplit(unlist(five_gff[utr]$Parent), ":")[[1]][1], ".", fixed = T)[[1]][1]
  if (runValue(strand(five_gff[utr]) == "-")) {
    related_features <- grep(id, five_gff$ID)
    utr_end <- max(end(five_gff[related_features[which(five_gff[related_features]$type == "five_prime_UTR")]]))

    original_features <- grep(id, original_gff$ID)

    end(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA" , "pseudogene", "pseudogenic_transcript"))]]) <- utr_end

  } else {
    related_features <- grep(id, five_gff$ID)
    utr_start <- min(start(five_gff[related_features[which(five_gff[related_features]$type == "five_prime_UTR")]]))

    original_features <- grep(id, original_gff$ID)

    start(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA" , "pseudogene", "pseudogenic_transcript"))]]) <- utr_start

  }
}

print(paste("Exporting the new gff file"))

export(original_gff, args[4], format="GFF3")
#export(original_gff, "test.gff", format="GFF3")
print("...all done!")

