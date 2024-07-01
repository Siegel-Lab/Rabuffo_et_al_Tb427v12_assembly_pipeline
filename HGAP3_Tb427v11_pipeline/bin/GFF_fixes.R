# To create a valid gff, gene and mRNA features need to be updated to include the UTRs

library(rtracklayer)

# Load UTRme output 'best score' and the original genome gff, which was used to run UTRme

args <- commandArgs(trailingOnly = TRUE)
five_gff <- import(args[1])
three_gff <- import(args[2])
original_gff <- import(args[3])

levels(original_gff$source) <- unique(c(levels(three_gff$source), levels(five_gff$source), levels(original_gff$source)))
levels(original_gff$type) <- unique(c(levels(three_gff$type), levels(five_gff$type), levels(original_gff$type)))

three_primes <- which(three_gff$type == "three_prime_UTR")
id <- ""
print(paste("Starting the gff correction for", length(three_primes), "3'UTRs...please wait"))
for (utr in three_primes) {
  # get the base gene_id
  if (id == strsplit(strsplit(unlist(three_gff[utr]$Parent), ":")[[1]][1], ".", fixed = T)[[1]][1]) next
  id <- strsplit(strsplit(unlist(three_gff[utr]$Parent), ":")[[1]][1], ".", fixed = T)[[1]][1]
  if (runValue(strand(three_gff[utr]) == "+")) {
    related_features <- grep(id, three_gff$ID)
    utr_end <- max(end(three_gff[related_features[which(three_gff[related_features]$type == "three_prime_UTR")]]))
    longest_UTR <- which(end(three_gff[related_features]) == utr_end)
    if (length(longest_UTR) > 1) next
    original_features <- grep(id, original_gff$ID)

#    end(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA"))]]) <- utr_end
    end(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA", "pseudogene", "pseudogenic_transcript"))]]) <- utr_end
    old_three <- grep("three_prime_UTR", original_gff[original_features]$type)
    if (length(old_three) > 0) {
      original_gff[original_features[old_three]] <- three_gff[longest_UTR]
    } else {
      pos <- original_features[length(original_features)]+1
      tmp <- original_gff[pos:length(original_gff)]
      original_gff[pos] <- three_gff[related_features[longest_UTR]]
      original_gff <- c(original_gff[1:pos], tmp)
    }
  } else {
    related_features <- grep(id, three_gff$ID)
    utr_start <- min(start(three_gff[related_features[which(three_gff[related_features]$type == "three_prime_UTR")]]))
    longest_UTR <- which(start(three_gff[related_features]) == utr_start)
    if (length(longest_UTR) > 1) next
    original_features <- grep(id, original_gff$ID)

#    start(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA"))]]) <- utr_start
    start(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA", "pseudogene", "pseudogenic_transcript"))]]) <- utr_start
    old_three <- grep("three_prime_UTR", original_gff[original_features]$type)
    if (length(old_three) > 0) {
      original_gff[original_features[old_three]] <- three_gff[longest_UTR]
    } else {
      pos <- original_features[length(original_features)]+1
      tmp <- original_gff[pos:length(original_gff)]
      original_gff[pos] <- three_gff[related_features[longest_UTR]]
      original_gff <- c(original_gff[1:pos], tmp)
    }
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
    longest_UTR <- which(end(five_gff[related_features]) == utr_end)
    if (length(longest_UTR) > 1) next
    original_features <- grep(id, original_gff$ID)
    
#    end(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA"))]]) <- utr_end
    end(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA", "pseudogene", "pseudogenic_transcript"))]]) <- utr_end
    old_five <- grep("five_prime_UTR", original_gff[original_features]$type)
    if (length(old_five) > 0) {
      original_gff[original_features[old_five]] <- five_gff[longest_UTR]
    } else {
      pos <- original_features[length(original_features)]+1
      tmp <- original_gff[pos:length(original_gff)]
      original_gff[pos] <- five_gff[related_features[longest_UTR]]
      original_gff <- c(original_gff[1:pos], tmp)
    }
  } else {
    related_features <- grep(id, five_gff$ID)
    utr_start <- min(start(five_gff[related_features[which(five_gff[related_features]$type == "five_prime_UTR")]]))
    longest_UTR <- which(start(five_gff[related_features]) == utr_start)
    if (length(longest_UTR) > 1) next
    original_features <- grep(id, original_gff$ID)
    
#    start(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA"))]]) <- utr_start
    start(original_gff[original_features[which(original_gff[original_features]$type %in% c("gene" , "mRNA", "pseudogene", "pseudogenic_transcript"))]]) <- utr_start
    old_five <- grep("five_prime_UTR", original_gff[original_features]$type)
    if (length(old_five) > 0) {
      original_gff[original_features[old_five]] <- five_gff[longest_UTR]
    } else {
      pos <- original_features[length(original_features)]+1
      tmp <- original_gff[pos:length(original_gff)]
      original_gff[pos] <- five_gff[related_features[longest_UTR]]
      original_gff <- c(original_gff[1:pos], tmp)
    }
  }
}

print(paste("Exporting the new gff file"))

export(original_gff, args[4])

print("...all done!")

