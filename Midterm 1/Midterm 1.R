#loading biostrings and msa to help with this lab.

library(Biostrings)
library(msa)
#set the working directory
setwd("C:/Users/domin/OneDrive/Documents/GitHub/Bioinformatics")
#Import DNA sequences from the FASTA file
dna_sequences <- readDNAStringSet("sequences.fasta")
#checking how many sequences were imported
length(dna_sequences)
#perform a multiple sequence alignment
alignment <- msa(dna_sequences)
#check aligned sequences 
alignment
#sequences imported and aligned
aligned_sequences <- DNAStringSet(alignment) 
#generate a consensus matrix to measure how good alignment is
consensus <- consensusMatrix(aligned_sequences)
consensus
#sequences are aligned well
#calculate consensus sequence
#remove gaps, was giving warning when I tried it another way
aligned_no_gaps <- gsub("-", "", as.character(aligned_sequences))
aligned_no_gaps <- DNAStringSet(aligned_no_gaps)
#try to calculate consensus sequence again
cons_seq <- consensusString(aligned_no_gaps)
cons_seq
#this worked much better
#calculate GC content 
gc_per_seq <- letterFrequency(aligned_no_gaps, letters = c("G", "C"), as.prob = TRUE)
gc_content <- rowSums(gc_per_seq) * 100
gc_content
#convert alignment to a character matrix 
alignment_matrix <- as.matrix(alignment)
dim(alignment_matrix)
#identify polymorphic position
polymorphic_cols <- apply(alignment_matrix, 2, function(col) length(unique(col)) > 1)
which(polymorphic_cols)
#identify different sequences
variable_positions <- alignment_matrix[, polymorphic_cols]
variable_positions
#homo sapiens 6 seems to have mutations such as SNPs.
#exporting to BLAST to figure out question 6
writeXStringSet(aligned_sequences, filepath = "my_sequences.fasta")
# i got this: Select seq LC121775.1	Homo sapiens hbb gene for beta globin, partial cds, note: HbLimassol Cd8(AAG>AAC) and this was the accession number LC121775.1
# i already figure out before that it was homo sapiens 6 that was the most different
#extracting homo sapiens 6
alignment_matrix <- as.matrix(alignment)
rownames(alignment_matrix)
most_diff_seq <- alignment_matrix["Homo_sapiens_6", ]
most_diff_DNA <- DNAString(paste(most_diff_seq, collapse = ""))
most_diff_DNA <- DNAString(most_diff_seq_str)
#translate DNA into protein
protein_seq <- translate(most_diff_DNA)
#create FASTA file 
writeXStringSet(AAStringSet(protein_seq), filepath = paste0(protein_seq, "protein.fasta"))
#accession number for number  is KAI2558340.1
#the gene is associated with Sickle Cell Anemia for question 9