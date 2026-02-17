if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

install.packages("seqinr")

library(Biostrings)
library(seqinr)

BiocManager::install("msa")

library(msa)

seq_1 <- readDNAStringSet("Hs10.fasta")
seq_2 <- readDNAStringSet("Hs1.fasta")
seq_3 <- readDNAStringSet("Hs128.fasta")
seqs_4 <- readDNAStringSet("Hs7.fasta")
seq_5 <- readDNAStringSet("Hs104.fasta")
sequences <- c(seq_1, seq_2, seq_3, seqs_4, seq_5)
sequences
alignment <- msa(sequences, method = "Muscle")

print(alignment, show="complete")
alignment_matrix <- as.matrix(alignment)

gc_counts <- apply(alignment_matrix, 1, function(seq) sum(seq == "G" | seq == "C"))


seq_length <- apply(alignment_matrix, 1, length)


gc_content <- gc_counts / seq_length * 100
gc_content
library(seqinr)

seqs_seqinr <- msaConvert(alignment, type = "seqinr::alignment")

dist_matrix <- dist.alignment(seqs_seqinr, "identity")  # percent identity
dist_matrix

as.matrix(dist_matrix)

library(Biostrings)

dna_seq <- sequences[[1]]       
protein_seq <- translate(dna_seq)
protein_seq

install.packages("phangorn")


library(phangorn)
Alignment_phyDat <- msaConvert(alignment, type="phangorn::phyDat")
write.phyDat(Alignment_phyDat, file="alignment.fasta", format="fasta")
getwd()
list.files()




