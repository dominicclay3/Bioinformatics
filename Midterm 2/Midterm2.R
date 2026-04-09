#load in biostrings
library(Biostrings)
#set working directory
setwd("C:/Users/domin/OneDrive/Documents/GitHub/Bioinformatics/Midterm 2")
#load in file
alignment <- readDNAStringSet("metazoa_alignment.gene.fasta")
#view names of sequences
names(alignment)
#extract homo sapiens sequence
hs_index <- grep("Homo_sapiens", names(alignment))
hs_sequence <- alignment[hs_index]
#remove gaps from alignment
hs_sequence_clean <- gsub("-", "", as.character(hs_sequence))
hs_sequence_clean <- DNAString(hs_sequence_clean)
#translate to protein
hs_protein <- translate(hs_sequence_clean)
hs_protein
hs_protein_set <- AAStringSet(hs_protein)
names(hs_protein_set) <- "Homo_sapiens_protein"
#save to fasta file
writeXStringSet(hs_protein_set, "Homo_sapiens_protein.fasta")
#Install some packages
if (!require("biomaRt")) install.packages("biomaRt")
if (!require("ggplot2")) install.packages("ggplot2")
library(biomaRt)
library(ggplot2)
#Connect Uniprot mapping to Ensembl

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#get GO terms 
go_results <- getBM(attributes = c('uniprotswissprot', 'go_id', 'name_1006', 'namespace_1003'),
                    filters = 'uniprotswissprot',
                    values = 'P54098',
                    mart = mart)
#Clean the data
go_results <- go_results[go_results$go_id != "", ]
#extract one term from each sub ontology
bp_term <- go_results[go_results$namespace_1003 == "biological_process", "name_1006"][1]
mf_term <- go_results[go_results$namespace_1003 == "molecular_function", "name_1006"][1]
cc_term <- go_results[go_results$namespace_1003 == "cellular_component", "name_1006"][1]
#print terms to console
cat("Biological Process:", bp_term, "\n")
cat("Molecular Function:", mf_term, "\n")
cat("Cellular Component:", cc_term, "\n")
#plot results
ggplot(go_results, aes(x = namespace_1003, fill = namespace_1003)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "GO Sub-Ontology Distribution for POLG",
       subtitle = "UniProt Accession: P54098",
       x = "Sub-Ontology Category",
       y = "Number of Associated GO Terms") +
  theme_minimal()
#save plot
ggsave("POLG_GO_Analysis.png", width = 8, height = 6, dpi = 300)
