if(!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}

library(Biostrings)

alignment <- readDNAStringSet("metazoa_alignment.gene.fasta")
names(alignment)

hs_index <- grep("Homo_sapiens", names(alignment))
hs_sequence <- alignment[hs_index]

hs_sequence_clean <- gsub("-", "", as.character(hs_sequence))
hs_sequence_clean <- DNAString(hs_sequence_clean)

hs_protein <- translate(hs_sequence_clean)

hs_protein_set <- AAStringSet(hs_protein)
names(hs_protein_set) <- "Homo_sapiens_protein"
writeXStringSet(hs_protein_set, "Homo_sapiens_protein.fasta")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("biomaRt", "ggplot2"))

library(biomaRt)
library(ggplot2)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

go_results <- getBM(attributes = c('uniprotswissprot', 'go_id', 'name_1006', 'namespace_1003'),
                    filters = 'uniprotswissprot',
                    values = 'P54098',
                    mart = mart)

go_results <- go_results[go_results$go_id != "", ]

bp_term <- go_results[go_results$namespace_1003 == "biological_process", "name_1006"][1]
mf_term <- go_results[go_results$namespace_1003 == "molecular_function", "name_1006"][1]
cc_term <- go_results[go_results$namespace_1003 == "cellular_component", "name_1006"][1]

print(bp_term)
print(mf_term)
print(cc_term)

ggplot(go_results, aes(x = namespace_1003, fill = namespace_1003)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "GO Sub-Ontology Distribution for POLG",
       subtitle = "UniProt Accession: P54098",
       x = "Sub-Ontology Category",
       y = "Number of Associated GO Terms") +
  theme_minimal()

ggsave("POLG_GO_Plot.png", width = 8, height = 6, dpi = 300)