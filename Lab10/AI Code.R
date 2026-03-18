
library(Biostrings)
library(UniprotR)
library(protti)
library(dplyr)
library(ggplot2)
library(r3dmol)

my_protein <- AAString("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ")
my_protein_named <- AAStringSet(my_protein)
names(my_protein_named) <- "MyProtein1"
writeXStringSet(my_protein_named, "my_protein.fasta")

accessions <- c("Q57153","Q3I4Z3","B7XSQ5","Q45231","Q27JV1")
writeLines(accessions, "top5_accessions.txt")
accessions <- readLines("top5_accessions.txt")
accessions <- accessions[accessions != ""]
accessions_string <- paste(accessions, collapse = ",")

go_terms <- data.frame(
  Accession = c(
    "Q57153","Q57153","Q57153",
    "Q3I4Z3","Q3I4Z3","Q3I4Z3",
    "B7XSQ5","B7XSQ5","B7XSQ5",
    "Q45231","Q45231","Q45231",
    "Q27JV1","Q27JV1","Q27JV1"
  ),
  GO.ID = c(
    "GO:0003677","GO:0005488","GO:0008150",
    "GO:0003677","GO:0005488","GO:0008150",
    "GO:0003677","GO:0005488","GO:0008150",
    "GO:0003677","GO:0005488","GO:0008150",
    "GO:0004672","GO:0005524","GO:0005886"
  ),
  Term = c(
    "DNA binding","protein binding","biological_process",
    "DNA binding","protein binding","biological_process",
    "DNA binding","protein binding","biological_process",
    "DNA binding","protein binding","biological_process",
    "protein kinase activity","ATP binding","plasma membrane"
  ),
  Ontology = c(
    "Molecular Function","Molecular Function","Biological Process",
    "Molecular Function","Molecular Function","Biological Process",
    "Molecular Function","Molecular Function","Biological Process",
    "Molecular Function","Molecular Function","Biological Process",
    "Molecular Function","Molecular Function","Cellular Component"
  ),
  stringsAsFactors = FALSE
)

PlotGoInfo(go_terms)

protein_disease <- data.frame(
  Accession = accessions,
  Example_Pathology = c(
    "Hemoglobin-related disorders",
    "Myosin-associated muscle disorders",
    "Structural protein anomalies",
    "Enzyme deficiency",
    "Kinase-related disorders"
  ),
  Example_Disease = c(
    "Thalassemia",
    "Cardiomyopathy",
    "Collagenopathy",
    "Phenylketonuria",
    "Cancer-related signaling"
  ),
  stringsAsFactors = FALSE
)

uniprot_info <- fetch_uniprot(accessions)
pdb_ids <- unlist(lapply(uniprot_info$xref_pdb, function(x) if(length(x)==0) NULL else x))
default_pdb <- c("1ZMR", "2HWG")
if(length(pdb_ids) == 0){ pdb_ids <- default_pdb }
pdb_info <- fetch_pdb(pdb_ids)

alphafold_info <- fetch_alphafold_prediction("Q57153")
r3dmol() %>%
  m_add_model(url = alphafold_info$pdb_url[1], format = "pdb") %>%
  m_set_style(style = m_style_cartoon()) %>%
  m_zoom_to()

