#set working directory
setwd("C:/Users/domin/OneDrive/Documents/GitHub/Bioinformatics/Lab10")
#package installation
install.packages("UniprotR")
install.packages("protti")
#load packages
library(UniprotR)
library(protti)
install.packages("r3dmol")
BiocManager::install("GenomicAlignments")
library(Biostrings)
library(r3dmol)
# turn dna data to protein 
dna <- readDNAStringSet("Hs7.fasta")
protein <- translate(dna)
protein
#write protein to fasta
writeXStringSet(protein, "protein.fasta")
#write accessions
accessions <- c("Q57153", "Q3I4Z3", "B7XSQ5", "Q45231", "Q27JV1")
print(accessions)
#Get gene ontology
go_terms <- GetProteinGOInfo(accessions)
#plot results
PlotGoInfo(go_terms)
#Save Go Plot to Working Directory
PlotGOAll(GOObj = go_terms, 
          Top = 10, 
          directorypath = getwd(), 
          width = 8, 
          height = 5)
#Get pathology and disease info
pathology_list <- list()
disease_list <- list()
for(acc in accessions)
  pathology <- GetPathology_Biotech(acc)
pathology_list[[acc]] <- pathology
diseases <- Get.diseases(acc)
disease_list[[acc]] <- diseases   
print("Pathology/Biotech info for Q57153:")
print(pathology_list[["Q57153"]])  
#Access structural information
uniprot_info <- fetch_uniprot(accessions)
head(uniprot_info)
niprot_info$xref_pdb[[1]]  
#pull available structural information
default_pdb <- c("1ZMR", "2HWG")
if(length(pdb_ids) == 0){
  pdb_ids <- default_pdb
}
pdb_info <- fetch_pdb(pdb_ids)
head(pdb_info)

