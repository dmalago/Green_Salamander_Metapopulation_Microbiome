# Spatial Variation of Skin-associated Microbiota in a Green Salamander Metapopulation
This repository contains all data and code related to the publication, "Spatial Variation of Skin-associated Microbiota in a Green Salamander Metapopulation"

Fastq files are housed in NCBI SRA under accession number PRJNA1176676


To perform the analyses at the Genus taxonomic level, add the following line of code where appropriate.

data_rarified<-tax_glom(data_rarified, taxrank="Genus", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))


To perform the analyses at the ASV level with all environmental microbiota removed, add the following lines of code where appropriate.

prune_negatives = function(physeq, negs, samps) {
  negs.n1 = prune_taxa(taxa_sums(negs)>=1, negs) 
  samps.n1 = prune_taxa(taxa_sums(samps)>=1, samps) 
  allTaxa <- names(sort(taxa_sums(physeq),TRUE))
  negtaxa <- names(sort(taxa_sums(negs.n1),TRUE))
  taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
  return(prune_taxa(taxa.noneg,samps.n1))
}

#Defining the samples that are blanks
neg.cts = subset_samples(data_rarified, Sample_Type == "Crevice")

#Defining the samples that are not
samples = subset_samples(data_rarified, Sample_Type != "Crevice")

#Creating a new Phyloseq object with only bacteria not found in the environment
data_rarified = prune_negatives(data_rarified,neg.cts,samples)


