################## Green Salamander Metapopulation Study Script ######################
### Written and directed by Daniel Malagon. Also written by executive producer: Sharon Bewick


########Loading required libraries, setting working directory, loading data#######

#remotes::install_github("KarstensLab/microshades")
#remotes::install_github("mikemc/speedyseq")

setwd("/Users/Daniel/Desktop/GS_microbiomes Objective 3")
library(lme4)
library("phyloseq")
library("ggplot2")  #Used for Graphics
library("VennDiagram") #Used for Venn Diagram
library(venneuler) #Used for Venn Diagram
library(gplots) #Used for Graphics
library(limma)
library(vegan)
library("metagMisc")
library("qiime2R")
library(devtools)
library(car)
library(metagenomeSeq)
library(tidyr)
library(dplyr)
library('stringr')
library('microshades')
library('speedyseq')
library('forcats')
library('cowplot')
library('knitr')
library('ggplot2')
library('PERMANOVA')
library('pairwiseAdonis')
library('ggsignif')
library('tidyverse')
library('ggpubr')
library('betapart')
library('stats')
library(ape)
library('geodist')
library(picante)
library(dplyr)
library(MASS)
library(lmerTest)
library("car")
library(cAIC4)
library("corrplot")
library(patchwork)
library("devtools")
library(data.table)
library(tidyverse)
library(DT)

#Importing data
SVs<-read_qza("ASV/table.qza")
SVs$type
taxonomy<-read_qza("ASV/taxonomy.qza")

#combine different objects into 1 phyloseq object
data<-qza_to_phyloseq(
  features="ASV/table.qza",
  "ASV/taxonomy.qza",
  tree="rooted-tree.qza",
  
  metadata = "metadata.txt")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7 
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")

#Removing ASVs found in sequenced blanks
prune_negatives = function(physeq, negs, samps) {
  negs.n1 = prune_taxa(taxa_sums(negs)>=1, negs) 
  samps.n1 = prune_taxa(taxa_sums(samps)>=1, samps) 
  allTaxa <- names(sort(taxa_sums(physeq),TRUE))
  negtaxa <- names(sort(taxa_sums(negs.n1),TRUE))
  taxa.noneg <- allTaxa[!(allTaxa %in% negtaxa)]
  return(prune_taxa(taxa.noneg,samps.n1))
}

#Defining the samples that are blanks
neg.cts = subset_samples(data, Sample_Type == "Blank")

#Defining the samples that are not
samples = subset_samples(data, Sample_Type != "Blank")

#Creating a new Phyloseq object with only bacteria not found in the blanks
data = prune_negatives(data,neg.cts,samples)

#Creating a new Phyloseq object removing non-site salamanders
data = subset_samples(data,is.na(Site)==FALSE)

#Creating a new Phyloseq object removing non-site salamanders
data = subset_samples(data, Bdfix==0)

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#remove taxa that have chloroplast and mitochondria at level 5, 6,7 
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")

#Rarefy your data, here I am rarefying to 5000 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=5000, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 


set.seed(1)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Find a variety of dissimilarity indices

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the euclidean, jaccard, bray, UniFrac and weighted-UniFrac indices
euclidean<-vegdist(t(otu_table(data_rarified)),method='euclidean',upper=TRUE,diag=TRUE,binary=FALSE)
jaccard<-vegdist(t(otu_table(data_rarified)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(data_rarified)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(data_rarified,weighted=FALSE)
wunifrac<-UniFrac(data_rarified,weighted=TRUE)

#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
euclidean_matrix<-as.matrix(euclidean,labels=TRUE)
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test for distances to within-site versus between-site salamander-controls

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Function to calculate a pseudo F statistic

pseudoF <- function(phylo_object,dist='Euclidean'){
  
  #Make an ordered list of sites
  site<-sample_data(phylo_object)$site
  #Make an ordered list of salamander vs. crevice
  type<-sample_data(phylo_object)$Sample_Type
  #Make a list of the unique sites
  unique_sites<-unique(site)
  
  #Find the distance matrix
  if (dist=='Euclidean'){
    euclidean<-vegdist(t(otu_table(phylo_object)),method='euclidean',upper=TRUE,diag=TRUE,binary=FALSE)
    focal_matrix<-as.matrix(euclidean,labels=TRUE)
  }else if (dist=='Jaccard'){
    jaccard<-vegdist(t(otu_table(phylo_object)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
    focal_matrix<-as.matrix(jaccard,labels=TRUE)
  }else if (dist=='Bray'){
    bray<-vegdist(t(otu_table(phylo_object)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
    focal_matrix<-as.matrix(bray,labels=TRUE)
  }else if (dist=='UniFrac'){
    unifrac<-UniFrac(phylo_object,weighted=FALSE)
    focal_matrix<-as.matrix(unifrac,labels=TRUE)
  }else{
    wunifrac<-UniFrac(phylo_object,weighted=TRUE)
    focal_matrix<-as.matrix(wunifrac,labels=TRUE)
  }
  
  
  matrix_list<-list()
  matrix_list_out<-list()
  type_list<-list()
  type_list_out<-list()
  for (k in 1:length(unique_sites)){
    temp<-which(site==unique_sites[k])
    tempout<-which(site!=unique_sites[k])
    matrix_list<-list.append(matrix_list,focal_matrix[temp,temp])
    matrix_list_out<-list.append(matrix_list_out,focal_matrix[tempout,temp])
    type_list<-list.append(type_list,type[temp])
    type_list_out<-list.append(type_list_out,type[tempout])
  }
  
  summer_in<-0
  summer_out<-0
  for (k in 1:length(unique_sites)){
    temp<-matrix_list[[k]]
    temp_out<-matrix_list_out[[k]]
    temp1<-which(type_list[[k]]=='Crevice')
    temp1_out<-which(type_list_out[[k]]=='Crevice')
    temp2<-which(type_list[[k]]=='Salamander')
    t1<-temp[temp1,temp2]
    summer_in<-summer_in+sum(t1^2)
    t2<-temp_out[temp1_out,temp2]
    summer_out<-summer_out+sum(t2^2)
  }
  
  temp1<-which(type=='Crevice')
  temp2<-which(type=='Salamander')
  fm<-focal_matrix[temp1,temp2]
  
  numerator<-summer_in/(length(unique_sites)-1)
  denominator<-summer_out/(length(temp1)-length(unique_sites))
  pseudoF<-numerator/denominator
  return(pseudoF)
}


data_rarified2<-prune_samples(sample_names(data_rarified)[which(is.na(sample_data(data_rarified)$site)==FALSE)],data_rarified)
data_rarified2<- prune_taxa(taxa_sums(data_rarified2) > 0, data_rarified2) 


true_data<-pseudoF(data_rarified2)

test_data<-c()
for (k in 1:1000){
  crevices<-which(sample_data(data_rarified2)$Sample_Type=='Crevice')
  crevice_sites<-sample_data(data_rarified2)$site[crevices]
  temp_rarified<-data_rarified2
  sample_data(temp_rarified)$site[crevices]<-crevice_sites[sample(1:length(crevice_sites),size=length(crevice_sites),replace=FALSE)]
  test_data<-c(test_data,pseudoF(temp_rarified))
}

stest_data<-sort(test_data)
stest_data[1:50]
p-value<-max(which(true_data-stest_data>0)/1000)

true_data<-pseudoF(data_rarified2,dist='Jaccard')

test_data<-c()
for (k in 1:1000){
  crevices<-which(sample_data(data_rarified2)$Sample_Type=='Crevice')
  crevice_sites<-sample_data(data_rarified2)$site[crevices]
  temp_rarified<-data_rarified2
  sample_data(temp_rarified)$site[crevices]<-crevice_sites[sample(1:length(crevice_sites),size=length(crevice_sites),replace=FALSE)]
  test_data<-c(test_data,pseudoF(temp_rarified,dist='Jaccard'))
}

stest_data<-sort(test_data)
stest_data[1:50]
pvalue<-max(which(true_data-stest_data>0)/1000)



true_data<-pseudoF(data_rarified2,dist='Bray')

test_data<-c()
for (k in 1:1000){
  crevices<-which(sample_data(data_rarified2)$Sample_Type=='Crevice')
  crevice_sites<-sample_data(data_rarified2)$site[crevices]
  temp_rarified<-data_rarified2
  sample_data(temp_rarified)$site[crevices]<-crevice_sites[sample(1:length(crevice_sites),size=length(crevice_sites),replace=FALSE)]
  test_data<-c(test_data,pseudoF(temp_rarified,dist='Bray'))
}

stest_data<-sort(test_data)
stest_data[1:50]
pvalue<-max(which(true_data-stest_data>0)/1000)


true_data<-pseudoF(data_rarified2,dist='UniFrac')

test_data<-c()
for (k in 1:100){
  crevices<-which(sample_data(data_rarified2)$Sample_Type=='Crevice')
  crevice_sites<-sample_data(data_rarified2)$site[crevices]
  temp_rarified<-data_rarified2
  sample_data(temp_rarified)$site[crevices]<-crevice_sites[sample(1:length(crevice_sites),size=length(crevice_sites),replace=FALSE)]
  test_data<-c(test_data,pseudoF(temp_rarified,dist='UniFrac'))
}

stest_data<-sort(test_data)
stest_data[1:5]
pvalue<-max(which(true_data-stest_data>0)/100)

true_data<-pseudoF(data_rarified2,dist='wUniFrac')

test_data<-c()
for (k in 1:100){
  crevices<-which(sample_data(data_rarified2)$Sample_Type=='Crevice')
  crevice_sites<-sample_data(data_rarified2)$site[crevices]
  temp_rarified<-data_rarified2
  sample_data(temp_rarified)$site[crevices]<-crevice_sites[sample(1:length(crevice_sites),size=length(crevice_sites),replace=FALSE)]
  test_data<-c(test_data,pseudoF(temp_rarified,dist='wUniFrac'))
}

stest_data<-sort(test_data)
stest_data[1:5]
pvalue<-max(which(true_data-stest_data>0)/100)



