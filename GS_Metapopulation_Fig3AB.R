################## Green Salamander Metapopulation Study Script ######################
### Written and directed by Daniel Malagon. Also written by executive producer: Sharon Bewick


########Loading required libraries, setting working directory, loading data#######

#remotes::install_github("KarstensLab/microshades")
#remotes::install_github("mikemc/speedyseq")

setwd("~/Desktop/GS_microbiomes Objective 3/Daniel_Green_Metapopulation")
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
library(rlist)

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


data_rarified2<-prune_samples(sample_names(data_rarified)[which(is.na(sample_data(data_rarified)$site)==FALSE)],data_rarified)
data_rarified2<- prune_taxa(taxa_sums(data_rarified2) > 0, data_rarified2) 




#Make an ordered list of sites
site<-sample_data(data_rarified2)$site
#Make an ordered list of salamander vs. crevice
type<-sample_data(data_rarified2)$Sample_Type
#Make a list of the unique sites
unique_sites<-unique(site)

hitall<-which(rowSums(otu_table(data_rarified2)[,which(sample_data(data_rarified2)$Sample_Type=='Crevice')])>0)
hitnone<-which(rowSums(otu_table(data_rarified2)[,which(sample_data(data_rarified2)$Sample_Type=='Crevice')])==0)


matrix_list<-list()
omatrix_list<-list()
type_list<-list()
otype_list<-list()
for (k in 1:length(unique_sites)){
  temp<-which(site==unique_sites[k])
  matrix_list<-list.append(matrix_list,otu_table(data_rarified2)[,temp])
  type_list<-list.append(type_list,type[temp])
  otemp<-which(site!=unique_sites[k])
  omatrix_list<-list.append(omatrix_list,otu_table(data_rarified2)[,otemp])
  otype_list<-list.append(otype_list,type[otemp])
}

list_my<-c()
elist_my<-c()
list_any<-c()
list_mysal<-c()
list_sal<-c()
sitename<-c()
samplename<-c()
for (k in 1:length(unique_sites)){
  
  #Find the taxa that are shared with the local crevices
  hit<-which(rowSums(matrix_list[[k]][,which(type_list[[k]]=='Crevice')])>0)
  hit1a<-intersect(hit,which(rowSums(omatrix_list[[k]][,which(otype_list[[k]]=='Crevice')])==0))
  hit1b<-intersect(hit,which(rowSums(omatrix_list[[k]][,which(otype_list[[k]]=='Crevice')])!=0))
  #Find the fraction of reads that are shared with the local crevice
  likemycrevicea<-colSums(matrix_list[[k]][hit1a,which(type_list[[k]]!='Crevice')])/5000
  likemycreviceb<-colSums(matrix_list[[k]][hit1b,which(type_list[[k]]!='Crevice')])/5000
  #Find the fraction of reads that are shared with any crevices and subtract the number shared with the local crevices
  likeanycrevice<-(colSums(matrix_list[[k]][hitall,which(type_list[[k]]!='Crevice')]))/5000-likemycrevicea-likemycreviceb
  
  #Find the taxa that are not on any crevices
  onlocsal<-which(rowSums(matrix_list[[k]][,which(type_list[[k]]=='Salamander')])>0)
  hit2<-setdiff(onlocsal,hitall)
  hit2<-intersect(hit2,which(rowSums(omatrix_list[[k]])==0))
  likemysalamanders<-colSums(matrix_list[[k]][hit2,which(type_list[[k]]!='Crevice')])/5000
  onlysals<-(colSums(matrix_list[[k]][,which(type_list[[k]]!='Crevice')]))/5000-likemycrevicea-likemycreviceb-likeanycrevice-likemysalamanders
  
  onlocsal<-which(rowSums(matrix_list[[k]][,which(type_list[[k]]=='Salamander')])>0)
  
  elist_my<-c(elist_my,likemycrevicea) 
  list_my<-c(list_my,likemycreviceb) 
  list_any<-c(list_any,likeanycrevice) 
  list_mysal<-c(list_mysal,likemysalamanders)
  list_sal<-c(list_sal,onlysals) 
  sitename<-c(sitename,rep(unique_sites[k],length(onlysals)))
  samplename<-c(samplename,colnames(matrix_list[[k]][,which(type_list[[k]]!='Crevice')]))
}

df<-data.frame(c(elist_my,list_my,list_any,list_mysal,list_sal),c(rep('local crevice only',length(list_my)),rep('local & other crevice',length(list_my)),rep('other crevice only',length(list_any)),rep('local salamanders only',length(list_mysal)),rep('multiple salamander populations',length(list_sal))),rep(sitename,5),rep(samplename,5))
colnames(df)<-c('abundance','location','site','sample')

ggplot(df, aes(fill=factor(location,levels=c('local crevice only','local & other crevice','other crevice only','local salamanders only','multiple salamander populations')), y=abundance, x=sample)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values = c("grey89", "grey70","grey50","darkseagreen1", "green"))+theme_bw()+theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border=element_blank())+
  theme(legend.title=element_blank(), axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 70, vjust = 1.5, hjust=1),axis.title=element_text(size=18),axis.title.x = element_text(vjust = 7.))




mean(elist_my+list_my+list_any)
median(elist_my+list_my+list_any)


list_my<-c()
elist_my<-c()
list_any<-c()
list_mysal<-c()
list_sal<-c()
sitename<-c()
samplename<-c()
for (k in 1:length(unique_sites)){
  
  #Find the taxa that are shared with the local crevices
  hit<-which(rowSums(matrix_list[[k]][,which(type_list[[k]]=='Crevice')])>0)
  hit1a<-intersect(hit,which(rowSums(omatrix_list[[k]][,which(otype_list[[k]]=='Crevice')])==0))
  hit1b<-intersect(hit,which(rowSums(omatrix_list[[k]][,which(otype_list[[k]]=='Crevice')])!=0))
  #Find the fraction of reads that are shared with the local crevice
  likemycrevicea<-colSums(sign(matrix_list[[k]])[hit1a,which(type_list[[k]]!='Crevice')])/colSums(sign(matrix_list[[k]])[,which(type_list[[k]]!='Crevice')])
  likemycreviceb<-colSums(sign(matrix_list[[k]])[hit1b,which(type_list[[k]]!='Crevice')])/colSums(sign(matrix_list[[k]])[,which(type_list[[k]]!='Crevice')])
  #Find the fraction of reads that are shared with any crevices and subtract the number shared with the local crevices
  likeanycrevice<-(colSums(sign(matrix_list[[k]])[hitall,which(type_list[[k]]!='Crevice')]))/colSums(sign(matrix_list[[k]])[,which(type_list[[k]]!='Crevice')])-likemycrevicea-likemycreviceb
  
  #Find the taxa that are not on any crevices
  onlocsal<-which(rowSums(matrix_list[[k]][,which(type_list[[k]]=='Salamander')])>0)
  hit2<-setdiff(onlocsal,hitall)
  hit2<-intersect(hit2,which(rowSums(omatrix_list[[k]])==0))
  likemysalamanders<-colSums(sign(matrix_list[[k]])[hit2,which(type_list[[k]]!='Crevice')])/colSums(sign(matrix_list[[k]])[,which(type_list[[k]]!='Crevice')])
  onlysals<-(colSums(sign(matrix_list[[k]])[,which(type_list[[k]]!='Crevice')]))/colSums(sign(matrix_list[[k]])[,which(type_list[[k]]!='Crevice')])-likemycrevicea-likemycreviceb-likeanycrevice-likemysalamanders
  
  
  elist_my<-c(elist_my,likemycrevicea) 
  list_my<-c(list_my,likemycreviceb) 
  list_any<-c(list_any,likeanycrevice) 
  list_mysal<-c(list_mysal,likemysalamanders)
  list_sal<-c(list_sal,onlysals) 
  sitename<-c(sitename,rep(unique_sites[k],length(onlysals)))
  samplename<-c(samplename,colnames(matrix_list[[k]][,which(type_list[[k]]!='Crevice')]))
}

df<-data.frame(c(elist_my,list_my,list_any,list_mysal,list_sal),c(rep('local crevice only',length(list_my)),rep('local & other crevice',length(list_my)),rep('other crevice only',length(list_any)),rep('local salamanders only',length(list_mysal)),rep('multiple salamander populations',length(list_sal))),rep(sitename,5),rep(samplename,5))
colnames(df)<-c('taxa','location','site','sample')

ggplot(df, aes(fill=factor(location,levels=c('local crevice only','local & other crevice','other crevice only','local salamanders only','multiple salamander populations')), y=taxa, x=sample)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values = c("grey89", "grey70","grey50","darkseagreen1","green"))+theme_bw()+theme(panel.background = element_blank(),panel.grid = element_blank(),panel.border=element_blank())+
  theme(legend.title=element_blank(), axis.ticks.x=element_blank(),axis.text.x = element_text(angle = 70, vjust = 1.5, hjust=1),axis.title=element_text(size=18),axis.title.x = element_text(vjust = 7.))


mean(elist_my+list_my+list_any)
median(elist_my+list_my+list_any)




