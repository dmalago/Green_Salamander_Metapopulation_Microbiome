################## Green Salamander Metapopulation Study Script Figure 6 ######################
### Written by Daniel Malagon. Also written by director and executive producer: Sharon Bewick


########Loading required libraries, setting working directory, loading data#######

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
library(MASS)
library(lmerTest)
library(cAIC4)
library("corrplot")
library(patchwork)
library(data.table)
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
data = subset_samples(data,is.na(site_all)==FALSE)

#Creating a new Phyloseq object removing Bd positive salamanders
data = subset_samples(data, Bdfix==0)



#Rarefy your data, here I am rarefying to 5000 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=5000, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
Crevicesdata<-subset_samples(data_rarified, Sample_Type == "Crevice")

data_rarified=subset_samples(data_rarified, Repeat == "0")

####Pull out only the salamander samples (no crevices)
data_rarified<-subset_samples(data_rarified, Sample_Type == "Salamander")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Calculate a variety of different ALPHA diversity metrics####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Find the richness for each sample
richness<-colSums(sign(otu_table(data_rarified)))

#Find the shannon diversity for each sample
shannon<-diversity(otu_table(t(data_rarified)),index='shannon')

#Find Faith's pd for each sample
faithpd<-pd(t(otu_table(data_rarified)), phy_tree(data_rarified), include.root=TRUE)
faiths<-faithpd[,1]

type<-sample_data(data_rarified)$Site

#Put all of your different diversity metrics into a dataframe
diversity_df<-data.frame(type,richness,shannon,faiths)
#Switch the ordering of the types (to the order you want them presented in your plots)
diversity_df$type<-factor(diversity_df$type,levels=c('HW1','BB','DNR','1250','1251','1292','TR2','3688','1477'))

#Put all site characteristics into a dataframe, Z-transforming
site<-sample_data(data_rarified)$Site
cc<-scale(sample_data(data_rarified)$Crevice_Count)
pa<-scale(sample_data(data_rarified)$Pop_Abundance)
os<-scale(sample_data(data_rarified)$Outcrop_Size)
pdc<-scale(sample_data(data_rarified)$Pop_Density_Crevices)
pdos<-scale(sample_data(data_rarified)$Pop_Density_Outcrop_Size)
regression_df<-data.frame(diversity_df,site,cc,os)
regression_df<-regression_df %>% na.omit()
regression_dfs<-data.frame(diversity_df,site,cc,pa,os,pdc,pdos)
regression_dfs<-regression_dfs %>% na.omit()

#Find the median values of site characteristics and alpha-diversities
median_richness<-regression_df %>% group_by(site) %>% summarise(Mean=median(richness))
median_shannon<-regression_df %>% group_by(site) %>% summarise(Mean=median(shannon))
median_faiths<-regression_df %>% group_by(site) %>% summarise(Mean=median(faiths))
median_richnesss<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(richness))
median_shannons<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(shannon))
median_faithss<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(faiths))
median_cc<-regression_df %>% group_by(site) %>% summarise(Mean=median(cc))
median_os<-regression_df %>% group_by(site) %>% summarise(Mean=median(os))
median_ccs<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(cc))
median_oss<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(os))
median_pa<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(pa))
median_pdc<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(pdc))
median_pdos<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(pdos))

#Make a dataframe of median values for regression
median_regression_df<-data.frame(median_richness$Mean,median_shannon$Mean,median_faiths$Mean,median_cc$Mean,median_os$Mean)
median_regression_dfs<-data.frame(median_richnesss$Mean,median_shannons$Mean,median_faithss$Mean,median_ccs$Mean,median_oss$Mean,median_pa$Mean,median_pdc$Mean,median_pdos$Mean)
colnames(median_regression_df)<-c('richness','shannon','faiths','cc','os')
colnames(median_regression_dfs)<-c('richness','shannon','faiths','cc','os','pa','pdc','pdos')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Run the regressions####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Richness regressions
lm_richness_cc<-lm(richness~cc,data=median_regression_df)
summary(lm_richness_cc)
lm_richness_os<-lm(richness~os,data=median_regression_df)
summary(lm_richness_os)
lm_richness_pa<-lm(richness~pa,data=median_regression_dfs)
summary(lm_richness_pa)
lm_richness_pdc<-lm(richness~pdc,data=median_regression_dfs)
summary(lm_richness_pdc)
lm_richness_pdos<-lm(richness~pdos,data=median_regression_dfs)
summary(lm_richness_pdos)

#Shannon regressions
lm_shannon_cc<-lm(shannon~cc,data=median_regression_df)
summary(lm_shannon_cc)
lm_shannon_os<-lm(shannon~os,data=median_regression_df)
summary(lm_shannon_os)
lm_shannon_pa<-lm(shannon~pa,data=median_regression_dfs)
summary(lm_shannon_pa)
lm_shannon_pdc<-lm(shannon~pdc,data=median_regression_dfs)
summary(lm_shannon_pdc)
lm_shannon_pdos<-lm(shannon~pdos,data=median_regression_dfs)
summary(lm_shannon_pdos)

#Faiths regressions
lm_faiths_cc<-lm(faiths~cc,data=median_regression_df)
summary(lm_faiths_cc)
lm_faiths_os<-lm(faiths~os,data=median_regression_df)
summary(lm_faiths_os)
lm_faiths_pa<-lm(faiths~pa,data=median_regression_dfs)
summary(lm_faiths_pa)
lm_faiths_pdc<-lm(faiths~pdc,data=median_regression_dfs)
summary(lm_faiths_pdc)
lm_faiths_pdos<-lm(faiths~pdos,data=median_regression_dfs)
summary(lm_faiths_pdos)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the regressions####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot richness regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')+ylab("Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2 <- ggplot(median_regression_df, aes(x=os, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')+ylab("Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')+ylab("Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')+ylab("Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')+ylab("Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

grid_layout_richness <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_richness

pdf(file = "Figure_6A-F.pdf",   
    width = 10, # The width of the plot in inches
    height = 4) # The height of the plot in inches
grid_layout_richness
dev.off()



#Plot shannon regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2 <- ggplot(median_regression_df, aes(x=os, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

grid_layout_shannon <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_shannon

#Plot faiths regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2 <- ggplot(median_regression_df, aes(x=os, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

grid_layout_faiths <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_faiths
