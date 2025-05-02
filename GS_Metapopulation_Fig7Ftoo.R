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
library(iNEXT.3D)

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
data_rarified<-prune_samples(sample_names(data_rarified)[sample_data(data_rarified)$Site!='HW1'],data_rarified)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Calculate a variety of different ALPHA diversity metrics####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


poplist<-unique(sample_data(data_rarified)$Site)
poplist<-poplist[poplist!='HW1']
pops<-list()
for (k in 1:length(poplist)){
  
  pops<-list.append(pops,sign(otu_table(prune_samples(sample_names(data_rarified)[which(sample_data(data_rarified)$Site==poplist[k])],data_rarified))))
  
}

names(pops)<-poplist

gtemp<-iNEXT3D(pops,diversity='TD',q=0,datatype='incidence_raw')

ff<-gtemp$TDiNextEst$size_based
gg<-ff[round(ff$mT)==10,]

#df<-data.frame(gtemp$TDAsyEst$Assemblage,gtemp$TDAsyEst$TD_asy,gtemp$TDAsyEst$qTD)
df<-data.frame(gg$Assemblage,gg$qTD)
#colnames(df)<-c('site','rich','metric')
colnames(df)<-c('site','rich')
#df_rich<-df[df$metric=='Species richness',]
#df_shannon<-df[df$metric=='Shannon diversity',]


#Put all site characteristics into a dataframe, Z-transforming
site<-sample_data(data_rarified)$Site
cc<-scale(sample_data(data_rarified)$Crevice_Count)
pa<-scale(sample_data(data_rarified)$Pop_Abundance)
os<-scale(sample_data(data_rarified)$Outcrop_Size)
pdc<-scale(sample_data(data_rarified)$Pop_Density_Crevices)
pdos<-scale(sample_data(data_rarified)$Pop_Density_Outcrop_Size)
regression_df<-data.frame(site,cc,os)
regression_df<-regression_df %>% na.omit()
regression_dfs<-data.frame(site,cc,pa,os,pdc,pdos)
regression_dfs<-regression_dfs %>% na.omit()

#Find the median values of site characteristics and alpha-diversities
median_cc<-regression_df %>% group_by(site) %>% summarise(Mean=median(cc))
median_os<-regression_df %>% group_by(site) %>% summarise(Mean=median(os))
median_ccs<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(cc))
median_oss<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(os))
median_pa<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(pa))
median_pdc<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(pdc))
median_pdos<-regression_dfs %>% group_by(site) %>% summarise(Mean=median(pdos))

sitelist<-median_cc$site
richer<-c()
similarity<-c()
for (k in 1:length(sitelist)){
  temp<-which(gg$Assemblage==sitelist[k])
  temp2<-which(poplist==sitelist[k])
  richer<-c(richer,gg$qTD[temp])
  similarity<-c(similarity,1-mean(as.matrix(vegdist(t(pops[[temp2]])))))
}

sitelist<-median_pa$site
richers<-c()
similaritys<-c()
for (k in 1:length(sitelist)){
  temp<-which(gg$Assemblage==sitelist[k])
  temp2<-which(poplist==sitelist[k])
  richers<-c(richers,gg$qTD[temp])
  similaritys<-c(similaritys,1-mean(as.matrix(vegdist(t(pops[[temp2]])))))
  
}


#Make a dataframe of median values for regression
median_regression_df<-data.frame(richer,similarity,median_cc$Mean,median_os$Mean)


median_regression_dfs<-data.frame(richers,similaritys,median_ccs$Mean,median_oss$Mean,median_pa$Mean,median_pdc$Mean,median_pdos$Mean)
colnames(median_regression_df)<-c('richness','similarity','cc','os')
colnames(median_regression_dfs)<-c('richness','similarity','cc','os','pa','pdc','pdos')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Run the regressions####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Richness regressions
lm_richness_cc<-lm(richness~cc,data=median_regression_dfs)
summary(lm_richness_cc)
lm_richness_os<-lm(richness~os,data=median_regression_dfs)
summary(lm_richness_os)
lm_richness_pa<-lm(richness~pa,data=median_regression_dfs)
summary(lm_richness_pa)
lm_richness_pdc<-lm(richness~pdc,data=median_regression_dfs)
summary(lm_richness_pdc)
lm_richness_pdos<-lm(richness~pdos,data=median_regression_dfs)
summary(lm_richness_pdos)

lm_similarity_cc<-lm(similarity~cc,data=median_regression_dfs)
summary(lm_similarity_cc)
lm_similarity_os<-lm(similarity~os,data=median_regression_dfs)
summary(lm_similarity_os)
lm_similarity_pa<-lm(similarity~pa,data=median_regression_dfs)
summary(lm_similarity_pa)
lm_similarity_pdc<-lm(similarity~pdc,data=median_regression_dfs)
summary(lm_similarity_pdc)
lm_similarity_pdos<-lm(similarity~pdos,data=median_regression_dfs)
summary(lm_similarity_pdos)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the regressions####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot richness regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')+ylab("Gamma Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2 <- ggplot(median_regression_df, aes(x=os, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')+ylab("Gamma Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')+ylab("Gamma Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')+ylab("Gamma Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=richness)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')+ylab("Gamma Richness")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

grid_layout_richness <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_richness


#Plot richness regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=similarity)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')+ylab("Similarity")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2 <- ggplot(median_regression_df, aes(x=os, y=similarity)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')+ylab("Similarity")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=similarity)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')+ylab("Similarity")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=similarity)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')+ylab("Similarity")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=similarity)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')+ylab("Similarity")+theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

grid_layout_richness <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_richness


pdf(file = "Figure_6A-F.pdf",   
    width = 10, # The width of the plot in inches
    height = 4) # The height of the plot in inches
grid_layout_richness
dev.off()



#Plot shannon regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')
p2 <- ggplot(median_regression_df, aes(x=os, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=shannon)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')

grid_layout_shannon <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_shannon

#Plot faiths regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')
p2 <- ggplot(median_regression_df, aes(x=os, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=faiths)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')

grid_layout_faiths <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_faiths
