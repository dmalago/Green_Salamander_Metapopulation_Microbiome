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

####Remove HW1 because there are not enough animals at that site
data_rarified<-subset_samples(data_rarified, Site != "HW1")

type_full<-sample_data(data_rarified)$Site
types<-unique(type_full)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Bootstrap beta diversity####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set.seed(1)

#The number from each treatment group you want to pool to assess core richness for each bootstrap
sample_size<-5

#The number of bootstrap samples you want to do
bootstrap_no<-7

#The treatment group lables
types<-unique(type_full)

#Define a matrix that will keep track of each treatment group core richness for each bootstrap sample
boot_tracker_turnover_jaccard<-data.frame()
boot_tracker_nestedness_jaccard<-data.frame()
boot_tracker_overall_jaccard<-data.frame()

#Define a matrix that will keep track of each treatment group core richness for each bootstrap sample
boot_tracker_turnover_jaccard_phy<-data.frame()
boot_tracker_nestedness_jaccard_phy<-data.frame()
boot_tracker_overall_jaccard_phy<-data.frame()


#boot_tracker_turnover_bray<-data.frame()
#boot_tracker_nestedness_bray<-data.frame()
boot_tracker_overall_bray<-data.frame()

#For each bootstrap sample...
for (j in 1:bootstrap_no){
  print(j)
  #Make a list of the individual animals you're going to randomly choose
  picklist<-c()
  pickgroups<-matrix(0,ncol=sample_size,nrow=length(types))
  for (k in 1:length(types)){
    #Randomly select the appropriate number of animals from each group
    pick<-sample(which(sample_data(data_rarified)$Site==types[k]),sample_size,replace=FALSE)
    #Add those animals to your list
    picklist<-c(picklist,pick)
    pickgroups[k,]<-pick
    
  }
  
  #Make a phyloseq object of only the individual animals you selected for this bootstrap sample
  ZB_physeq<-prune_samples(sample_names(data_rarified)[picklist],data_rarified)
  
  #Rarefy the OTU table specifically for this bootstrap sample
  rZOTU<-rarefy_even_depth(ZB_physeq,verbose=FALSE,rngseed = sample(1000,1))
  
  
  
  
  #Initialize a vector of counts of core microbes for each treatment group
  turnover_jaccard<-c()
  nestedness_jaccard<-c()
  overall_jaccard<-c()
  turnover_jaccard_phy<-c()
  nestedness_jaccard_phy<-c()
  overall_jaccard_phy<-c()
  # turnover_bray<-c()
  # nestedness_bray<-c()
  overall_bray<-c()
  
  #For each treatment group...
  for (k in 1:length(types)){
    #Pull out the microbiomes from that treatment group and put them in their own phyloseq object
    rZOTUtemp<-prune_samples(sample_names(data_rarified)[pickgroups[k,]],rZOTU)
    #Sum up the presences of each microbial taxon in the subsetted phyloseq object... count how many microbial taxa are present on as many or more individuals than the defined core fraction
    temp_community_matrix<-otu_table(rZOTUtemp)
    temp_tree<-phy_tree(rZOTUtemp)
    class(temp_community_matrix) <- "matrix"
    temp_betapart_jaccard<-beta.multi(sign(t(temp_community_matrix)),index.family='jaccard')
    temp_betapart_jaccard_phy<-phylo.beta.multi(sign(t(temp_community_matrix)),temp_tree,index.family='jaccard')
    temp_betapart_bray<-beta.multi.abund(t(temp_community_matrix),index.family='bray')
    
    turnover_jaccard<-c(turnover_jaccard,temp_betapart_jaccard$beta.JTU)
    nestedness_jaccard<-c(nestedness_jaccard,temp_betapart_jaccard$beta.JNE)
    overall_jaccard<-c(overall_jaccard,temp_betapart_jaccard$beta.JAC)
    turnover_jaccard_phy<-c(turnover_jaccard_phy,temp_betapart_jaccard_phy$phylo.beta.JTU)
    nestedness_jaccard_phy<-c(nestedness_jaccard_phy,temp_betapart_jaccard_phy$phylo.beta.JNE)
    overall_jaccard_phy<-c(overall_jaccard_phy,temp_betapart_jaccard_phy$phylo.beta.JAC)
    
    # turnover_bray<-c(turnover_bray,temp_betapart_bray$beta.BRAY.BAL)
    # nestedness_bray<-c(nestedness_bray,temp_betapart_bray$beta.BRAY.GRA)
    overall_bray<-c(overall_bray,temp_betapart_bray$beta.BRAY)
  }
  
  #Record the core counts for each treatment group from this particular bootstrap
  boot_tracker_turnover_jaccard<-rbind(boot_tracker_turnover_jaccard,turnover_jaccard)
  boot_tracker_nestedness_jaccard<-rbind(boot_tracker_nestedness_jaccard,nestedness_jaccard)
  boot_tracker_overall_jaccard<-rbind(boot_tracker_overall_jaccard,overall_jaccard)
  boot_tracker_turnover_jaccard_phy<-rbind(boot_tracker_turnover_jaccard_phy,turnover_jaccard_phy)
  boot_tracker_nestedness_jaccard_phy<-rbind(boot_tracker_nestedness_jaccard_phy,nestedness_jaccard_phy)
  boot_tracker_overall_jaccard_phy<-rbind(boot_tracker_overall_jaccard_phy,overall_jaccard_phy)
  
  # boot_tracker_turnover_bray<-rbind(boot_tracker_turnover_bray,turnover_bray)
  # boot_tracker_nestedness_bray<-rbind(boot_tracker_nestedness_bray,nestedness_bray)
  boot_tracker_overall_bray<-rbind(boot_tracker_overall_bray,overall_bray)
  
}

#Name the bootstrap columns based on the treatment groups
colnames(boot_tracker_turnover_jaccard)<-types
#colnames(boot_tracker_turnover_bray)<-types
colnames(boot_tracker_nestedness_jaccard)<-types
# colnames(boot_tracker_nestedness_bray)<-types
colnames(boot_tracker_overall_jaccard)<-types
colnames(boot_tracker_turnover_jaccard_phy)<-types
#colnames(boot_tracker_turnover_bray)<-types
colnames(boot_tracker_nestedness_jaccard_phy)<-types
# colnames(boot_tracker_nestedness_bray)<-types
colnames(boot_tracker_overall_jaccard_phy)<-types
colnames(boot_tracker_overall_bray)<-types

#Turn your bootstraps into two lists...
#One with the treatment group for each measurement
type<-c()
#And one with the core richness for each measurement
jaccard_turnover<-c()
jaccard_turnover_phy<-c()
#bray_turnover<-c()
jaccard_nestedness<-c()
jaccard_nestedness_phy<-c()
#bray_nestedness<-c()
jaccard_overall<-c()
jaccard_overall_phy<-c()
bray_overall<-c()

for (k in 1:length(types)){
  #A list of the treatment group for each pooled bootstrap
  type<-c(type,rep(types[k],bootstrap_no))
  #A list of the core richness for each pooled bootstrap
  jaccard_turnover<-c(jaccard_turnover,boot_tracker_turnover_jaccard[,k])
  jaccard_nestedness<-c(jaccard_nestedness,boot_tracker_nestedness_jaccard[,k])
  jaccard_overall<-c(jaccard_overall,boot_tracker_overall_jaccard[,k])
  jaccard_turnover_phy<-c(jaccard_turnover_phy,boot_tracker_turnover_jaccard_phy[,k])
  jaccard_nestedness_phy<-c(jaccard_nestedness_phy,boot_tracker_nestedness_jaccard_phy[,k])
  jaccard_overall_phy<-c(jaccard_overall_phy,boot_tracker_overall_jaccard_phy[,k])
  
  #  bray_turnover<-c(bray_turnover,boot_tracker_turnover_bray[,k])
  #  bray_nestedness<-c(bray_nestedness,boot_tracker_nestedness_bray[,k])
  bray_overall<-c(bray_overall,boot_tracker_overall_bray[,k])
}

#Put your richness and treatment group data into a dataframe
diversity_df<-data.frame(type,jaccard_turnover,jaccard_nestedness,jaccard_overall,bray_overall,jaccard_turnover_phy,jaccard_nestedness_phy,jaccard_overall_phy)#,bray_turnover,bray_nestedness)
#Switch the ordering of the types (to the order you want them presented in your plots)
diversity_df$type<-factor(diversity_df$type,levels=c('BB','DNR','1250','1251','1292','TR2','3688','1477'))

site<-c('BB','DNR','1250','1251','1292','TR2','3688','1477')

#Put all site characteristics into a dataframe
sites<-sample_data(data_rarified)$Site
ccs<-sample_data(data_rarified)$Crevice_Count
pas<-sample_data(data_rarified)$Pop_Abundance
oss<-sample_data(data_rarified)$Outcrop_Size
pdcs<-sample_data(data_rarified)$Pop_Density_Crevices
pdoss<-sample_data(data_rarified)$Pop_Density_Outcrop_Size

cc<-c()
pa<-c()
os<-c()
pdc<-c()
pdos<-c()
for (k in 1:length(site)){
  temp<-which(sites == site[k])
  cc<-c(cc,ccs[temp[1]])
  pa<-c(pa,pas[temp[1]])
  os<-c(os,oss[temp[1]])
  pdc<-c(pdc,pdcs[temp[1]])
  pdos<-c(pdos,pdoss[temp[1]])
}

#Find the median values of site characteristics and alpha-diversities####
median_jaccard<-diversity_df %>% group_by(type) %>% summarise(Mean=median(jaccard_overall))
median_jaccard_phy<-diversity_df %>% group_by(type) %>% summarise(Mean=median(jaccard_overall_phy))
median_bray<-diversity_df %>% group_by(type) %>% summarise(Mean=median(bray_overall))



#Make a dataframe of median values for regression
median_regression_df<-data.frame(median_jaccard$Mean,median_jaccard_phy$Mean,median_bray$Mean,cc,os)
median_regression_dfs<-data.frame(median_jaccard$Mean,median_jaccard_phy$Mean,median_bray$Mean,cc,os,pa,pdc,pdos)
median_regression_dfs<-median_regression_dfs %>% na.omit()

colnames(median_regression_df)<-c('jaccard','unifrac','bray','cc','os')
colnames(median_regression_dfs)<-c('jaccard','unifrac','bray','cc','os','pa','pdc','pdos')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Run the regressions####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Richness regressions
lm_jaccard_cc<-lm(jaccard~cc,data=median_regression_df)
summary(lm_jaccard_cc)
lm_jaccard_os<-lm(jaccard~os,data=median_regression_df)
summary(lm_jaccard_os)
lm_jaccard_pa<-lm(jaccard~pa,data=median_regression_dfs)
summary(lm_jaccard_pa)
lm_jaccard_pdc<-lm(jaccard~pdc,data=median_regression_dfs)
summary(lm_jaccard_pdc)
lm_jaccard_pdos<-lm(jaccard~pdos,data=median_regression_dfs)
summary(lm_jaccard_pdos)

#unifrac regressions
lm_unifrac_cc<-lm(unifrac~cc,data=median_regression_df)
summary(lm_unifrac_cc)
lm_unifrac_os<-lm(unifrac~os,data=median_regression_df)
summary(lm_unifrac_os)
lm_unifrac_pa<-lm(unifrac~pa,data=median_regression_dfs)
summary(lm_unifrac_pa)
lm_unifrac_pdc<-lm(unifrac~pdc,data=median_regression_dfs)
summary(lm_unifrac_pdc)
lm_unifrac_pdos<-lm(unifrac~pdos,data=median_regression_dfs)
summary(lm_unifrac_pdos)

#bray regressions
lm_bray_cc<-lm(bray~cc,data=median_regression_df)
summary(lm_bray_cc)
lm_bray_os<-lm(bray~os,data=median_regression_df)
summary(lm_bray_os)
lm_bray_pa<-lm(bray~pa,data=median_regression_dfs)
summary(lm_bray_pa)
lm_bray_pdc<-lm(bray~pdc,data=median_regression_dfs)
summary(lm_bray_pdc)
lm_bray_pdos<-lm(bray~pdos,data=median_regression_dfs)
summary(lm_bray_pdos)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the regressions####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot jaccard regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=jaccard)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2 <- ggplot(median_regression_df, aes(x=os, y=jaccard)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=jaccard)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=jaccard)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=jaccard)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

grid_layout_jaccard <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_jaccard

pdf(file = "Figure_6F-J.pdf",   
    width = 10, # The width of the plot in inches
    height = 4) # The height of the plot in inches
grid_layout_jaccard
dev.off()

#Plot unifrac regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=unifrac)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')
p2 <- ggplot(median_regression_df, aes(x=os, y=unifrac)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=unifrac)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=unifrac)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=unifrac)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')

grid_layout_unifrac <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_unifrac

#Plot bray regressions
p1 <- ggplot(median_regression_df, aes(x=cc, y=bray)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('crevice count')
p2 <- ggplot(median_regression_df, aes(x=os, y=bray)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('outcrop size')
p3 <- ggplot(median_regression_dfs, aes(x=pa, y=bray)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('pop. abundance')
p4 <- ggplot(median_regression_dfs, aes(x=pdc, y=bray)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (crevice)')
p5 <- ggplot(median_regression_dfs, aes(x=pdos, y=bray)) + geom_point()+geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,color='red')+xlab('density (outcrop)')

grid_layout_bray <- plot_grid(p1, p2, p3,p4,p5, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4,1/4,1/4), rel_widths = c(1/4, 1/4, 1/4,1/4,1/4))
grid_layout_bray
