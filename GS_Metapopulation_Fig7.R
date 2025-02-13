################## Green Salamander Metapopulation Study Script Figure 7 ######################
### Written by Daniel Malagon. Also written by director and executive producer: Sharon Bewick

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


########Loading required libraries, setting working directory, loading data#######
setwd("~/Desktop/GS_microbiomes Objective 3/Daniel_Green_Metapopulation")
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
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=5000, replace=F, verbose = TRUE,trimOTUs = 0)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
Crevicesdata<-subset_samples(data_rarified, Sample_Type == "Crevice")

data_rarified=subset_samples(data_rarified, Repeat == "0")


####consider only salamander microbiomes
data_rarified<-subset_samples(data_rarified, Sample_Type == "Salamander")

      
      
      

setwd("~/Desktop/GS_microbiomes Objective 3/Daniel_Green_Metapopulation/Antifungal")

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
library(MASS)
library(lmerTest)
library("car")
library(cAIC4)
library("corrplot")
library(patchwork)
library(data.table)
library(tidyverse)
library(DT)
library("readr")
library("magrittr")
library("dplyr")
library("purrr")
library("stringr")




#Preparing the Antifungal ASV table using code from https://github.com/AmphiBac/AmphiBac-Database/ 
Antifungal_table = qiime2R::read_qza("Antifungal_Matches_99/clustered_table.qza")

Antifungal_table = as.data.frame(Antifungal_table$data)

TotalAntiFungal = Antifungal_table %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

TotalAntiFungal = as.data.frame(TotalAntiFungal)

TotalAntiFungal = cbind(SampleID = rownames(TotalAntiFungal), TotalAntiFungal)

rownames(TotalAntiFungal) <- NULL



rardepth <- as.numeric(5000)

TotalAntiFungal_data = TotalAntiFungal %>%
  mutate(Propor_TotalAntiFungal = TotalAntiFungal/rardepth)




Metadata = readr::read_tsv("metadata.tsv")
names(Metadata)[names(Metadata) == "sample-id"] <- "SampleID"

Metadata_Antifungal = TotalAntiFungal_data %>%
  left_join(Metadata, by = "SampleID")

Antifungal_tablePA = qiime2R::read_qza("Antifungal_Matches_99/PresAbs_clustered_table.qza")

Antifungal_tablePA = as.data.frame(Antifungal_tablePA$data)

AntiFungal_Richness = Antifungal_tablePA %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

AntiFungal_Richness = as.data.frame(AntiFungal_Richness)

AntiFungal_Richness = cbind(SampleID = rownames(AntiFungal_Richness), AntiFungal_Richness)

rownames(AntiFungal_Richness) <- NULL

Metadata_Antifungal_Richness = Metadata_Antifungal %>%
  left_join(AntiFungal_Richness, by = "SampleID")



#readr::write_tsv(Metadata_Antifungal_Richness,"Metadata_Antifungal_Predictions.txt")

#Turn the above text file into a csv file named "Antifungal_data.csv"
Antifungal_data<-read.csv("Antifungal_data.csv")
##### Analyses ####

#Antifungal richness across populations
a<-ggplot(Antifungal_data, aes(x = reorder(Site, AntiFungal_Richness), y = AntiFungal_Richness, fill = Site)) + 
  geom_violin(width = 1) + scale_fill_manual( values =  c("1250"= "red", "1251" = "blue",  "1292" = "forestgreen", "1477" ="purple", "3688" = "darkorange","BB" = "yellow", "DNR" = "sienna", "HW1" ="pink",  "TR2" = "gray27"))+
  stat_summary(fun="mean", color = "black", size= 1.25)+
  labs( 
    title = "",  
    x = "Population", 
    y = "Antifungal ASV Richness" 
  ) + 
  geom_point(position = position_jitter(seed = 1, width = 0.25, ), color = "grey") +
  theme_classic()+
  theme(text = element_text(size=30))+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


Richness_site_kruskal_wallis<-kruskal.test(AntiFungal_Richness~Site, data = Antifungal_data)
Richness_site_kruskal_wallis

pwt<-pairwise.wilcox.test(Antifungal_data$AntiFungal_Richness, Antifungal_data$Site,
                          p.adjust.method = "BH")
pwt

pwt2<-pairwise.wilcox.test(Antifungal_data$AntiFungal_Richness, Antifungal_data$Site,
                          p.adjust.method = "none")
pwt2


b<-ggplot(Antifungal_data, aes(x = reorder(Site, Propor_TotalAntiFungal), y = Propor_TotalAntiFungal, fill = Site)) + 
  geom_violin(width = 1) + scale_fill_manual( values =  c("1250"= "red", "1251" = "blue",  "1292" = "forestgreen", "1477" ="purple", "3688" = "darkorange","BB" = "yellow", "DNR" = "sienna", "HW1" ="pink",  "TR2" = "gray27"))+
  stat_summary(fun="mean", color = "black", size= 1.25)+
  labs( 
    title = "",  
    x = "Population", 
    y = "Proportion of Microbiota Reads Mapping Antifungal" 
  ) + 
  geom_point(position = position_jitter(seed = 1, width = 0.25, ), color = "grey") +
  theme_classic()+
  theme(text = element_text(size=30))+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

Propor_site_kruskal_wallis<-kruskal.test(Propor_TotalAntiFungal~Site, data = Antifungal_data)
Propor_site_kruskal_wallis

pwt<-pairwise.wilcox.test(Antifungal_data$AntiFungal_Richness, Antifungal_data$Site,
                          p.adjust.method = "BH")
pwt



#Make your plot

pdf(file = "Figure_7.pdf",   
    width = 10, # The width of the plot in inches
    height = 6) # The height of the plot in inches
a+b
dev.off()







#Antifungal Mantel tests#####

antifungaldata<-qza_to_phyloseq(
  features="Antifungal_Matches_99/clustered_table.qza", metadata = "metadata.tsv")



lonlat<-data.frame(sample_data(data_rarified)$Long,sample_data(data_rarified)$Lat)
colnames(lonlat)=c('longitude','latitutde')

dist_matrix<-geodist(lonlat,measure="geodesic")
rownames(dist_matrix)<-sample_names(data_rarified)
colnames(dist_matrix)<-sample_names(data_rarified)


#Find  jaccard, bray, UniFrac and weighted-UniFrac indices
jaccard<-vegdist(t(otu_table(antifungaldata)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(antifungaldata)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)

#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)


#The mantel tests
mantel_test_jaccard  <- mantel(jaccard_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_jaccard

mantel_test_bray  <- mantel(bray_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_bray


#Set Distance classes for correlogram analyses
ddd<-c()
ddd<-c(ddd,0)
ddd<-c(ddd,100)
ddd<-c(ddd,1000)
ddd<-c(ddd,5000)
ddd<-c(ddd,10000)
ddd<-c(ddd,15000)
ddd<-c(ddd,20000)
ddd<-c(ddd,25000)
ddd<-c(ddd,45000)



#Mantel Correlograms
mantel_correlog_jaccard<-mantel.correlog(jaccard_matrix,dist_matrix,nperm=9999,cutoff=FALSE, break.pts=ddd)
plot(mantel_correlog_jaccard, alpha=0.05)




mantel_correlog_bray<-mantel.correlog(bray_matrix,dist_matrix,nperm=9999,cutoff=FALSE, break.pts=ddd)

plot(mantel_correlog_bray, alpha=0.05)






