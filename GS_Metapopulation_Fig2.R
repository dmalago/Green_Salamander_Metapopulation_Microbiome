################## Spatial Variation of Skin-associated Microbiota in a Green Salamander Metapopulation Figure 2 Script ######################

#####Loading required libraries, setting working directory#####
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


#####Importing data####
SVs<-read_qza("ASV/table.qza")
SVs$type
taxonomy<-read_qza("ASV/taxonomy.qza")

#Combine different objects into 1 phyloseq object
data<-qza_to_phyloseq(
  features="ASV/table.qza",
  "ASV/taxonomy.qza",
  tree="rooted-tree.qza",
  
  metadata = "metadata.txt")

#Remove uncharacterized phyla
data <- subset_taxa(data, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove phyla that aren't bacterial
data <- subset_taxa(data, Phylum!="Chytridiomycota"& Phylum !="Ascomycota" & Phylum != "Chlorophyta" & Phylum != "Apicomplexa" & Phylum != "Ciliophora"  & Phylum != "Phragmoplastophyta" & Phylum != "Zoopagomycota")

#Remove taxa that have chloroplast and mitochondria at level 5, 6,7 
data<-subset_taxa(data, Family != "Chloroplast")
data<-subset_taxa(data, Genus != "Mitochondria")

#Creating function for removing ASVs found in sequenced blanks, this comes from creater of Phyloseq package
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

#Creating a new Phyloseq object removing Bd positive salamanders, Bdfix column includes individuals positive for Bd and/or displaying clinical signs of chytridiomycosis
data = subset_samples(data, Bdfix==0)



#Rarefy your data, here I am rarefying to 5000 sequences per sample
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=5000, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying, removing samples from previously swabbed individuals, and creating a crevice microbiota dataset
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
Crevicesdata<-subset_samples(data_rarified, Sample_Type == "Crevice")
data_rarified=subset_samples(data_rarified, Repeat == "0")


write_rds(data_rarified, file ="data_rarified.rds")

data_rarified <- merge_phyloseq(data_rarified,Crevicesdata)

#Creating Taxonomic Bar Plots####
#Find the total abundance of each phylum
summed_phyla<-rowMeans(data.frame(otu_table(tax_glom(data_rarified,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list<-data.frame(tax_table(tax_glom(data_rarified,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant... we will explicitly show the four most abundant
sorted_phyla_list<-phyla_list[order(-summed_phyla)]

#Prepare the OTU table in the correct format for the plotting program (this must be done with level 6 data including genera!)
mdf_prep <- prep_mdf(data_rarified)
#Tell the function to plot the five most abundant phyla (you could pick something different... )
color_objs_GP <- create_color_dfs(mdf_prep,selected_groups = c("Proteobacteria", "Actinobacteriota", "Acidobacteriota","Planctomycetota","Bacteroidota"),cvd=TRUE)
#Extract the OTU table and color choices (again, putting things in the right format for the function)
mdf_GP <- color_objs_GP$mdf
cdf_GP <- color_objs_GP$cdf
cdf_GP <- color_reassign(cdf_GP,group_assignment = c("Proteobacteria", "Actinobacteriota", "Acidobacteriota","Planctomycetota","Bacteroidota"),color_assignment = c("micro_cvd_green", "micro_cvd_orange","micro_cvd_blue","micro_cvd_purple","micro_cvd_turquoise"))

#Define the legend
GP_legend <-custom_legend(mdf_GP, cdf_GP,legend_text_size = 14)

#Define the plot that you will be making
plot <- plot_microshades(mdf_GP, cdf_GP)
#Define all of the formatting aspects of the plot
plot_diff <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 6)) +
  facet_grid(~fct_relevel(Site,'HW1','BB','DNR','1250','1251','1292','TR2','3688','1477','Crevice'), scale="free_x", space = "free_x") +
  theme(axis.text.x = element_text(size= 6)) +
  theme(plot.margin = margin(6,20,6,6))+theme(text = element_text(size=25))+ylab("Relative Abundance")+theme(strip.background = element_rect(fill = "lightgrey"),strip.text = element_text(size = 8)) +  
  theme(panel.spacing.x = unit(.1, "cm")) 

#Make the plot
plot_grid(plot_diff, GP_legend,  rel_widths = c(.8, .2))


pdf(file = "Figure_2.pdf",   
    width = 18, # The width of the plot in inches
    height = 8) # The height of the plot in inches

plot_grid(plot_diff, GP_legend,  rel_widths = c(.8, .2))
dev.off()


######Subset Data and Create Data Files#####
####This will divide your phyloseq object into 2, one with all salamander samples and one with environment
Salamanders<-subset_samples(data_rarified, Sample_Type == "Salamander")
Environment<-subset_samples(data_rarified, Sample_Type != "Salamander")
Crevicesdata<-subset_samples(data_rarified, Sample_Type == "Crevice")



#Find the most abundant phyla on salamanders and crevices
summed_phyla_salamanders<-rowMeans(data.frame(otu_table(tax_glom(Salamanders,taxrank="Phylum"))))
#Find the names of the different phyla
phyla_list_salamanders<-data.frame(tax_table(tax_glom(Salamanders,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundant... we will explicitly show the four most abundant
sorted_phyla_list_salamanders<-phyla_list_salamanders[order(-summed_phyla_salamanders)][1:5]
sorted_phyla_percentages_salamanders<-summed_phyla_salamanders[order(-summed_phyla_salamanders)][1:5]*100/5000
data.frame(sorted_phyla_list_salamanders,sorted_phyla_percentages_salamanders)

summed_phyla_crevices<-rowMeans(data.frame(otu_table(tax_glom(Crevicesdata,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list_crevices<-data.frame(tax_table(tax_glom(Crevicesdata,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_phyla_list_crevices<-phyla_list_crevices[order(-summed_phyla_crevices)][1:5]
sorted_phyla_percentages_crevices<-summed_phyla_crevices[order(-summed_phyla_crevices)][1:5]*100/5000
data.frame(sorted_phyla_list_crevices,sorted_phyla_percentages_crevices)



#Find the most abundant genera on salamanders and crevices
summed_genus_salamanders<-rowMeans(data.frame(otu_table(tax_glom(Salamanders,taxrank="Genus"))))
#Find the names of the different phylum
genus_list_salamanders<-data.frame(tax_table(tax_glom(Salamanders,taxrank="Genus")))[,6]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_genus_list_salamanders<-genus_list_salamanders[order(-summed_genus_salamanders)][1:5]
sorted_genus_percentages_salamanders<-summed_genus_salamanders[order(-summed_genus_salamanders)][1:5]*100/5000
data.frame(sorted_genus_list_salamanders,sorted_genus_percentages_salamanders)

summed_genus_crevices<-rowMeans(data.frame(otu_table(tax_glom(Crevicesdata,taxrank="Genus"))))
#Find the names of the different phyla
genus_list_crevices<-data.frame(tax_table(tax_glom(Crevicesdata,taxrank="Genus")))[,6]
family_list_crevices<-data.frame(tax_table(tax_glom(Crevicesdata,taxrank="Genus")))[,5]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_genus_list_crevices<-genus_list_crevices[order(-summed_genus_crevices)][1:5]
sorted_family_list_crevices<-family_list_crevices[order(-summed_genus_crevices)][1:5]
sorted_genus_percentages_crevices<-summed_genus_crevices[order(-summed_genus_crevices)][1:5]*100/5000
data.frame(sorted_genus_list_crevices,sorted_genus_percentages_crevices)

