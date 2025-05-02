################## Green Salamander Metapopulation Study Script Figure 3 ######################
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

data_rarified <- merge_phyloseq(data_rarified,Crevicesdata)



#####Figure 3C##### ######

#ANCOM Crevice vs Salamander Microbiomes
data_rarified<-merge_phyloseq(data_rarified,Crevicesdata)
unload("phyloseq")

#ANCOMBC2 will not run with ASVs as is,it requires something taxonomic, so I am going to change the species names of all my ASVs to the actual ASVs.

df<-c()
df<-as.data.frame(tax_table(data_rarified))
df$Species<-c()

library(tibble)
df <- tibble::rownames_to_column(df, "Species")
rownames(df)<-df$Species
df2 <- df[, c("Kingdom", "Phylum" , "Class",   "Order",   "Family",  "Genus", "Species" )]
data_rarified@tax_table<-tax_table(as.matrix(df2))
rownames(tax_table(data_rarified))<-df2$Species

library(ANCOMBC)
Type_ANCOM<-ancombc2(data_rarified, group = "Sample_Type", p_adj_method ="holm", fix_formula = "Sample_Type", tax_level = "Species",verbose = TRUE)

#saveRDS(Type_ANCOM,file = "Type_ANCOM.rds")

output<-readRDS(file = "Type_ANCOM.rds")

library(data.table)
library(tidyverse)
library(DT)

res_prim = output$res

setDT(res_prim, keep.rownames = TRUE)[]

names(res_prim)[1] <- "taxon"
df_type = res_prim 

COLORS = c("Positive LFC" = "green", "Negative LFC" = "gray")

df_fig_type = df_type %>%
  dplyr::filter(diff_Sample_TypeSalamander == "TRUE") %>%
  dplyr::arrange(desc(lfc_Sample_TypeSalamander)) %>%
  dplyr::mutate(direct = ifelse(lfc_Sample_TypeSalamander > 0, "Positive LFC", "Negative LFC"),
                color = "black")




df_fig_type$taxon = factor(df_fig_type$taxon, levels = df_fig_type$taxon)

ASV_to_Genera<-df_fig_type$taxon

genera_collapsed<-c("Methylobacterium-Methylorubrum", "Xanthobacteraceae", "Luteolibacter", "Nordella", "Pseudonocardia", "uncultured_Sandaracinaceae", "Rudaea", "Paenibacillus", "Dokdonella", "uncultured_Ktedonobacteraceae", "Nocardia", "Mesorhizobium", "Chitinophaga", "Rhodanobacter", "Reyranella", "Nitrososphaeraceae", "Candidatus_Nitrocosmicus")

df_fig_type$taxon = factor(genera_collapsed, levels = genera_collapsed)

df_fig_type$direct = factor(df_fig_type$direct, 
                            levels = c("Positive LFC", "Negative LFC"))


fig_type = df_fig_type %>%
  ggplot(aes(x = taxon, y = lfc_Sample_TypeSalamander, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_Sample_TypeSalamander - se_Sample_TypeSalamander, ymax = lfc_Sample_TypeSalamander + se_Sample_TypeSalamander), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = "Taxon", y = "Log fold change", 
       title = "") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,
                                   color = df_fig_type$color))+theme(text = element_text(size=50))+  scale_fill_manual(values = COLORS, name = NULL)

fig_type




pdf(file = "ANCOM Crevice v Salamander.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
fig_type
dev.off()






####Figure 3D #####

######Indicator Species Analyses #####
library('indicspecies')
library('data.table')

#Indicator for Salamander vs Environment


set.seed(1)

type_environment<-sample_data(data_rarified)$Sample_Type

#indval_environment = multipatt(t(otu_table(data_rarified)), type_environment, control = how(nperm=50000), duleg=TRUE)

#saveRDS(indval_environment, file = "indval_environment_final.rds")
indval_environment<-readRDS(file = "indval_environment_final.rds")
summary(indval_environment,indvalcomp=TRUE)
#extract table of stats
indisp.sign<-as.data.table(indval_environment$sign, keep.rownames=TRUE)
abundant_taxa<-rownames(otu_table(data_rarified))[which(apply(otu_table(data_rarified),MARGIN=1,FUN=max)/colSums(otu_table(data_rarified))[1]>0.01)]
indisp.sign<-indisp.sign[which(indisp.sign$rn %in% abundant_taxa),]
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
df_environment<-data.frame(indisp.sign[p.value.bh<=0.05, ])



id_names<-c()
genus_names<-c()
rn<-c()

#for (k in 1:length(taxonomy[["data"]]$Taxon)){
  namer<-strsplit(strsplit(taxonomy[["data"]]$Taxon[k],'g__')[[1]][2],';')[[1]][1]
  if (is.na(namer)){
    namer<-paste('undescribed', strsplit(strsplit(taxonomy[["data"]]$Taxon[k],'f__')[[1]][2],';')[[1]][1])
  }
  if (namer == 'uncultured'){
    namer2<-strsplit(strsplit(taxonomy[["data"]]$Taxon[k],'f__')[[1]][2],';')[[1]][1]
    if (namer2 == 'uncultured'){
      namer<-paste('uncultured', strsplit(strsplit(taxonomy[["data"]]$Taxon[k],'o__')[[1]][2],';')[[1]][1])
    }
    else{
      namer<-paste('uncultured', strsplit(strsplit(taxonomy[["data"]]$Taxon[k],'f__')[[1]][2],';')[[1]][1])
    }
  }
  if (namer == '1174-901-12'){
    namer<-paste('unclassified', strsplit(strsplit(taxonomy[["data"]]$Taxon[k],'f__')[[1]][2],';')[[1]][1])
  }
  if (namer == '1921-2'){
    namer<-paste('unclassified', strsplit(strsplit(taxonomy[["data"]]$Taxon[k],'f__')[[1]][2],';')[[1]][1])
  }
  if (namer == '67-14'){
    namer<-paste('unclassified', strsplit(strsplit(taxonomy[["data"]]$Taxon[k],'f__')[[1]][2],';')[[1]][1])
  }
  if (namer == 'JG30-KF-AS9'){
    namer<-paste('unclassified', strsplit(strsplit(taxonomy[["data"]]$Taxon[k],'f__')[[1]][2],';')[[1]][1])
  }
  
genus_names<-c(genus_names,namer)id_names<-c(id_names,taxonomy[["data"]]$Feature.ID[k])
taxa_info<-data.frame(id_names, genus_names)


colnames(taxa_info)[1] <- "rn"
df_environment<-inner_join(df_environment,taxa_info, by = "rn")



#write.csv(df_environment,'Environment_ASV_indicators.csv')
# Excel formula for coding Crevice v Salamander in column called Location=IF(E2=1,"Crevice","Salamander")
df_environment<-read.csv("Environment_ASV_indicators.csv")


df_environment<-subset(df_environment, p.value.bh<=0.05)

colnames(df_environment)[colnames(df_environment) == 'P.value'] <- 'P_value'
colnames(df_environment)[colnames(df_environment) == 'stat'] <- 'Indicator_value'



length(which(df_environment$Location =="Salamander")) #7 ASVs are indicative of a salamander microbiome
length(which(df_environment$Location =="Crevice")) #27 ASVs are indicative of a crevice microbiome



pdf(file = "Env_v_Sal_Indicator.pdf", width = 15, height = 10)



d<-ggplot(df_environment, aes(Location, genus_names)) +
  geom_point(aes(color = p.value.bh, size = Indicator_value)) +
  scale_color_gradient(low = "red", high = "gold") +
  scale_x_discrete(position = "top") +
  scale_size(range = c(2, 10)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size=10),legend.title = element_text(size=14), legend.key.width = unit(1, 'cm'))+
  theme(text = element_text(size=25))+ labs(title = "")
d
dev.off()



