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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate a variety of different ALPHA diversity metrics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Find the richness for each sample
richness<-colSums(sign(otu_table(data_rarified)))

#Find the shannon diversity for each sample
shannon<-diversity(otu_table(t(data_rarified)),index='shannon')

#Find Faith's pd for each sample
faithpd<-pd(t(otu_table(data_rarified)), phy_tree(data_rarified), include.root=TRUE)
faiths<-faithpd[,1]

type<-sample_data(data_rarified)$Sample_Type

#Put all of your different diversity metrics into a dataframe
diversity_df<-data.frame(type,richness,shannon,faiths)
#Switch the ordering of the types (to the order you want them presented in your plots)
diversity_df$type<-factor(diversity_df$type,levels=c('Salamander','Crevice'))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test whether there are differences in ALPHA diversity between groups

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Richness
w_richness<-wilcox.test(richness ~ type, data = diversity_df)
w_richness

#shannon Entropy
w_shannon<-wilcox.test(shannon ~ type, data = diversity_df)
w_shannon

#Faith's PD
w_faiths<-wilcox.test(faiths ~ type, data = diversity_df)
w_faiths

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Make Violin plots of the various diversity metrics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


##########################Richness Violin Plot#########################################################################################################################################

#Define violin plot
p_richness <- ggplot(diversity_df, aes(x=type, y=richness)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_richness<-p_richness+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))
#Choose the size of font for the legend title and lables
p_richness<-p_richness+theme(legend.title = element_text(size = 25))+theme(legend.text = element_text(size = 20))
#Choose the violin colors for each group
p_richness<-p_richness+scale_fill_manual(values=c("green", "darkgrey"))+labs(x= "Location", y = "Richness")
#Add boxplots inside the violins
p_richness<-p_richness+geom_boxplot(aes(fill=type),width=0.2)+geom_point(position = position_jitter(seed = 1, width = 0.25, )) +theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#Add the p-value for the Wilcox test somewhere on your figure (you may have to change the x and y positions of this label)
#p_richness<-p_richness+annotate("text", size=6,x=1.5, y=1100, label= "Wilcox Test")
#p_richness<-p_richness+annotate("text", size=6,x=1.5, y=1025, label= paste("p = ",round(w_richness$p.value,3)))


#Make your plot
p_richness
pdf(file = "Figure_3A.pdf",
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
p_richness
dev.off()

##########################Shannon Diversity Violin Plot#########################################################################################################################################

#Define violin plot
p_shannon <- ggplot(diversity_df, aes(x=type, y=shannon)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_shannon<-p_shannon+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))
#Choose the size of font for the legend title and lables
p_shannon<-p_shannon+theme(legend.title = element_text(size = 25))+theme(legend.text = element_text(size = 20))
#Choose the violin colors for each group
p_shannon<-p_shannon+scale_fill_manual(values=c("green", "darkgrey"))
#Add boxplots inside the violins
p_shannon<-p_shannon+geom_boxplot(aes(fill=type),width=0.2)+geom_point(position = position_jitter(seed = 1, width = 0.25, )) +labs(x= "Location", y = "shannon Diversity")+
  theme_classic()
#Add the p-value for the Wilcox test somewhere on your figure (you may have to change the x and y positions of this label)
#p_shannon<-p_shannon+annotate("text", size=6,x=1.5, y=7, label= "Wilcox Test")
#p_shannon<-p_shannon+annotate("text", size=6,x=1.5, y=6.5, label= paste("p = ",round(w_shannon$p.value,3)))


#Make your plot
p_shannon

#####Faith's PD Violin Plot#########################################################################################################################################

#Define violin plot
p_faiths <- ggplot(diversity_df, aes(x=type, y=faiths)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_faiths<-p_faiths+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))
#Choose the size of font for the legend title and lables
p_faiths<-p_faiths+theme(legend.title = element_text(size = 25))+theme(legend.text = element_text(size = 20))
#Choose the violin colors for each group
p_faiths<-p_faiths+scale_fill_manual(values=c("green", "darkgrey"))
#Add boxplots inside the violins
p_faiths<-p_faiths+geom_boxplot(aes(fill=type),width=0.2)+geom_point(position = position_jitter(seed = 1, width = 0.25, )) +labs(x= "Location", y = "Faith's PD")+
  theme_classic()
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
#p_faiths<-p_faiths+annotate("text", size=6,x=1.5, y=65, label= "Wilcox Test")
#p_faiths<-p_faiths+annotate("text", size=6,x=1.5, y=60, label= paste("p = ",round(w_faiths$p.value,3)))


#Make your plot
p_faiths


#Remove legends
p_richness<-p_richness + guides(fill="none")
p_shannon<-p_shannon + guides(fill="none")
p_faiths<-p_faiths+ guides(fill="none")
#####SI Figure####
pdf(file = "SI_Alpha_Crevice_v_Salamander_ASV.pdf",
    width = 10, # The width of the plot in inches
    height = 6) # The height of the plot in inches

p_richness+p_shannon+p_faiths

dev.off()

####Figure 3B####

#                 Set your color scheme and the type of each sample in a vector

colvec<-c("green", "darkgrey")
type<-sample_data(data_rarified)$Sample_Type
types<-c('Salamander','Crevice')
#Convert names of groups to numbers
no_type<-rep(1,length(type))
for (k in 1:length(types)){
  no_type[which(type==types[k])]<-k
}

set.seed(1)


#                 Find a variety of dissimilarity indices

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



#                 Perform 2D PCA


#Run PCA analysis
pca <- prcomp(t(otu_table(data_rarified)))

#Find variation explained by axes
variance <- (pca$sdev)^2
varPercent <- variance/sum(variance) * 100

#Find loadings on axes
loadings <- pca$rotation

#Find coordinates for plotting (identical to PCoA on Euclidean distances)
scores <- pca$x
scores_2<-scores[,1:2]

pcoa_euclidean<-pcoa(euclidean)$vectors[,1:2]
pcoa_jaccard<-pcoa(jaccard)$vectors[,1:2]
pcoa_bray<-pcoa(bray)$vectors[,1:2]
pcoa_unifrac<-pcoa(unifrac)$vectors[,1:2]
pcoa_wunifrac<-pcoa(wunifrac)$vectors[,1:2]


#           Plot the PCA

par(mar = c(5, 5,5, 5))
plot(pcoa_euclidean,type='n',cex.lab=1.75,xlab='PCA1',ylab='PCA2')
ordihull(pcoa_euclidean,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_euclidean, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_euclidean[which(!(rownames(pcoa_euclidean) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")

#Find the centroid (based on mean) of each treatment group in the PCoA plot
xmeans_euclidean<-c()    #Value of the centroid along the first PCoA axis
ymeans_euclidean<-c()    #Value of the centroid along the second PCoA axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_euclidean<-c(xmeans_euclidean,mean(pcoa_euclidean[pts,1]))
  ymeans_euclidean<-c(ymeans_euclidean,mean(pcoa_euclidean[pts,2]))
}

#Plot the centroid for each group overtop of your PCoA scatterplot
for (k in 1:length(types)){
  points(xmeans_euclidean[k],ymeans_euclidean[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_euclidean[k],ymeans_euclidean[k],col='black',pch=1,cex=2)
}



#           Plot the PCoA Jaccard

#Plot the Jaccard PCoA 
plot(pcoa_jaccard,type='n',cex.lab=1.75,xlab='PCoA1 (3.1%)',ylab='PCoA2 (2.8%)')
ordihull(pcoa_jaccard,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_jaccard, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_jaccard[which(!(rownames(pcoa_jaccard) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")

#Find the centroid (based on mean) of each treatment group in the PCoA plot
xmeans_jaccard<-c()    #Value of the centroid along the first PCoA axis
ymeans_jaccard<-c()    #Value of the centroid along the second PCoA axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_jaccard<-c(xmeans_jaccard,mean(pcoa_jaccard[pts,1]))
  ymeans_jaccard<-c(ymeans_jaccard,mean(pcoa_jaccard[pts,2]))
}

#Plot the centroid for each group overtop of your PCOA scatterplot
for (k in 1:length(types)){
  points(xmeans_jaccard[k],ymeans_jaccard[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_jaccard[k],ymeans_jaccard[k],col='black',pch=1,cex=2)
}



pdf(file = "Figure_3B.pdf",
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

#Plot the Jaccard PCoA 
plot(pcoa_jaccard,type='n',cex.lab=1,xlab='PCoA1 (3.1%)',ylab='PCoA2 (2.8%)')
ordihull(pcoa_jaccard,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_jaccard, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_jaccard[which(!(rownames(pcoa_jaccard) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")
#legend("bottomright", title="Location",
            # c("Salamander","Crevice"), fill=c("green","gray"), horiz=TRUE, cex=.75)

#Find the centroid (based on mean) of each treatment group in the PCoA plot
xmeans_jaccard<-c()    #Value of the centroid along the first PCoA axis
ymeans_jaccard<-c()    #Value of the centroid along the second PCoA axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_jaccard<-c(xmeans_jaccard,mean(pcoa_jaccard[pts,1]))
  ymeans_jaccard<-c(ymeans_jaccard,mean(pcoa_jaccard[pts,2]))
}

#Plot the centroid for each group overtop of your PCOA scatterplot
for (k in 1:length(types)){
  points(xmeans_jaccard[k],ymeans_jaccard[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_jaccard[k],ymeans_jaccard[k],col='black',pch=1,cex=2)
}
dev.off()

#Find variation explained by axes, look at the $values portion of pcoa_jaccard when run with not just the first 2 axes. the relative_eig gives the %variation explained by each axis.


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCoA Bray-Curtis

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the Bray PCoA 
pdf(file = "SI_Bray_Curtis_Crevice_v_Salamander.pdf",  
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

plot(pcoa_bray,type='n',cex.lab=1.5,xlab='PCoA1 (8.6%)',ylab='PCoA2 (4.9%')
ordihull(pcoa_bray,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_bray, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_bray[which(!(rownames(pcoa_bray) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")
#legend("bottomright", title="Location",
      # c("Salamander","Crevice"), fill=c("green","gray"), horiz=TRUE, cex=.75)


#Find the centroid (based on mean) of each treatment group in the PCoA plot
xmeans_bray<-c()    #Value of the centroid along the first PCoA axis
ymeans_bray<-c()    #Value of the centroid along the second PCoA axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_bray<-c(xmeans_bray,mean(pcoa_bray[pts,1]))
  ymeans_bray<-c(ymeans_bray,mean(pcoa_bray[pts,2]))
}

#Plot the centroid for each group overtop of your PCoA scatterplot
for (k in 1:length(types)){
  points(xmeans_bray[k],ymeans_bray[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_bray[k],ymeans_bray[k],col='black',pch=1,cex=2)
}


while (!is.null(dev.list()))  dev.off()

#           Plot the PCoA UniFrac

pdf(file = "SI_Unifrac_Crevice_v_Salamander.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

#Plot the Unifrac PCoA 
plot(pcoa_unifrac,type='n',cex.lab=1.5,xlab='PCoA1 (5.3%)',ylab='PCoA2 (2.6%)')
ordihull(pcoa_unifrac,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_unifrac, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_unifrac[which(!(rownames(pcoa_unifrac) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")
legend("bottomright", title="Location",
       c("Salamander","Crevice"), fill=c("green","gray"), horiz=TRUE, cex=.75)

#Find the centroid (based on mean) of each treatment group in the PCoA plot
xmeans_unifrac<-c()    #Value of the centroid along the first PCoA axis
ymeans_unifrac<-c()    #Value of the centroid along the second PCoA axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_unifrac<-c(xmeans_unifrac,mean(pcoa_unifrac[pts,1]))
  ymeans_unifrac<-c(ymeans_unifrac,mean(pcoa_unifrac[pts,2]))
}

#Plot the centroid for each group overtop of your PCoA scatterplot
for (k in 1:length(types)){
  points(xmeans_unifrac[k],ymeans_unifrac[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_unifrac[k],ymeans_unifrac[k],col='black',pch=1,cex=2)
}
dev.off()

#           Plot the PCoA weighted-UniFrac



pdf(file = "SI_W_Unifrac_Crevice_v_Salamander.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches
#Plot the weighted Unifrac PCoA 
plot(pcoa_wunifrac,type='n',cex.lab=1.5,xlab='PCoA1 (49%)',ylab='PCoA2 (7.7%)')
ordihull(pcoa_wunifrac,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_wunifrac, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_wunifrac[which(!(rownames(pcoa_wunifrac) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")
legend("bottomright", title="Location",
       c("Salamander","Crevice"), fill=c("green","gray"), horiz=TRUE, cex=.75)


#Find the centroid (based on mean) of each treatment group in the PCoA plot
xmeans_wunifrac<-c()    #Value of the centroid along the first PCoA axis
ymeans_wunifrac<-c()    #Value of the centroid along the second PCoA axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_wunifrac<-c(xmeans_wunifrac,mean(pcoa_wunifrac[pts,1]))
  ymeans_wunifrac<-c(ymeans_wunifrac,mean(pcoa_wunifrac[pts,2]))
}

#Plot the centroid for each group overtop of your PCoA scatterplot
for (k in 1:length(types)){
  points(xmeans_wunifrac[k],ymeans_wunifrac[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_wunifrac[k],ymeans_wunifrac[k],col='black',pch=1,cex=2)
}
dev.off()


#           Perform Jaccard PERMANOVA


#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_jaccard<-list(Data=t(otu_table(data_rarified)),D=jaccard_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_jaccard=PERMANOVA(inputpermanova_jaccard, factor(type),nperm=1000)





#                 Perform PERMANOVA Tests on Bray Curtis Distances



#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_bray<-list(Data=t(otu_table(data_rarified)),D=bray_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_bray=PERMANOVA(inputpermanova_bray, factor(type),nperm=1000)



#                 Perform PERMANOVA Tests on UniFrac Distances



#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_unifrac<-list(Data=t(otu_table(data_rarified)),D=unifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_unifrac=PERMANOVA(inputpermanova_unifrac, factor(type),nperm=1000)



#                 Perform PERMANOVA Tests on weighted UniFrac Distances



#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_wunifrac<-list(Data=t(otu_table(data_rarified)),D=wunifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_wunifrac=PERMANOVA(inputpermanova_wunifrac, factor(type),nperm=1000)



#                 Summarize PERMANOVA Tests 

Fstats<-c(permanova_by_group_jaccard$Initial$Global[5],permanova_by_group_bray$Initial$Global[5],permanova_by_group_unifrac$Initial$Global[5],permanova_by_group_wunifrac$Initial$Global[5])
pvalues<-c(permanova_by_group_jaccard$pvalue,permanova_by_group_bray$pvalue,permanova_by_group_unifrac$pvalue,permanova_by_group_wunifrac$pvalue)

permanova_stats<-data.frame(Fstats,pvalues)
rownames(permanova_stats)<-c('Jaccard','Bray-Curtis','Unifrac','weighted Unifrac')


#                 Test for homogeneity of variances



bd<-betadisper(euclidean,as.factor(type))
dispersion_anova_euclidean<-anova(bd)

bd<-betadisper(jaccard,as.factor(type))
dispersion_anova_jaccard<-anova(bd)

bd<-betadisper(bray,as.factor(type))
dispersion_anova_bray<-anova(bd)

bd<-betadisper(unifrac,as.factor(type))
dispersion_anova_unifrac<-anova(bd)

bd<-betadisper(wunifrac,as.factor(type))
dispersion_anova_wunifrac<-anova(bd)

Fs<-c(dispersion_anova_euclidean$`F value`[1],dispersion_anova_jaccard$`F value`[1],dispersion_anova_bray$`F value`[1],dispersion_anova_unifrac$`F value`[1],dispersion_anova_wunifrac$`F value`[1])
ps<-c(dispersion_anova_euclidean$`Pr(>F)`[1],dispersion_anova_jaccard$`Pr(>F)`[1],dispersion_anova_bray$`Pr(>F)`[1],dispersion_anova_unifrac$`Pr(>F)`[1],dispersion_anova_wunifrac$`Pr(>F)`[1])

dispersion_stats<-data.frame(Fs,ps)
rownames(dispersion_stats)<-c('Euclidean','Jaccard','Bray-Curtis','Unifrac','weighted Unifrac')

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



