################## Green Salamander Metapopulation Study Script Figure 5 ######################
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

type<-sample_data(data_rarified)$Site
year<-sample_data(data_rarified)$Year

#Put all of your different diversity metrics into a dataframe
diversity_df<-data.frame(type,year,richness,shannon,faiths)
#Switch the ordering of the types (to the order you want them presented in your plots)
diversity_df$type<-factor(diversity_df$type,levels=c('HW1','BB','DNR','1250','1251','1292','TR2','3688','1477'))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test whether there are differences in ALPHA diversity between groups

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Richness
kw_richness<-kruskal.test(richness ~ type, data = diversity_df)
if (kw_richness$p.value<0.05){
  print('Alpha Richness')
  pw_richness<-pairwise.wilcox.test(diversity_df$richness, diversity_df$type,p.adjust.method = "BH")
  print(pw_richness)
  pw_richness_nocorrect<-pairwise.wilcox.test(diversity_df$richness, diversity_df$type,p.adjust.method = "none")
}

#Shannon Entropy
kw_shannon<-kruskal.test(shannon ~ type, data = diversity_df)
if (kw_shannon$p.value<0.05){
  print('Alpha Shannon')
  pw_shannon<-pairwise.wilcox.test(diversity_df$shannon, diversity_df$type,p.adjust.method = "BH")
  print(pw_shannon)
}

#Faith's PD
kw_faiths<-kruskal.test(faiths ~ type, data = diversity_df)
if (kw_faiths$p.value<0.05){
  print('Alpha Faiths')
  pw_faiths<-pairwise.wilcox.test(diversity_df$faiths, diversity_df$type,p.adjust.method = "BH")
  pw_faiths_nocorrect<-pairwise.wilcox.test(diversity_df$faiths, diversity_df$type,p.adjust.method = "none")
  print(pw_faiths)
}


ff<-glm(richness~type+as.factor(year),data=diversity_df)
gg<-stepAIC(ff)
summary(gg)

ff<-glm(shannon~type+as.factor(year),data=diversity_df)
gg<-stepAIC(ff)
summary(gg)

ff<-glm(richness~type+as.factor(year),data=diversity_df)
gg<-stepAIC(ff)
summary(gg)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Make Violin plots of the various diversity metrics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Do you want to show p-values for non-significant pairwise comparisons? (yes or no)
p_value_ns<-'no'

##########################Richness Violin Plot#########################################################################################################################################

#Define violin plot
p_richness <- ggplot(diversity_df, aes(x=type, y=richness)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_richness<-p_richness+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_richness<-p_richness+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_richness<-p_richness+scale_fill_manual(values=c('pink','yellow','brown','red','blue','green','grey','orange',"purple"))
#Add boxplots inside the violins
p_richness<-p_richness+geom_boxplot(aes(fill=type),width=0.2)+geom_point(position = position_jitter(seed = 1, width = 0.1, )) +theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
#p_richness<-p_richness+annotate("text", size=6,x=2, y=100, label= "Kruskal-Wallis")
#p_richness<-p_richness+annotate("text", size=6,x=2, y=70, label= paste("p = ",round(kw_richness$p.value,4)))
p_richness<-p_richness + labs(x = "Population")
p_richness<-p_richness + labs(fill = "Population")+labs(y = "ASV Richness")+theme(text = element_text(size=25))


 

#Make your plot
pdf(file = "Figure_4B.pdf",   
    width = 10, # The width of the plot in inches
    height = 8) # The height of the plot in inches
p_richness
dev.off()



##########################Shannon Diversity Violin Plot#########################################################################################################################################

#Define violin plot
p_shannon <- ggplot(diversity_df, aes(x=type, y=shannon)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_shannon<-p_shannon+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_shannon<-p_shannon+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_shannon<-p_shannon+scale_fill_manual(values=c('pink','yellow','brown','red','blue','green','grey','orange',"purple"))
#Add boxplots inside the violins
p_shannon<-p_shannon+geom_boxplot(aes(fill=type),width=0.2)+geom_point(position = position_jitter(seed = 1, width = 0.1, )) +
  theme_classic()+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_shannon<-p_shannon + labs(fill = "Population")+labs(x = "Population")
p_shannon<-p_shannon +labs(y = "Shannon Diversity")+theme(text = element_text(size=25))

#Make your plot
p_shannon

##########################Faith's PD Violin Plot#########################################################################################################################################

#Define violin plot
p_faiths <- ggplot(diversity_df, aes(x=type, y=faiths)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_faiths<-p_faiths+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 15))
#Choose the size of font for the legend title and lables
p_faiths<-p_faiths+theme(legend.title = element_text(size = 20))+theme(legend.text = element_text(size = 15))
#Choose the violin colors for each group
p_faiths<-p_faiths+scale_fill_manual(values=c('pink','yellow','brown','red','blue','green','grey','orange',"purple"))
#Add boxplots inside the violins
p_faiths<-p_faiths+geom_boxplot(aes(fill=type),width=0.2)+geom_point(position = position_jitter(seed = 1, width = 0.1, )) + labs(fill = "Population")+labs(x = "Population")+
  theme_classic()+labs(y = "Faith's PD")+theme(text = element_text(size=25))+theme_linedraw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+theme(text = element_text(size=25))
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)


#Make your plot
p_faiths


########SI Figure 5, ########
pdf(file = "SI_Figure_5.pdf", width = 27, height = 8)
p_richness + p_shannon + p_faiths
dev.off()

######Figure 4A #######
####Only consider the salamanders (not the crevices)

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




data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=5000, replace=F, verbose = TRUE)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
Crevicesdata<-subset_samples(data_rarified, Sample_Type == "Crevice")

data_rarified=subset_samples(data_rarified, Repeat == "0")
data_rarified<-subset_samples(data_rarified, Sample_Type == "Salamander")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Set your color scheme and the type of each sample in a vector

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colvec<-c('pink','yellow','brown','red','blue','green','grey','orange',"purple")
type<-sample_data(data_rarified)$Site
types<-c('HW1','BB','DNR','1250','1251','1292','TR2','3688','1477')
#Convert names of groups to numbers
no_type<-rep(1,length(type))
for (k in 1:length(types)){
  no_type[which(type==types[k])]<-k
}

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

#                 Perform 2D PCA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCoA Jaccard

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make your plot
pdf(file = "Figure_4C.pdf",   
    width = 9, # The width of the plot in inches
    height = 7) # The height of the plot in inches

#Plot the Jaccard PCoA 
par(mar=c(4.1, 4.1, 4.1, 4.1), xpd=TRUE)
plot(pcoa_jaccard,type='n',cex.lab=1.5,xlab='PCoA1 (3.7%)',ylab='PCoA2 (2.9%)')
ordihull(pcoa_jaccard,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_jaccard, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_jaccard[which(!(rownames(pcoa_jaccard) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")
legend("topright", title="Population",
       c('HW1','BB','DNR','1250','1251','1292','TR2','3688','1477'),inset=c(-.7,0), fill=c(colvec), cex=1.25,horiz= FALSE, ncol = 1)


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

#Plot the centroid for each group overtop of your PCoA scatterplot
for (k in 1:length(types)){
  points(xmeans_jaccard[k],ymeans_jaccard[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_jaccard[k],ymeans_jaccard[k],col='black',pch=1,cex=2)
}

dev.off()
#Find variation explained by axes, look at the $values portion of pcoa_jaccard when run with not just the first 2 axes. the relative_eig gives the %variation explained by each axis.

#pcoa_jaccard<-pcoa(jaccard)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCoA Bray-Curtis

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the Bray PCoA 
plot(pcoa_bray,type='n',cex.lab=1.75,xlab='PCoA1',ylab='PCoA2')
ordihull(pcoa_bray,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_bray, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_bray[which(!(rownames(pcoa_bray) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")


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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCoA UniFrac

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the Unifrac PCoA 
plot(pcoa_unifrac,type='n',cex.lab=1.75,xlab='PCoA1',ylab='PCoA2')
ordihull(pcoa_unifrac,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_unifrac, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_unifrac[which(!(rownames(pcoa_unifrac) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCoA weighted-UniFrac

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the weighted Unifrac PCoA 
plot(pcoa_wunifrac,type='n',cex.lab=1.75,xlab='PCoA1',ylab='PCoA2')
ordihull(pcoa_wunifrac,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_wunifrac, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_wunifrac[which(!(rownames(pcoa_wunifrac) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")


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

Year<-sample_data(data_rarified)$Year
Site<-sample_data(data_rarified)$Site

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Perform Euclidean PERMANOVA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_euclidean<-list(Data=t(otu_table(data_rarified)),D=euclidean_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_euclidean=PERMANOVA(inputpermanova_euclidean, factor(type),nperm=1000)

test<-adonis2(euclidean_matrix~Year+Site,by='margin')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Perform Jaccard PERMANOVA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_jaccard<-list(Data=t(otu_table(data_rarified)),D=jaccard_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_jaccard=PERMANOVA(inputpermanova_jaccard, factor(type),nperm=1000)

test<-adonis2(jaccard_matrix~Year+Site,by='margin')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Bray Curtis Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_bray<-list(Data=t(otu_table(data_rarified)),D=bray_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_bray=PERMANOVA(inputpermanova_bray, factor(type),nperm=1000)

test<-adonis2(bray_matrix~Year+Site,by='margin')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on UniFrac Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_unifrac<-list(Data=t(otu_table(data_rarified)),D=unifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_unifrac=PERMANOVA(inputpermanova_unifrac, factor(type),nperm=1000)

test<-adonis2(unifrac_matrix~Year+Site,by='margin')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on weighted UniFrac Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_wunifrac<-list(Data=t(otu_table(data_rarified)),D=wunifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_wunifrac=PERMANOVA(inputpermanova_wunifrac, factor(type),nperm=1000)

test<-adonis2(wunifrac_matrix~Year+Site,by='margin')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Summarize PERMANOVA Tests 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fstats<-c(permanova_by_group_euclidean$Initial$Global[5],permanova_by_group_jaccard$Initial$Global[5],permanova_by_group_bray$Initial$Global[5],permanova_by_group_unifrac$Initial$Global[5],permanova_by_group_wunifrac$Initial$Global[5])
pvalues<-c(permanova_by_group_euclidean$pvalue,permanova_by_group_jaccard$pvalue,permanova_by_group_bray$pvalue,permanova_by_group_unifrac$pvalue,permanova_by_group_wunifrac$pvalue)

permanova_stats<-data.frame(Fstats,pvalues)
rownames(permanova_stats)<-c('Euclidean','Jaccard','Bray-Curtis','Unifrac','weighted Unifrac')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test for homogeneity of variances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



####SI Figure 11A#####

pdf(file = "Figure_11A.pdf",   
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

#Plot the Jaccard PCoA 
plot(pcoa_jaccard,type='n',cex.lab=1.25,xlab='PCoA1 (3.7%)',ylab='PCoA2 (2.9%)')
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

#Plot the centroid for each group overtop of your PCoA scatterplot
for (k in 1:length(types)){
  points(xmeans_jaccard[k],ymeans_jaccard[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_jaccard[k],ymeans_jaccard[k],col='black',pch=1,cex=2)
}

dev.off()




####SI Figure 11B#####

pdf(file = "Figure_11B.pdf",   
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

plot(pcoa_bray,type='n',cex.lab=1.25,xlab='PCoA1 (10%)',ylab='PCoA2 (5%)')
ordihull(pcoa_bray,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_bray, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_bray[which(!(rownames(pcoa_bray) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")


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

dev.off()




####SI Figure 11C#####

pdf(file = "Figure_11C.pdf",   
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

#           Plot the PCoA UniFrac

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the Unifrac PCoA 
plot(pcoa_unifrac,type='n',cex.lab=1.25,xlab='PCoA1 (5.6%)',ylab='PCoA2 (3.1%)')
ordihull(pcoa_unifrac,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_unifrac, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_unifrac[which(!(rownames(pcoa_unifrac) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")

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










####SI Figure 11D#####

pdf(file = "Figure_11D.pdf",   
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

#           Plot the PCoA weighted-UniFrac

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the weighted Unifrac PCoA 
plot(pcoa_wunifrac,type='n',cex.lab=1.25,xlab='PCoA1 (30.6%)',ylab='PCoA2 (13.6%)')
ordihull(pcoa_wunifrac,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_wunifrac, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_wunifrac[which(!(rownames(pcoa_wunifrac) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")


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






######Indicator Species Analyses #####
library('indicspecies')
library('data.table')

#Indicator for Salamander vs Environment


set.seed(1)

type_pop<-sample_data(data_rarified)$Site

#indval_pop = multipatt(t(otu_table(data_rarified)), type_pop, control = how(nperm=50000), duleg=TRUE)

#saveRDS(indval_pop, file = "indval_pop_final.rds")
indval_pop<-readRDS(file = "indval_pop_final.rds")
summary(indval_pop,indvalcomp=TRUE)
#extract table of stats
indisp.sign<-as.data.table(indval_pop$sign, keep.rownames=TRUE)
abundant_taxa<-rownames(otu_table(data_rarified))[which(apply(otu_table(data_rarified),MARGIN=1,FUN=max)/colSums(otu_table(data_rarified))[1]>0.01)]
indisp.sign<-indisp.sign[which(indisp.sign$rn %in% abundant_taxa),]
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
df_pop<-data.frame(indisp.sign[p.value.bh<=0.05, ])



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
genus_names<-c(genus_names,namer)
id_names<-c(id_names,taxonomy[["data"]]$Feature.ID[k])

taxa_info<-data.frame(id_names, genus_names)


colnames(taxa_info)[1] <- "rn"
df_environment<-inner_join(df_environment,taxa_info, by = "rn")



#write.csv(df_environment,'Environment_ASV_indicators.csv')
# Excel formula for coding Crevice v Salamander in column called Location=IF(E2=1,"Crevice","Salamander")
df_environment<-read.csv("Environment_ASV_indicators.csv")


df_environment<-subset(df_environment, p.value.bh<=0.05)

colnames(df_environment)[colnames(df_environment) == 'P.value'] <- 'P_value'
colnames(df_environment)[colnames(df_environment) == 'stat'] <- 'Indicator_value'

#Figure 4C
poplist<-unique(type)
poplist<-poplist[poplist!='HW1']
pops<-list()
for (k in 1:length(poplist)){
  
  pops<-list.append(pops,sign(otu_table(prune_samples(sample_names(data_rarified)[which(sample_data(data_rarified)$Site==poplist[k])],data_rarified))))
  
}

names(pops)<-poplist

gtemp<-iNEXT3D(pops,diversity='TD',q=0,datatype='incidence_raw')



ggiNEXT3D(gtemp,type = 1, facet.var = "Order.q",color.var='Assemblage')+scale_shape_manual(values=c(19,19,19,19,19,19,19,19))+
  labs(title="",
       x ="Salamanders", y = "ASV Gamma Diversity")+
  scale_fill_manual(values=c('red','blue','green','purple','orange','yellow','brown','grey'))+
  scale_color_manual(values=c('red','blue','green','purple','orange','yellow','brown','grey'))+theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.position = "none"  # This line removes the legend
  )





