################## Green Salamander Metapopulation Study Script Figure 5 ######################

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
data_rarified = rarefy_even_depth(data, rngseed=1, sample.size=5000, replace=F, verbose = TRUE,trimOTUs = 0)

#Removing any taxa which aren't present following rarefying
data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
Crevicesdata<-subset_samples(data_rarified, Sample_Type == "Crevice")

data_rarified=subset_samples(data_rarified, Repeat == "0")


####consider only salamander microbiomes
data_rarified<-subset_samples(data_rarified, Sample_Type == "Salamander")


##### Salamander Mantel Tests#####

lonlat<-data.frame(sample_data(data_rarified)$Long,sample_data(data_rarified)$Lat)
colnames(lonlat)=c('longitude','latitutde')

dist_matrix<-geodist(lonlat,measure="geodesic")
rownames(dist_matrix)<-sample_names(data_rarified)
colnames(dist_matrix)<-sample_names(data_rarified)



#Find  jaccard, bray, UniFrac and weighted-UniFrac indices
jaccard<-vegdist(t(otu_table(data_rarified)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(data_rarified)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(data_rarified,weighted=FALSE)
wunifrac<-UniFrac(data_rarified,weighted=TRUE)

#Perform PcoA
pcoa_jaccard<-pcoa(jaccard)$vectors[,1:2]
pcoa_bray<-pcoa(bray)$vectors[,1:2]
pcoa_unifrac<-pcoa(unifrac)$vectors[,1:2]
pcoa_wunifrac<-pcoa(wunifrac)$vectors[,1:2]



#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)

#The mantel tests
mantel_test_jaccard  <- mantel(jaccard_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_jaccard

mantel_test_bray  <- mantel(bray_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_bray

mantel_test_unifrac  <- mantel(unifrac_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_unifrac

mantel_test_wunifrac  <- mantel(wunifrac_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_wunifrac






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

jac_df<-as.data.frame(mantel_correlog_jaccard$mantel.res)
jac_df$Sig<-c("1","1","1","1","0","0","1","0")

color_palette <- c("white", "black")

######Figure 5A#####
pdf(file = "Fig_5A.pdf", width = 7, height = 7)
ggplot(jac_df, aes(x = class.index, y = Mantel.cor, fill = Sig, group=1)) +
  geom_line()+
  geom_point(size = 7, shape = 21) +  # Point plot with fill based on Petal.Width
  scale_fill_manual(values = color_palette) +  # Manual fill colors
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
        fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1)+theme_classic()+theme(text = element_text(size=25))


dev.off()



mantel_correlog_bray<-mantel.correlog(bray_matrix,dist_matrix,nperm=9999,cutoff=FALSE, break.pts=ddd)

plot(mantel_correlog_bray, alpha=0.05)


mantel_correlog_unifrac<-mantel.correlog(unifrac_matrix,dist_matrix,nperm=9999,cutoff=FALSE,break.pts=ddd)
plot(mantel_correlog_unifrac, alpha=0.05)


mantel_correlog_wunifrac<-mantel.correlog(wunifrac_matrix,dist_matrix,cutoff=FALSE,nperm=9999, break.pts=ddd)

plot(mantel_correlog_wunifrac, alpha=0.05)

pdf(file = "Mantel_Salamander_Geographic.pdf", width = 7, height = 7)
par(mfrow = c(2, 2))
plot(mantel_correlog_jaccard, alpha=0.05); mtext("Jaccard", side = 3, line = 0, cex = 2.2)
plot(mantel_correlog_bray, alpha=0.05); mtext("Bray-Curtis", side = 3, line = 0, cex = 2.2)
plot(mantel_correlog_unifrac, alpha=0.05); mtext("Unifrac", side = 3, line = 0, cex = 2.2)
plot(mantel_correlog_wunifrac, alpha=0.05); mtext("Weighted Unifrac", side = 3, line = 0, cex = 2.2)
dev.off()

#### Salamander Mantel with Genetic Distances ####

genetic_sites<-c('DNR','1251','1250','3688','1477','HW1')
genetics_rarified<-prune_samples(sample_data(data_rarified)$Site %in% genetic_sites,data_rarified)
Fst_dist<-read.csv("Fst_6.csv",row.names=1,header=TRUE)
colnames(Fst_dist)<-c('DNR','1251','1250','3688','HW1','1477')
rownames(Fst_dist)<-c('DNR','1251','1250','3688','HW1','1477')
Msize<-length(sample_names(genetics_rarified))
M<-matrix(ncol=Msize,nrow=Msize)
for (j in 1:Msize){
  for (k in j:Msize){
    
    site1<-sample_data(genetics_rarified)$Site[j]
    site2<-sample_data(genetics_rarified)$Site[k]
    find1<-which(colnames(Fst_dist)==site1)
    find2<-which(rownames(Fst_dist)==site2)
    
    M[j,k]<-Fst_dist[find1,find2]
    M[k,j]<-M[j,k]
    
  }
}
colnames(M)<-c(sample_names(genetics_rarified))
rownames(M)<-c(sample_names(genetics_rarified))

M<-as.data.frame(M)

#Find the euclidean, jaccard, bray, UniFrac and weighted-UniFrac indices
jaccard<-vegdist(t(otu_table(genetics_rarified)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(genetics_rarified)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(genetics_rarified,weighted=FALSE)
wunifrac<-UniFrac(genetics_rarified,weighted=TRUE)

#Perform PcoA
pcoa_jaccard<-pcoa(jaccard)$vectors[,1:2]
pcoa_bray<-pcoa(bray)$vectors[,1:2]
pcoa_unifrac<-pcoa(unifrac)$vectors[,1:2]
pcoa_wunifrac<-pcoa(wunifrac)$vectors[,1:2]



#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)



#write.csv(M, "M.csv")

hello<-read.csv("M.csv")

df <- hello[ -c(1) ]


#I manually changed all the negative values in the df dataframe to 0. Mantel test doesn't like the negative values.
mantel_test_jaccard  <- mantel(jaccard_matrix, df, method = "spearman", permutations = 9999, na.rm = TRUE)

mantel_test_bray  <- mantel(bray_matrix, df, method = "spearman", permutations = 9999, na.rm = TRUE)


mantel_test_unifrac  <- mantel(unifrac_matrix, df, method = "spearman", permutations = 9999, na.rm = TRUE)

mantel_test_wunifrac  <- mantel(wunifrac_matrix, df, method = "spearman", permutations = 9999, na.rm = TRUE)





gen <-c()
gen<-c(gen,0)
gen<-c(gen,.01)
gen<-c(gen,.05)
gen<-c(gen,.1)
gen<-c(gen,.2)


mantel_correlog_jaccard<-mantel.correlog(jaccard_matrix,df,nperm=9999,cutoff=FALSE, break.pts=gen)
plot(mantel_correlog_jaccard, alpha=0.05)


jac_df<-as.data.frame(mantel_correlog_jaccard$mantel.res)
jac_df$Sig<-c("1","1","1","1")

color_palette <- c("black", "white")

######Figure 5B#####
pdf(file = "Fig_5B.pdf", width = 7, height = 7)
ggplot(jac_df, aes(x = class.index, y = Mantel.cor, fill = Sig, group=1)) +
  geom_line()+
  geom_point(size = 7, shape = 21) +  # Point plot with fill based on Petal.Width
  scale_fill_manual(values = color_palette) +  # Manual fill colors
  labs(title = "",
       x = "Distance Class (Fst)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1)+theme_classic()+theme(text = element_text(size=25))


dev.off()

mantel_correlog_bray<-mantel.correlog(bray_matrix,df,nperm=9999,cutoff=FALSE, break.pts=gen)

plot(mantel_correlog_bray, alpha=0.05)


mantel_correlog_unifrac<-mantel.correlog(unifrac_matrix,df,nperm=9999,cutoff=FALSE,break.pts=gen)
plot(mantel_correlog_unifrac, alpha=0.05)

mantel_correlog_wunifrac<-mantel.correlog(wunifrac_matrix,df,cutoff=FALSE,nperm=9999, break.pts=gen)

plot(mantel_correlog_wunifrac, alpha=0.05)

pdf(file = "Mantel_Salamander_Genetic.pdf", width = 7, height = 7)
par(mfrow = c(2, 2))
plot(mantel_correlog_jaccard, alpha=0.05); mtext("Jaccard", side = 3, line = 0, cex = 2.2)
plot(mantel_correlog_bray, alpha=0.05); mtext("Bray-Curtis", side = 3, line = 0, cex = 2.2)
plot(mantel_correlog_unifrac, alpha=0.05); mtext("Unifrac", side = 3, line = 0, cex = 2.2)
plot(mantel_correlog_wunifrac, alpha=0.05); mtext("Weighted Unifrac", side = 3, line = 0, cex = 2.2)
dev.off()

#####Partial Mantel Tests



lonlat<-data.frame(sample_data(genetics_rarified)$Long,sample_data(genetics_rarified)$Lat)
colnames(lonlat)=c('longitude','latitutde')

dist_matrix<-geodist(lonlat,measure="geodesic")
rownames(dist_matrix)<-sample_names(genetics_rarified)
colnames(dist_matrix)<-sample_names(genetics_rarified)


#Find the euclidean, jaccard, bray, UniFrac and weighted-UniFrac indices
jaccard<-vegdist(t(otu_table(genetics_rarified)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(genetics_rarified)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(genetics_rarified,weighted=FALSE)
wunifrac<-UniFrac(genetics_rarified,weighted=TRUE)

#Perform PcoA
pcoa_jaccard<-pcoa(jaccard)$vectors[,1:2]
pcoa_bray<-pcoa(bray)$vectors[,1:2]
pcoa_unifrac<-pcoa(unifrac)$vectors[,1:2]
pcoa_wunifrac<-pcoa(wunifrac)$vectors[,1:2]



#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)


mantel.partial(jaccard_matrix, scale(dist_matrix), scale(df), method="spearman", permutations=9999)
mantel.partial(jaccard_matrix, scale(df), scale(dist_matrix),  method="spearman", permutations=9999)

#####Crevice Mantel Tests#####

Crevicesdata<-subset_samples(Crevicesdata, Lat != "NA")
lonlat<-data.frame(sample_data(Crevicesdata)$Long,sample_data(Crevicesdata)$Lat)
colnames(lonlat)=c('longitude','latitutde')

dist_matrix<-geodist(lonlat,measure="geodesic")
rownames(dist_matrix)<-sample_names(Crevicesdata)
colnames(dist_matrix)<-sample_names(Crevicesdata)



#Find the euclidean, jaccard, bray, UniFrac and weighted-UniFrac indices
jaccard<-vegdist(t(otu_table(Crevicesdata)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(Crevicesdata)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(Crevicesdata,weighted=FALSE)
wunifrac<-UniFrac(Crevicesdata,weighted=TRUE)

#Perform PcoA
pcoa_jaccard<-pcoa(jaccard)$vectors[,1:2]
pcoa_bray<-pcoa(bray)$vectors[,1:2]
pcoa_unifrac<-pcoa(unifrac)$vectors[,1:2]
pcoa_wunifrac<-pcoa(wunifrac)$vectors[,1:2]



#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)

#Running Mantel tests
mantel_test_jaccard  <- mantel(jaccard_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_jaccard

mantel_test_bray  <- mantel(bray_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_bray

mantel_test_unifrac  <- mantel(unifrac_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_unifrac 

mantel_test_wunifrac  <- mantel(wunifrac_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_wunifrac

#Creating Distance Classes
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



mantel_correlog_jaccard<-mantel.correlog(jaccard_matrix,dist_matrix,nperm=9999,cutoff=FALSE, break.pts=ddd)
plot(mantel_correlog_jaccard, alpha=0.05)

jac_df<-as.data.frame(mantel_correlog_jaccard$mantel.res)
jac_df$Sig<-c("0","0","0","0","0","0","0","0")

color_palette <- c("white", "black")

######Figure 5C#####
pdf(file = "Fig_5C.pdf", width = 7, height = 7)
ggplot(jac_df, aes(x = class.index, y = Mantel.cor, fill = Sig, group=1)) +
  geom_line()+
  geom_point(size = 7, shape = 21) +  # Point plot with fill based on Petal.Width
  scale_fill_manual(values = color_palette) +  # Manual fill colors
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1)+theme_classic()+theme(text = element_text(size=25))


dev.off()


mantel_correlog_bray<-mantel.correlog(bray_matrix,dist_matrix,nperm=9999,cutoff=FALSE, break.pts=ddd)

plot(mantel_correlog_bray, alpha=0.05)

mantel_correlog_unifrac<-mantel.correlog(unifrac_matrix,dist_matrix,nperm=9999,cutoff=FALSE,break.pts=ddd)
plot(mantel_correlog_unifrac, alpha=0.05)

mantel_correlog_wunifrac<-mantel.correlog(wunifrac_matrix,dist_matrix,cutoff=FALSE,nperm=9999, break.pts=ddd)

plot(mantel_correlog_wunifrac, alpha=0.05)

pdf(file = "Mantel_Crevice_Geographic.pdf", width = 7, height = 7)
par(mfrow = c(2, 2))
plot(mantel_correlog_jaccard, alpha=0.05); mtext("Jaccard", side = 3, line = 0, cex = 2.2)
plot(mantel_correlog_bray, alpha=0.05); mtext("Bray-Curtis", side = 3, line = 0, cex = 2.2)
plot(mantel_correlog_unifrac, alpha=0.05); mtext("Unifrac", side = 3, line = 0, cex = 2.2)
plot(mantel_correlog_wunifrac, alpha=0.05); mtext("Weighted Unifrac", side = 3, line = 0, cex = 2.2)
dev.off()
plot(mantel_correlog_jaccard, alpha=0.05)






#####Distance Decay ######

#Here I am turning my distance matrix and dissimilairy matrices into vectors 
lonlat<-data.frame(sample_data(data_rarified)$Long,sample_data(data_rarified)$Lat)
colnames(lonlat)=c('longitude','latitutde')

dist_matrix<-geodist(lonlat,measure="geodesic")
rownames(dist_matrix)<-sample_names(data_rarified)
colnames(dist_matrix)<-sample_names(data_rarified)



#Find  jaccard, bray, UniFrac and weighted-UniFrac indices
jaccard<-vegdist(t(otu_table(data_rarified)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(data_rarified)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(data_rarified,weighted=FALSE)
wunifrac<-UniFrac(data_rarified,weighted=TRUE)

#Perform PcoA
pcoa_jaccard<-pcoa(jaccard)$vectors[,1:2]
pcoa_bray<-pcoa(bray)$vectors[,1:2]
pcoa_unifrac<-pcoa(unifrac)$vectors[,1:2]
pcoa_wunifrac<-pcoa(wunifrac)$vectors[,1:2]



#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)


dist<-c()
jac<-c()
dist<-c(dist_matrix)
dist[dist==0] <- NA
jac<-c(jaccard_matrix)
jac[jac==0] <- NA
jac2<-log(1-jac)
jac_df<-data.frame(jac2,dist)
jac_df2<-na.omit(jac_df)
jac_df2 <- jac_df2[apply(jac_df2, 1, function(x) all(is.finite(x))), ]
plot(jac_df2$dist,exp(jac_df2$jac2),xlab="Distance (meters)", ylab="Jaccard Similarity", cex.lab =1.5);thing1<-lm(jac_df2$jac2 ~ jac_df2$dist, data = jac_df2);abline(thing1, col="red")
summary(thing1)


dist<-c()
bray<-c()
dist<-c(dist_matrix)
dist[dist==0] <- NA
bray<-c(bray_matrix)
bray[bray==0] <- NA
bray2<-log(1-bray)
bray_df<-data.frame(bray2,dist)
bray_df2<-na.omit(bray_df)
bray_df2 <- bray_df2[apply(bray_df2, 1, function(x) all(is.finite(x))), ]
plot(bray_df2$dist,exp(bray_df2$bray2),xlab="Distance (meters)", ylab="bray Similarity", cex.lab =1.5);thing1<-lm(bray_df2$bray2 ~ bray_df2$dist, data = bray_df2);abline(thing1, col="red")
summary(thing1)

dist<-c()
unifrac<-c()
dist<-c(dist_matrix)
dist[dist==0] <- NA
unifrac<-c(unifrac_matrix)
unifrac[unifrac==0] <- NA
unifrac2<-log(1-unifrac)
unifrac_df<-data.frame(unifrac2,dist)
unifrac_df2<-na.omit(unifrac_df)
unifrac_df2 <- unifrac_df2[apply(unifrac_df2, 1, function(x) all(is.finite(x))), ]
plot(unifrac_df2$dist,exp(unifrac_df2$unifrac2),xlab="Distance (meters)", ylab="Unifrac Similarity", cex.lab =1.5);thing1<-lm(unifrac_df2$unifrac2 ~ unifrac_df2$dist, data = unifrac_df2);abline(thing1, col="red")
summary(thing1)


dist<-c()
wunifrac<-c()
dist<-c(dist_matrix)
dist[dist==0] <- NA
wunifrac<-c(wunifrac_matrix)
wunifrac[wunifrac==0] <- NA
wunifrac2<-log(1-wunifrac)
wunifrac_df<-data.frame(wunifrac2,dist)
wunifrac_df2<-na.omit(wunifrac_df)
wunifrac_df2 <- wunifrac_df2[apply(wunifrac_df2, 1, function(x) all(is.finite(x))), ]
plot(wunifrac_df2$dist,exp(wunifrac_df2$wunifrac2),xlab="Distance (meters)", ylab="wunifrac Similarity", cex.lab =1.5);thing1<-lm(wunifrac_df2$wunifrac2 ~ wunifrac_df2$dist, data = wunifrac_df2);abline(thing1, col="red")
summary(thing1)




pdf(file = "Distance_Decay.pdf", width = 15, height = 12)
par(mfrow = c(2, 2))

plot(dist,jac2,xlab="Distance (meters)", ylab="Jaccard Similarity", cex.lab =1.5);thing1<-lm((jac2) ~ dist);abline(thing1, col="red")

plot(dist,bray2, xlab="Distance (meters)", ylab="Bray-Curtis Similarity", cex.lab =1.5);thing2<-lm(bray2 ~ dist); abline(thing2, col = "red")

plot(dist,unifrac2, xlab="Distance (meters)", ylab="Unifrac Similarity", cex.lab =1.5);thing3<-lm(unifrac2 ~ dist); abline(thing3, col = "red")

plot(dist,wunifrac2, xlab="Distance (meters)", ylab="Weighted Unifrac Similarity", cex.lab =1.5); thing4<-lm(wunifrac2 ~ dist); abline(thing4, col = "red")

dev.off()



#######Write biom file #######
# As a biom file

library(biomformat);packageVersion("biomformat")
## [1] ‘1.6.0’

otu<-as(otu_table(data_rarified),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"otu_biom.biom")









