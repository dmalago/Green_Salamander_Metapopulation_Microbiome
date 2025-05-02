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

#SO = salamander only, SE = shared with environment, EO = environment only, ES = shared with salamanders, E = all environment, S = all salamander 

group_list<-c('S','SO','SE','E','EO','ES')

lister_curve_jaccard<-list()
lister_points_jaccard<-list()
lister_curve_bray<-list()
lister_points_bray<-list()
lister_curve_unifrac<-list()
lister_points_unifrac<-list()
lister_curve_wunifrac<-list()
lister_points_wunifrac<-list()
glister_curve_jaccard<-list()
glister_points_jaccard<-list()
glister_curve_bray<-list()
glister_points_bray<-list()
glister_curve_unifrac<-list()
glister_points_unifrac<-list()
glister_curve_wunifrac<-list()
glister_points_wunifrac<-list()
lister_x_jaccard<-list()
lister_x_bray<-list()
lister_x_unifrac<-list()
lister_x_wunifrac<-list()
glister_x_jaccard<-list()
glister_x_bray<-list()
glister_x_unifrac<-list()
glister_x_wunifrac<-list()
lister_mantel_jaccard<-list()
lister_mantel_jaccard_sig<-list()
lister_mantel_jaccard_class<-list()
lister_mantel_bray<-list()
lister_mantel_bray_sig<-list()
lister_mantel_bray_class<-list()
lister_mantel_unifrac<-list()
lister_mantel_unifrac_sig<-list()
lister_mantel_unifrac_class<-list()
lister_mantel_wunifrac<-list()
lister_mantel_wunifrac_sig<-list()
lister_mantel_wunifrac_class<-list()
mantel_jaccard_test<-c()
mantel_bray_test<-c()
mantel_unifrac_test<-c()
mantel_wunifrac_test<-c()
glister_mantel_jaccard<-list()
glister_mantel_jaccard_sig<-list()
glister_mantel_jaccard_class<-list()
glister_mantel_bray<-list()
glister_mantel_bray_sig<-list()
glister_mantel_bray_class<-list()
glister_mantel_unifrac<-list()
glister_mantel_unifrac_sig<-list()
glister_mantel_unifrac_class<-list()
glister_mantel_wunifrac<-list()
glister_mantel_wunifrac_sig<-list()
glister_mantel_wunifrac_class<-list()
gmantel_jaccard_test<-c()
gmantel_bray_test<-c()
gmantel_unifrac_test<-c()
gmantel_wunifrac_test<-c()

for (list_k in 1:length(group_list)){
  
which_group=group_list[list_k]

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

Cphyloseq<-subset_samples(data, Sample_Type == "Crevice")
Sphyloseq<-subset_samples(data, Sample_Type == "Salamander")

if (which_group == 'SO'){
  
  data_rarified<-prune_taxa(taxa_sums(Cphyloseq)==0,data)
  
  #Rarefy your data, here I am rarefying to 5000 sequences per sample
  data_rarified = rarefy_even_depth(data_rarified, rngseed=1, sample.size=5000, replace=F, verbose = TRUE,trimOTUs = 0)
  
  #Removing any taxa which aren't present following rarefying
  data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
  #Remove repeat sampled animals
  data_rarified=subset_samples(data_rarified, Repeat == "0")
  ####consider only salamander microbiomes
  data_rarified<-subset_samples(data_rarified, Sample_Type == "Salamander")
  #Removing any taxa which aren't present following rarefying
  data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
  
} else if (which_group == 'SE'){
  
  data_rarified<-prune_taxa(taxa_sums(Cphyloseq)>0,data)
  
  #Rarefy your data, here I am rarefying to 5000 sequences per sample
  data_rarified = rarefy_even_depth(data_rarified, rngseed=1, sample.size=5000, replace=F, verbose = TRUE,trimOTUs = 0)
  
  #Removing any taxa which aren't present following rarefying
  data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
  #Remove repeat sampled animals
  data_rarified=subset_samples(data_rarified, Repeat == "0")
  ####consider only salamander microbiomes
  data_rarified<-subset_samples(data_rarified, Sample_Type == "Salamander")
  #Removing any taxa which aren't present following rarefying
  data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
  
} else if (which_group == 'E'){
  
  data_rarified<-Cphyloseq
  #Rarefy your data, here I am rarefying to 5000 sequences per sample
  data_rarified = rarefy_even_depth(data_rarified, rngseed=1, sample.size=5000, replace=F, verbose = TRUE,trimOTUs = 0)
  
  #Removing any taxa which aren't present following rarefying
  data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
  
} else if (which_group=='EO'){

  data_rarified<-prune_taxa(taxa_sums(Sphyloseq)==0,data)
  #Rarefy your data, here I am rarefying to 5000 sequences per sample
  data_rarified = rarefy_even_depth(data_rarified, rngseed=1, sample.size=5000, replace=F, verbose = TRUE,trimOTUs = 0)
  
  #Removing any taxa which aren't present following rarefying
  data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
  
  ####consider only crevice microbiomes
  data_rarified<-subset_samples(data_rarified, Sample_Type == "Crevice")
  #Removing any taxa which aren't present following rarefying
  data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
  
  } else if (which_group == 'ES'){
    
    data_rarified<-prune_taxa(taxa_sums(Sphyloseq)>0,data)
    #Rarefy your data, here I am rarefying to 5000 sequences per sample
    data_rarified = rarefy_even_depth(data_rarified, rngseed=1, sample.size=5000, replace=F, verbose = TRUE,trimOTUs = 0)
    
    #Removing any taxa which aren't present following rarefying
    data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
    
    ####consider only crevice microbiomes
    data_rarified<-subset_samples(data_rarified, Sample_Type == "Crevice")
    #Removing any taxa which aren't present following rarefying
    data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
    
    }else{
  
  data_rarified<-data
  #Rarefy your data, here I am rarefying to 5000 sequences per sample
  data_rarified = rarefy_even_depth(data_rarified, rngseed=1, sample.size=5000, replace=F, verbose = TRUE,trimOTUs = 0)
  
  #Removing any taxa which aren't present following rarefying
  data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
  #Remove repeat sampled animals
  data_rarified=subset_samples(data_rarified, Repeat == "0")
  ####consider only salamander microbiomes
  data_rarified<-subset_samples(data_rarified, Sample_Type == "Salamander")
  #Removing any taxa which aren't present following rarefying
  data_rarified<- prune_taxa(taxa_sums(data_rarified) > 0, data_rarified) 
  
}


##### Salamander Mantel Tests#####

lonlat<-data.frame(sample_data(data_rarified)$Long,sample_data(data_rarified)$Lat)
colnames(lonlat)=c('longitude','latitutde')

dist_matrix<-geodist(lonlat,measure="geodesic")
rownames(dist_matrix)<-sample_names(data_rarified)
colnames(dist_matrix)<-sample_names(data_rarified)

df<-data.frame(sample_data(data_rarified)$site,sample_data(data_rarified)$Year)
colnames(df)<-c('site','year')
yearxsite<-table(df)

#Find  jaccard, bray, UniFrac and weighted-UniFrac indices
jaccard<-vegdist(t(otu_table(data_rarified)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(data_rarified)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(data_rarified,weighted=FALSE)
wunifrac<-UniFrac(data_rarified,weighted=TRUE)

#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)




##### Salamander Distance Decay#####

dist_matrixU <- dist_matrix

#Find the vector from the upper triangle
dist_matrixU<-dist_matrixU[upper.tri(dist_matrixU)] 

##### Jaccard#####

#Convert to similarity
jaccard_matrixU <- 1-jaccard_matrix

#Find the vector from the upper triangle
jaccard_matrixU<-jaccard_matrixU[upper.tri(jaccard_matrixU)] 

#Remove any identical distances (should be very few, they blow up in log transform)
jdist<-dist_matrixU[jaccard_matrixU!=0]
lister_x_jaccard<-list.append(lister_x_jaccard,jdist)

jmet<-jaccard_matrixU[jaccard_matrixU!=0]

#Perform regression
jfit<-glm(log(jmet)~jdist)
summary(jfit)
print(paste('intercept:',round(exp(jfit$coefficients[1]),3)))
print(paste('slope:',jfit$coefficients[2]))

#Plot regression
x<-1:1:50000
y<-jfit$coefficients[[1]]+x*jfit$coefficients[[2]]
plot(jdist,jmet,cex.lab=1.5,xlab='Distance',ylab='Similarity',cex.axis=1.2,pch=16)
lines(x,exp(y),col='red')

lister_points_jaccard<-list.append(lister_points_jaccard,jmet)
lister_curve_jaccard<-list.append(lister_curve_jaccard,y)

##### Bray#####

#Convert to similarity
bray_matrixU <- 1-bray_matrix

#Find the vector from the upper triangle
bray_matrixU<-bray_matrixU[upper.tri(bray_matrixU)] 

#Remove any identical distances (should be very few, they blow up in log transform)
jdist<-dist_matrixU[bray_matrixU!=0]
lister_x_bray<-list.append(lister_x_bray,jdist)

jmet<-bray_matrixU[bray_matrixU!=0]

#Perform regression
jfit<-glm(log(jmet)~jdist)
summary(jfit)
print(paste('intercept:',round(exp(jfit$coefficients[1]),3)))
print(paste('slope:',jfit$coefficients[2]))

#Plot regression
x<-1:1:50000
y<-jfit$coefficients[[1]]+x*jfit$coefficients[[2]]
plot(jdist,jmet,cex.lab=1.5,xlab='Distance',ylab='Similarity',cex.axis=1.2,pch=16)
lines(x,exp(y),col='red')

lister_points_bray<-list.append(lister_points_bray,jmet)
lister_curve_bray<-list.append(lister_curve_bray,y)


##### UniFrac#####

#Convert to similarity
unifrac_matrixU <- 1-unifrac_matrix

#Find the vector from the upper triangle
unifrac_matrixU<-unifrac_matrixU[upper.tri(unifrac_matrixU)] 

#Remove any identical distances (should be very few, they blow up in log transform)
jdist<-dist_matrixU[unifrac_matrixU!=0]
lister_x_unifrac<-list.append(lister_x_unifrac,jdist)

jmet<-unifrac_matrixU[unifrac_matrixU!=0]

#Perform regression
jfit<-glm(log(jmet)~jdist)
summary(jfit)
print(paste('intercept:',round(exp(jfit$coefficients[1]),3)))
print(paste('slope:',jfit$coefficients[2]))

#Plot regression
x<-1:1:50000
y<-jfit$coefficients[[1]]+x*jfit$coefficients[[2]]
plot(jdist,jmet,cex.lab=1.5,xlab='Distance',ylab='Similarity',cex.axis=1.2,pch=16)
lines(x,exp(y),col='red')


lister_points_unifrac<-list.append(lister_points_unifrac,jmet)
lister_curve_unifrac<-list.append(lister_curve_unifrac,y)

##### weighted UniFrac#####

#Convert to similarity
wunifrac_matrixU <- 1-wunifrac_matrix

#Find the vector from the upper triangle
wunifrac_matrixU<-wunifrac_matrixU[upper.tri(wunifrac_matrixU)] 

#Remove any identical distances (should be very few, they blow up in log transform)
jdist<-dist_matrixU[wunifrac_matrixU!=0]
lister_x_wunifrac<-list.append(lister_x_wunifrac,jdist)

jmet<-wunifrac_matrixU[wunifrac_matrixU!=0]

#Perform regression
jfit<-glm(log(jmet)~jdist)
summary(jfit)
print(paste('intercept:',round(exp(jfit$coefficients[1]),3)))
print(paste('slope:',jfit$coefficients[2]))

#Plot regression
x<-1:1:50000
y<-jfit$coefficients[[1]]+x*jfit$coefficients[[2]]
plot(jdist,jmet,cex.lab=1.5,xlab='Distance',ylab='Similarity',cex.axis=1.2,pch=16)
lines(x,exp(y),col='red')

lister_points_wunifrac<-list.append(lister_points_wunifrac,jmet)
lister_curve_wunifrac<-list.append(lister_curve_wunifrac,y)




#The mantel tests
mantel_test_jaccard  <- mantel(jaccard_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_jaccard
mantel_jaccard_test<-c(mantel_jaccard_test,mantel_test_jaccard$signif)

mantel_test_bray  <- mantel(bray_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_bray
mantel_bray_test<-c(mantel_bray_test,mantel_test_bray$signif)

mantel_test_unifrac  <- mantel(unifrac_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_unifrac
mantel_unifrac_test<-c(mantel_unifrac_test,mantel_test_unifrac$signif)

mantel_test_wunifrac  <- mantel(wunifrac_matrix, dist_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_test_wunifrac
mantel_wunifrac_test<-c(mantel_wunifrac_test,mantel_test_wunifrac$signif)


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

lister_mantel_jaccard<-list.append(lister_mantel_jaccard,jac_df$Mantel.cor)
lister_mantel_jaccard_class<-list.append(lister_mantel_jaccard_class,jac_df$class.index)
lister_mantel_jaccard_sig<-list.append(lister_mantel_jaccard_sig,sign(-sign(mantel_correlog_jaccard$mantel.res[,5]-0.0499999)+1))

color_palette <- c("white", "black")

######Figure 5A#####
#pdf(file = "Fig_5A.pdf", width = 7, height = 7)
ggplot(jac_df, aes(x = class.index, y = Mantel.cor, fill = Sig, group=1)) +
  geom_line()+
  geom_point(size = 7, shape = 21) +  # Point plot with fill based on Petal.Width
  scale_fill_manual(values = color_palette) +  # Manual fill colors
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
        fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10))


#dev.off()



mantel_correlog_bray<-mantel.correlog(bray_matrix,dist_matrix,nperm=9999,cutoff=FALSE, break.pts=ddd)
bray_df<-as.data.frame(mantel_correlog_jaccard$mantel.res)

lister_mantel_bray<-list.append(lister_mantel_bray,bray_df$Mantel.cor)
lister_mantel_bray_class<-list.append(lister_mantel_bray_class,bray_df$class.index)
lister_mantel_bray_sig<-list.append(lister_mantel_bray_sig,sign(-sign(mantel_correlog_bray$mantel.res[,5]-0.0499999)+1))


mantel_correlog_unifrac<-mantel.correlog(unifrac_matrix,dist_matrix,nperm=9999,cutoff=FALSE,break.pts=ddd)
plot(mantel_correlog_unifrac, alpha=0.05)
unifrac_df<-as.data.frame(mantel_correlog_unifrac$mantel.res)

lister_mantel_unifrac<-list.append(lister_mantel_unifrac,unifrac_df$Mantel.cor)
lister_mantel_unifrac_class<-list.append(lister_mantel_unifrac_class,unifrac_df$class.index)
lister_mantel_unifrac_sig<-list.append(lister_mantel_unifrac_sig,sign(-sign(mantel_correlog_unifrac$mantel.res[,5]-0.0499999)+1))


mantel_correlog_wunifrac<-mantel.correlog(wunifrac_matrix,dist_matrix,cutoff=FALSE,nperm=9999, break.pts=ddd)
wunifrac_df<-as.data.frame(mantel_correlog_unifrac$mantel.res)

lister_mantel_wunifrac<-list.append(lister_mantel_wunifrac,wunifrac_df$Mantel.cor)
lister_mantel_wunifrac_class<-list.append(lister_mantel_wunifrac_class,wunifrac_df$class.index)
lister_mantel_wunifrac_sig<-list.append(lister_mantel_wunifrac_sig,sign(-sign(mantel_correlog_wunifrac$mantel.res[,5]-0.0499999)+1))


if (list_k<4){

#### Salamander Mantel with Genetic Distances ####

genetic_sites<-c('DNR','1251','1250','3688','1477','HW1')
genetics_rarified<-prune_samples(sample_data(data_rarified)$Site %in% genetic_sites,data_rarified)
Fst_dist<-read.csv("Fst_7.csv",row.names=1,header=TRUE)

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
    
    M[j,k]<-as.numeric(Fst_dist[find1,find2])
    M[k,j]<-M[j,k]
    
  }
}
colnames(M)<-c(sample_names(genetics_rarified))
rownames(M)<-c(sample_names(genetics_rarified))

#M<-as.data.frame(M)

#Find the euclidean, jaccard, bray, UniFrac and weighted-UniFrac indices
jaccard<-vegdist(t(otu_table(genetics_rarified)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(genetics_rarified)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(genetics_rarified,weighted=FALSE)
wunifrac<-UniFrac(genetics_rarified,weighted=TRUE)

#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)



#write.csv(M, "M.csv")

#hello<-read.csv("M.csv",row.names = 1)

df <- M#hello[ -c(1) ]
df[df<0]<-0

dist_matrixU <- as.matrix(df)
dist_matrixU<-dist_matrixU[upper.tri(dist_matrixU)] 
#Convert to similarity
jaccard_matrixU <- 1-jaccard_matrix

#Find the vector from the upper triangle
jaccard_matrixU<-jaccard_matrixU[upper.tri(jaccard_matrixU)] 

#Remove any identical distances (should be very few, they blow up in log transform)
jdist<-dist_matrixU[jaccard_matrixU!=0]
glister_x_jaccard<-list.append(glister_x_jaccard,jdist)

jmet<-jaccard_matrixU[jaccard_matrixU!=0]

#Perform regression
jfit<-glm(log(jmet)~jdist)
summary(jfit)
print(paste('intercept:',round(exp(jfit$coefficients[1]),3)))
print(paste('slope:',jfit$coefficients[2]))

#Plot regression
x<-c(0,0.001*(1:1:170))
y<-jfit$coefficients[[1]]+x*jfit$coefficients[[2]]
plot(jdist,jmet,cex.lab=1.5,xlab='Distance',ylab='Similarity',cex.axis=1.2,pch=16)
lines(x,exp(y),col='red')

glister_points_jaccard<-list.append(glister_points_jaccard,jmet)
glister_curve_jaccard<-list.append(glister_curve_jaccard,y)

##### Bray#####

#Convert to similarity
bray_matrixU <- 1-bray_matrix

#Find the vector from the upper triangle
bray_matrixU<-bray_matrixU[upper.tri(bray_matrixU)] 

#Remove any identical distances (should be very few, they blow up in log transform)
jdist<-dist_matrixU[bray_matrixU!=0]
glister_x_bray<-list.append(glister_x_bray,jdist)

jmet<-bray_matrixU[bray_matrixU!=0]

#Perform regression
jfit<-glm(log(jmet)~jdist)
summary(jfit)
print(paste('intercept:',round(exp(jfit$coefficients[1]),3)))
print(paste('slope:',jfit$coefficients[2]))

#Plot regression
x<-c(0,0.001*(1:1:170))
y<-jfit$coefficients[[1]]+x*jfit$coefficients[[2]]
plot(jdist,jmet,cex.lab=1.5,xlab='Distance',ylab='Similarity',cex.axis=1.2,pch=16)
lines(x,exp(y),col='red')

glister_points_bray<-list.append(glister_points_bray,jmet)
glister_curve_bray<-list.append(glister_curve_bray,y)


##### UniFrac#####

#Convert to similarity
unifrac_matrixU <- 1-unifrac_matrix

#Find the vector from the upper triangle
unifrac_matrixU<-unifrac_matrixU[upper.tri(unifrac_matrixU)] 

#Remove any identical distances (should be very few, they blow up in log transform)
jdist<-dist_matrixU[unifrac_matrixU!=0]
glister_x_unifrac<-list.append(glister_x_unifrac,jdist)

jmet<-unifrac_matrixU[unifrac_matrixU!=0]

#Perform regression
jfit<-glm(log(jmet)~jdist)
summary(jfit)
print(paste('intercept:',round(exp(jfit$coefficients[1]),3)))
print(paste('slope:',jfit$coefficients[2]))

#Plot regression
x<-c(0,0.001*(1:1:170))
y<-jfit$coefficients[[1]]+x*jfit$coefficients[[2]]
plot(jdist,jmet,cex.lab=1.5,xlab='Distance',ylab='Similarity',cex.axis=1.2,pch=16)
lines(x,exp(y),col='red')


glister_points_unifrac<-list.append(glister_points_unifrac,jmet)
glister_curve_unifrac<-list.append(glister_curve_unifrac,y)

##### weighted UniFrac#####

#Convert to similarity
wunifrac_matrixU <- 1-wunifrac_matrix

#Find the vector from the upper triangle
wunifrac_matrixU<-wunifrac_matrixU[upper.tri(wunifrac_matrixU)] 

#Remove any identical distances (should be very few, they blow up in log transform)
jdist<-dist_matrixU[wunifrac_matrixU!=0]
glister_x_wunifrac<-list.append(glister_x_wunifrac,jdist)

jmet<-wunifrac_matrixU[wunifrac_matrixU!=0]

#Perform regression
jfit<-glm(log(jmet)~jdist)
summary(jfit)
print(paste('intercept:',round(exp(jfit$coefficients[1]),3)))
print(paste('slope:',jfit$coefficients[2]))

#Plot regression
x<-c(0,0.001*(1:1:170))
y<-jfit$coefficients[[1]]+x*jfit$coefficients[[2]]
plot(jdist,jmet,cex.lab=1.5,xlab='Distance',ylab='Similarity',cex.axis=1.2,pch=16)
lines(x,exp(y),col='red')

glister_points_wunifrac<-list.append(glister_points_wunifrac,jmet)
glister_curve_wunifrac<-list.append(glister_curve_wunifrac,y)



#I manually changed all the negative values in the df dataframe to 0. Mantel test doesn't like the negative values.
mantel_test_jaccard  <- mantel(jaccard_matrix, df, method = "spearman", permutations = 9999, na.rm = TRUE)
gmantel_jaccard_test<-c(gmantel_jaccard_test,mantel_test_jaccard$signif)

mantel_test_bray  <- mantel(bray_matrix, df, method = "spearman", permutations = 9999, na.rm = TRUE)
gmantel_bray_test<-c(gmantel_bray_test,mantel_test_bray$signif)

mantel_test_unifrac  <- mantel(unifrac_matrix, df, method = "spearman", permutations = 9999, na.rm = TRUE)
gmantel_unifrac_test<-c(gmantel_unifrac_test,mantel_test_unifrac$signif)

mantel_test_wunifrac  <- mantel(wunifrac_matrix, df, method = "spearman", permutations = 9999, na.rm = TRUE)
gmantel_wunifrac_test<-c(gmantel_wunifrac_test,mantel_test_wunifrac$signif)





gen <-c()
gen<-c(gen,0)
gen<-c(gen,.01)
gen<-c(gen,.05)
gen<-c(gen,.1)
gen<-c(gen,.2)


mantel_correlog_jaccard<-mantel.correlog(jaccard_matrix,df,nperm=9999,cutoff=FALSE, break.pts=gen)
jac_df<-as.data.frame(mantel_correlog_jaccard$mantel.res)
glister_mantel_jaccard<-list.append(glister_mantel_jaccard,jac_df$Mantel.cor)
glister_mantel_jaccard_class<-list.append(glister_mantel_jaccard_class,jac_df$class.index)
glister_mantel_jaccard_sig<-list.append(glister_mantel_jaccard_sig,sign(-sign(mantel_correlog_jaccard$mantel.res[,5]-0.0499999)+1))

mantel_correlog_bray<-mantel.correlog(bray_matrix,df,nperm=9999,cutoff=FALSE, break.pts=gen)
bray_df<-as.data.frame(mantel_correlog_bray$mantel.res)
glister_mantel_bray<-list.append(glister_mantel_bray,bray_df$Mantel.cor)
glister_mantel_bray_class<-list.append(glister_mantel_bray_class,bray_df$class.index)
glister_mantel_bray_sig<-list.append(glister_mantel_bray_sig,sign(-sign(mantel_correlog_bray$mantel.res[,5]-0.0499999)+1))

mantel_correlog_unifrac<-mantel.correlog(unifrac_matrix,df,nperm=9999,cutoff=FALSE, break.pts=gen)
unifrac_df<-as.data.frame(mantel_correlog_unifrac$mantel.res)
glister_mantel_unifrac<-list.append(glister_mantel_unifrac,unifrac_df$Mantel.cor)
glister_mantel_unifrac_class<-list.append(glister_mantel_unifrac_class,unifrac_df$class.index)
glister_mantel_unifrac_sig<-list.append(glister_mantel_unifrac_sig,sign(-sign(mantel_correlog_unifrac$mantel.res[,5]-0.0499999)+1))

mantel_correlog_wunifrac<-mantel.correlog(wunifrac_matrix,df,nperm=9999,cutoff=FALSE, break.pts=gen)
wunifrac_df<-as.data.frame(mantel_correlog_wunifrac$mantel.res)
glister_mantel_wunifrac<-list.append(glister_mantel_wunifrac,wunifrac_df$Mantel.cor)
glister_mantel_wunifrac_class<-list.append(glister_mantel_wunifrac_class,wunifrac_df$class.index)
glister_mantel_wunifrac_sig<-list.append(glister_mantel_wunifrac_sig,sign(-sign(mantel_correlog_wunifrac$mantel.res[,5]-0.0499999)+1))

}
}


xv<-1:1:50000
a1<-ggplot()+
  geom_point(aes(x=lister_x_jaccard[[1]],y=(lister_points_jaccard[[1]])),size=2,fill='green',col='green',alpha=0.2)+xlab('Distance (m)')+ylab('Similarity')+ylim(0,0.3)+
  geom_point(aes(x=lister_x_jaccard[[2]],y=(lister_points_jaccard[[2]])),size=2,fill='greenyellow',col='greenyellow',alpha=0.2)+
  geom_point(aes(x=lister_x_jaccard[[3]],y=(lister_points_jaccard[[3]])),size=2,fill='darkolivegreen4',col='darkolivegreen4',alpha=0.2)+
  geom_line(aes(y=exp(lister_curve_jaccard[[1]]),x=xv),linewidth=0.5,col='green')+
  geom_line(aes(y=exp(lister_curve_jaccard[[2]]),x=xv),linewidth=0.5,col='greenyellow')+
  geom_line(aes(y=exp(lister_curve_jaccard[[3]]),x=xv),linewidth=0.5,col='darkolivegreen4')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

a2<-ggplot()+
  geom_point(aes(x=lister_x_jaccard[[4]],y=(lister_points_jaccard[[4]])),size=2,fill='navajowhite3',col='navajowhite3',alpha=0.2)+xlab('Distance (m)')+ylab('Similarity')+ylim(0,0.3)+
  geom_point(aes(x=lister_x_jaccard[[5]],y=(lister_points_jaccard[[5]])),size=2,fill='grey40',col='grey40',alpha=0.2)+
  geom_point(aes(x=lister_x_jaccard[[6]],y=(lister_points_jaccard[[6]])),size=2,fill='darkseagreen3',col='darkseagreen3',alpha=0.2)+
  geom_line(aes(y=exp(lister_curve_jaccard[[4]]),x=xv),linewidth=0.5,col='navajowhite3')+
  geom_line(aes(y=exp(lister_curve_jaccard[[5]]),x=xv),linewidth=0.5,col='grey40')+
  geom_line(aes(y=exp(lister_curve_jaccard[[6]]),x=xv),linewidth=0.5,col='darkseagreen3')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))



b1<-ggplot()+ylim(0,0.7)+
  geom_point(aes(x=lister_x_bray[[1]],y=(lister_points_bray[[1]])),size=2,fill='green',col='green',alpha=0.2)+xlab('Distance (m)')+ylab('Similarity')+
  geom_point(aes(x=lister_x_bray[[2]],y=(lister_points_bray[[2]])),size=2,fill='greenyellow',col='greenyellow',alpha=0.2)+
  geom_point(aes(x=lister_x_bray[[3]],y=(lister_points_bray[[3]])),size=2,fill='darkolivegreen4',col='darkolivegreen4',alpha=0.2)+
  geom_line(aes(y=exp(lister_curve_bray[[1]]),x=xv),linewidth=0.5,col='green')+
  geom_line(aes(y=exp(lister_curve_bray[[2]]),x=xv),linewidth=0.5,col='greenyellow')+
  geom_line(aes(y=exp(lister_curve_bray[[3]]),x=xv),linewidth=0.5,col='darkolivegreen4')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

b2<-ggplot()+ylim(0,0.7)+
  geom_point(aes(x=lister_x_bray[[4]],y=(lister_points_bray[[4]])),size=2,fill='navajowhite3',col='navajowhite3',alpha=0.2)+xlab('Distance (m)')+ylab('Similarity')+
  geom_point(aes(x=lister_x_bray[[5]],y=(lister_points_bray[[5]])),size=2,fill='grey40',col='grey40',alpha=0.2)+
  geom_point(aes(x=lister_x_bray[[6]],y=(lister_points_bray[[6]])),size=2,fill='darkseagreen3',col='darkseagreen3',alpha=0.2)+
  geom_line(aes(y=exp(lister_curve_bray[[4]]),x=xv),linewidth=0.5,col='navajowhite3')+
  geom_line(aes(y=exp(lister_curve_bray[[5]]),x=xv),linewidth=0.5,col='grey40')+
  geom_line(aes(y=exp(lister_curve_bray[[6]]),x=xv),linewidth=0.5,col='darkseagreen3')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


c1<-ggplot()+
  geom_point(aes(x=lister_x_unifrac[[1]],y=(lister_points_unifrac[[1]])),size=2,fill='green',col='green',alpha=0.2)+xlab('Distance (m)')+ylab('Similarity')+
  geom_point(aes(x=lister_x_unifrac[[2]],y=(lister_points_unifrac[[2]])),size=2,fill='greenyellow',col='greenyellow',alpha=0.2)+
  geom_point(aes(x=lister_x_unifrac[[3]],y=(lister_points_unifrac[[3]])),size=2,fill='darkolivegreen4',col='darkolivegreen4',alpha=0.2)+
  geom_line(aes(y=exp(lister_curve_unifrac[[1]]),x=xv),linewidth=0.5,col='green')+
  geom_line(aes(y=exp(lister_curve_unifrac[[2]]),x=xv),linewidth=0.5,col='greenyellow')+
  geom_line(aes(y=exp(lister_curve_unifrac[[3]]),x=xv),linewidth=0.5,col='darkolivegreen4')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

c2<-ggplot()+
  geom_point(aes(x=lister_x_unifrac[[4]],y=(lister_points_unifrac[[4]])),size=2,fill='navajowhite3',col='navajowhite3',alpha=0.2)+xlab('Distance (m)')+ylab('Similarity')+
  geom_point(aes(x=lister_x_unifrac[[5]],y=(lister_points_unifrac[[5]])),size=2,fill='grey40',col='grey40',alpha=0.2)+
  geom_point(aes(x=lister_x_unifrac[[6]],y=(lister_points_unifrac[[6]])),size=2,fill='darkseagreen3',col='darkseagreen3',alpha=0.2)+
  geom_line(aes(y=exp(lister_curve_unifrac[[4]]),x=xv),linewidth=0.5,col='navajowhite3')+
  geom_line(aes(y=exp(lister_curve_unifrac[[5]]),x=xv),linewidth=0.5,col='grey40')+
  geom_line(aes(y=exp(lister_curve_unifrac[[6]]),x=xv),linewidth=0.5,col='darkseagreen3')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


d1<-ggplot()+
  geom_point(aes(x=lister_x_wunifrac[[1]],y=(lister_points_wunifrac[[1]])),size=2,fill='green',col='green',alpha=0.2)+xlab('Distance (m)')+ylab('Similarity')+
  geom_point(aes(x=lister_x_wunifrac[[2]],y=(lister_points_wunifrac[[2]])),size=2,fill='greenyellow',col='greenyellow',alpha=0.2)+
  geom_point(aes(x=lister_x_wunifrac[[3]],y=(lister_points_wunifrac[[3]])),size=2,fill='darkolivegreen4',col='darkolivegreen4',alpha=0.2)+
  geom_line(aes(y=exp(lister_curve_wunifrac[[1]]),x=xv),linewidth=0.5,col='green')+
  geom_line(aes(y=exp(lister_curve_wunifrac[[2]]),x=xv),linewidth=0.5,col='greenyellow')+
  geom_line(aes(y=exp(lister_curve_wunifrac[[3]]),x=xv),linewidth=0.5,col='darkolivegreen4')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

d2<-ggplot()+
  geom_point(aes(x=lister_x_wunifrac[[4]],y=(lister_points_wunifrac[[4]])),size=2,fill='navajowhite3',col='navajowhite3',alpha=0.2)+xlab('Distance (m)')+ylab('Similarity')+
  geom_point(aes(x=lister_x_wunifrac[[5]],y=(lister_points_wunifrac[[5]])),size=2,fill='grey40',col='grey40',alpha=0.2)+
  geom_point(aes(x=lister_x_wunifrac[[6]],y=(lister_points_wunifrac[[6]])),size=2,fill='darkseagreen3',col='darkseagreen3',alpha=0.2)+
  geom_line(aes(y=exp(lister_curve_wunifrac[[4]]),x=xv),linewidth=0.5,col='navajowhite3')+
  geom_line(aes(y=exp(lister_curve_wunifrac[[5]]),x=xv),linewidth=0.5,col='grey40')+
  geom_line(aes(y=exp(lister_curve_wunifrac[[6]]),x=xv),linewidth=0.5,col='darkseagreen3')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


newcol<-c('green','greenyellow','darkolivegreen4','navajowhite3','grey40','darkseagreen3')


onezerolist_jaccard<-list()
for (k in 1:length(group_list)){
  temp<-as.vector(lister_mantel_jaccard_sig[[k]])
  temp[temp==1]<-newcol[k]
  temp[temp==0]<-'white'
  onezerolist_jaccard<-list.append(onezerolist_jaccard,temp)
}


a3<-ggplot()+ylim(-0.4,0.4)+
  geom_point(aes(x = lister_mantel_jaccard_class[[1]], y = lister_mantel_jaccard[[1]]),fill=onezerolist_jaccard[[1]],color='green',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_jaccard_class[[1]], y = lister_mantel_jaccard[[1]]),color='green')+
  geom_point(aes(x = lister_mantel_jaccard_class[[2]], y = lister_mantel_jaccard[[2]]),fill=onezerolist_jaccard[[2]],color='greenyellow',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_jaccard_class[[2]], y = lister_mantel_jaccard[[2]]),color='greenyellow')+
  geom_point(aes(x = lister_mantel_jaccard_class[[3]], y = lister_mantel_jaccard[[3]]),fill=onezerolist_jaccard[[3]],color='darkolivegreen4',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_jaccard_class[[3]], y = lister_mantel_jaccard[[3]]),color='darkolivegreen4')+
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


a4<-ggplot()+ylim(-0.4,0.4)+
  geom_point(aes(x = lister_mantel_jaccard_class[[4]], y = lister_mantel_jaccard[[4]]),fill=onezerolist_jaccard[[4]],color='navajowhite3',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_jaccard_class[[4]], y = lister_mantel_jaccard[[4]]),color='navajowhite3')+
  geom_point(aes(x = lister_mantel_jaccard_class[[5]], y = lister_mantel_jaccard[[5]]),fill=onezerolist_jaccard[[5]],color='grey40',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_jaccard_class[[5]], y = lister_mantel_jaccard[[5]]),color='grey40')+
  geom_point(aes(x = lister_mantel_jaccard_class[[6]], y = lister_mantel_jaccard[[6]]),fill=onezerolist_jaccard[[6]],color='darkseagreen3',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_jaccard_class[[6]], y = lister_mantel_jaccard[[6]]),color='darkseagreen3')+
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


onezerolist_bray<-list()
for (k in 1:length(group_list)){
  temp<-as.vector(lister_mantel_bray_sig[[k]])
  temp[temp==1]<-newcol[k]
  temp[temp==0]<-'white'
  onezerolist_bray<-list.append(onezerolist_bray,temp)
}


b3<-ggplot()+ylim(-0.4,0.3)+
  geom_point(aes(x = lister_mantel_bray_class[[1]], y = lister_mantel_bray[[1]]),fill=onezerolist_bray[[1]],color='green',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_bray_class[[1]], y = lister_mantel_bray[[1]]),color='green')+
  geom_point(aes(x = lister_mantel_bray_class[[2]], y = lister_mantel_bray[[2]]),fill=onezerolist_bray[[2]],color='greenyellow',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_bray_class[[2]], y = lister_mantel_bray[[2]]),color='greenyellow')+
  geom_point(aes(x = lister_mantel_bray_class[[3]], y = lister_mantel_bray[[3]]),fill=onezerolist_bray[[3]],color='darkolivegreen4',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_bray_class[[3]], y = lister_mantel_bray[[3]]),color='darkolivegreen4')+
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


b4<-ggplot()+ylim(-0.4,0.3)+
  geom_point(aes(x = lister_mantel_bray_class[[4]], y = lister_mantel_bray[[4]]),fill=onezerolist_bray[[4]],color='navajowhite3',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_bray_class[[4]], y = lister_mantel_bray[[4]]),color='navajowhite3')+
  geom_point(aes(x = lister_mantel_bray_class[[5]], y = lister_mantel_bray[[5]]),fill=onezerolist_bray[[5]],color='grey40',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_bray_class[[5]], y = lister_mantel_bray[[5]]),color='grey40')+
  geom_point(aes(x = lister_mantel_bray_class[[6]], y = lister_mantel_bray[[6]]),fill=onezerolist_bray[[6]],color='darkseagreen3',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_bray_class[[6]], y = lister_mantel_bray[[6]]),color='darkseagreen3')+
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


onezerolist_unifrac<-list()
for (k in 1:length(group_list)){
  temp<-as.vector(lister_mantel_unifrac_sig[[k]])
  temp[temp==1]<-newcol[k]
  temp[temp==0]<-'white'
  onezerolist_unifrac<-list.append(onezerolist_unifrac,temp)
}


c3<-ggplot()+ylim(-0.4,0.4)+
  geom_point(aes(x = lister_mantel_unifrac_class[[1]], y = lister_mantel_unifrac[[1]]),fill=onezerolist_unifrac[[1]],color='green',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_unifrac_class[[1]], y = lister_mantel_unifrac[[1]]),color='green')+
  geom_point(aes(x = lister_mantel_unifrac_class[[2]], y = lister_mantel_unifrac[[2]]),fill=onezerolist_unifrac[[2]],color='greenyellow',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_unifrac_class[[2]], y = lister_mantel_unifrac[[2]]),color='greenyellow')+
  geom_point(aes(x = lister_mantel_unifrac_class[[3]], y = lister_mantel_unifrac[[3]]),fill=onezerolist_unifrac[[3]],color='darkolivegreen4',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_unifrac_class[[3]], y = lister_mantel_unifrac[[3]]),color='darkolivegreen4')+
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


c4<-ggplot()+ylim(-0.4,0.4)+
  geom_point(aes(x = lister_mantel_unifrac_class[[4]], y = lister_mantel_unifrac[[4]]),fill=onezerolist_unifrac[[4]],color='navajowhite3',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_unifrac_class[[4]], y = lister_mantel_unifrac[[4]]),color='navajowhite3')+
  geom_point(aes(x = lister_mantel_unifrac_class[[5]], y = lister_mantel_unifrac[[5]]),fill=onezerolist_unifrac[[5]],color='grey40',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_unifrac_class[[5]], y = lister_mantel_unifrac[[5]]),color='grey40')+
  geom_point(aes(x = lister_mantel_unifrac_class[[6]], y = lister_mantel_unifrac[[6]]),fill=onezerolist_unifrac[[6]],color='darkseagreen3',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_unifrac_class[[6]], y = lister_mantel_unifrac[[6]]),color='darkseagreen3')+
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


onezerolist_wunifrac<-list()
for (k in 1:length(group_list)){
  temp<-as.vector(lister_mantel_wunifrac_sig[[k]])
  temp[temp==1]<-newcol[k]
  temp[temp==0]<-'white'
  onezerolist_wunifrac<-list.append(onezerolist_wunifrac,temp)
}


d3<-ggplot()+ylim(-0.4,0.4)+
  geom_point(aes(x = lister_mantel_wunifrac_class[[1]], y = lister_mantel_wunifrac[[1]]),fill=onezerolist_wunifrac[[1]],color='green',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_wunifrac_class[[1]], y = lister_mantel_wunifrac[[1]]),color='green')+
  geom_point(aes(x = lister_mantel_wunifrac_class[[2]], y = lister_mantel_wunifrac[[2]]),fill=onezerolist_wunifrac[[2]],color='greenyellow',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_wunifrac_class[[2]], y = lister_mantel_wunifrac[[2]]),color='greenyellow')+
  geom_point(aes(x = lister_mantel_wunifrac_class[[3]], y = lister_mantel_wunifrac[[3]]),fill=onezerolist_wunifrac[[3]],color='darkolivegreen4',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_wunifrac_class[[3]], y = lister_mantel_wunifrac[[3]]),color='darkolivegreen4')+
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


d4<-ggplot()+ylim(-0.4,0.4)+
  geom_point(aes(x = lister_mantel_wunifrac_class[[4]], y = lister_mantel_wunifrac[[4]]),fill=onezerolist_wunifrac[[4]],color='navajowhite3',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_wunifrac_class[[4]], y = lister_mantel_wunifrac[[4]]),color='navajowhite3')+
  geom_point(aes(x = lister_mantel_wunifrac_class[[5]], y = lister_mantel_wunifrac[[5]]),fill=onezerolist_wunifrac[[5]],color='grey40',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_wunifrac_class[[5]], y = lister_mantel_wunifrac[[5]]),color='grey40')+
  geom_point(aes(x = lister_mantel_wunifrac_class[[6]], y = lister_mantel_wunifrac[[6]]),fill=onezerolist_wunifrac[[6]],color='darkseagreen3',size=3,shape=21) +
  geom_line(aes(x = lister_mantel_wunifrac_class[[6]], y = lister_mantel_wunifrac[[6]]),color='darkseagreen3')+
  labs(title = "",
       x = "Distance Class (m)", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


xv2<-c(0,0.001*(1:1:170))
a5<-ggplot()+
  geom_point(aes(x=glister_x_jaccard[[1]],y=(glister_points_jaccard[[1]])),size=2,fill='green',col='green',alpha=0.2)+xlab('Genetic Distance')+ylab('Similarity')+ylim(0,0.3)+xlim(0,0.17)+
  geom_point(aes(x=glister_x_jaccard[[2]],y=(glister_points_jaccard[[2]])),size=2,fill='greenyellow',col='greenyellow',alpha=0.2)+
  geom_point(aes(x=glister_x_jaccard[[3]],y=(glister_points_jaccard[[3]])),size=2,fill='darkolivegreen4',col='darkolivegreen4',alpha=0.2)+
  geom_line(aes(y=exp(glister_curve_jaccard[[1]]),x=xv2),linewidth=0.5,col='green')+
  geom_line(aes(y=exp(glister_curve_jaccard[[2]]),x=xv2),linewidth=0.5,col='greenyellow')+
  geom_line(aes(y=exp(glister_curve_jaccard[[3]]),x=xv2),linewidth=0.5,col='darkolivegreen4')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

b5<-ggplot()+
  geom_point(aes(x=glister_x_bray[[1]],y=(glister_points_bray[[1]])),size=2,fill='green',col='green',alpha=0.2)+xlab('Genetic Distance')+ylab('Similarity')+ylim(0,0.7)+xlim(0,0.17)+
  geom_point(aes(x=glister_x_bray[[2]],y=(glister_points_bray[[2]])),size=2,fill='greenyellow',col='greenyellow',alpha=0.2)+
  geom_point(aes(x=glister_x_bray[[3]],y=(glister_points_bray[[3]])),size=2,fill='darkolivegreen4',col='darkolivegreen4',alpha=0.2)+
  geom_line(aes(y=exp(glister_curve_bray[[1]]),x=xv2),linewidth=0.5,col='green')+
  geom_line(aes(y=exp(glister_curve_bray[[2]]),x=xv2),linewidth=0.5,col='greenyellow')+
  geom_line(aes(y=exp(glister_curve_bray[[3]]),x=xv2),linewidth=0.5,col='darkolivegreen4')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

c5<-ggplot()+
  geom_point(aes(x=glister_x_unifrac[[1]],y=(glister_points_unifrac[[1]])),size=2,fill='green',col='green',alpha=0.2)+xlab('Genetic Distance')+ylab('Similarity')+ylim(0,0.6)+xlim(0,0.17)+
  geom_point(aes(x=glister_x_unifrac[[2]],y=(glister_points_unifrac[[2]])),size=2,fill='greenyellow',col='greenyellow',alpha=0.2)+
  geom_point(aes(x=glister_x_unifrac[[3]],y=(glister_points_unifrac[[3]])),size=2,fill='darkolivegreen4',col='darkolivegreen4',alpha=0.2)+
  geom_line(aes(y=exp(glister_curve_unifrac[[1]]),x=xv2),linewidth=0.5,col='green')+
  geom_line(aes(y=exp(glister_curve_unifrac[[2]]),x=xv2),linewidth=0.5,col='greenyellow')+
  geom_line(aes(y=exp(glister_curve_unifrac[[3]]),x=xv2),linewidth=0.5,col='darkolivegreen4')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

d5<-ggplot()+
  geom_point(aes(x=glister_x_wunifrac[[1]],y=(glister_points_wunifrac[[1]])),size=2,fill='green',col='green',alpha=0.2)+xlab('Genetic Distance')+ylab('Similarity')+xlim(0,0.17)+
  geom_point(aes(x=glister_x_wunifrac[[2]],y=(glister_points_wunifrac[[2]])),size=2,fill='greenyellow',col='greenyellow',alpha=0.2)+
  geom_point(aes(x=glister_x_wunifrac[[3]],y=(glister_points_wunifrac[[3]])),size=2,fill='darkolivegreen4',col='darkolivegreen4',alpha=0.2)+
  geom_line(aes(y=exp(glister_curve_wunifrac[[1]]),x=xv2),linewidth=0.5,col='green')+
  geom_line(aes(y=exp(glister_curve_wunifrac[[2]]),x=xv2),linewidth=0.5,col='greenyellow')+
  geom_line(aes(y=exp(glister_curve_wunifrac[[3]]),x=xv2),linewidth=0.5,col='darkolivegreen4')+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

onezerolist_jaccard<-list()
for (k in 1:3){
  temp<-as.vector(glister_mantel_jaccard_sig[[k]])
  temp[temp==1]<-newcol[k]
  temp[temp==0]<-'white'
  onezerolist_jaccard<-list.append(onezerolist_jaccard,temp)
}

a6<-ggplot()+ylim(-0.3,0.5)+
  geom_point(aes(x = glister_mantel_jaccard_class[[1]], y = glister_mantel_jaccard[[1]]),fill=onezerolist_jaccard[[1]],color='green',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_jaccard_class[[1]], y = glister_mantel_jaccard[[1]]),color='green')+
  geom_point(aes(x = glister_mantel_jaccard_class[[2]], y = glister_mantel_jaccard[[2]]),fill=onezerolist_jaccard[[2]],color='greenyellow',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_jaccard_class[[2]], y = glister_mantel_jaccard[[2]]),color='greenyellow')+
  geom_point(aes(x = glister_mantel_jaccard_class[[3]], y = glister_mantel_jaccard[[3]]),fill=onezerolist_jaccard[[3]],color='darkolivegreen4',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_jaccard_class[[3]], y = glister_mantel_jaccard[[3]]),color='darkolivegreen4')+
  labs(title = "",
       x = "Genetic Distance Class", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


onezerolist_bray<-list()
for (k in 1:3){
  temp<-as.vector(glister_mantel_bray_sig[[k]])
  temp[temp==1]<-newcol[k]
  temp[temp==0]<-'white'
  onezerolist_bray<-list.append(onezerolist_bray,temp)
}

b6<-ggplot()+ylim(-0.3,0.5)+
  geom_point(aes(x = glister_mantel_bray_class[[1]], y = glister_mantel_bray[[1]]),fill=onezerolist_bray[[1]],color='green',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_bray_class[[1]], y = glister_mantel_bray[[1]]),color='green')+
  geom_point(aes(x = glister_mantel_bray_class[[2]], y = glister_mantel_bray[[2]]),fill=onezerolist_bray[[2]],color='greenyellow',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_bray_class[[2]], y = glister_mantel_bray[[2]]),color='greenyellow')+
  geom_point(aes(x = glister_mantel_bray_class[[3]], y = glister_mantel_bray[[3]]),fill=onezerolist_bray[[3]],color='darkolivegreen4',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_bray_class[[3]], y = glister_mantel_bray[[3]]),color='darkolivegreen4')+
  labs(title = "",
       x = "Genetic Distance Class", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

onezerolist_unifrac<-list()
for (k in 1:3){
  temp<-as.vector(glister_mantel_unifrac_sig[[k]])
  temp[temp==1]<-newcol[k]
  temp[temp==0]<-'white'
  onezerolist_unifrac<-list.append(onezerolist_unifrac,temp)
}
c6<-ggplot()+ylim(-0.3,0.5)+
  geom_point(aes(x = glister_mantel_unifrac_class[[1]], y = glister_mantel_unifrac[[1]]),fill=onezerolist_unifrac[[1]],color='green',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_unifrac_class[[1]], y = glister_mantel_unifrac[[1]]),color='green')+
  geom_point(aes(x = glister_mantel_unifrac_class[[2]], y = glister_mantel_unifrac[[2]]),fill=onezerolist_unifrac[[2]],color='greenyellow',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_unifrac_class[[2]], y = glister_mantel_unifrac[[2]]),color='greenyellow')+
  geom_point(aes(x = glister_mantel_unifrac_class[[3]], y = glister_mantel_unifrac[[3]]),fill=onezerolist_unifrac[[3]],color='darkolivegreen4',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_unifrac_class[[3]], y = glister_mantel_unifrac[[3]]),color='darkolivegreen4')+
  labs(title = "",
       x = "Genetic Distance Class", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


onezerolist_wunifrac<-list()
for (k in 1:3){
  temp<-as.vector(glister_mantel_wunifrac_sig[[k]])
  temp[temp==1]<-newcol[k]
  temp[temp==0]<-'white'
  onezerolist_wunifrac<-list.append(onezerolist_wunifrac,temp)
}
d6<-ggplot()+ylim(-0.3,0.5)+
  geom_point(aes(x = glister_mantel_wunifrac_class[[1]], y = glister_mantel_wunifrac[[1]]),fill=onezerolist_wunifrac[[1]],color='green',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_wunifrac_class[[1]], y = glister_mantel_wunifrac[[1]]),color='green')+
  geom_point(aes(x = glister_mantel_wunifrac_class[[2]], y = glister_mantel_wunifrac[[2]]),fill=onezerolist_wunifrac[[2]],color='greenyellow',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_wunifrac_class[[2]], y = glister_mantel_wunifrac[[2]]),color='greenyellow')+
  geom_point(aes(x = glister_mantel_wunifrac_class[[3]], y = glister_mantel_wunifrac[[3]]),fill=onezerolist_wunifrac[[3]],color='darkolivegreen4',size=3,shape=21) +
  geom_line(aes(x = glister_mantel_wunifrac_class[[3]], y = glister_mantel_wunifrac[[3]]),color='darkolivegreen4')+
  labs(title = "",
       x = "Genetic Distance Class", y = "Mantel Correlation",
       fill = "Significance") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.5)+theme_classic()+theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))



library("gridExtra")
library("grid")
library("cowplot")


grid_layout <- plot_grid(a1, a5, a2, a3, a6, a4, align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4, 1/4))
grid_layout

grid_layout <- plot_grid(b1, b5, b2, b3, b6, b4, align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4, 1/4))
grid_layout

grid_layout <- plot_grid(c1, c5, c2, c3, c6, c4, align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4, 1/4))
grid_layout

grid_layout <- plot_grid(d1, d5, d2, d3, d6, d4, align = "v", nrow = 2, rel_heights = c(1/4, 1/4, 1/4, 1/4), rel_widths = c(1/4, 1/4, 1/4, 1/4))
grid_layout


group<-c(rep('full salamander microbiota',10),rep('salamander exclusive microbiota',10),rep('salamander shared microbiota',10))
xgroup<-1:1:30
ygroup<-1:1:30
dl<-data.frame(group,xgroup,ygroup)
ggplot()+
  geom_point(data=dl, aes(x=xgroup,y=ygroup,color=group))+scale_color_manual(values=c('green','greenyellow','darkolivegreen4'))+theme_set(theme_bw() + theme(legend.key=element_blank())) 


group<-c(rep('full crevice microbiota',10),rep('crevice exclusive microbiota',10),rep('crevice shared microbiota',10))
xgroup<-1:1:30
ygroup<-1:1:30
dl<-data.frame(group,xgroup,ygroup)
dl$group <- factor(dl$group, levels = c('full crevice microbiota','crevice exclusive microbiota','crevice shared microbiota'))
ggplot()+
  geom_point(data=dl, aes(x=xgroup,y=ygroup,color=group))+scale_color_manual(values=c('navajowhite3','grey40','darkseagreen3'))+theme_set(theme_bw() + theme(legend.key=element_blank())) 



