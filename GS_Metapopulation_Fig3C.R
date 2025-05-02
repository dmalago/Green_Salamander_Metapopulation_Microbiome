#Adapted from Adam Burns - 2/10/2015
#From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
#Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE.
#spp: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
#pool: A community table for defining source community (optional; Default=NULL).
#taxon: A table listing the taxonomic calls for each otu, with OTU ids as row names and taxonomic classifications as columns.
#If stats=TRUE the function will return fitting statistics.
#If stats=FALSE the function will return a table of observed and predicted values for each otu.

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
library(minpack.lm)
library(Hmisc)
library(stats4)

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

Cphyloseq<-subset_samples(data, Sample_Type == "Crevice")
Sphyloseq<-subset_samples(data, Sample_Type == "Salamander")


pool<-NULL
pps<-t(otu_table(Sphyloseq))
spp<-t(otu_table(Cphyloseq))
taxon<-rownames(otu_table(Cphyloseq))

options(warn=-1)

#Calculate the number of individuals per community
Nfit <- mean(apply(spp, 1, sum))

#Calculate the average relative abundance of each taxa across communities
if(is.null(pool)){
  p.m <- apply(spp, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/Nfit
} else {
  p.m <- apply(pool, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/Nfit
}

#Calculate the occurrence frequency of each taxa across communities
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]

#Combine
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]

#pps<-pps[,which(colnames(pps) %in% names(p))]
N <- mean(apply(pps, 1, sum))

#Calculate the limit of detection
d = 1/N

##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
m.ci <- confint(m.fit, 'm', level=0.95)

##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
sncm.LL <- function(m, sigma){
  R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
  R = dnorm(R, 0, sigma)
  -sum(log(R))
}
m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))

##Calculate Akaike's Information Criterion (AIC)
aic.fit <- AIC(m.mle, k=2)
bic.fit <- BIC(m.mle)

##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))

pred.ci <- binconf(freq.pred*nrow(pps), nrow(pps), alpha=0.01, method="wilson", return.df=TRUE)
pred.ci5 <- binconf(freq.pred*nrow(pps), nrow(pps), alpha=0.05, method="wilson", return.df=TRUE)





pps2<-pps[,which(colnames(pps) %in% names(p))]



#Calculate the occurrence frequency of each taxa across communities
pps2.bi <- 1*(pps2>0)
freq <- apply(pps2.bi, 2, mean)

#Combine
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]  

dataf<-data.frame(p,freq,freq.pred,pred.ci[,2:3],pred.ci5[,2:3])
colnames(dataf)<-c('p','freq','freq.pred','pred.lwr','pred.upr','pred.lwr5','pred.upr5')


col<-c()  
for (k in 1:length(p)){
  if (dataf$freq[k]>dataf$pred.upr[k]){
    col<-c(col,'red')
  }else if (dataf$freq[k]<dataf$pred.lwr[k]){
    col<-c(col,'blue')
  }else{
    col<-c(col,'black')
  }
}


ggplot() +
  geom_point(data=dataf, aes(x=log10(p), y=freq), pch=21, fill=col, col=col,alpha=.2,size=1.5)+
  geom_line(color='black', data=dataf, size=1, aes(y=freq.pred, x=log10(p)),linewidth=0.5, alpha=.25) +
  geom_line(color='black', lty='twodash', size=1, data=dataf, aes(y=pred.upr, x=log10(p)),linewidth=0.5, alpha=.25)+
  geom_line(color='black', lty='twodash', size=1, data=dataf, aes(y=pred.lwr, x=log10(p)),linewidth=0.5, alpha=.25)+
  geom_line(color='black', lty='dotted', size=1, data=dataf, aes(y=pred.upr5, x=log10(p)),linewidth=0.5, alpha=.25)+
  geom_line(color='black', lty='dotted', size=1, data=dataf, aes(y=pred.lwr5, x=log10(p)),linewidth=0.5, alpha=.25)+
  labs(x="log10(mean relative abundance in crevices)", y="Occupancy on salamanders")+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))


pps3<-pps[,which(!(colnames(pps) %in% names(p)))]
pps3.bi <- 1*(pps3>0)
freq3 <- apply(pps3.bi, 2, mean)
freq3 <- freq3[freq3 != 0]
ggplot() +
  geom_point(aes(x=rep(0,length(freq3)), y=freq3), pch=21, fill='red', col='red',alpha=.2,size=1.5)+xlim(-0.01,0.01)+ylim(0,1)+scale_x_discrete(limits=c(-1,0,1))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",linewidth=0.01),panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

