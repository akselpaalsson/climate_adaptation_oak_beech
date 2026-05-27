

library(LEA)
library(adegenet)
library(viridis)

library(vcfR)

rm(list=ls())

#### Set variables ####
labsp <- c("Quercus spp.")
cols <- c("#E69F00")
Ks <- 10
opt.K <- c(4)

setwd("/path/to/working_directory")


#### Set the paths ####
getwd()
dir.path <- "../res/"
dat.path <- "../data/"


#### Set the working directory ####
setwd(dir.path)
getwd()


#### Set the species ####
sp <- 1 # 

#### Convert the genotype file from ped to lfmm format (might have to be done in R if R studio crashes!) ####
### The ped file was created from a vcf file
gc() #free memory
LEA::ped2lfmm(input.file=paste(dat.path,"populations_no277.snps_minDP3_meanDP10_maxDP80_NA0.90_reduced_HWE0.0001.MAF0.01.thin1kb.snps.ped",sep=""),output.file=paste(dat.path, "populations_no277.snps_minDP3_meanDP10_maxDP80_NA0.90_reduced_HWE0.0001.MAF0.01.thin1kb.snps.lfmm",sep=""), force=T)

#run sNMF (might have to be done in R if R studio crashes!)
project.missing <- snmf(input.file=paste(dat.path, "populations_no277.snps_minDP3_meanDP10_maxDP80_NA0.90_reduced_HWE0.0001.MAF0.01.thin1kb.snps.lfmm",sep=""),
                        K=1:10, entropy=T, repetitions=10, CPU=8, project="new") # http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf



#### Inspect cross-entropy criterion across K values ####
proj.snmf <- load.snmfProject(paste(dat.path,"populations_no277.snps_minDP3_meanDP10_maxDP80_NA0.90_reduced_HWE0.0001.MAF0.01.thin1kb.snps.snmfProject",sep=""))
pdf(paste(dir.path, "snmf/SNMF_CrossEntropy.pdf", sep=""), width=8.24*0.5, height=8.24*0.5)
plot(proj.snmf, cex=1.2, col="blue", pch=19)
dev.off()
#define best K
K <- 5
##use above as additional diagnostic


#### Impute data for missing genotypes using ancestry coefficients #####run sNMF (might have to be done in R if R studio crashes!)
#best <- which.min(cross.entropy(proj.snmf, K=K))
#impute(proj.snmf, paste(dat.path,"populations_minDP3_meanDP10_maxDP80_NA0.95_reduced_HWE0.0001.MAF0.05.snps.lfmm",sep=""), method="mode", K=opt.K[sp], run=best)

#### Visualise barplot of the Q-matrix across K values ####
# load the metadata file for sample names
sample.info <- read.table(paste(dat.path, "Common-Ring.sample.info.txt",sep=""), sep="\t", header=T)
head(sample.info)
sample.vcf <- read.table(paste(dat.path, "samples_no77_NA.90_HWE0.0001.MAf0.01_think1kb.txt",sep=""), sep="\t", header=F)
colnames(sample.vcf)[1] <- "Sample"
head(sample.vcf)
sample.info.reduced <- merge(sample.info, sample.vcf, by.x="SampleID", by.y="Sample", sort=F)
dim(sample.info.reduced)
head(sample.info.reduced)
head(sample.info.reduced[1:7])


#head(dapc4$posterior)
#K4_ap <- cbind(sample.info.reduced[1:7],dapc4$posterior)
#head(K4_ap)
#sort
#K4_ap <- K4_ap[order(K4_ap$Prov2, K4_ap$Site),]
#head(K4_ap)

# finally we loop the SNMF from K=2 to K=10 to show the development of cluster membership (takes some time...)

gc() #clear memory
par(mfrow=c(3,1))
for (clust in 2:10) 
{
  best <- which.min(cross.entropy(proj.snmf, K=clust))
  ap_best <- read.table(paste(dat.path,"populations_no277.snps_minDP3_meanDP10_maxDP80_NA0.90_reduced_HWE0.0001.MAF0.01.thin1kb.snps.snmf/K", clust,"/run",best,"/populations_no277.snps_minDP3_meanDP10_maxDP80_NA0.90_reduced_HWE0.0001.MAF0.01.thin1kb.snps_r",best,".",clust,".Q", sep=""))
  assign("ap", cbind(sample.info.reduced[1:7],ap_best))
  ap <- ap[order(ap$Prov2, ap$Site),]
  #for the comparison with DAPC, we also write a file with assignment probabilities
  write.table(ap, paste("snmf/ap_DAPC_K", clust,".txt", sep=""), sep="\t", quote=F, row.names=F)
  ind.nb <- as.data.frame(table(ap$Prov2)); colnames(ind.nb) <- c("Population","Individuals")
  axis.position <-c()
  assign("ap_1", as.matrix(t(ap)))
  barplot(ap_1[8:(8+clust-1),], main=paste("Posterior assignment probailities K=", clust, sep=""), col=viridis(clust), space=0,border=NA, xaxt='n', las=3)
  for (j in 1:NROW(ind.nb)) 
  {
    axis.position[j+1] <- sum(axis.position[j], ind.nb[j,"Individuals"])
    abline(v=axis.position, col="black", lty=1, lwd=1)
    text.position <- NULL
    text.position <- axis.position[j+1] - (ind.nb[j,"Individuals"] / 2)
    text(text.position, -0.2, ind.nb[j,"Population"], srt=90, cex=1.0, xpd=T, las=1)
  }   
}

#####



#### Set the paths ####
getwd()
dir.path <- "../res/"
dat.path <- "../data/"

setwd("/path/to/working_directory")
assignments<-read.table("populations_no277.snps_minDP3_meanDP10_maxDP80_NA0.90_reduced_HWE0.0001.MAF0.01.thin1kb.snps_r1.3.Q",  sep="", header=FALSE)

#### Set the working directory ####
setwd(dir.path)
getwd()

sample.info <- read.table(paste(dat.path, "Common-Ring.sample.info.txt",sep=""), sep="\t", header=T)
head(sample.info)
sample.vcf <- read.table(paste(dat.path, "samples_no77_NA.90_HWE0.0001.MAf0.01_think1kb.txt",sep=""), sep="\t", header=F)
colnames(sample.vcf)[1] <- "Sample"
head(sample.vcf)
sample.info.reduced <- merge(sample.info, sample.vcf, by.x="SampleID", by.y="Sample", sort=F)
dim(sample.info.reduced)
head(sample.info.reduced)
head(sample.info.reduced[1:7])
colnames(assignments)[1] <- "other1"
colnames(assignments)[2] <- "pet"
colnames(assignments)[3] <- "other2"
dev.off()
hist(assignments$other1)
hist(assignments$pet)
hist(assignments$other2)

samples_and_assignments<-cbind(sample.info.reduced,assignments)
petrea_only <- subset(samples_and_assignments, pet >.5)
head(assignments)
#this keeps 489 samples

sampleID_column <- as.data.frame(petrea_only$SampleID); colnames(sampleID_column)[1] <- "SampleID"

write.table(sampleID_column, file="petrea_0.5.txt", row.names = FALSE)
