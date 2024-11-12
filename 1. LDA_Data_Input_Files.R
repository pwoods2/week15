## This code is to format the input files for to run LDA for Week 15 of this study.
## It 38 samples from 5 different genotypes: ACKG (n=7), ACK (n=8), HACK (n=8), HACKG (n=8) and WT (n=7). 
## They can also be combined as ACK/G (n=15), HACK/G (n=16) and WT (n=7)


## Import & Cleanup Normalized counts
#counts <- read.csv("../Mouse_Exfoliome_Normalization/Exfoliome-TMM_Normalization.csv")
counts <- read.csv("Extra_Files_Needed/Exfoliome-TMM_Normalization.csv")
library(tibble)
counts <- counts %>% remove_rownames %>% column_to_rownames(var = "X")

## Function to make data noisy...most likely it needs it and it saves a lot of time to do it first
noisy.data <- function(x){
  set.seed(123)
  noise.df <- as.data.frame(matrix(replicate(nrow(x)*ncol(x),rnorm(1,0,0.00001)), ncol = ncol(x), byrow = TRUE))
  data.w.noise <- x + noise.df
  return(data.w.noise)  
}


## Import & Cleanup Metadata
library(openxlsx)
#meta <- as.data.frame(read.xlsx("../../../Mouse_Longitudinal_Study_Metadata-01242023DM.xlsx", sheet = "Joined_Metadata"))
meta <- as.data.frame(read.xlsx("Extra_Files_Needed/Mouse_Longitudinal_Study_Metadata-01242023DM.xlsx", sheet = "Joined_Metadata"))
library(dplyr)
meta <- filter(meta,Timepoint == c("W15"))
meta <- meta[,-c(1:6,8:9)]

## Order metadata by Genotype
meta <- meta[order(factor(meta$Genotype, levels = c('WT','ACKG','ACK','HACKG','HACK'))),]

## Import gene list
#colon.genes <- read.xlsx("/mnt/matrix/roo/Genelists/Full_Gene_List-updated-10232023.xlsx")
colon.genes <- read.xlsx("Extra_Files_Needed/Full_Gene_List-updated-10232023.xlsx")
colon.genes <- filter(colon.genes,Category == c("colon.genes"))
colon.genes <- colon.genes[,-c(1:4,6:11)]
colon.genes <- unique(colon.genes)

## Subset counts data
colon.gene.counts <- counts[rownames(counts) %in% colon.genes,]

## Create Key for annotation
numbering <- 0:(nrow(colon.gene.counts)-1)
numbering <- as.data.frame(numbering)
colon.genes <- cbind(numbering,rownames(colon.gene.counts))
names(colon.genes) <- c("Number","Ensembl.ID")

## Get names for Mouse Genes for Key
#mouse.genes <- read.csv("/mnt/matrix/roo/refs/GRCm38.94-mouse/Mouse_biotypes_annotation.csv")[-1]
mouse.genes <- read.csv("Extra_Files_Needed/Mouse_biotypes_annotation.csv")[-1]
mouse.genes <- mouse.genes[,-c(3:4)]
names(mouse.genes) <- c("Ensembl.ID","Gene.Name")

colon.genes <- left_join(colon.genes,mouse.genes)

write.csv(colon.genes,"Annotation_Files/LDA-Gene-Annotation_Colon-Genes.csv",row.names = FALSE)


## Tidy Data
## Colon genes
colon.genes.num <- paste(nrow(colon.gene.counts))
colon.gene.counts <- as.data.frame(t(colon.gene.counts))
colon.gene.counts <- left_join(meta,rownames_to_column(colon.gene.counts),by=c("Exfoliome.ID"="rowname"))
labels <- colon.gene.counts[,1:3]


#################
## Subset Data ##
#################
WT.ACK <- colon.gene.counts %>% filter(Genotype=='WT'|Genotype=='ACK') %>% column_to_rownames(var = "Exfoliome.ID")
WT.ACK <- WT.ACK[,-c(1:2)]
WT.ACK <- noisy.data(WT.ACK)

WT.ACKG <- colon.gene.counts %>% filter(Genotype=='WT'|Genotype=='ACKG') %>% column_to_rownames(var = "Exfoliome.ID")
WT.ACKG <- WT.ACKG[,-c(1:2)]
WT.ACKG <- noisy.data(WT.ACKG)

WT.HACK <- colon.gene.counts %>% filter(Genotype=='WT'|Genotype=='HACK') %>% column_to_rownames(var = "Exfoliome.ID")
WT.HACK <- WT.HACK[,-c(1:2)]
WT.HACK <- noisy.data(WT.HACK)

WT.HACKG <- colon.gene.counts %>% filter(Genotype=='WT'|Genotype=='HACKG') %>% column_to_rownames(var = "Exfoliome.ID")
WT.HACKG <- WT.HACKG[,-c(1:2)]
WT.HACKG <- noisy.data(WT.HACKG)

HACK.G.WT <- colon.gene.counts %>% filter(Genotype=='WT'|Grouped.Genotypes=='HACK/G') %>% column_to_rownames(var = "Exfoliome.ID")
HACK.G.WT <- HACK.G.WT[,-c(1:2)]
HACK.G.WT <- noisy.data(HACK.G.WT)

ACK.G.WT <- colon.gene.counts %>% filter(Genotype=='WT'|Grouped.Genotypes=='ACK/G') %>% column_to_rownames(var = "Exfoliome.ID")
ACK.G.WT <- ACK.G.WT[,-c(1:2)]
ACK.G.WT <- noisy.data(ACK.G.WT)


########################
## Write Output files ##
########################
## File path to folder to save all input files
fp <- c("../LDA/Input_Files")

## Set the count to match the dataset
genes.count <- colon.genes.num

## ACK-WT File
## File name
filename <- paste(fp,"WT_ACK.txt",sep = "/")
## Adds the numbers at the top of input file indicating the number of samples.
sample.count <- paste(length(which(colon.gene.counts$Genotype=="WT")),",",length(which(colon.gene.counts$Genotype=="ACK")),sep = "")
header <- paste(sample.count,genes.count, sep = "\n")
#Save header & data frame to file
writeLines(header,filename)
write.table(WT.ACK,filename,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE, eol = "\n")

## ACKG-WT File
filename <- paste(fp,"WT_ACKG.txt",sep = "/")
sample.count <- paste(length(which(colon.gene.counts$Genotype=='WT')),",",length(which(colon.gene.counts$Genotype=='ACKG')),sep = "")
header <- paste(sample.count,genes.count, sep = "\n")
writeLines(header,filename)
write.table(WT.ACKG,filename,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE, eol = "\n")

## HACK-WT File
filename <- paste(fp,"WT_HACK.txt",sep = "/")
sample.count <- paste(length(which(colon.gene.counts$Genotype=='WT')),",",length(which(colon.gene.counts$Genotype=='HACK')),sep = "")
header <- paste(sample.count,genes.count, sep = "\n")
writeLines(header,filename)
write.table(WT.HACK,filename,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE, eol = "\n")

## HACKG-WT File
filename <- paste(fp,"WT_HACKG.txt",sep = "/")
sample.count <- paste(length(which(colon.gene.counts$Genotype=='WT')),",",length(which(colon.gene.counts$Genotype=='HACKG')),sep = "")
header <- paste(sample.count,genes.count, sep = "\n")
writeLines(header,filename)
write.table(WT.HACKG,filename,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE, eol = "\n")

## ACK & ACKG-WT File
filename <- paste(fp,"WT_ACKG-ACK.txt",sep = "/")
sample.count <- paste(length(which(colon.gene.counts$Genotype=='WT')),",",length(which(colon.gene.counts$Grouped.Genotypes=='ACK/G')),sep = "")
header <- paste(sample.count,genes.count, sep = "\n")
writeLines(header,filename)
write.table(ACK.G.WT,filename,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE, eol = "\n")

## HACK & HACKG-WT File
filename <- paste(fp,"WT_HACKG-HACK.txt",sep = "/")
sample.count <- paste(length(which(colon.gene.counts$Genotype=='WT')),",",length(which(colon.gene.counts$Grouped.Genotypes=='HACK/G')),sep = "")
header <- paste(sample.count,genes.count, sep = "\n")
writeLines(header,filename)
write.table(HACK.G.WT,filename,append=TRUE,quote=FALSE,row.names = FALSE,col.names = FALSE, eol = "\n")
