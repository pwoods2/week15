## This code is to annotate and perform delta calculations for LDA 
setwd("/Users/pwoods/week15August2024")

library(tidyr) # separate function
library(dplyr)
library(stringr)

##Import & Format Gene List. I only use the columns for the gene index and gene name and drop the ensembl ID for this
gene.key <- read.csv("/Users/pwoods/week15August2024/Annotation_Files/LDA-Gene-Key-03202024.csv")[,c(1,3)] %>%
  setNames(.,c("Gene.Number", "Gene.Name")) # added this to avoid conflict with the Numbering column in gene.key

##Import LDA Results as output by the LDA code. It imports the columns with the gene index and bresub only as written.
feature.1 <- read.csv("Output-Files/1-Feature-Output/1-Feature-WT-HACK-HACKG-Output.txt",sep = "\t")[,-c(2:4)]
feature.2 <- read.csv("Output-Files/2-Feature-Output/2-Feature-WT-HACK-HACKG-Output.txt",sep = "\t")[,-c(2:4)]
feature.3 <- read.csv("Output-Files/3-Feature-Output/3-Feature-WT-HACK-HACKG-Output.txt",sep = "\t")[,-c(2:4)]

## 1 feature annotation
feature.1.annotation <- feature.1 %>% 
  mutate_at("gene.index", str_replace, "\\[", "") %>% ## Removes first bracket from gene index
  mutate_at("gene.index", str_replace, "\\]", "") %>% ## Removes second bracket from gene index
  mutate_all(.,as.numeric) %>%  ## Makes all data numeric
  setNames(.,c("Gene.Number","bresub")) %>% ## Sets column names 
  left_join(.,gene.key) ## Join gene index with key gene index to add gene names


## 2 feature annotation
feature.2.annotation <- feature.2 %>% mutate_at("gene.index", str_replace, "\\[", "") %>% ## Removes first bracket from gene index
  mutate_at("gene.index", str_replace, "\\]", "") %>% ## Removes second bracket from gene index
  separate(gene.index,c("Gene.A.Number","Gene.B.Number")) %>% ## Separates the gene index into new columns
  mutate_all(.,as.numeric) %>% ## Makes all data numeric
  left_join(.,feature.1.annotation, by=c("Gene.A.Number"="Gene.Number")) %>% ## Join gene index with key gene index to add gene names and bresub for gene A
  left_join(.,feature.1.annotation, by=c("Gene.B.Number"="Gene.Number")) %>% ## Join gene index with key gene index to add gene names and bresub for gene B
  replace(is.na(.),max(feature.1.annotation$bresub)) %>% ## Fills any NA values (genes with no bresub from one feature) with the maximum error from 1 feature
  mutate(Delta.Gene.A=bresub.y-bresub.x) %>% ## Creates column Delta.Gene.A & calculates the bresub improvement from 1 to 2 feature classification for gene A
  mutate(Delta.Gene.B=bresub-bresub.x) %>% ## Creates column Delta.Gene.B & calculates the bresub improvement from 1 to 2 feature classification for gene B
  mutate(Smallest.Delta=pmin(Delta.Gene.A,Delta.Gene.B)) %>% ## Creates column for Smallest.Delta and selects the smallest amount of improvement from all delta columns
  select(.,-c(bresub.y,bresub)) %>% ## Remove bresub values from 1 feature used for calculations
  setNames(.,c("Gene.A.Number","Gene.B.Number","bresub","Gene.A.Name","Gene.B.Name","Delta-Gene.A","Delta-Gene.B","Min.Delta")) ## Set column names

## Create Pair-Index on 2 feature annotation helping annotate 3-feature delta calculations
feature.pair.index <- feature.2.annotation %>% 
  mutate(pair.index=str_c(Gene.A.Number,'-',Gene.B.Number)) %>% ## Creates a "pair.index" column in two feature annotation
  select(.,c("pair.index","bresub")) ## Selects only the columns pair.index and bresub to keep

## 3 feature annotation
feature.3.annotation <- feature.3 %>% 
  mutate_at("gene.index", str_replace, "\\[", "") %>% ## Removes first bracket from gene index
  mutate_at("gene.index", str_replace, "\\]", "") %>% ## Removes second bracket from gene index
  separate(gene.index,c("Gene.A.Number","Gene.B.Number","Gene.C.Number")) %>% ## Separates the gene index into new columns
  mutate_all(.,as.numeric) %>% ## Makes all data numeric
  left_join(.,feature.1.annotation,by=c("Gene.A.Number"="Gene.Number")) %>% ## Join gene A index with 1 feature annotation to add gene name and bresub for gene A
  left_join(.,feature.1.annotation,by=c("Gene.B.Number"="Gene.Number")) %>% ## Join gene B index with 1 feature annotation to add gene name and bresub for gene B
  left_join(.,feature.1.annotation,by=c("Gene.C.Number"="Gene.Number")) %>% ## Join gene C index with 1 feature annotation to add gene name and bresub for gene C
  mutate(Delta.Gene.A=bresub.y-bresub.x) %>% ## Creates column Delta.Gene.A & calculates the bresub improvement from 1 to 2 feature classification for gene A
  mutate(Delta.Gene.B=bresub.x.x-bresub.x) %>% ## Creates column Delta.Gene.B & calculates the bresub improvement from 1 to 2 feature classification for gene B
  mutate(Delta.Gene.C=bresub.y.y-bresub.x) %>%  ## Creates column Delta.Gene.C & calculates the bresub improvement from 1 to 2 feature classification for gene C
  mutate(pair.index.AB=str_c(Gene.A.Number,'-',Gene.B.Number)) %>% ## Creates a pair index for Genes A & B
  mutate(pair.index.AC=str_c(Gene.A.Number,'-',Gene.C.Number)) %>% ## Creates a pair index for Genes A & C
  mutate(pair.index.BC=str_c(Gene.B.Number,'-',Gene.C.Number)) %>% ## Creates a pair index for Genes B & C
  left_join(.,feature.pair.index,by=c("pair.index.AB"="pair.index")) %>% ## Join pair index AB with pair index to add and bresub for pair AB
  left_join(.,feature.pair.index,by=c("pair.index.AC"="pair.index")) %>% ## Join pair index AC with pair index to add and bresub for pair AC
  left_join(.,feature.pair.index,by=c("pair.index.BC"="pair.index")) %>% ## Join pair index BC with pair index to add and bresub for pair BC
  mutate(Delta.Pair.AB=bresub.x.x.x-bresub.x) %>% ## Creates column Delta.Gene.AB & calculates the bresub improvement from 2 to 3 feature classification for gene pair AB
  mutate(Delta.Pair.AC=bresub.y.y.y-bresub.x) %>% ## Creates column Delta.Gene.AC & calculates the bresub improvement from 2 to 3 feature classification for gene pair AC
  mutate(Delta.Pair.BC=bresub-bresub.x) %>% ## Creates column Delta.Gene.BC & calculates the bresub improvement from 2 to 3 feature classification for gene pair BC
  replace(is.na(.),max(feature.2.annotation$bresub)) %>% ## Fills any NA values (genes with no bresub from two feature) with the maximum error from 2 feature
  mutate(Min.Delta=pmin(Delta.Gene.A,Delta.Gene.B,Delta.Gene.C,Delta.Pair.AB,Delta.Pair.AC,Delta.Pair.BC)) %>% ## Creates column for Smallest.Delta and selects the smallest amount of improvement from all delta columns
  select(.,-c(bresub.y,bresub.x.x,bresub.y.y,pair.index.AB,pair.index.AC,pair.index.BC,bresub.x.x.x,bresub.y.y.y,bresub)) %>% ## Remove pair.indexes & bresub values from 1 & 2 feature used for calculations
  setNames(.,c("Gene.A.Number","Gene.B.Number","Gene.C.Number","bresub","Gene.A.Name","Gene.B.Name","Gene.C.Name","Delta-Gene.A","Delta-Gene.B","Delta-Gene.C",
               "Delta-Pair.AB","Delta-Pair.AC","Delta-Pair.BC","Min.Delta")) ## Set column names

## Save 1 feature annotated results
write.csv(feature.1.annotation, file="Calculation-Output/1-Feature-WT-HACK-HACKG-Calculation.csv",row.names=FALSE)

## Save 2 feature annotated results
write.csv(feature.2.annotation, file="Calculation-Output/2-Feature-WT-HACK-HACKG-Calculation.csv",row.names=FALSE)

## Save 3 feature annotated results
write.csv(feature.3.annotation, file="Calculation-Output/3-Feature-WT-HACK-HACKG-Calculation.csv",row.names=FALSE)

