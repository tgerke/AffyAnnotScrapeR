#######################################################################################
# load and clean affymetrix transcript cluster annotation for Human Gene 1.0 ST array

#######################################################################################
# set directory, load libraries

#setwd("~/Documents/affy") 
#obviously localy machine dependent for now (can't post the Affy data)

library(Biobase)
library(data.table)

options(stringsAsFactors=FALSE)

#######################################################################################
# load transcript cluster annotation, downloaded from Affymetrix.com 01/10/2017

#windows
annot <- fread("M:/lab/Lab_Gerke/ukb/AXIOM/Axiom_UKB_WCSG.na35.annot.csv", 
               sep=",", header=TRUE, skip=19, showProgress=FALSE)
#mac
annot <- fread("/Volumes/Lab_Gerke/ukb/AXIOM/Axiom_UKB_WCSG.na35.annot.csv", 
               sep=",", header=TRUE, skip=19, showProgress=FALSE)

#######################################################################################
# create lists of data frames which contain gene level annotation 

# x will be annot[,"Associated Gene"]
getgeneannot <- function(x) {
   list1 <- sapply(x, function(y) strsplit(y, " /// "), USE.NAMES=FALSE)
   list1 <- lapply(list1, function(x) {
      strsplit(x, " // ")
   })
   list2 <- lapply(list1, function(x) {
      dat <- as.data.frame(matrix(unlist(x), nrow=length(x), byrow=TRUE), stringsAsFactors=FALSE)
      if (ncol(dat)==1) {dat <- data.frame(t(rep("---", 7)))}
      names(dat) <- c("transcript_accession", "SNP_gene_relationship", "distance", "unigene_cluster_id", "gene","ncbi_gene_name","genbank_description")
      return(dat)
   })
}
geneannot <- sapply(annot[,"Associated Gene"], getgeneannot, USE.NAMES=FALSE)
names(geneannot) <- annot[,"Probe Set ID"]

