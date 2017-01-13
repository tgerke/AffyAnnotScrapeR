#######################################################################################
# load and clean affymetrix probeset annotation for Human Exon 1.0 ST array

#######################################################################################
# set directory, load libraries

setwd("~/Documents/affy") 
#obviously localy machine dependent for now (can't post the Affy data)

library(Biobase)
library(data.table)

options(stringsAsFactors=FALSE)

#######################################################################################
# load transcript cluster annotation, downloaded from Affymetrix.com 01/10/2017

annot <- fread("HuEx-1_0-st-v2_download011317/HuEx-1_0-st-v2.na36.hg19.probeset.csv", 
               sep=",", header=TRUE, skip=22, showProgress=FALSE)

#######################################################################################
# create lists of data frames which contain mRNA and gene level annotation 

# x will be annot$mrna_assignment
getmrnaannot <- function(x) {
   list1 <- sapply(x, function(y) strsplit(y, " /// "), USE.NAMES=FALSE)
   list1 <- lapply(list1, function(x) {
      # no TRUNCATED fields in the HuExon file
      #anytrunc <- grep("TRUNCATED", x)
      #if (length(anytrunc)>0) {x <- x[-anytrunc]}
      strsplit(x, " // ")
   })
   list2 <- lapply(list1, function(x) {
      dat <- as.data.frame(matrix(unlist(x), nrow=length(x), byrow=TRUE), stringsAsFactors=FALSE)
      if (ncol(dat)==1) {dat <- data.frame(t(rep("---", 6)))}
      names(dat) <- c("accession", "assignment_seqname", "assignment_score", "direct_probes",
                      "possible_probes", "assignment_xhyb")
      return(dat)
   })
}
mrnaannot <- sapply(annot$mrna_assignment, getmrnaannot, USE.NAMES=FALSE)
names(mrnaannot) <- annot$probeset_id

temp <- lapply(list2)

# x will be annot$gene_assignment
getgeneannot <- function(x) {
   list1 <- sapply(x, function(y) strsplit(y, " /// "), USE.NAMES=FALSE)
   list1 <- lapply(list1, function(x) {
      # no TRUNCATED fields in the HuExon file
      #anytrunc <- grep("TRUNCATED", x)
      #if (length(anytrunc)>0) {x <- x[-anytrunc]}
      strsplit(x, " // ")
   })
   list2 <- lapply(list1, function(x) {
      dat <- as.data.frame(matrix(unlist(x), nrow=length(x), byrow=TRUE), stringsAsFactors=FALSE)
      if (ncol(dat)==1) {dat <- data.frame(t(rep("---", 2)))}
      names(dat) <- c("accession", "gene_symbol")
      return(dat)
   })
}
geneannot <- sapply(annot$gene_assignment, getgeneannot, USE.NAMES=FALSE)
names(geneannot) <- annot$probeset_id

