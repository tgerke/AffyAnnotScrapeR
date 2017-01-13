#######################################################################################
# load and clean affymetrix transcript cluster annotation for Human Gene 1.0 ST array

#######################################################################################
# set directory, load libraries

setwd("~/Documents/affy") 
#obviously localy machine dependent for now (can't post the Affy data)

library(Biobase)
library(data.table)

options(stringsAsFactors=FALSE)

#######################################################################################
# load transcript cluster annotation, downloaded from Affymetrix.com 01/10/2017

annot <- fread("HuGene-1_0-st-v1_download011017/HuGene-1_0-st-v1.na36.hg19.transcript.csv", 
               sep=",", header=TRUE, skip=23, showProgress=FALSE)

#######################################################################################
# create lists of data frames which contain mRNA and gene level annotation 

# x will be annot$mrna_assignment
getmrnaannot <- function(x) {
   list1 <- sapply(x, function(y) strsplit(y, " /// "), USE.NAMES=FALSE)
   list1 <- lapply(list1, function(x) {
      anytrunc <- grep("TRUNCATED", x)
      if (length(anytrunc)>0) {x <- x[-anytrunc]}
      strsplit(x, " // ")
   })
   list2 <- lapply(list1, function(x) {
      dat <- as.data.frame(matrix(unlist(x), nrow=length(x), byrow=TRUE), stringsAsFactors=FALSE)
      if (ncol(dat)==1) {dat <- data.frame(t(rep("---", 9)))}
      names(dat) <- c("accession", "source_name", "description", "assignment_seqname", "assignment_score",
                      "assignment_coverage", "direct_probes", "possible_probes", "assignment_xhyb")
      return(dat)
   })
}
mrnaannot <- sapply(annot$mrna_assignment, getmrnaannot, USE.NAMES=FALSE)
names(mrnaannot) <- annot$transcript_cluster_id

# x will be annot$gene_assignment
getgeneannot <- function(x) {
   list1 <- sapply(x, function(y) strsplit(y, " /// "), USE.NAMES=FALSE)
   list1 <- lapply(list1, function(x) {
      anytrunc <- grep("TRUNCATED", x)
      if (length(anytrunc)>0) {x <- x[-anytrunc]}
      strsplit(x, " // ")
   })
   list2 <- lapply(list1, function(x) {
      dat <- as.data.frame(matrix(unlist(x), nrow=length(x), byrow=TRUE), stringsAsFactors=FALSE)
      if (ncol(dat)==1) {dat <- data.frame(t(rep("---", 5)))}
      names(dat) <- c("accession", "gene_symbol", "gene_title", "cytoband", "entrez_gene_id")
      return(dat)
   })
}
geneannot <- sapply(annot$gene_assignment, getgeneannot, USE.NAMES=FALSE)
names(geneannot) <- annot$transcript_cluster_id

annot <- annot[,c("mrna_assignment", "gene_assignment"):=NULL]
out <- list(annot, mrnaannot, geneannot)

save(out, file="HuGene10STtranscript.RData")
