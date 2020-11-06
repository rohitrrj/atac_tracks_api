library(Gviz)
library(rtracklayer)
library(trackViewer)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(httr)
require(jsonlite)

getGeneIDsFromTxDb_RRJ <- function(gr, txdb){
  stopifnot(is(gr, "GRanges"))
  stopifnot(length(gr)>0)
  stopifnot(is(txdb, "TxDb"))
  if(length(gr)>1){
    warning("The length of gr is greater than 1. Only first genomic location will be used.")
    gr <- gr[1]
  }
  genes <- genes(txdb, columns="gene_id", single.strand.genes.only=FALSE)
  genes <- subsetByOverlaps(genes, gr)
  return(names(genes))
}
