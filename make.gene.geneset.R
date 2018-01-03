#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
source('figure.util.R')
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(qvalue)
library(GSEABase)
len <- function(...) length(geneIds(...))


cammel.tab.file <- 'tables/genes/full-0/significant_genes.txt.gz'
tot.tab.file <- 'tables/genes/full-0/total_genes.txt.gz'

cammel.tab <- read_tsv(cammel.tab.file)
tot.tab <- read_tsv(tot.tab.file)

geneset.file <- 'genesets/c2.cp.kegg.v6.0.symbols.gmt'

geneset <- getGmt(geneset.file)

geneset.genes <- unique(do.call(c, lapply(geneset, geneIds)))
tot.genes <- GeneSet(unique(tot.tab$hgnc), setName = 'Tot')
geneset.tot <- GeneSet(geneset.genes, setName = 'geneset') & tot.genes
cammel.genes <- GeneSet(unique(cammel.tab$hgnc), setName = 'CaMMEL') & geneset.tot
n.cammel <- len(cammel.genes)
n.tot <- len(geneset.tot)

## x, q: vector of quantiles representing the number of white balls
##       drawn without replacement from an urn which contains both
##       black and white balls.
##
##    m: the number of white balls in the urn.
##    n: the number of black balls in the urn.
##    k: the number of balls drawn from the urn.
##    p: probability, it must be between 0 and 1.
hg.p <- function(gs) {
    q <- len(gs & cammel.genes)
    m <- len(gs & tot.genes)
    n <- n.tot - m
    k <- n.cammel
    pv <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
    data.frame(set.name = setName(gs), set.size = m, overlap = q, p.val = pv)
}

out <- do.call(rbind, lapply(geneset, hg.p))
out <- data.frame(out, q.val = p.adjust(out$p.val, method = 'fdr')) %>% arrange(p.val)

dir.create('genesets/result', recursive = TRUE)
write_tsv(out, path = 'genesets/result/kegg.txt')

out.pandoc <- ret <- pandoc.table.return(out, row.names = FALSE, style = 'simple',
                                         split.tables = 200, digits = 2)

cat(out.pandoc, file = 'genesets/result/kegg.md')



