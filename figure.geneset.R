#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

tab.file <- argv[1] # e.g., 'tables/bootstrap_gene_significant.txt.gz'
tab.tot.file <- argv[2] # e.g., 'tables/bootstrap_gene.txt.gz'
gsc.file <- argv[3] # e.g., 'genesets/c2.cp.reactome.v6.0.symbols.gmt'
out.file <- argv[4]

options(stringsAsFators = FALSE)

library(dplyr)
library(readr)
library(piano)
library(cba)
library(proxy)

source('util.R')
source('figure.util.R')

gsc <- loadGSC(gsc.file)
tab <- read_tsv(tab.file, col_names = TRUE)

## p-value cutoff
q.cutoff <- 0.20

## basic test with piano
stat.tab <- read_tsv(tab.tot.file, col_names = TRUE)
.stat <- stat.tab %>% select(hgnc, lodds) %>%
    group_by(hgnc) %>%
    slice(which.max(lodds))

gsc.stat <- .stat$lodds
names(gsc.stat) <- .stat$hgnc

gsa <- runGSA(gsc.stat, gsc = gsc, geneSetStat = 'mean', nPerm = 5e4,
              gsSizeLim = c(10, 1000))

gsa.tab <- GSAsummaryTable(gsa, save = FALSE)

gsa.tab <- data.frame(gs = gsa.tab[, 'Name'],
                      p.val = gsa.tab[, 'p (dist.dir.up)'],
                      q.val = gsa.tab[, 'p adj (dist.dir.up)'])

significant.path <- gsa.tab %>% dplyr::filter(q.val < q.cutoff) %>%
    dplyr::select(gs)

if(nrow(significant.path) < 1) q()

rm.header <- function(s) paste(strsplit(s, '[_]')[[1]][-1], collapse = '_')

significant.path <- as.character(significant.path$gs)
significant.path <- as.vector(sapply(significant.path, rm.header))

################################################################
## gene -> genesets
gsc.tab <- lapply(1:length(gsc$gsc), function(j) {
    gsc.name <- names(gsc$gsc)[j]
    data.frame(hgnc = gsc$gsc[[j]], gs = gsc.name)
})

gsc.tab <- do.call(rbind, gsc.tab)

## assign membership and focus on significant ones
tab.annot <- tab %>% dplyr::select(hgnc, theta, theta.se) %>%
    right_join(gsc.tab, by = 'hgnc') %>% na.omit() %>%
    mutate(path = sapply(gs, rm.header)) %>%
    dplyr::filter(path %in% significant.path)

if(nrow(tab.annot) < 1) q()

## sort rows and columns -- by membership matrix
genes <- unique(tab.annot$hgnc)
pathways <- unique(tab.annot$path)

tab.annot <- tab.annot %>%
    mutate(gene.idx = match(hgnc, genes),
           path.idx = match(path, pathways))

M <- matrix(0, nrow = length(pathways), ncol = length(genes))
tab.idx <- tab.annot %>% dplyr::select(path.idx, gene.idx) %>% as.matrix()
M[tab.idx] <- 1

optimal.row.order <- function(mat) {
    dist.fn <- function(x, y) 1 - cor(x, y, method = 'spearman')
    dd <- proxy::dist(mat, method = dist.fn)
    dd[!is.finite(dd)] <- 0
    hh <- hclust(dd)
    oo <- cba::order.optimal(dd, merge = hh$merge)
    return(oo$order)
}

if(nrow(M) > 2 && ncol(M) > 1) {
    path.order <- optimal.row.order(M)
} else {
    path.order <- 1:nrow(M)
}

if(ncol(M) > 2 && nrow(M) > 1) {
    gene.order <- optimal.row.order(t(M))
} else {
    gene.order <- 1:ncol(M)
}

tab.sort <- tab.annot
tab.sort$path <- factor(tab.sort$path, pathways[path.order])
tab.sort$hgnc <- factor(tab.sort$hgnc, genes[gene.order])

plt <- ggplot(tab.sort) + theme_bw() +
    geom_tile(aes(x = hgnc, y = path, fill = sign(theta)),
              size = .5, color = 'gray20') +
    scale_fill_gradient2(low = 'blue', high = 'yellow', guide = FALSE) +
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 90,
              hjust = 1, vjust = .5),
          axis.title = element_blank())

w <- length(genes) * .1 + 1 + max(sapply(pathways, nchar)) * .05
h <- length(pathways) * .1 + 1

ggsave(filename = out.file, plot = plt, width = w, height = h)
       
