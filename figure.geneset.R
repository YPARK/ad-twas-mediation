#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

tab.file <- argv[1] # e.g., tab.file <- 'tables/bootstrap_gene_significant.txt.gz'
tab.tot.file <- argv[2] # e.g., tab.tot.file <- 'tables/bootstrap_gene.txt.gz'
gsc.file <- argv[3] # e.g., gsc.file <- 'genesets/c2.cp.reactome.v6.0.symbols.gmt'
out.file <- argv[4]

if(file.exists(out.file)) q()

options(stringsAsFators = FALSE)

library(dplyr)
library(readr)
library(piano)
library(cba)
library(proxy)
library(ggplot2)

source('util.R')

rm.header <- function(s) paste(strsplit(s, '[_]')[[1]][-1], collapse = '_')

gsc <- loadGSC(gsc.file)
tab <- read_tsv(tab.file, col_names = TRUE)
stat.tab <- read_tsv(tab.tot.file, col_names = TRUE)

################################################################
## gene -> genesets
gsc.tab <- lapply(1:length(gsc$gsc), function(j) {
    gsc.name <- names(gsc$gsc)[j]
    data.frame(hgnc = gsc$gsc[[j]], gs = gsc.name)
})

gsc.tab <- do.call(rbind, gsc.tab)

## q-value cutoff
q.cutoff <- 0.25

.stat <- stat.tab %>% select(hgnc, lodds) %>%
    group_by(hgnc) %>%
    slice(which.max(lodds))

gsc.stat <- .stat$lodds
names(gsc.stat) <- .stat$hgnc

gsa <- runGSA(geneLevelStats = gsc.stat, gsc = gsc, geneSetStat = 'median',
              nPerm = 1e5, gsSizeLim = c(10, 500), verbose = TRUE)              

gsa.tab <- GSAsummaryTable(gsa, save = FALSE)

gsa.tab <- data.frame(gs = gsa.tab[, 'Name'],
                      p.val = gsa.tab[, 'p (dist.dir.up)'],
                      q.val = gsa.tab[, 'p adj (dist.dir.up)'],
                      p.val.mix = gsa.tab[, 'p (mix.dir.up)'],
                      q.val.mix = gsa.tab[, 'p adj (mix.dir.up)']) %>%
                          mutate(gs = sapply(gs, rm.header))

significant.path <- gsa.tab %>%
    dplyr::filter(q.val < q.cutoff | q.val.mix < q.cutoff) %>%
    dplyr::select(gs)

## show top 20 pathways
pdf(file = gsub(out.file, pattern = '.pdf', replacement = '-top20.pdf'), width = 6, height = 8)
GSAheatmap(gsa, cutoff=20, adjusted=FALSE, ncharLabel=50, cellnote='pvalue', columnnames='full',
           colorkey=TRUE, colorgrad=NULL, cex=NULL)
dev.off()

## save the list of p-values
write.tab.named(gsa.tab %>% arrange(p.val), file = gzfile(out.file %&&% '.txt.gz'))

if(nrow(significant.path) < 1) {
    log.msg('No significant pathway\n\n')
    q()
}

significant.path <- unique(as.character(significant.path$gs))

################################################################
## assign membership and focus on significant ones
tab.annot <- tab %>% dplyr::select(hgnc, theta, theta.se) %>%
    right_join(gsc.tab, by = 'hgnc') %>% na.omit() %>%
    mutate(path = sapply(gs, rm.header)) %>%
    dplyr::filter(path %in% significant.path)

if(nrow(tab.annot) < 1) {
    log.msg('No overlap with significant genes\n\n')
    q()
}

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

w <- length(genes) * .15 + 1 + max(sapply(pathways, nchar)) * .05
h <- length(pathways) * .1 + 1

ggsave(filename = out.file, plot = plt, width = w, height = h)
       
log.msg('Saved %s\n\n', out.file)
