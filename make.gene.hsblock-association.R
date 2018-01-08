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


mem.file <- 'network/hsblock/functional-slim/final.dict.gz'
pair.file <- 'network/functional-slim.pairs'
geneset.file <- 'genesets/c2.cp.kegg.v6.0.symbols.gmt'

module.map.file <- 'network/hsblock/functional-slim/module_map.txt'
module.hg.file <- 'network/hsblock/functional-slim/module_kegg.txt.gz'
module.fig.file <- 'network/hsblock/functional-slim/module_kegg.pdf'
mem.out.file <- 'network/hsblock/functional-slim/module.dict'
cam.out.file <- 'network/hsblock/functional-slim/cammel.dict'

pair.tab <- read_tsv(pair.file, col_names = c('gene1', 'gene2'))
cammel.tab.file <- 'tables/genes/full-0/significant_genes.txt.gz'
tot.tab.file <- 'tables/genes/full-0/total_genes.txt.gz'
cammel.tab <- read_tsv(cammel.tab.file)
tot.tab <- read_tsv(tot.tab.file)



## look at modules containg at least one of signficant genes
mem.tab <- read_tsv(mem.file, col_names = c('hgnc', 'module'))

modules <- cammel.tab %>% left_join(mem.tab) %>%
    dplyr::select(module) %>% na.omit() %>% unique()

mem.tab <- mem.tab %>% right_join(modules) %>% na.omit()

################################################################
## test geneset enrichment
tot.genes <- GeneSet(unique(tot.tab$hgnc), setName = 'Tot')

geneset <- getGmt(geneset.file)
geneset.genes <- unique(do.call(c, lapply(geneset, geneIds)))

geneset.tot <- GeneSet(geneset.genes, setName = 'geneset') & tot.genes

n.tot <- len(geneset.tot)

.mk.module.gs <- function(m) GeneSet(unique(subset(mem.tab, module == m)$hgnc), setName = as.character(m))

modules.gs <- lapply(modules$module, .mk.module.gs)

## x, q: vector of quantiles representing the number of white balls
##       drawn without replacement from an urn which contains both
##       black and white balls.
##
##    m: the number of white balls in the urn.
##    n: the number of black balls in the urn.
##    k: the number of balls drawn from the urn.
##    p: probability, it must be between 0 and 1.
hg.p <- function(gs, module.genes) {
    n.module <- len(module.genes)
    q <- len(gs & module.genes)
    m <- len(gs & tot.genes)
    n <- n.tot - m
    k <- n.module
    pv <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
    data.frame(set.name = setName(gs),
               module.name = setName(module.genes),
               set.size = m,
               module.size = k,
               overlap = q,
               p.val = pv)
}

scan.gs <- function(mm) {
    ret <- do.call(rbind, lapply(geneset, hg.p, module.genes = mm))
    print(ret %>% arrange(p.val) %>% head(3), 3)
    return(ret)
}

hg.tab <- do.call(rbind, lapply(modules.gs, scan.gs) )

################################################################
## sort modules and rename them

hg.pair <- hg.tab %>%
    dplyr::rename(row = module.name, col = set.name) %>%
    dplyr::filter(overlap > 0) %>%
    dplyr::mutate(weight = overlap / set.size) %>%
    dplyr::select(row, col, weight)

hg.order <- order.pair(hg.pair)

module.map <- data.frame(module = as.integer(hg.order$rows)) %>%
    dplyr::mutate(module.sort = 1:n()) %>%
    dplyr::mutate(module.name = module)

cammel.annot <- cammel.tab %>% left_join(mem.tab) %>%
    left_join(module.map)

size.tab <- cammel.annot %>% dplyr::group_by(module.sort) %>% na.omit() %>%
    dplyr::summarize(n.genes = n())


plt.size <- gg.plot(size.tab) +
    geom_bar(aes(x = module.sort, y = n.genes), stat = 'identity')

write_tsv(module.map %>% left_join(size.tab) %>% dplyr::select(module, module.sort, n.genes),
          path = module.map.file)

out.tab <- hg.tab %>%
    mutate(module.name = as.integer(module.name)) %>%
    left_join(module.map) %>% left_join(size.tab) %>%
    dplyr::select(-module.name)

write_tsv(out.tab, path = gzfile(module.hg.file))

## show only FWER < 0.05 pathways
hg.filter <- hg.tab %>% dplyr::filter(p.adjust(p.val) < 0.05) %>%
    dplyr::select(module.name, set.name, p.val) %>% unique()

p.cutoff <- max(hg.filter$p.val)

hg.tab.ordered <-
    out.tab %>%
    dplyr::filter(set.name %in% hg.filter$set.name) %>%
    dplyr::filter(module %in% hg.filter$module.name) %>%
    dplyr::filter(n.genes >= 2, overlap / module.size >= .2)

hg.tab.ordered$set.name <- factor(hg.tab.ordered$set.name,
                                  hg.order$cols,
                                  sapply(hg.order$cols, gsub, pattern = 'KEGG_', replacement = ''))

.aes <- aes(x = as.factor(module.sort),
            y = set.name,
            fill = overlap/module.size)

plt.pval <-
    gg.plot() +
    geom_tile(data = hg.tab.ordered, .aes, color = 'black') +
    geom_point(data = hg.tab.ordered %>% dplyr::filter(p.val <= p.cutoff),
               .aes, pch = 4, color = '#0055FF', size = .5) +     
    scale_fill_continuous('overlap', low = 'white', high = 'red') +
    scale_x_discrete(position = 'top') +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 0),
          axis.title = element_blank())

ggsave(filename = module.fig.file, plot = plt.pval, width = 10, height = 10, units = 'in',
       useDingbats = FALSE)

## output vertex with module labels
mem.out.tab <- mem.tab %>%
    left_join(out.tab %>% dplyr::select(module, module.sort) %>% unique()) %>%
    dplyr::select(-module)

write_tsv(mem.out.tab, path = mem.out.file)

## output cammel genes matched
cam.out.tab <- cammel.tab %>% left_join(mem.out.tab) %>%
    na.omit()

write_tsv(cam.out.tab, path = cam.out.file)
