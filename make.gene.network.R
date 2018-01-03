#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
source('figure.util.R')
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

cammel.tab.file <- 'tables/genes/full-0/significant_genes.txt.gz'
tot.tab.file <- 'tables/genes/full-0/total_genes.txt.gz'
net.file <- 'network/BIOGRID-ORGANISM-Homo_sapiens-3.4.155.tab2.txt.gz'

deg.fig.file <- 'network/degree.pdf'

cammel.tab <- read_tsv(cammel.tab.file)
tot.tab <- read_tsv(tot.tab.file)
net.tab <- read_tsv(net.file)

## overall PPI network
net.overall <- net.tab %>%
    rename(pair.type = 'Experimental System Type') %>%
    filter(pair.type == 'physical') %>%
    select(starts_with('Official')) %>%
    rename(gene.1 = 'Official Symbol Interactor A',
           gene.2 = 'Official Symbol Interactor B') %>%
    filter(gene.1 != gene.2) %>%
    filter(gene.1 %in% tot.tab$hgnc, gene.2 %in% tot.tab$hgnc) %>%
    unique()

net.dup <- net.overall %>%
    mutate(temp = gene.2, gene.2 = gene.1, gene.1 = temp) %>%
    select(-temp) %>% rbind(net.overall) %>%
    unique()

degree.tab <- net.dup %>% group_by(gene.1) %>%
    summarize(d = n()) %>%
    rename(gene = gene.1)


get.deg.hist <- function(deg.tab) {
    ret <- deg.tab %>% group_by(d) %>%
        summarize(freq = n())
    ntot <- sum(ret$freq)
    ret <- ret %>% mutate(pr = freq / ntot)
    return(ret)
}

degree.cammel <- degree.tab %>% filter(gene %in% cammel.tab$hgnc)

degree.hist <- get.deg.hist(degree.tab)
degree.hist.cammel <- get.deg.hist(degree.cammel)

deg.df <- rbind(degree.hist %>% mutate(network = 'overall'),
                degree.hist.cammel %>% mutate(network = 'CaMMEL'))

## q1. degree distributions differ?
plt.deg <- gg.plot(deg.df) +
    geom_point(aes(x = d, y = pr, shape = network, color = network)) +
    scale_x_log10(breaks = c(1, 5, 10, 50, 100, 500, 1000)) +
    scale_y_log10(breaks = c(0.1, 0.01, 0.001, 1e-4)) +
    xlab('Degree') + ylab('Empirical Probability') +
    scale_shape_manual(values = c(19, 3)) +
    theme(legend.position = c(.9,.9), legend.justification = c(1,1))

ggsave(filename = deg.fig.file, plot = plt.deg, units = 'in', width = 4, height = 4,
       useDingbats = FALSE)

## q2. more interaction between the selected genes?
gene.in.out <- function(selected) {

    .selected <- selected %>% mutate(idx = 1:n())
    .net.pairs <- left_join(.selected, net.dup)
    .net.pairs.in <- .net.pairs %>% filter(gene.2 %in% .selected$gene.1)
    .net.pairs.out <- .net.pairs %>% filter(!(gene.2 %in% .selected$gene.1))

    deg.comp <- .net.pairs %>% group_by(idx, gene.1) %>%
        summarize(d.tot = n()) %>%
            left_join(.net.pairs.in %>% group_by(idx,gene.1) %>% summarize(d.in = n()) %>% as.data.frame()) %>%
                mutate(d.in = ifelse(is.na(d.in), 0, d.in)) %>%
                    left_join(.net.pairs.out %>% group_by(idx,gene.1) %>% summarize(d.out= n()) %>% as.data.frame()) %>%
                        mutate(d.out = ifelse(is.na(d.out), 0, d.out)) %>%
                            mutate(freq = d.in / (d.in + d.out)) %>%
                                as.data.frame()

    return(deg.comp)
}

selected <- cammel.tab %>%
    filter(hgnc %in% net.dup$gene.1) %>%
    select(hgnc) %>% unique() %>%
    rename(gene.1 = hgnc)

## how to select genes 
deg.matched <-
    degree.cammel %>%
    left_join(degree.tab %>% rename(null.gene = gene), by = 'd') %>%
    as.data.frame()

################################################################
n.perm <- 1000
n.genes <- nrow(obs.stat)
r.idx <- rep(1:n.perm, n.genes)

selected.null <- deg.matched %>% group_by(gene) %>%
    sample_n(size = n.perm, replace = TRUE) %>% as.data.frame() %>%
    rename(gene.1 = null.gene) %>%
    cbind(r = r.idx)

null.gene.in.out <- function(.r) {
    .null.stat <- gene.in.out(selected.null %>% filter(r == .r))
    .null.stat %>%
        summarize(d.in = mean(d.in), d.out = mean(d.out)) %>%
            mutate(r = .r) %>%
                as.data.frame()
}

null.stat <- do.call(rbind, lapply(1:n.perm, null.gene.in.out))

obs.stat <- gene.in.out(selected) %>%
    summarize(d.in = mean(d.in), d.out = mean(d.out))


stat.tab <- rbind(obs.stat %>% mutate(r = 'obs'), null.stat %>% mutate(r = 'r.' %&&% r))
write_tsv(stat.tab, path = 'network/degree_stat.txt')
## p-value 0.76

################################################################
net.pairs <- net.dup %>%
    filter(gene.1 %in% cammel.tab$hgnc, gene.2 %in% cammel.tab$hgnc) %>%
    unique()

temp <- net.dup %>%
    filter(gene.1 %in% cammel.tab$hgnc | gene.2 %in% cammel.tab$hgnc) %>%
    unique()

temp <- union(temp$gene.1, temp$gene.2)

net.pairs.relaxed <- net.dup %>%
    filter(gene.1 %in% temp, gene.2 %in% temp) %>%
    unique()

################################################################
write_tsv(net.pairs, path = 'network/ppi-strict.pairs', col_names = FALSE)
write_tsv(net.pairs.relaxed, path = 'network/ppi.pairs', col_names = FALSE)
