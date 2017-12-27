#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
library(dplyr)
library(tidyr)
library(readr)

gene.tab.file <- 'tables/genes/full-0/significant_genes.txt.gz'
net.file <- 'network/BIOGRID-ORGANISM-Homo_sapiens-3.4.155.tab2.txt.gz'

gene.tab <- read_tsv(gene.tab.file)

net.tab <- read_tsv(net.file)

net.pairs <- net.tab %>%
    rename(pair.type = 'Experimental System Type') %>%
    filter(pair.type == 'physical') %>%
    select(starts_with('Official')) %>%
    rename(gene.1 = 'Official Symbol Interactor A',
           gene.2 = 'Official Symbol Interactor B') %>%
    filter(gene.1 != gene.2) %>%
    filter(gene.1 %in% gene.tab$hgnc, gene.2 %in% gene.tab$hgnc) %>%
    unique()

net.pairs.relaxed <- net.tab %>%
    rename(pair.type = 'Experimental System Type') %>%
    filter(pair.type == 'physical') %>%
    select(starts_with('Official')) %>%
    rename(gene.1 = 'Official Symbol Interactor A',
           gene.2 = 'Official Symbol Interactor B') %>%
    filter(gene.1 != gene.2) %>%
    filter(gene.1 %in% gene.tab$hgnc | gene.2 %in% gene.tab$hgnc) %>%
    unique()

write_tsv(net.pairs, path = 'network/ppi-strict.pairs', col_names = TRUE)
write_tsv(net.pairs.relaxed, path = 'network/ppi.pairs', col_names = TRUE)
