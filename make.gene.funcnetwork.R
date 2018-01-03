#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
source('figure.util.R')
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

cammel.tab.file <- 'tables/genes/full-0/significant_genes.txt.gz'
net.file <- 'network/FIsInGene_022717_with_annotations.txt'

cammel.tab <- read_tsv(cammel.tab.file)
net.tab <- read_tsv(net.file)

temp <- net.tab %>%
    filter(Gene1 %in% cammel.tab$hgnc | Gene2 %in% cammel.tab$hgnc) %>%
    filter(Gene1 != Gene2)

temp <- union(temp$Gene1, temp$Gene2)

net.pairs <- net.tab %>% 
    filter(Gene1 %in% temp, Gene2 %in% temp)

net.dup <- net.pairs %>%
    filter(Direction == '-') %>%
    mutate(temp = Gene2, Gene2 = Gene1, Gene1 = temp) %>%
    select(-temp) %>%
    unique()

write_tsv(rbind(net.pairs, net.dup), path = 'network/functional.pairs', col_names = FALSE)

write_tsv(net.pairs %>% select(Gene1, Gene2) %>% unique(),
          path = 'network/functional-slim.pairs', col_names = FALSE)

