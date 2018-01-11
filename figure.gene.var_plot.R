#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
source('figure.util.R')
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)

cammel.tab.file <- 'tables/genes/full-0/significant_genes.txt.gz'
tot.tab.file <- 'tables/genes/full-0/total_genes.txt.gz'
cammel.tab <- read_tsv(cammel.tab.file)
tot.tab <- read_tsv(tot.tab.file)

dir.create('figures/genes/', recursive = TRUE, showWarnings = FALSE)

## accumulate total amount of genetic variance
ld.tab <- tot.tab %>% dplyr::group_by(chr, ld.lb, ld.ub) %>%
    slice(which.max(v.med)) %>%
    arrange(chr, ld.lb) %>%
    dplyr::select(chr, ld.lb, ld.ub, ld.gwas.loc, ld.gwas.p,
                  hgnc, v.med, v.med.tot, v.dir, v.resid) %>%
    arrange(desc(v.med.tot)) %>%
    as.data.frame()

ld.cammel.tab <- ld.tab %>%
    right_join(cammel.tab %>% dplyr::select(chr, ld.lb, ld.ub) %>% unique()) %>%
    arrange(desc(v.med.tot)) %>%
    na.omit()

chr.tab <- ld.tab %>%
    group_by(chr) %>%
    summarize(len = sum(ld.ub - ld.lb),
              v.med.tot = sum(v.med.tot),
              v.dir = sum(v.dir),
              v.resid = sum(v.resid))
    
sum.tab <- ld.tab %>%
    summarize(v.med = sum(v.med.tot), v.dir = sum(v.dir), v.resid = sum(v.resid))

v.tot <- sum(sum.tab)

## summary for each chromosome

p <- gg.plot(chr.tab, aes(x = len / 1e6, y = v.med.tot / v.tot)) +
    geom_point() +
    geom_text_repel(aes(label = 'chr' %&&% chr %&&% '_' %&&% signif(v.med.tot/v.tot, 3))) +
    xlab('included genomic length (Mb)') +
    ylab('prop. var. explained by gene expr.')

ggsave(filename = 'figures/genes/pve_chr_scatter.pdf', plot = p, width = 3, height = 3, unit = 'in')

p <- gg.plot(chr.tab %>% mutate(cum.v.med = cumsum(v.med.tot)/v.tot)) +
    theme(panel.grid = element_blank()) +
    geom_bar(aes(x = chr, y = v.med.tot/v.tot), stat = 'identity', color = 'orange', fill = 'orange') +
    geom_text_repel(aes(x = chr, y = cum.v.med, label = signif(cum.v.med, 2)), size = 3) +
    geom_line(aes(x = chr, y = cum.v.med)) +
    geom_point(aes(x = chr, y = cum.v.med), pch = 3) +
    scale_x_continuous(breaks = 1:22) +
    ylab('cumulative prop. var. explained by gene expr.') +
    xlab('chromosomes')

ggsave(filename = 'figures/genes/pve_chr_cum.pdf', plot = p, width = 3, height = 3, unit = 'in')

## summary for LD blocks
.df <- data.frame(v.med = ld.tab$v.med.tot / v.tot,
                  v.med.cum = cumsum(ld.tab$v.med.tot) / v.tot,
                  gene = ld.tab$hgnc) %>%
    mutate(idx = 1:n())
                  
.df.cammel <- data.frame(v.med = ld.cammel.tab$v.med.tot / v.tot,
                         v.med.cum = cumsum(ld.cammel.tab$v.med.tot) / v.tot,
                         gene = ld.cammel.tab$hgnc) %>%
    mutate(idx = 1:n())

p0 <- gg.plot(.df) +
    geom_bar(aes(x = idx, y = v.med), stat = 'identity', color = 'orange', fill = 'orange') +
    ylab('prop. var. explained by gene expr.') +
    theme(axis.title.x = element_blank()) +
    scale_x_log10()

ggsave(filename = 'figures/genes/pve_ld.pdf', plot = p0, width = 3, height = 3) 

p1 <- gg.plot(.df) +
    geom_line(aes(x = idx, y = v.med.cum)) +
    geom_line(data = .df.cammel, aes(x = idx, y = v.med.cum), color = 'red') +
    ylab('cumulative') +
    xlab('independent genomic regions (sorted)') +
    scale_x_log10(breaks = c(1, 10, 100, nrow(.df.cammel)))

ggsave(filename = 'figures/genes/pve_ld_cum.pdf', plot = p1, width = 3, height = 3) 

## summary for genes
gene.tab <- tot.tab %>% group_by(hgnc) %>% slice(which.max(PVE)) %>%
    as.data.frame() %>%
    arrange(desc(v.med)) %>%
    mutate(v.med.cum = cumsum(v.med) / v.tot) %>%
    mutate(v.med = v.med / v.tot) %>%
    dplyr::select(chr, tss, tes, strand, hgnc, v.med, v.med.cum) %>%
    mutate(idx = 1:n())

cammel.gene.tab <- gene.tab %>% dplyr::filter(hgnc %in% cammel.tab$hgnc) %>%
    mutate(idx = 1:n())

v.rng <- range(cammel.gene.tab$v.med.cum)

p <- gg.plot(gene.tab, aes(x = idx, y = v.med.cum)) +
    geom_line(color = 'gray20') +
    geom_line(data = cammel.gene.tab, color = 'red') +
    scale_x_log10(breaks = c(1, 10, 100, nrow(cammel.gene.tab))) +
    scale_y_continuous(breaks = round(v.rng, 4)) +
    ylab('cumulative prop. var. explained') +
    xlab('ranked genes')

ggsave(filename = 'figures/genes/pve_gene_cum.pdf', plot = p, width = 3, height = 3, unit = 'in')

