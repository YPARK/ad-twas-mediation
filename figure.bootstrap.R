#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)
source('util.R')
source('mediation.R')
library(readr)
library(zqtl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(latex2exp)
source('figure.util.R')

gene.tab.file <- 'tables/bootstrap_gene.txt.gz'
gene.tab <- read.table(gene.tab.file, header = TRUE) %>% na.omit()
gene.tab.sig.file <- 'tables/bootstrap_gene_significant.txt.gz'
gene.tab.sig <- read.table(gene.tab.sig.file, header = TRUE) %>% na.omit()

## Figure 1. global effect sizes
draw.effect.block <- function(.df) {

    .df.sig <- .df %>% filter(p.val.marginal < 1e-6, p.val.direct < 1e-6, lodds > 2,
                              abs(max.gwas.z) > 3)

    plt <- ggplot() + theme_bw() +
        geom_linerange(data = .df %>% filter(lodds > 0),
                       aes(x = (tss+tes)/2e6, ymin = theta - 2 * theta.se, ymax = theta + 2* theta.se),
                       color = 'gray40', alpha = 0.5)

    plt <- plt +
        geom_point(data = .df, aes(x = (tss+tes)/2e6, y = theta, size = 1/(1+exp(-lodds))),
                   alpha = 0.8, color = 'gray40', show.legend = FALSE)

    plt <- plt +
        geom_point(data = .df.sig, aes(x = (tss+tes)/2e6, y = theta),
                   color = 'red', pch = 4, size = 1) +
                       facet_grid(.~chr, scales = 'free', space = 'free') +
                           scale_size_continuous(limits = c(0, 1), range = c(0, 1))
    
    plt <- plt +
        geom_text_repel(data = .df.sig %>% filter(theta > 0), aes(x = (tss+tes)/2e6, y = theta+2*theta.se, label = hgnc), nudge_y = .1, size = 2, segment.color = 'green', segment.alpha = 0.5)

    plt <- plt +
        geom_text_repel(data = .df.sig %>% filter(theta < 0), aes(x = (tss+tes)/2e6, y = theta-2*theta.se, label = hgnc), nudge_y = -.1, size = 2, segment.color = 'green', segment.alpha = 0.5) +
            ylab('Mediation effect') + xlab('Genomic location (Mb)')
    
    return(plt)
}

p.list <- list(draw.effect.block(gene.tab %>% filter(chr <= 3)),
               draw.effect.block(gene.tab %>% filter(chr >= 4, chr <= 7)),
               draw.effect.block(gene.tab %>% filter(chr >= 7, chr <= 11)),
               draw.effect.block(gene.tab %>% filter(chr >= 12, chr <= 16)),
               draw.effect.block(gene.tab %>% filter(chr >= 17)))

out.dir <- 'figures/genes/'
out.files <- out.dir %&&% 'global_' %&&% 1:length(p.list) %&&% '.pdf'

sapply(1:length(p.list), function(j) ggsave(plot = p.list[[j]], filename = out.files[j],
                                            useDingbats = FALSE, width = 10, height = 3))

## Figure 2. bootstrap results chr by chr

draw.boot.chr <- function(.chr) {

    .df.marg.file <- 'bootstrap/marginal_IGAP_rosmap_eqtl_hs-lm_' %&&% .chr %&&% '.mediation.gz'
    .df.direct.file <- 'bootstrap/direct_IGAP_rosmap_eqtl_hs-lm_' %&&% .chr %&&% '.mediation.gz'

    .cols <- c('ensg', 'theta', 'theta.se', 'lodds',
               'emp.p', 'fd', 'nboot',
               'lodds.mean', 'lodds.se', 'cauchy.location', 'cauchy.scale',
               't.m', 't.s', 'cauchy.p', 't.p', 'norm.p',
               'hgnc', 'tss', 'tes', 'strand',
               'max.gwas.theta', 'max.gwas.z',
               'chr', 'ld.1', 'ld.2', 'n.snps')
    
    .df.marg <- read.table(.df.marg.file, col.names = .cols)
    .df.direct <- read.table(.df.direct.file, col.names = .cols)

    .sigmoid <- function(x) 1/(1+exp(-x))

    .aes <- aes(x = (tss+tes)/2e6,
                ymin = .sigmoid(lodds.mean - 4 * lodds.se),
                ymax = .sigmoid(lodds.mean + 4 * lodds.se))

    plt <- ggplot() + theme_bw() + geom_ribbon(data = .df.direct, .aes, fill = '#99FF99', alpha = 0.3)
        
    plt <- plt + geom_ribbon(data = .df.marg, .aes, fill = '#FF9999', alpha = 0.3)

    plt <- plt + 
        geom_point(data = .df.marg, aes(x = (tss+tes)/2e6, y = .sigmoid(lodds)),
                   color = 'gray40', show.legend = FALSE, size = .5)

    .df.sig <- .df.marg %>% filter(lodds > 2, abs(max.gwas.z) > 3)

    if(nrow(.df.sig) > 0) {

        plt <- plt +
            geom_point(data = .df.sig, aes(x = (tss+tes)/2e6, y = .sigmoid(lodds)),
                       color = 'red', pch = 4, size = 1)
        
        plt <- plt +
            geom_text_repel(data = .df.sig, aes(x = (tss+tes)/2e6, y = .sigmoid(lodds), label = hgnc), nudge_y = .1, size = 2, segment.color = 'green', segment.alpha = 0.5)
        
    }    

    plt <- plt + xlab('Genomic location (Mb)') + ylab('Mediation PIP')

    file.name <- 'figures/genes/bootstrap_chr' %&&% .chr %&&% '.pdf'
    n <- nrow(.df.marg)

    ggsave(filename = file.name, plot = plt, useDingbats = FALSE, width = n/300 + .5, height = 3)
}

for(chr in 1:22) {
    draw.boot.chr(chr)
}
