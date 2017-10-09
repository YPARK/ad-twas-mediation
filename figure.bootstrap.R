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
draw.effect.chr <- function(.chr) {

    .df <- gene.tab %>% filter(chr == .chr)
    .df.sig <- gene.tab.sig %>% filter(chr == .chr, lodds > 2)

    plt <- ggplot() + theme_bw() +
        geom_linerange(data = .df %>% filter(lodds > 0),
                       aes(x = (tss+tes)/2, ymin = theta - 2 * theta.se, ymax = theta + 2* theta.se),
                       color = 'gray40', alpha = 0.5)

    plt <- plt +
        geom_point(data = .df, aes(x = (tss+tes)/2, y = theta, size = 1/(1+exp(-lodds))),
                   alpha = 0.8, color = 'gray40', show.legend = FALSE)

    plt <- plt +
        geom_point(data = .df.sig, aes(x = (tss+tes)/2, y = theta),
                   color = 'red', pch = 4, size = 1) +
                       scale_size_continuous(limits = c(0, 1), range = c(0, 1))
    
    plt <- plt +
        geom_text_repel(data = .df.sig %>% filter(theta > 0), aes(x = (tss+tes)/2, y = theta+2*theta.se, label = hgnc),
                        nudge_y = .1, size = 4, segment.color = 'green', segment.alpha = 0.5)

    plt <- plt +
        geom_text_repel(data = .df.sig %>% filter(theta < 0), aes(x = (tss+tes)/2, y = theta-2*theta.se, label = hgnc),
                        nudge_y = -.1, size = 4, segment.color = 'green', segment.alpha = 0.5) +
            ylab('Mediation effect') + xlab('genomic location')
    
    file.name <- 'figures/genes/effect_chr' %&&% .chr %&&% '.pdf'

    blk.width <- 1e6
    genome.range <- range(c(.df$tss, .df$tes))
    genome.length <- genome.range[2] - genome.range[1]
    genome.w <- genome.length / blk.width * 0.05 + 1

    .gwas.mb <- function() {
        function(x) format(x/1e6, big.mark=',')
    }

    genome.x.scale <- scale_x_continuous(limits = genome.range + c(-1000, 1000), expand = c(0, 0),
                                         labels = .gwas.mb())

    plt <- plt + genome.x.scale

    ggsave(filename = file.name, plot = plt, useDingbats = FALSE, limitsize = FALSE, width = genome.w + .5, height = 3)
}

for(chr in 1:22) {
    draw.effect.chr(chr)
}

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
    
    .df.marg <- read_tsv(.df.marg.file, col_names = .cols) %>%
        filter(ensg %in% gene.tab$ensg)

    .df.direct <- read_tsv(.df.direct.file, col_names = .cols) %>%
        filter(ensg %in% gene.tab$ensg)

    .sigmoid <- function(x) 1/(1+exp(-x))

    .aes <- aes(x = (tss+tes)/2,
                ymin = .sigmoid(lodds.mean - 4 * lodds.se),
                ymax = .sigmoid(lodds.mean + 4 * lodds.se),
                y = .sigmoid(lodds.mean + 4 * lodds.se))

    plt <- ggplot() + theme_bw() +
        geom_ribbon(data = .df.direct, .aes, fill = '#FFFF99', alpha = 0.3)

    plt <- plt +
        geom_line(data = .df.direct, .aes, color = '#99FF99', size = .5)
    
##    plt <- plt + geom_ribbon(data = .df.marg, .aes, fill = '#FF9999', alpha = 0.3)

    .df <- gene.tab %>% filter(chr == .chr)
    .df.sig <- gene.tab.sig %>% filter(chr == .chr, lodds > 2)

    plt <- plt +
        geom_line(data = .df.marg, .aes, color = '#FF9999', size = .5)

    plt <- plt + 
        geom_point(data = .df, aes(x = (tss+tes)/2, y = .sigmoid(lodds)),
                   color = 'gray40', show.legend = FALSE, size = .5)

    if(nrow(.df.sig) > 0) {

        plt <- plt +
            geom_point(data = .df.sig, aes(x = (tss+tes)/2, y = .sigmoid(lodds)),
                       color = 'red', pch = 4, size = 1)
        
        plt <- plt +
            geom_text_repel(data = .df.sig, aes(x = (tss+tes)/2, y = .sigmoid(lodds), label = hgnc),
                            size = 4, segment.color = 'green', segment.alpha = 0.5)        
    }    

    plt <- plt + xlab('Genomic location (Mb)') + ylab('Mediation PIP')

    file.name <- 'figures/genes/bootstrap_chr' %&&% .chr %&&% '.pdf'

    blk.width <- 1e6
    genome.range <- range(c(.df$tss, .df$tes))
    genome.length <- genome.range[2] - genome.range[1]
    genome.w <- genome.length / blk.width * 0.05 + 1

    .gwas.mb <- function() {
        function(x) format(x/1e6, big.mark=',')
    }

    genome.x.scale <- scale_x_continuous(limits = genome.range + c(-1000, 1000), expand = c(0, 0),
                                         labels = .gwas.mb())

    plt <- plt + genome.x.scale

    ggsave(filename = file.name, plot = plt, useDingbats = FALSE, limitsize = FALSE, width = genome.w + .5, height = 3)
}

for(chr in 1:22) {
    draw.boot.chr(chr)
}
