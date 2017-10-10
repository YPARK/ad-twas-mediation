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

gene.tab.strict.file <- 'tables/bootstrap_gene_strict.txt.gz'
gene.tab.strict <- read.table(gene.tab.strict.file, header = TRUE) %>% na.omit()

gene.tab <- gene.tab %>% group_by(hgnc) %>%
    slice(which.max(lodds))

gene.tab.sig <- gene.tab.sig %>% group_by(hgnc) %>%
    slice(which.max(lodds))

gene.tab.strict <- gene.tab.strict %>% group_by(hgnc) %>%
    slice(which.max(lodds))

## Figure 0. just combine all
.genome.mb <- function() {
    function(x) format(x/1e6, big.mark=',')
}

get.genome.scale <- function(df, blk.width = 1e6) {

    genome.range <- range(c(df$tss, df$tes))
    genome.length <- genome.range[2] - genome.range[1]
    genome.w <- genome.length / blk.width * 0.05 + 1

    genome.x.scale <- scale_x_continuous(limits = genome.range + c(-1000, 1000), expand = c(0, 0),
                                         labels = .genome.mb())

    return(genome.x.scale)
}

x.scale <- get.genome.scale(gene.tab, blk.width = 1e6)

.thm <- theme(plot.background = element_blank(),
              panel.background = element_blank(),
              strip.background = element_blank())
              
plt <- ggplot(gene.tab) + theme_classic() + .thm +
    theme(axis.text.x = element_text(size = 6, angle = 30)) +
    xlab('genomic location (mb)') +
    ylab('mediation effect') +
    scale_x_continuous(labels = .genome.mb())

.aes.lr <- aes(x = (tss + tes)/2,
               ymin = theta - 2 * theta.se,
               ymax = theta + 2 * theta.se)

plt <- plt + geom_linerange(data = gene.tab.strict, .aes.lr, alpha = 0.5)

.aes <- aes(x = (tss + tes)/2, y = theta, color = (chr %% 2 == 0),
            alpha = 1/(1+exp(-lodds)), size = 1/(1+exp(-lodds)))

plt <- plt + geom_point(.aes, show.legend = FALSE) +
    scale_color_manual(values = c('#3333FF', 'gray40'))
plt <- plt +
    facet_grid(.~chr, space = 'free', scales = 'free')

plt <- plt + scale_size_continuous(limits = c(0, 1), range = c(0, 2))

.aes.lab <- aes(x = (tss + tes)/2, y = theta, label = hgnc)

plt <- plt + geom_text_repel(data = gene.tab.strict, .aes.lab, direction = 'y',
                      size = 4, segment.color = 'green', color = 'red')

ggsave(filename = 'figures/mediation_gene_effect.pdf', plot = plt, width = 16, height = 8,
       useDingbats = FALSE, limitsize = FALSE)

################################################################
## best GWAS vs log-odds
.df <- gene.tab %>%
    mutate(gwas.p = pmin(l10.p.two(abs(max.gwas.z)),20)) %>%
    mutate(twas.p = pmin(l10.p.two(abs(nwas.z)),20))

.df.strict <- .df %>% filter(ensg %in% gene.tab.strict$ensg)

.lodds.cutoff <- min(.df.strict$lodds)

.aes.gwas.lodds <- aes(x = gwas.p, y = lodds)

.grad <- scale_fill_gradientn(colors = c('gray', 'black'), trans = 'log10')

p1 <- ggplot(.df, .aes.gwas.lodds) + theme_bw() + .thm +
    geom_hex(bins = 50) + xlab('-log10 GWAS') + ylab('logit PIP') +
    geom_hline(yintercept = .lodds.cutoff, color = 'green', lty = 2) +
    .grad 

p1 <- p1 + geom_point(data = .df.strict, aes(x = gwas.p, y = lodds), pch = 3, color = 'red') +
    geom_text_repel(data = .df.strict, aes(x = gwas.p, y = lodds, label = hgnc),
                    size = 4, direction = 'x', segment.color = 'green', segment.alpha = 0.5)

.aes.gwas.twas <- aes(x = gwas.p, y = twas.p)

.df.twas <- .df %>% filter(lodds < -2, twas.p > -log10(2.5e-6))

p2 <-
    ggplot(.df, aes(x = twas.p, y = lodds)) + theme_bw() + .thm +
    geom_hex(bins = 50) + .grad +
    geom_hline(yintercept = .lodds.cutoff, color = 'green', lty = 2) +
    geom_vline(xintercept = -log10(2.5e-6), color = 'blue', lty = 2)

p2 <- p2 +
    geom_point(data = .df.strict, aes(x = twas.p, y = lodds), pch = 3, color = 'red') +
    geom_point(data = .df.twas, aes(x = twas.p, y = lodds), pch = 1, color = 'blue')

p2 <- p2 +
    geom_text_repel(data = .df.strict, aes(x = twas.p, y = lodds, label = hgnc),
                    size = 4, direction = 'x', segment.color = 'green', segment.alpha = 0.5) +
    geom_text_repel(data = .df.twas, aes(x = twas.p, y = lodds, label = hgnc),
                    size = 2, nudge_x = 1, nudge_y = -1, segment.color = 'gray80',
                    color = 'blue', segment.alpha = 0.5)

p2 <- p2 + xlab('-log10 TWAS') + ylab('logit PIP')

pdf(file = 'figures/mediation_gene_twas.pdf', width = 16, height = 8, useDingbats = FALSE)
print(grid.hcat(list(p1, p2)))
dev.off()



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

    genome.x.scale <- get.genome.scale(.df)

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

    genome.x.scale <- get.genome.scale(.df)

    plt <- plt + genome.x.scale

    ggsave(filename = file.name, plot = plt, useDingbats = FALSE, limitsize = FALSE, width = genome.w + .5, height = 3)
}

for(chr in 1:22) {
    draw.boot.chr(chr)
}
