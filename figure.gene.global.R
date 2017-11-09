#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
source('mediation.R')
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(latex2exp)
source('figure.util.R')
library(scales)

draw.plots <- function(qtl.data) {

    gene.tab.file <- 'tables/genes/' %&&% qtl.data %&&% '/total_genes.txt.gz'
    gene.tab.sig.file <- 'tables/genes/' %&&% qtl.data %&&% '/significant_genes.txt.gz'

    dir.create('figures/genes/', recursive = TRUE, showWarnings = FALSE)

    gene.tab <- read_tsv(gene.tab.file, col_names = TRUE) %>% na.omit() %>%
        group_by(hgnc) %>% slice(which.max(lodds))

    gene.tab.sig <- read_tsv(gene.tab.sig.file, col_names = TRUE) %>% na.omit() %>%
        group_by(hgnc) %>% slice(which.max(lodds))

    logit.trans <- trans_new(".logit", function(x) log(x) - log(1-x), function(y) 1/(1+exp(-y)))

    ## GWAS vs PIP vs PVE
    df <- gene.tab %>%
        mutate(ld.gwas.p = pmin(-log10(ld.gwas.p), 20)) %>%
            mutate(qtl.p = pmin(-log10(2*pnorm(abs(best.qtl.z), lower.tail = FALSE)), 20)) %>%
                mutate(pip = 1/(1+exp(-lodds))) %>%
                    mutate(gene.loc = ifelse(strand == '+', tss, tes)) %>%
                        as.data.frame()
    
    df.sig <- df %>% filter(hgnc %in% gene.tab.sig$hgnc)

    df.strict <- df.sig %>% filter(PVE >= 5e-2, ld.gwas.p > 4)

    scale.grad <- scale_fill_gradientn(colors = c('#FFFFFF', '#FFAA00'), trans = 'log10')

    gwas.pip.plot <- function(transpose = FALSE) {

        .pip.cutoff <- min(df.sig$pip)
        
        if(transpose) {
            aes.gwas <- aes(y = ld.gwas.p, x = pip, label = hgnc)
        } else {
            aes.gwas <- aes(x = ld.gwas.p, y = pip, label = hgnc)
        }

        ret <- gg.plot(df, aes.gwas) + geom_hex(bins = 40, color = 'gray') + scale.grad 

        if(transpose) {
            ret <- ret +
                geom_vline(xintercept = .pip.cutoff, color = 'red', lty = 2, size = 1) +
                    ylab('best -log10 GWAS p-value in LD') + xlab('posterior probability of mediation')

            ret <- ret + scale_x_continuous(breaks = c(0.01, 0.5, 0.99), trans = logit.trans)
        } else {
            ret <- ret +
                geom_hline(yintercept = .pip.cutoff, color = 'red', lty = 2, size = 1) +
                    xlab('best -log10 GWAS p-value in LD') + ylab('posterior probability of mediation')

            ret <- ret + scale_y_continuous(breaks = c(0.01, 0.5, 0.99), trans = logit.trans)
        }

        ret <- ret + geom_point(data = df.strict, color = 'black', fill = 'red', pch = 21) +
            geom_text_repel(data = df.strict, size = 4, segment.color = 'green', segment.alpha = 0.5)

        return(ret)
    }

    qtl.pip.plot <- function(transpose = FALSE) {

        .pip.cutoff <- min(df.sig$pip)
        
        if(transpose) {
            aes.qtl <- aes(y = qtl.p, x = pip, label = hgnc)
        } else {
            aes.qtl <- aes(x = qtl.p, y = pip, label = hgnc)
        }

        ret <- gg.plot(df, aes.qtl) + geom_hex(bins = 40, color = 'gray') + scale.grad 

        if(transpose) {
            ret <- ret +
                geom_vline(xintercept = .pip.cutoff, color = 'red', lty = 2, size = 1) +
                    ylab('best -log10 QTL p-value') + xlab('posterior probability of mediation')

            ret <- ret + scale_x_continuous(breaks = c(0.01, 0.5, 0.99), trans = logit.trans)
        } else {
            ret <- ret +
                geom_hline(yintercept = .pip.cutoff, color = 'red', lty = 2, size = 1) +
                    xlab('best -log10 QTL p-value') + ylab('posterior probability of mediation')

            ret <- ret + scale_y_continuous(breaks = c(0.01, 0.5, 0.99), trans = logit.trans)
        }

        ret <- ret + geom_point(data = df.strict, color = 'black', fill = 'red', pch = 21) +
            geom_text_repel(data = df.strict, size = 4, segment.color = 'green', segment.alpha = 0.5)

        return(ret)
    }

    pve.pip.plot <- function() {

        ret <- gg.plot(df, aes(x = pip, y = PVE)) +
            geom_hex(bins = 40, color = 'gray') + scale.grad

        ret <- ret + geom_point(data = df.strict, aes(x = pip, y = PVE),
                                pch = 21, fill = 'red', color = 'black')
        ret <- ret +
            geom_text_repel(data = df.strict, aes(x = pip, y = PVE, label = hgnc),
                            size = 4, segment.color = 'green',
                            segment.alpha = 0.5, direction = 'y')

        ret <- ret + xlab('posterior probability of mediation') +
            ylab('proportion of variance explained by mediation')

        ret <- ret + scale_x_continuous(limits = c(1e-3, 10), breaks = c(1e-2, .5, .99), trans = logit.trans) +
            scale_y_continuous(limits = c(1e-10, 10), breaks = c(1e-2, .5, .99), trans = logit.trans)

        ret <- ret + theme(legend.position = c(1, 0), legend.justification = c(1, 0))

        return(ret)
    }

    gwas.pve.plot <- function(transpose = FALSE) {

        if(transpose) {
            aes.gwas <- aes(y = ld.gwas.p, x = PVE, label = hgnc)
        } else {
            aes.gwas <- aes(x = ld.gwas.p, y = PVE, label = hgnc)
        }

        ret <- gg.plot(df, aes.gwas) +
            geom_hex(bins = 40, color = 'gray') + scale.grad

        ret <- ret + geom_point(data = df.strict, pch = 21, color = 'black', fill = 'red') +
            geom_text_repel(data = df.strict,
                            size = 4, segment.color = 'green',
                            segment.alpha = 0.5)

        if(transpose) {
            ret <- ret + xlab('proportion of variance explained by mediation') +
                ylab('best -log10 GWAS p-value in LD')

            ret <- ret + scale_x_continuous(limits = c(1e-8, 10),
                                            breaks = c(1e-4, 1e-2, 1e-1, 0.99), trans = logit.trans)
        } else {
            ret <- ret + ylab('proportion of variance explained by mediation') +
                xlab('best -log10 GWAS p-value in LD')

            ret <- ret + scale_y_continuous(limits = c(1e-8, 10),
                                            breaks = c(1e-4, 1e-2, 1e-1, 0.99), trans = logit.trans)
        }

        return(ret)
    }

    gwas.qtl.plot <- function(transpose = FALSE) {

        if(transpose) {
            aes.gwas <- aes(y = gwas.p, x = qtl.p, label = hgnc)
        } else {
            aes.gwas <- aes(x = gwas.p, y = qtl.p, label = hgnc)
        }

        ret <- gg.plot(df, aes.gwas) +
            geom_hex(bins = 40, color = 'gray') + scale.grad

        ret <- ret + geom_point(data = df.strict, pch = 21, color = 'black', fill = 'red') +
            geom_text_repel(data = df.strict,
                            size = 4, segment.color = 'green',
                            segment.alpha = 0.5)

        if(transpose) {
            ret <- ret + xlab('best -log10 QTL p-value') +
                ylab('best -log10 GWAS p-value')
        } else {
            ret <- ret + ylab('best -log10 QTL p-value') +
                xlab('best -log10 GWAS p-value')
        }

        return(ret)
    }

    qtl.pve.plot <- function(transpose = FALSE) {

        if(transpose) {
            aes.qtl <- aes(y = qtl.p, x = PVE, label = hgnc)
        } else {
            aes.qtl <- aes(x = qtl.p, y = PVE, label = hgnc)
        }

        ret <- gg.plot(df, aes.qtl) +
            geom_hex(bins = 40, color = 'gray') + scale.grad

        ret <- ret + geom_point(data = df.strict, pch = 21, color = 'black', fill = 'red') +
            geom_text_repel(data = df.strict,
                            size = 4, segment.color = 'green',
                            segment.alpha = 0.5)

        if(transpose) {
            ret <- ret + xlab('proportion of variance explained by mediation') +
                ylab('best -log10 QTL p-value')

            ret <- ret + scale_x_continuous(limits = c(1e-8, 10),
                                            breaks = c(1e-4, 1e-2, 1e-1, 0.99), trans = logit.trans)
        } else {
            ret <- ret + ylab('proportion of variance explained by mediation') +
                xlab('best -log10 QTL p-value')

            ret <- ret + scale_y_continuous(limits = c(1e-8, 10),
                                            breaks = c(1e-4, 1e-2, 1e-1, 0.99), trans = logit.trans)
        }

        return(ret)
    }


    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_gwas_pip.pdf',
           plot = gwas.pip.plot(), width = 6, height = 4, units = 'in')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_qtl_pip.pdf',
           plot = qtl.pip.plot(), width = 6, height = 4, units = 'in')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_pip_gwas.pdf',
           plot = gwas.pip.plot(TRUE), width = 6, height = 4, units = 'in')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_pve_pip.pdf',
           plot = pve.pip.plot(), width = 6, height = 4, units = 'in')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_gwas_pve.pdf',
           plot = gwas.pve.plot(), width = 6, height = 4, units = 'in')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_gwas_qtl.pdf',
           plot = gwas.qtl.plot(), width = 6, height = 4, units = 'in')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_qtl_pve.pdf',
           plot = qtl.pve.plot(), width = 6, height = 4, units = 'in')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_pve_gwas.pdf',
           plot = gwas.pve.plot(TRUE), width = 6, height = 4, units = 'in')


    ################################################################
    ## Global effect sizes
    ## color genes by GWAS significance

    .genome.mb <- function() {
        function(x) round(x/1e6)
    }

    .thm <- theme(legend.direction = 'horizontal', legend.position = 'bottom',
                  panel.spacing.x = unit(0.05, 'lines'), panel.grid.minor.x = element_blank(),
                  panel.border = element_blank(), strip.text = element_text(size = 8),
                  axis.text = element_text(size = 8),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.title = element_text(size = 12))
    
    ret <-
        gg.plot(df %>% arrange(gwas.p),
                aes(x = gene.loc, y = theta, color = pmin(ld.gwas.p, 8), size = PVE, label = hgnc))

    ret <- ret +
        geom_linerange(data = df.strict,
                       aes(ymin = theta - 2*theta.se, ymax = theta + 2*theta.se), size = 1.5)

    ret <- ret + geom_point() +
        scale_size_continuous('proportion of variance explained', range = c(0, 5), breaks = c(1e-2, 1e-1, .99))

    ret <- ret +
        scale_color_gradientn('-log10 GWAS p-value', colors = c('gray80', 'gray40', 'purple', 'purple'))

    ret <- ret +
        facet_grid(.~chr, scales = 'free', space = 'free') +
            scale_x_continuous(labels = .genome.mb(), breaks = c(0, 5e7, 1e8, 1.5e8, 2e8)) + .thm

    ret <- ret + xlab('genomic location (mb)') + ylab('mediation effect')

    df.strict.neg <- df.strict %>% filter(theta < 0)

    df.strict.pos <- df.strict %>% filter(theta > 0)

    ret <- ret + geom_text_repel(data = df.strict.neg, 
                                 aes(y = theta), color = 'blue', segment.alpha = 0.5,
                                 nudge_y = -.2, size = 3, segment.size = .5, segment.color = 'green')

    ret <- ret + geom_point(data = df.strict.neg, pch = 22, size = 3, color = 'green')

    ret <- ret + geom_text_repel(data = df.strict.pos,
                                 aes(y = theta), color = 'red', segment.alpha = 0.5,
                                 nudge_y = .2, size = 3, segment.size = .5, segment.color = 'orange')

    ret <- ret + geom_point(data = df.strict.pos, pch = 22, size = 3, color = 'orange')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_mediation_effect.pdf',
           plot = ret, width = 8, height = 4, units = 'in',
           useDingbats = FALSE, limitsize = FALSE)
}

qtl.data.vec <- c('full-' %&&% c(5, 3, 0), 'hic-' %&&% c(5, 3, 0))

for(qtl.data in qtl.data.vec) {
    draw.plots(qtl.data)
}
