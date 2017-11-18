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

draw.plots <- function(qtl.data, p.cutoff = 1e-4) {

    gene.tab.file <- 'tables/genes/' %&&% qtl.data %&&% '/total_genes.txt.gz'

    dir.create('figures/genes/', recursive = TRUE, showWarnings = FALSE)

    gene.tab <- read_tsv(gene.tab.file, col_names = TRUE) %>% na.omit() %>%
        group_by(hgnc) %>% slice(which.max(lodds))

    logit.trans <- trans_new(".logit", function(x) log(x) - log(1-x), function(y) 1/(1+exp(-y)))

    ## GWAS vs PIP / PVE
    df <- gene.tab %>%
        mutate(ld.gwas.p = pmin(-log10(ld.gwas.p), 20)) %>%
            mutate(qtl.p = pmin(-log10(2*pnorm(abs(best.qtl.z), lower.tail = FALSE)), 20)) %>%
                mutate(pip = 1/(1+exp(-lodds))) %>%
                    mutate(gene.loc = ifelse(strand == '+', tss, tes)) %>%
                        as.data.frame()
    
    df.strict <- df %>% filter(PVE >= 5e-2, qval < 1e-4, ld.gwas.p > 4)

    scale.grad <- scale_fill_gradientn(colors = c('white', 'gray20', 'gray20'), trans = 'sqrt')
    
    df.pve.hist <- df %>% mutate(pve.logit = round((log(PVE) - log(1 - PVE)) * 2) / 2 ) %>%
        group_by(pve.logit) %>% summarize(count = n()) %>%
            mutate(PVE = 1/(1+exp(-pve.logit )))

    df.pve.hist.sig <- df %>% filter(qval < 1e-4) %>%
        mutate(pve.logit = round((log(PVE) - log(1 - PVE)) * 2) / 2 ) %>%
            group_by(pve.logit) %>% summarize(count = n()) %>%
                mutate(PVE = 1/(1+exp(-pve.logit ))) %>%
                    arrange(desc(PVE))
    
    avg.pve.sig <- df %>% filter(qval < 1e-4) %>% summarize(pve = mean(PVE))
    avg.pve.strict <- df %>% filter(qval < 1e-4, ld.gwas.p > 4) %>% summarize(pve = mean(PVE))

    df.pve.hist.strict <- df %>% filter(qval < 1e-4, ld.gwas.p > 4) %>%
        mutate(pve.logit = round((log(PVE) - log(1 - PVE)) * 2) / 2 ) %>%
            group_by(pve.logit) %>% summarize(count = n()) %>%
                mutate(PVE = 1/(1+exp(-pve.logit ))) %>%
                    arrange(desc(PVE))

    gwas.pve.plot <- function() {

        brk <- c(1e-8, 1e-4, 1e-3, 1e-2, 1e-1, 0.99, round(avg.pve.sig$pve, 2), 5e-2)                 

        pve.logit.scale <- scale_y_continuous(limits = c(1e-8, 10),
                                              breaks = brk,
                                              trans = logit.trans)
        
        p1 <- gg.plot(df.pve.hist, aes(y = PVE, yend = PVE, x = 0, xend = count))

        ## p1 <- p1 +
        ##     geom_rect(aes(xmin = 0, xmax = max(df.pve.hist$count), ymin = 5e-2, ymax = .99),
        ##               fill = '#FFFFDD')
        ## p1 <- p1 +
        ##     geom_rect(aes(xmin = 0, xmax = max(df.pve.hist$count), ymin = 1e-8, ymax = 5e-2),
        ##           fill = '#F0F0F0')

        p1 <- p1 +
            geom_hline(yintercept = 5e-2, color = 'red', lty = 'solid', size = .5)

        p1 <- p1 + geom_hline(yintercept = avg.pve.sig$pve, color = 'green')

        p1 <- p1 +
            geom_segment(size = 3, color = 'gray40') +
                scale_x_reverse() + xlab('count') +
                    pve.logit.scale +
                        ylab('proportion of variance explained by mediation')
        
        pve.aes <- aes(y = PVE, yend = PVE, x = 0, xend = count)

        p1 <- p1 + geom_segment(data = df.pve.hist.sig, pve.aes, size = 1.5, color = 'green')
        p1 <- p1 + geom_segment(data = df.pve.hist.strict, pve.aes, size = 1, color = 'blue') +
            geom_path(data = df.pve.hist.strict, aes(y = PVE, x = count),
                      size = 1, color = 'blue')

        p2 <- gg.plot(df, aes(x = ld.gwas.p, y = PVE, label = hgnc))

        p2 <- p2 + geom_hline(yintercept = 5e-2, color = 'red', lty = 'solid', size = .5)
            
        ## p2 <- p2 + geom_vline(xintercept = 4, color = 'green', lty = 'solid', size = .5) +
        ## p2 <- p2 +geom_vline(xintercept = -log10(5e-8), color = 'green', lty = 'dashed', size = .5)
        ## p2 <- p2 + geom_rect(aes(xmin = 0, xmax = 4, ymin = 1e-8, ymax = .99),
        ##                      fill = '#F0F0F0')
        ## p2 <- p2 + geom_rect(aes(xmin = 4, xmax = -log10(5e-8), ymin = 1e-8, ymax = .99),
        ##                      fill = '#FFFFDD')
        ## p2 <- p2 + geom_rect(aes(xmin = -log10(5e-8), xmax = 20, ymin = 1e-8, ymax = .99),
        ##                      fill = '#FFFFCC')

        p2 <- p2 +
            geom_hex(bins = 50, color = 'gray') + scale.grad +
                xlab('best -log10 GWAS p-value in LD') +
                    pve.logit.scale +
                        theme(axis.title.y = element_blank(),
                              axis.text.y = element_blank())
        
        p2 <- p2 +
            geom_point(data = df.strict, pch = 21, color = 'black', fill = 'red') +
                geom_text_repel(data = df.strict, size = 4, segment.color = '#FFDDDD')
        
        ret <- grid.hcat(list(p1, p2), widths = c(2, 3.5))
        return(ret)
    }

    gwas.pip.plot <- function() {

        .pip.cutoff <- min((df %>% filter(qval < 1e-4))$pip)
        aes.gwas <- aes(x = ld.gwas.p, y = pip, label = hgnc)
        ret <- gg.plot(df, aes.gwas) + geom_hex(bins = 40, color = 'gray') + scale.grad 
        ret <- ret +
            geom_hline(yintercept = .pip.cutoff, color = 'red', size = 1) +
                xlab('best -log10 GWAS p-value in LD') + ylab('posterior probability of mediation')
        
        ret <- ret + scale_y_continuous(breaks = c(0.01, 0.5, 0.99), trans = logit.trans)
        ret <- ret + geom_point(data = df.strict, color = 'black', fill = 'red', pch = 21) +
            geom_text_repel(data = df.strict, size = 4, segment.color = 'green', segment.alpha = 0.5)

        return(ret)
    }


    ## QTL vs pip
    qtl.pip.plot <- function() {

        .pip.cutoff <- min((df %>% filter(qval < 1e-4))$pip)
        aes.qtl <- aes(x = qtl.p, y = pip, label = hgnc)

        ret <- gg.plot(df, aes.qtl) + geom_hex(bins = 40, color = 'gray') + scale.grad 
        ret <- ret +
            geom_hline(yintercept = .pip.cutoff, color = 'red', lty = 2, size = 1) +
                xlab('best -log10 eQTL p-value') + ylab('posterior probability of mediation')
        
        ret <- ret + scale_y_continuous(breaks = c(0.01, 0.5, 0.99), trans = logit.trans)
        ret <- ret + geom_point(data = df.strict, color = 'black', fill = 'red', pch = 21) +
            geom_text_repel(data = df.strict, size = 4, segment.color = 'green', segment.alpha = 0.5)

        return(ret)
    }

    gwas.qtl.plot <- function(transpose = FALSE, ld.level = TRUE) {

        if(transpose) {
            if(ld.level) {
                aes.gwas <- aes(y = ld.gwas.p, x = qtl.p, label = hgnc)
            } else {
                aes.gwas <- aes(y = gwas.p, x = qtl.p, label = hgnc)
            }
        } else {
            if(ld.level) {
                aes.gwas <- aes(x = ld.gwas.p, y = qtl.p, label = hgnc)
            } else {
                aes.gwas <- aes(x = gwas.p, y = qtl.p, label = hgnc)
            }
        }

        ret <- gg.plot(df, aes.gwas) +
            geom_hex(bins = 40, color = 'gray') + scale.grad

        ret <- ret + geom_point(data = df.strict, pch = 21, color = 'black', fill = 'red') +
            geom_text_repel(data = df.strict,
                            size = 4, segment.color = 'green',
                            segment.alpha = 0.5)

        if(transpose) {
            ret <- ret + xlab('best -log10 eQTL p-value in cis') +
                ylab('best -log10 GWAS p-value' %&&% ifelse(ld.level, ' in LD', ''))
        } else {
            ret <- ret + ylab('best -log10 eQTL p-value in cis') +
                xlab('best -log10 GWAS p-value' %&&% ifelse(ld.level, ' in LD', ''))
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

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_gwas_pve.pdf',
           plot = gwas.pve.plot(), width = 6, height = 4, units = 'in')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_gwas_qtl.pdf',
           plot = gwas.qtl.plot(), width = 6, height = 4, units = 'in')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_qtl_pve.pdf',
           plot = qtl.pve.plot(), width = 6, height = 4, units = 'in')

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
        gg.plot(df %>% filter(ld.gwas.p > 4) %>% arrange(gwas.p),
                aes(x = gene.loc, y = theta, label = hgnc))

    ret <- ret +
        geom_errorbar(data = df %>% filter(ld.gwas.p > 4, pip > .5),
                      aes(ymin = theta - 2*theta.se, ymax = theta + 2*theta.se),
                      color = 'gray40')

    ret <- ret + geom_point(aes(fill = pmin(ld.gwas.p, 8), size = PVE), pch = 21) +
        scale_size_continuous('proportion of variance explained', range = c(0, 5), breaks = c(1e-2, 1e-1, .99))

    ret <- ret +
        scale_fill_gradientn('-log10 GWAS p-value', colors = c('gray80', 'gray40', 'purple'))

    ret <- ret +
        facet_grid(.~chr, scales = 'free', space = 'free') +
            scale_x_continuous(labels = .genome.mb(), breaks = c(0, 5e7, 1e8, 1.5e8, 2e8)) + .thm

    ret <- ret + xlab('genomic location (mb)') + ylab('mediation effect')

    df.strict.neg <- df.strict %>% filter(theta < 0)

    df.strict.pos <- df.strict %>% filter(theta > 0)

    ret <- ret + geom_text_repel(data = df.strict.neg, 
                                 aes(y = theta), color = 'blue', segment.alpha = 0.5,
                                 nudge_y = -.2, size = 3, segment.size = .5, segment.color = 'green')

    ret <- ret + geom_text_repel(data = df.strict.pos,
                                 aes(y = theta), color = 'red', segment.alpha = 0.5,
                                 nudge_y = .2, size = 3, segment.size = .5, segment.color = 'orange')

    ggsave(filename = 'figures/genes/Fig_global_' %&&% qtl.data %&&% '_mediation_effect.pdf',
           plot = ret, width = 8, height = 4, units = 'in',
           limitsize = FALSE)
}

draw.plots('full-0')
