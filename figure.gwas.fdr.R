#!/usr/bin/env Rscript

library(qvalue)
library(dplyr)

source('util.R')
source('mediation.R')
source('figure.util.R')
library(scales)

gwas.files <- 'IGAP/chr' %&&% 1:22 %&&% '.txt.gz'

gwas.tab <- do.call(rbind, lapply(gwas.files, read.gwas.tab))

p.val.tab <- gwas.tab %>% mutate(pval = 2 * pnorm(abs(gwas.z),lower.tail = FALSE)) %>%
    select(rs, pval)

q.obj <- qvalue(p.val.tab$pval)


fdr.thres <- c(1e-4, 1e-2, seq(.1, 1, .1))
pval.cutoff <- sapply(fdr.thres, function(tt) max(q.obj$pvalues[q.obj$qvalues < tt]))

plt.df <- data.frame(fdr.thres, pval.cutoff)

l10.trans <- trans_new('.l10', function(x) -log10(x), function(y) 10^(-y))

plt <- gg.plot(plt.df, aes(x = fdr.thres, y = pval.cutoff)) +
    geom_point() +
    geom_line() +
        geom_text(aes(x = fdr.thres + .01, y = pval.cutoff, label = signif(pval.cutoff, 2))) +
        scale_y_continuous(breaks = c(1e-4, 1e-3, 1e-2, 1e-1), trans = l10.trans) +
        scale_x_continuous(breaks = fdr.thres) +
        xlab('controlled FDR') + ylab('p-value cutoff')

ggsave(filename = 'figures/IGAP_FDR.pdf', plot = plt, width = 4, height = 4, units = 'in')
