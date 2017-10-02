#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)
source('util.R')
library(dplyr)
library(qvalue)

dir.create('tables', recursive = TRUE, showWarnings = FALSE)

read.mediation <- function(in.hdr, in.cols) {
    .files <- in.hdr %&&% 1:22 %&&% '.mediation.gz'
    .dat <- lapply(.files, read.table, col.names = in.cols)
    return(do.call(rbind, .dat))
}

expr.cols <- c('ensg', 'theta', 'theta.se', 'lodds',
               'emp.p', 'fd', 'nboot',
               'lodds.mean', 'lodds.se', 'cauchy.location', 'cauchy.scale',
               't.m', 't.s', 'cauchy.p', 't.p', 'norm.p',
               'hgnc', 'tss', 'tes', 'strand',
               'max.gwas.theta', 'max.gwas.z',
               'chr', 'ld.1', 'ld.2', 'n.snps')

data.direct.tab <- read.mediation('bootstrap/direct_IGAP_rosmap_eqtl_hs-lm_', expr.cols)
data.marginal.tab <- read.mediation('bootstrap/marginal_IGAP_rosmap_eqtl_hs-lm_', expr.cols)


## generate null lodds statistics
combined.p.val <- function(.df) {
    n <- nrow(.df)
    lo.boot <- sweep(sweep(.rnorm(n, 100), 1, .df$lodds.se, `*`), 1, .df$lodds.mean, `+`)
    lo.boot <- as.vector(lo.boot)

    ret <- .df %>% select(chr, ld.1, ld.2, ensg, hgnc, tss, tes, theta, theta.se, lodds, max.gwas.z)
    p.val <- empPvals(ret$lodds, lo.boot)
    q.val <- qvalue(p.val, fdr.level = 0.05, pi0.method = 'bootstrap')
    ret <- data.frame(ret, p.val = p.val, q.val = q.val$qvalues)
    return(ret)
}

marginal.df <- combined.p.val(data.marginal.tab)
direct.df <- combined.p.val(data.direct.tab)

stat.cols <- c('chr', 'ld.1', 'ld.2', 'ensg', 'hgnc', 'tss', 'tes', 'theta', 'theta.se', 'lodds', 'max.gwas.z')

stat.tab <- marginal.df %>% rename(p.val.marginal = p.val, q.val.marginal = q.val) %>%
    left_join(direct.df %>% rename(p.val.direct = p.val, q.val.direct = q.val), by = stat.cols)

out.tab <- stat.tab %>%
    filter(p.val.marginal < 1e-6, p.val.direct < 1e-6) %>%
    filter(abs(theta/theta.se) > abs(qnorm(1e-6/2)), abs(max.gwas.z) > 3) %>%
    arrange(chr, tss)

write.tab.named(stat.tab, file = gzfile('tables/bootstrap_gene.txt.gz'))
write.tab.named(out.tab, file = gzfile('tables/bootstrap_gene_significant.txt.gz'))

