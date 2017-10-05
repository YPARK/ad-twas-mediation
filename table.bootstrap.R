#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)
source('util.R')
library(dplyr)
library(qvalue)
library(readr)

dir.create('tables', recursive = TRUE, showWarnings = FALSE)

read.mediation <- function(in.hdr, in.cols) {
    .files <- in.hdr %&&% 1:22 %&&% '.mediation.gz'
    .dat <- lapply(.files, read_tsv, col_names = in.cols)
    return(do.call(rbind, .dat))
}

herit.genes.tab <- read_tsv('heritable.genes.qvalue.txt.gz', col_names = TRUE)
h.genes <- herit.genes.tab %>% filter(q.val < 0.1) %>% dplyr::select(ENSG)
h.genes <- h.genes$ENSG

stwas.tab <- read_tsv('stwas-rosmap.fdr.txt.gz', col_names = TRUE) %>%
    rename(stwas.p.val = p.val, stwas.q.val = q.val)

expr.cols <- c('ensg', 'theta', 'theta.se', 'lodds',
               'emp.p', 'fd', 'nboot',
               'lodds.mean', 'lodds.se', 'cauchy.location', 'cauchy.scale',
               't.m', 't.s', 'cauchy.p', 't.p', 'norm.p',
               'hgnc', 'tss', 'tes', 'strand',
               'max.gwas.theta', 'max.gwas.z',
               'chr', 'ld.1', 'ld.2', 'n.snps')

## select only heritable genes
data.direct.tab <- read.mediation('bootstrap/direct_IGAP_rosmap_eqtl_hs-lm_', expr.cols) %>%
    filter(ensg %in% h.genes)

data.marginal.tab <- read.mediation('bootstrap/marginal_IGAP_rosmap_eqtl_hs-lm_', expr.cols) %>%
    filter(ensg %in% h.genes)


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

stat.tab <- marginal.df %>%
    rename(p.val.marginal = p.val, q.val.marginal = q.val)

stat.tab <- stat.tab %>%
    left_join(direct.df %>% rename(p.val.direct = p.val, q.val.direct = q.val),
              by = stat.cols)

stat.tab <- stat.tab %>% left_join(stwas.tab, by = 'ensg')

out.tab <- stat.tab %>%
    filter(p.val.marginal < 2.5e-6, p.val.direct < 2.5e-6) %>%
    arrange(chr, tss)

out.strict.tab <- stat.tab %>%
    filter(p.val.marginal < 2.5e-6, p.val.direct < 2.5e-6, lodds > 0) %>%
    filter(abs(max.gwas.z) > abs(qnorm(1e-4/2))) %>%
    arrange(chr, tss)

write.tab.named(stat.tab, file = gzfile('tables/bootstrap_gene.txt.gz'))

write.tab.named(out.tab, file = gzfile('tables/bootstrap_gene_significant.txt.gz'))

write.tab.named(out.strict.tab, file = gzfile('tables/bootstrap_gene_strict.txt.gz'))

