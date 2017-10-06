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

read.nwas <- function(in.hdr, in.cols) {
    .files <- in.hdr %&&% 1:22 %&&% '.nwas.gz'
    .dat <- lapply(.files, read_tsv, col_names = in.cols)
    return(do.call(rbind, .dat))
}

herit.genes.tab <- read_tsv('heritable.genes.qvalue.txt.gz', col_names = TRUE)
h.genes <- herit.genes.tab %>% filter(q.val < 0.1) %>% dplyr::select(ENSG)
h.genes <- h.genes$ENSG

expr.cols <- c('ensg', 'theta', 'theta.se', 'lodds',
               'emp.p', 'fd', 'nboot',
               'lodds.mean', 'lodds.se', 'cauchy.location', 'cauchy.scale',
               't.m', 't.s', 'cauchy.p', 't.p', 'norm.p',
               'hgnc', 'tss', 'tes', 'strand',
               'max.gwas.theta', 'max.gwas.z',
               'chr', 'ld.lb', 'ld.ub', 'n.snps')

## select only heritable genes
data.direct.tab <- read.mediation('bootstrap/direct_IGAP_rosmap_eqtl_hs-lm_', expr.cols) %>%
    filter(ensg %in% h.genes)

data.marginal.tab <- read.mediation('bootstrap/marginal_IGAP_rosmap_eqtl_hs-lm_', expr.cols) %>%
    filter(ensg %in% h.genes)

nwas.cols <- c('ensg', 'nwas.z', 'chr', 'ld.lb', 'ld.ub', 'n.snps')

nwas.tab <- read.nwas('nwas/IGAP_rosmap_eqtl_hs-lm_', nwas.cols) %>%
    select(-n.snps) %>%
    mutate(p.val.nwas = 2 * pnorm(abs(nwas.z), lower.tail = FALSE))

## generate null lodds statistics
combined.p.val <- function(.df) {
    n <- nrow(.df)
    lo.boot <- sweep(sweep(.rnorm(n, 100), 1, .df$lodds.se, `*`), 1, .df$lodds.mean, `+`)
    lo.boot <- as.vector(lo.boot)

    ret <- .df %>% select(chr, ld.lb, ld.ub, ensg, hgnc, tss, tes, theta, theta.se, lodds, max.gwas.z)
    p.val <- empPvals(ret$lodds, lo.boot)
    q.val <- qvalue(p.val, fdr.level = 0.05, pi0.method = 'bootstrap')
    ret <- data.frame(ret, p.val = p.val, q.val = q.val$qvalues)
    return(ret)
}

marginal.df <- combined.p.val(data.marginal.tab) %>%
    rename(p.val.marginal = p.val, q.val.marginal = q.val)

direct.df <- combined.p.val(data.direct.tab) %>%
    rename(p.val.direct = p.val, q.val.direct = q.val)    

.temp <- direct.df %>% dplyr::select(chr, ld.lb, ld.ub, ensg, p.val.direct, q.val.direct)

stat.tab <- marginal.df %>%
    left_join(.temp, by = c('chr', 'ld.lb', 'ld.ub', 'ensg')) %>%
    left_join(nwas.tab, by = c('chr', 'ld.lb', 'ld.ub', 'ensg'))

out.tab <- stat.tab %>%
    filter(p.val.marginal < 2.5e-6, p.val.direct < 2.5e-6, lodds > 0, p.val.nwas < 2.5e-6,
           sign(theta) == sign(nwas.z)) %>%
    arrange(chr, tss)

.temp <- out.tab %>% select(chr, ld.lb, ld.ub) %>% unique()

out.ld.tab <- stat.tab %>%
    right_join(.temp, by = c('chr', 'ld.lb', 'ld.ub'))

out.strict.tab <- out.tab %>%
    filter(abs(max.gwas.z) > abs(qnorm(1e-4/2))) %>%
    arrange(chr, tss)

.temp <- out.strict.tab %>% select(chr, ld.lb, ld.ub) %>% unique()

out.strict.ld.tab <- stat.tab %>%
    right_join(.temp, by = c('chr', 'ld.lb', 'ld.ub'))

write.tab.named(stat.tab, file = gzfile('tables/bootstrap_gene.txt.gz'))

write.tab.named(out.tab, file = gzfile('tables/bootstrap_gene_significant.txt.gz'))

write.tab.named(out.ld.tab, file = gzfile('tables/bootstrap_ld_significant.txt.gz'))

write.tab.named(out.strict.tab, file = gzfile('tables/bootstrap_gene_strict.txt.gz'))

write.tab.named(out.strict.ld.tab, file = gzfile('tables/bootstrap_ld_strict.txt.gz'))
