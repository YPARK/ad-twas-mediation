#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
library(dplyr)
library(qvalue)
library(readr)

dir.create('tables', recursive = TRUE, showWarnings = FALSE)

read.chr <- function(in.hdr, in.cols, file.name) {
    .files <- in.hdr %&&% 1:22 %&&% '.' %&&% file.name %&&% '.gz'
    .dat <- lapply(.files, read_tsv, col_names = in.cols)
    return(do.call(rbind, .dat))
}

read.boot <- function(in.hdr) {
    in.cols <- c('ensg', 'theta', 'theta.se', 'lodds',
                 'emp.p', 'fd', 'nboot',
                 'lodds.mean', 'lodds.se', 'cauchy.location', 'cauchy.scale',
                 't.m', 't.s', 'cauchy.p', 't.p', 'norm.p',
                 'hgnc', 'tss', 'tes', 'strand',
                 'max.gwas.theta', 'max.gwas.z',
                 'chr', 'ld.lb', 'ld.ub', 'n.snps')
    read.chr(in.hdr, in.cols, 'mediation')
}

read.nwas <- function(in.hdr) {
    in.cols <- c('ensg', 'nwas.z', 'chr', 'ld.lb', 'ld.ub', 'n.snps')
    read.chr(in.hdr, in.cols, 'nwas')
}

read.pve <- function(in.hdr) {
    in.cols <- c('chr', 'ld.lb', 'ld.ub', 'med.id', 'theta', 'theta.se', 'lodds', 'hgnc',
                 'tss', 'tes', 'strand', 'best.gwas.z', 'best.gwas.rs', 'best.gwas.loc',
                 'best.qtl.z', 'best.qtl.rs', 'best.qtl.loc',
                 'v.med', 'v.med.tot', 'v.dir', 'v.resid')    
    read.chr(in.hdr, in.cols, 'pve')
}

herit.genes.tab <- read_tsv('heritable.genes.qvalue.txt.gz', col_names = TRUE)
h.genes <- herit.genes.tab %>% filter(q.val < 0.1) %>% dplyr::select(ENSG)
h.genes <- h.genes$ENSG

## Bootstrap p-value distributions
boot.direct.tab <- read.boot('bootstrap/direct_IGAP_rosmap_eqtl_hs-lm_') %>%
    filter(ensg %in% h.genes)

boot.marginal.tab <- read.boot('bootstrap/marginal_IGAP_rosmap_eqtl_hs-lm_') %>%
    filter(ensg %in% h.genes)

null.lodds <- function(.df) {
    n <- nrow(.df)
    lo.boot <- sweep(sweep(.rnorm(n, 100), 1, .df$lodds.se, `*`), 1, .df$lodds.mean, `+`)
    lo.boot <- as.vector(lo.boot)
    return(lo.boot)
}

null.direct <- null.lodds(boot.direct.tab)
null.marginal <- null.lodds(boot.marginal.tab)

## read PVE calculation and add p-values
pve.tab <- read.pve('pve/IGAP_rosmap_eqtl_hs-lm_') %>%
    rename(ensg = med.id) %>%
    filter(ensg %in% h.genes)

pval.direct <- empPvals(pve.tab$lodds, null.direct)
pval.marginal <- empPvals(pve.tab$lodds, null.marginal)

pve.tab <- cbind(pve.tab,
                 p.val.direct = pval.direct,
                 p.val.marginal = pval.marginal)

## read TWAS
nwas.tab <- read.nwas('nwas/IGAP_rosmap_eqtl_hs-lm_') %>%
    select(-n.snps) %>%
    mutate(p.val.nwas = 2 * pnorm(abs(nwas.z), lower.tail = FALSE))


## Overall statistics
stat.tab <- pve.tab %>%
    left_join(nwas.tab, by = c('chr', 'ld.lb', 'ld.ub', 'ensg'))

## significant mediation
out.tab <- stat.tab %>%
    filter(p.val.marginal < 2.5e-6, p.val.direct < 2.5e-6, lodds > 0) %>%
    arrange(chr, tss)

.temp <- out.tab %>% select(chr, ld.lb, ld.ub) %>% unique()

out.ld.tab <- stat.tab %>%
    right_join(.temp, by = c('chr', 'ld.lb', 'ld.ub'))

out.strict.tab <- out.tab %>%
    filter(abs(best.gwas.z) > abs(qnorm(1e-4/2)),
           abs(best.qtl.z) > abs(qnorm(1e-4/2))) %>%
    arrange(chr, tss)

.temp <- out.strict.tab %>% select(chr, ld.lb, ld.ub) %>% unique()

out.strict.ld.tab <- stat.tab %>%
    right_join(.temp, by = c('chr', 'ld.lb', 'ld.ub'))


## significant TWAS but no mediation
out.twas.tab <- stat.tab %>%
    filter(p.val.nwas < 2.5e-6, lodds < -2)

.temp <- out.twas.tab %>% select(chr, ld.lb, ld.ub) %>% unique()

out.twas.ld.tab <- stat.tab %>%
    right_join(.temp, by = c('chr', 'ld.lb', 'ld.ub'))

################################################################
write.tab.named(stat.tab, file = gzfile('tables/genes.txt.gz'))
write.tab.named(out.tab, file = gzfile('tables/genes_significant.txt.gz'))
write.tab.named(out.ld.tab, file = gzfile('tables/genes_ld_significant.txt.gz'))
write.tab.named(out.strict.tab, file = gzfile('tables/genes_strict.txt.gz'))
write.tab.named(out.strict.ld.tab, file = gzfile('tables/genes_ld_strict.txt.gz'))
write.tab.named(out.twas.tab, file = gzfile('tables/genes_twas.txt.gz'))
write.tab.named(out.twas.ld.tab, file = gzfile('tables/genes_ld_twas.txt.gz'))
