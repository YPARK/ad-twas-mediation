#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 4) q()

qtl.file <- argv[1] # e.g., qtl.file <- 'eqtl/21/b10.qtl-raw-y0.gz'
snp.file <- argv[2] # e.g., snp.file <- 'eqtl/21/b10.snps.gz'
gwas.file <- argv[3] # e.g., gwas.file <- 'data/IGAP/chr21.ld_bed.gz'
out.file <- argv[4] # e.g., 'temp.ucsc_bed.gz'

if(file.exists(out.file)) q()

dir.create(dirname(out.file), recursive = TRUE)

library(dplyr)
options(stringsAsFactors = FALSE, scipen=999)

.read.tab <- function(.file, .cols) read.table(.file, col.names = .cols, sep = '\t')

snp.cols <- c('chr', 'rs', '.', 'snp.loc', 'eqtl.a1', 'eqtl.a2')
qtl.cols <- c('snp.loc', 'null.gene', 'eqtl.theta', 'eqtl.z')
gwas.cols <- c('chr', 'snp.loc.1', 'snp.loc', 'rs',
               'gwas.a1', 'gwas.a2', 'gwas.z', 'gwas.theta', 'gwas.se', 'ld')

snp.tab <- .read.tab(snp.file, snp.cols)
qtl.tab <- .read.tab(qtl.file, qtl.cols)
gwas.tab <- .read.tab(gwas.file, gwas.cols) %>%
    mutate(chr = sapply(chr, gsub, pattern = 'chr', replacement = '')) %>%
        mutate(chr = as.integer(chr)) %>%
            select(-snp.loc.1)

out <- qtl.tab %>% left_join(snp.tab, by = c('snp.loc')) %>%
    left_join(gwas.tab, by = c('snp.loc', 'rs', 'chr')) %>%
        na.omit()

out <- out %>%
    filter(((gwas.a1 == eqtl.a1) & (gwas.a2 == eqtl.a2)) | ((gwas.a1 == eqtl.a2) & (gwas.a2 == eqtl.a1))) %>%
        mutate(gwas.z.flip = if_else(eqtl.a1 != gwas.a1, -gwas.z, gwas.z)) %>%
            mutate(gwas.theta.flip = if_else(eqtl.a1 != gwas.a1, -gwas.theta, gwas.theta))

out <- out %>% mutate(snp.loc.1 = snp.loc - 1) %>%
    select(chr, snp.loc.1, snp.loc, rs,
           eqtl.a1, eqtl.a2, eqtl.theta, eqtl.z,
           gwas.theta.flip, gwas.se, gwas.z.flip,
           ld, null.gene)

write.table(out, file = gzfile(out.file),
            col.names = FALSE, row.names = FALSE, quote = FALSE,
            sep = '\t')
