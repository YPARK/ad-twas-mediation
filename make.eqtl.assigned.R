#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 5) q()

qtl.file <- argv[1] # e.g., qtl.file <- 'hic-eqtl/CO/1/b1.qtl-hs-lm.gz'
snp.file <- argv[2] # e.g., snp.file <- 'hic-eqtl/CO/1/b1.snps.gz'
gene.file <- argv[3] # e.g., gene.file <- 'hic-eqtl/CO/1/b1.genes.gz'
gwas.file <- argv[4] # e.g., gwas.file <- 'data/IGAP/chr1.ld_bed.gz'
out.file <- argv[5] # e.g., out.file <- 'temp.ucsc_bed.gz'

if(file.exists(out.file)) {
    q()
}

dir.create(dirname(out.file), recursive = TRUE)

options(stringsAsFactors = FALSE, scipen=999)
library(readr)
library(dplyr)
source('util.R')

.files <- c(qtl.file, snp.file, gene.file, gwas.file)
if(!all(sapply(.files, file.exists))) {
    if(file.exists(gene.file)) {
        system('printf "" | gzip > ' %&&% out.file)
        log.msg('no QTL connection\n')
        q()
    }
    log.msg('incomplete list of input files:\n\n%s\n', paste(.files, collpse='\n'))
    q()
}

.read.tab <- function(.file, .cols) {
    read_tsv(.file, col_names = .cols)
}

snp.cols <- c('chr', 'rs', '.', 'snp.loc', 'eqtl.a1', 'eqtl.a2')
gene.cols <- c('chr', 'tss', 'tes', 'strand', 'ensg', 'hgnc')
qtl.cols <- c('snp.loc', 'ensg', 'eqtl.theta', 'eqtl.z')
gwas.cols <- c('chr', 'snp.loc.1', 'snp.loc', 'rs',
               'gwas.a1', 'gwas.a2', 'gwas.z', 'gwas.theta', 'gwas.se', 'ld')

qtl.tab <- .read.tab(qtl.file, qtl.cols)

if(nrow(qtl.tab) == 0) {
    system('printf "" | gzip > ' %&&% out.file)
}

snp.tab <- .read.tab(snp.file, snp.cols)
gene.tab <- .read.tab(gene.file, gene.cols)
gwas.tab <- .read.tab(gwas.file, gwas.cols) %>%
    mutate(chr = sapply(chr, gsub, pattern = 'chr', replacement = '')) %>%
        mutate(chr = as.integer(chr)) %>%
            select(-snp.loc.1)

out <- qtl.tab %>% left_join(snp.tab, by = c('snp.loc')) %>%
    left_join(gene.tab, by = c('chr', 'ensg')) %>%
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
           ld, ensg, hgnc, tss, tes, strand)

write.table(out, file = gzfile(out.file),
            col.names = FALSE, row.names = FALSE, quote = FALSE,
            sep = '\t')
