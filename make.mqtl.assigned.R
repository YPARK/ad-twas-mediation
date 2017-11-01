#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 5) q()

qtl.file <- argv[1] # e.g., qtl.file <- 'mqtl/chr1/b1.qtl-hs-lm.gz'
snp.file <- argv[2] # e.g., snp.file <- 'mqtl/chr1/b1.snps.gz'
probe.file <- argv[3] # e.g., probe.file <- 'mqtl/chr1/b1.probes.gz'
gwas.file <- argv[4] # e.g., gwas.file <- 'data/IGAP/chr1.ld_bed.gz'
out.file <- argv[5] # e.g., out.file <- 'temp.ucsc_bed.gz'

options(stringsAsFactors = FALSE, scipen=999)
library(readr)
library(dplyr)
source('util.R')

if(file.exists(out.file)) {
    log.msg('File : %s exists\n', out.file)
    q()
}

dir.create(dirname(out.file), recursive = TRUE)

.files <- c(qtl.file, snp.file, probe.file, gwas.file)
if(!all(sapply(.files, file.exists))) {
    if(file.exists(probe.file)) {
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

snp.cols <- c('chr', 'rs', '.', 'snp.loc', 'mqtl.a1', 'mqtl.a2')
probe.cols <- c('cg', 'chr', 'cg.loc')
qtl.cols <- c('snp.loc', 'cg', 'mqtl.theta', 'mqtl.z')
gwas.cols <- c('chr', 'snp.loc.1', 'snp.loc', 'rs',
               'gwas.a1', 'gwas.a2', 'gwas.z', 'gwas.theta', 'gwas.se', 'ld')

qtl.tab <- .read.tab(qtl.file, qtl.cols)

if(nrow(qtl.tab) == 0) {
    system('printf "" | gzip > ' %&&% out.file)
}

snp.tab <- .read.tab(snp.file, snp.cols)
probe.tab <- .read.tab(probe.file, probe.cols)
gwas.tab <- .read.tab(gwas.file, gwas.cols) %>%
    mutate(chr = sapply(chr, gsub, pattern = 'chr', replacement = '')) %>%
        mutate(chr = as.integer(chr)) %>%
            select(-snp.loc.1)

out <- qtl.tab %>% left_join(snp.tab, by = c('snp.loc')) %>%
    left_join(probe.tab, by = c('chr', 'cg')) %>%
        left_join(gwas.tab, by = c('snp.loc', 'rs', 'chr')) %>%
            na.omit()

out <- out %>%
    filter(((gwas.a1 == mqtl.a1) & (gwas.a2 == mqtl.a2)) | ((gwas.a1 == mqtl.a2) & (gwas.a2 == mqtl.a1))) %>%
        mutate(gwas.z.flip = if_else(mqtl.a1 != gwas.a1, -gwas.z, gwas.z)) %>%
            mutate(gwas.theta.flip = if_else(mqtl.a1 != gwas.a1, -gwas.theta, gwas.theta))

out <- out %>% mutate(snp.loc.1 = snp.loc - 1) %>%
    select(chr, snp.loc.1, snp.loc, rs,
           mqtl.a1, mqtl.a2, mqtl.theta, mqtl.z,
           gwas.theta.flip, gwas.se, gwas.z.flip,
           ld, cg, cg.loc)

write.table(out, file = gzfile(out.file),
            col.names = FALSE, row.names = FALSE, quote = FALSE,
            sep = '\t')
