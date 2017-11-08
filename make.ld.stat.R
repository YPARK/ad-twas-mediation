#!/usr/bin/env Rscript

## Take LD statistics
source('util.R')
library(readr)
library(dplyr)
library(tidyr)

take.best.gwas <- function(chr) {

    ld.file = 'stat/IGAP/ld/' %&&% chr %&&% '.ld.gz'
    gwas.file = 'IGAP/chr' %&&% chr %&&% '.txt.gz'

    ld.cols <- c('chr', 'ld.lb', 'ld.ub', 'ld.overlap', 'ld.size')

    gwas.cols <- c('chr', 'snp.loc.1', 'snp.loc', 'rs', 'gwas.a1', 'gwas.a2',
                   'gwas.z', 'gwas.theta', 'gwas.se')

    ld.tab <- read_tsv(ld.file, col_names = ld.cols)

    gwas.tab <- read_tsv(gwas.file, col_names = gwas.cols) %>%
        mutate(chr = gsub(chr, pattern = 'chr', replacement = ''))

    .best.gwas <- function(ld.idx) {

        ld.lb <- ld.tab[ld.idx, ]$ld.lb
        ld.ub <- ld.tab[ld.idx, ]$ld.ub

        ret <- gwas.tab %>% filter(snp.loc.1 >= ld.lb, snp.loc <= ld.ub) %>%
            slice(which.max(abs(gwas.z))) %>%
                mutate(ld.lb, ld.ub) %>%
                    select(chr, ld.lb, ld.ub, gwas.z, rs, snp.loc)
        return(ret)
    }

    n.ld <- nrow(ld.tab)

    best.gwas.tab <- do.call(rbind, lapply(1:n.ld, .best.gwas))
    return(best.gwas.tab)
}

best.gwas.tab <- do.call(rbind, lapply(1:22, take.best.gwas)) %>%
    mutate(gwas.p = 2 * pnorm(abs(gwas.z), lower.tail = FALSE))

write_tsv(best.gwas.tab, path = gzfile('stat/IGAP/ld.summary.txt.gz'))

