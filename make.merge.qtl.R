#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

# e.g., argv = c(TRUE, 'out.file', 'stat/IGAP/hic-data/CO/hs-lm/1.eqtl_bed.gz', 'stat/IGAP/hic-data/HC/hs-lm/1.eqtl_bed.gz')

if(length(argv) < 4) {
    q()
}

n.arg <- length(argv)

is.eqtl <- as.logical(argv[1])
out.file <- argv[2]
in.files <- argv[-(1:2)]

dir.create(dirname(out.file), recursive = TRUE)

library(dplyr)
library(readr)

eqtl.cols <- c('chr', 'snp.loc.1', 'snp.loc', 'rs',
               'qtl.a1', 'qtl.a2', 'qtl.theta', 'qtl.z',
               'gwas.theta', 'gwas.se', 'gwas.z',
               'ld', 'med.id', 'hgnc', 'tss', 'tes', 'strand')

mqtl.cols <- c('chr', 'snp.loc.1', 'snp.loc', 'rs',
               'qtl.a1', 'qtl.a2', 'qtl.theta', 'qtl.z',
               'gwas.theta', 'gwas.se', 'gwas.z',
               'ld', 'med.id', 'cg.loc')

.read.tsv <- function(...) {
    if (is.eqtl) { ret <- read_tsv(..., col_names = eqtl.cols) }
    else { ret <- read_tsv(..., col_names = mqtl.cols) }
    return(ret)
}

data.tab <- do.call(rbind, lapply(in.files, .read.tsv))

ret <- data.tab %>%
    group_by(chr, snp.loc.1, snp.loc, rs, qtl.a1, qtl.a2, med.id) %>%
        slice(which.max(abs(qtl.z)))

write_tsv(ret, path = gzfile(out.file), col_names = FALSE)

