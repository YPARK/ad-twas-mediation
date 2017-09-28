#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

ld.file <- argv[1]   # e.g., ld.file = 'stat/IGAP/ld/19.ld.gz'
plink.hdr <- argv[2] # e.g., plink.hdr = 'geno/rosmap1709-chr19'

out.file <- argv[3]

dir.create(dirname(out.file), recursive = TRUE)

source('util.R')
source('mediation.R')
options(stringsAsFactors = FALSE)

if(file.exists(out.file)) {
    log.msg('File exists : %s\n\n\n', out.file)
    q()
}

cutoff <- 1000

## temporary directory
temp.dir <- system('mktemp -d ' %&&% out.file %&&% '-temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

ld.tab <- read.table(ld.file, col.names = c('chr', 'ld.lb', 'ld.ub', 'n.snp.ld', 'n.qtl.ld'))
out.tab <- NULL

for(ld.idx in 1:nrow(ld.tab)) {

    ld.info <- ld.tab[ld.idx, ]
    x.bim <- subset.plink.snps(ld.info, temp.dir, plink.hdr)
    n.snps <- nrow(x.bim)
    log.msg('num SNPs = %d\n\n', n.snps)

    if(nrow(x.bim) >= cutoff)
        out.tab <- rbind(out.tab, data.frame(ld.info))
}

write.tab(out.tab, gzfile(out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished\n')
