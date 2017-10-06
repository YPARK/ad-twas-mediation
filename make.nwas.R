#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) q()

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/1.ld.gz'
sum.file <- argv[2]                     # e.g., sum.file = 'stat/IGAP/data/hs-lm/1.eqtl_bed.gz'
plink.hdr <- argv[3]                    # e.g., plink.hdr = 'geno/rosmap1709-chr1'
ld.idx <- as.integer(argv[4])           # e.g., ld.idx = 37
is.eqtl <- as.logical(argv[5])          # e.g., is.eqtl = TRUE
out.file <- argv[6]                     # e.g., out.file = 'temp.gz'

dir.create(dirname(out.file), recursive = TRUE)

source('util.R')
options(stringsAsFactors = FALSE)

if(file.exists(out.file)) {
    log.msg('File exists : %s\n\n\n', out.file)
    q()
}

library(zqtl)
library(dplyr)
library(methods)
library(broom)
source('mediation.R')
n.snp.cutoff <- 100

## temporary directory
temp.dir <- system('mktemp -d ' %&&% out.file %&&% '-temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

ld.tab <- read.table(ld.file, col.names = c('chr', 'ld.lb', 'ld.ub', 'n.snp.ld', 'n.qtl.ld'))
ld.info <- ld.tab[ld.idx, ]

plink <- subset.plink(ld.info, temp.dir, plink.hdr)
x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))

sum.stat.out <- extract.sum.stat(ld.info, sum.file, x.bim, temp.dir, is.eqtl)

sum.stat <- sum.stat.out$sum.stat
mediators <- sum.stat.out$mediators

## Only work on sufficiently large LD blocks
n.snps <- sum.stat %>% select(snp.loc) %>% unique() %>% nrow()
if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n\n\n', n.snps, n.snp.cutoff)
    system('echo | gzip > ' %&&% out.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

zqtl.data <- make.zqtl.data(plink, sum.stat, mediators)


svd.out <- take.ld.svd(zqtl.data$X, options = list(do.stdize = TRUE, eigen.tol = 1e-2))
W.t <- sweep(svd.out$V.t, 1, svd.out$D, `/`)
## V.t <- sweep(svd.out$V.t, 1, svd.out$D, `*`)
## LD.inv <- t(W.t) %*% W.t
## LD <- t(V.t) %*% V.t

gwas.z <- zqtl.data$gwas.theta / zqtl.data$gwas.se

take.twas <- function(g) {

    qtl.z <- (zqtl.data$qtl.theta %c% g) / (zqtl.data$qtl.se %c% g)
    qtl.z[is.na(qtl.z)] <- 0

    ## This is too memory-intensive
    ## qtl.z.poly <- LD.inv %*% qtl.z
    ## num <- t(qtl.z.poly) %*% gwas.z
    ## denom <- t(qtl.z.poly) %*% LD %*% qtl.z.poly

    num <- t(W.t %*% gwas.z) %*% (W.t %*% qtl.z)
    denom <- t(W.t %*% qtl.z) %*% (W.t %*% qtl.z)
    log.msg('TWAS finished [%d / %d]\n', g, n.med)

    return(signif(as.numeric(num/sqrt(denom + 1e-16)), 4))
}

n.med <- ncol(zqtl.data$qtl.theta)
out.tab <- data.frame(med.id = mediators$med.id,
                      sapply(1:n.med, take.twas),
                      ld = paste(ld.info[1:3], collapse = '\t'),
                      ld.size = n.snps)

if(is.eqtl) names(out.tab)[2] <- 'TWAS'
if(!is.eqtl) names(out.tab)[2] <- 'MWAS'

write.tab(out.tab, file = gzfile(out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished zQTL\n\n\n')
