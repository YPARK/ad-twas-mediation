#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 8) q()

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/1.ld.gz'
sum.file <- argv[2]                     # e.g., sum.file = 'stat/IGAP/data/hs-lm/1.mqtl_bed.gz'
plink.hdr <- argv[3]                    # e.g., plink.hdr = 'geno/rosmap1709-chr1'
ld.idx <- as.integer(argv[4])           # e.g., ld.idx = 15
gwas.sample.size <- as.numeric(argv[5]) # e.g., gwas.sample.size = 74000
qtl.sample.size <- as.numeric(argv[6])  # e.g., qtl.sample.size = 350
is.eqtl <- as.logical(argv[7])          # e.g., is.eqtl = FALSE
out.hdr <- argv[8]                      # e.g., out.hdr = 'temp'

dir.create(dirname(out.hdr), recursive = TRUE)

source('util.R')
options(stringsAsFactors = FALSE)

z.out.file <- out.hdr %&&% '.mediation.gz'
snp.out.file <- out.hdr %&&% '.qtl.gz'

if(file.exists(z.out.file)) {
    log.msg('File exists : %s\n\n\n', z.out.file)
    q()
}

library(zqtl)
library(dplyr)
library(methods)
library(broom)
source('mediation.R')

n.snp.cutoff <- 100

## temporary directory
temp.dir <- system('mktemp -d ' %&&% out.hdr %&&% '-temp.XXXX',
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
    log.msg('This LD block is too small : %d SNPs < %d\n', n.snps, n.snp.cutoff)
    q()
}

zqtl.data <- make.zqtl.data(plink, sum.stat, mediators)

## Just run the zQTL w/o bootstrapping, but finemap
vb.opt <- list(pi = -0, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 5000, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2, nsample = 10, print.interv = 100,
               nboot = 0, med.finemap = TRUE, med.lodds.cutoff = 0,
               weight.y = TRUE, weight.m = TRUE)

z.out <- fit.med.zqtl(zqtl.data$gwas.theta, zqtl.data$gwas.se,
                      zqtl.data$qtl.theta, zqtl.data$qtl.se,
                      X = zqtl.data$X, n = gwas.sample.size,
                      n.med = qtl.sample.size, options = vb.opt)

phase.1 <- melt.effect(z.out$param.mediated, mediators$med.id, 1) %>%
    rename(med.id = Var1) %>%
        mutate(theta.se = theta.var) %>%
            dplyr::select(med.id, theta, theta.se, lodds, -theta.var, -Var2)

if('param.mediated' %in% names(z.out$finemap)) {

    .med <- mediators$med.id[z.out$finemap$mediators]

    phase.2 <- melt.effect(z.out$finemap$param.mediated, 2, .med) %>%
        rename(med.id = Var2, theta.2 = theta, lodds.2 = lodds) %>%
            mutate(theta.se.2 = sqrt(theta.var)) %>%
                dplyr::select(med.id, theta.2, theta.se.2, lodds.2, -Var1, -theta.var)

    out.tab <- phase.1 %>% left_join(phase.2, by = 'med.id') %>%
        left_join(mediators, by = 'med.id')

    snp.tab <- melt.effect(z.out$finemap$param.qtl, zqtl.data$snps$rs, .med) %>%
        filter(lodds > -0) %>%
            rename(rs = Var1, med.id = Var2) %>%
                mutate(theta.se = sqrt(theta.var)) %>%
                    select(med.id, rs, theta, theta.se, lodds) %>%
                        left_join(zqtl.data$snps, by = 'rs')

    if(nrow(snp.tab) > 0) {
        write.tab(snp.tab, file = gzfile(snp.out.file))
    }

} else {
    out.tab <- phase.1 %>%
        mutate(theta.2 = NA, theta.se.2 = NA, lodds.2 = NA) %>%
            left_join(mediators, by = 'med.id')
}

out.tab <- out.tab %>% mutate(ld = paste(ld.info[1:3], collapse = '\t'), ld.size = n.snps)

write.tab(out.tab, file = gzfile(z.out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished zQTL\n\n\n')
