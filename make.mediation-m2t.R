#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 9) q()

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/1.ld.gz'
mqtl.stat.file <- argv[2]               # e.g., mqtl.stat.file = 'stat/IGAP/data/hs-lm/1.mqtl_bed.gz'
eqtl.stat.file <- argv[3]               # e.g., eqtl.stat.file = 'stat/IGAP/data/hs-lm/1.eqtl_bed.gz'
plink.hdr <- argv[4]                    # e.g., plink.hdr = 'geno/rosmap1709-chr1'
ld.idx <- as.integer(argv[5])           # e.g., ld.idx = 15
gwas.sample.size <- as.numeric(argv[6]) # e.g., gwas.sample.size = 74000
mqtl.sample.size <- as.numeric(argv[7]) # e.g., mqtl.sample.size = 598
eqtl.sample.size <- as.numeric(argv[8]) # e.g., eqtl.sample.size = 356
out.hdr <- argv[9]                      # e.g., out.hdr = 'temp'

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

eqtl.stat <- extract.sum.stat(ld.info, eqtl.stat.file, x.bim, temp.dir, is.eqtl = TRUE)

eqtl.sum.stat <- eqtl.stat$sum.stat
eqtl.mediators <- eqtl.stat$mediators

## Only work on sufficiently large LD blocks
n.snps <- eqtl.sum.stat %>% select(snp.loc) %>% unique() %>% nrow()
if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n', n.snps, n.snp.cutoff)
    q()
}

mqtl.stat <- extract.sum.stat(ld.info, mqtl.stat.file, x.bim, temp.dir, is.eqtl = FALSE)

mqtl.sum.stat <- mqtl.stat$sum.stat
mqtl.mediators <- mqtl.stat$mediators

zqtl.data.eqtl <- make.zqtl.data(plink, eqtl.sum.stat, eqtl.mediators)
.zqtl.data.mqtl <- make.zqtl.data(plink, mqtl.sum.stat, mqtl.mediators)
m2t <- match(.zqtl.data.mqtl$snps$rs, zqtl.data.eqtl$snps$rs)
.zqtl.data.mqtl[['X']] <- NULL
zqtl.data.mqtl <- lapply(.zqtl.data.mqtl, function(xx) xx %r% m2t)

## Just run the zQTL w/o bootstrapping, but finemap
vb.opt <- list(pi = -0, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 5000, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2, nsample = 10, print.interv = 100,
               nboot = 0, med.finemap = TRUE, med.lodds.cutoff = 0)

z.out <- fit.med.zqtl(zqtl.data.eqtl$qtl.theta,
                      zqtl.data.eqtl$qtl.se,
                      zqtl.data.mqtl$qtl.theta,
                      zqtl.data.mqtl$qtl.se,
                      X = zqtl.data.eqtl$X,
                      n = eqtl.sample.size,
                      n.med = mqtl.sample.size,
                      options = vb.opt)

phase.1 <- melt.effect(z.out$param.mediated, mqtl.mediators$med.id, eqtl.mediators$med.id) %>%
    rename(cg = Var1, ensg = Var2) %>%
        mutate(theta.se = sqrt(theta.var)) %>%
            dplyr::select(ensg, cg, theta, theta.se, lodds)

.temp.ensg <- eqtl.mediators %>% rename(ensg = med.id) %>% dplyr::select(-gwas.theta,-gwas.z)
.temp.cg <- mqtl.mediators %>% rename(cg = med.id) %>% dplyr::select(-gwas.theta,-gwas.z)

if('param.mediated' %in% names(z.out$finemap)) {

    .cg <- mqtl.mediators$med.id[z.out$finemap$mediators]
    .ensg <- eqtl.mediators$med.id

    phase.2 <- melt.effect(z.out$finemap$param.mediated, .ensg, .cg) %>%
        rename(ensg = Var1, cg = Var2, theta.2 = theta, lodds.2 = lodds) %>%
            mutate(theta.se.2 = sqrt(theta.var)) %>%
                dplyr::select(ensg, cg, theta.2, theta.se.2, lodds.2)

    out.tab <- phase.1 %>% left_join(phase.2, by = c('ensg', 'cg')) %>%
        left_join(.temp.ensg, by = 'ensg') %>%
            left_join(.temp.cg, by = 'cg')

    snp.tab <- melt.effect(z.out$finemap$param.qtl, zqtl.data.mqtl$snps$rs, .cg) %>%
        filter(lodds > -0) %>%
            rename(rs = Var1, cg = Var2) %>%
                mutate(theta.se = sqrt(theta.var)) %>%
                    select(cg, rs, theta, theta.se, lodds)

    if(nrow(snp.tab) > 0) {
        write.tab(snp.tab, file = gzfile(snp.out.file))
    }

} else {
    out.tab <- phase.1 %>%
        mutate(theta.2 = NA, theta.se.2 = NA, lodds.2 = NA) %>%
        left_join(.temp.ensg, by = 'ensg') %>%
            left_join(.temp.cg, by = 'cg')
}

out.tab <- out.tab %>% mutate(ld = paste(ld.info[1:3], collapse = '\t'), ld.size = n.snps)

write.tab(out.tab, file = gzfile(z.out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished M -> T\n\n\n')
