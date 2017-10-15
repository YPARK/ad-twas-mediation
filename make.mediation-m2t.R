#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) q()

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/4.ld.gz'
mqtl.stat.file <- argv[2]               # e.g., mqtl.stat.file = 'stat/IGAP/data/hs-lm/4.mqtl_bed.gz'
eqtl.stat.file <- argv[3]               # e.g., eqtl.stat.file = 'stat/IGAP/data/hs-lm/4.eqtl_bed.gz'
plink.hdr <- argv[4]                    # e.g., plink.hdr = 'geno/rosmap1709-chr4'
ld.idx <- as.integer(argv[5])           # e.g., ld.idx = 68
out.hdr <- argv[6]                      # e.g., out.hdr = 'temp'

## Fit two models:
##
## 1. G -> M + T -> AD
## 2. G -> M -> T
##
## Save
## 1. mediation effects
## 2. bootstrapped samples
##

dir.create(dirname(out.hdr), recursive = TRUE)

source('util.R')
options(stringsAsFactors = FALSE)

joint.gene.file <- out.hdr %&&% '.gene-mediation.gz'
joint.cpg.file <- out.hdr %&&% '.cpg-mediation.gz'

joint.null.file <- out.hdr %&&% '.joint-null.gz'
m2t.out.file <- out.hdr %&&% '.mediation.gz'
m2t.null.file <- out.hdr %&&% '.m2t-null.gz'

files <- c(joint.gene.file, joint.cpg.file, m2t.out.file, joint.null.file, m2t.null.file)

if(all(sapply(files, file.exists))) {
    log.msg('Files exists : %s\n\n\n', paste(files, collapse = ', '))
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

## best QTL for each gene
best.eqtl.tab <- find.argmax.snps(eqtl.sum.stat)

## Only work on sufficiently large LD blocks
n.snps <- eqtl.sum.stat %>% select(snp.loc) %>% unique() %>% nrow()
if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n', n.snps, n.snp.cutoff)
    system('printf "" | gzip > ' %&&% m2t.out.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

mqtl.stat <- extract.sum.stat(ld.info, mqtl.stat.file, x.bim, temp.dir, is.eqtl = FALSE)
mqtl.sum.stat <- mqtl.stat$sum.stat
mqtl.mediators <- mqtl.stat$mediators

n.snps <- mqtl.sum.stat %>% select(snp.loc) %>% unique() %>% nrow()
if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n', n.snps, n.snp.cutoff)
    system('printf "" | gzip > ' %&&% m2t.out.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

## best QTL for each CpG
best.mqtl.tab <- find.argmax.snps(mqtl.sum.stat)

zqtl.data.eqtl <- make.zqtl.data(plink, eqtl.sum.stat, eqtl.mediators)
zqtl.data.mqtl <- make.zqtl.data(plink, mqtl.sum.stat, mqtl.mediators)

.snps.eqtl <- zqtl.data.eqtl$snps %>% dplyr::select(rs) %>% mutate(eqtl.pos = 1:n())
.snps.mqtl <- zqtl.data.mqtl$snps %>% dplyr::select(rs) %>% mutate(mqtl.pos = 1:n())
.matched <- .snps.eqtl %>% left_join(.snps.mqtl, by = 'rs') %>% na.omit()

zqtl.data.mqtl$X <- zqtl.data.mqtl$X %c% .matched$mqtl.pos
zqtl.data.mqtl$qtl.theta <- zqtl.data.mqtl$qtl.theta %r% .matched$mqtl.pos
zqtl.data.mqtl$qtl.se <- zqtl.data.mqtl$qtl.se %r% .matched$mqtl.pos
zqtl.data.mqtl$gwas.theta <- zqtl.data.mqtl$gwas.theta %r% .matched$mqtl.pos
zqtl.data.mqtl$gwas.se <- zqtl.data.mqtl$gwas.se %r% .matched$mqtl.pos

zqtl.data.eqtl$X <- zqtl.data.eqtl$X %c% .matched$eqtl.pos
zqtl.data.eqtl$qtl.theta <- zqtl.data.eqtl$qtl.theta %r% .matched$eqtl.pos
zqtl.data.eqtl$qtl.se <- zqtl.data.eqtl$qtl.se %r% .matched$eqtl.pos
zqtl.data.eqtl$gwas.theta <- zqtl.data.eqtl$gwas.theta %r% .matched$eqtl.pos
zqtl.data.eqtl$gwas.se <- zqtl.data.eqtl$gwas.se %r% .matched$eqtl.pos

################################################################
## [eQTL, mQTL]
joint.qtl.theta <- cbind(zqtl.data.eqtl$qtl.theta, zqtl.data.mqtl$qtl.theta)
joint.qtl.se <- cbind(zqtl.data.eqtl$qtl.se, zqtl.data.mqtl$qtl.se)
joint.mediators <- c(eqtl.mediators$med.id, mqtl.mediators$med.id)

log.msg('Have all data ready:\n%s mediators\n\n', length(joint.mediators))

## 1. joint model with short number of bootstrap steps (store them all)
vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 5000, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2, nsample = 10, print.interv = 100,
               nboot = 10, bootstrap.method = 2, med.finemap = FALSE)

joint.out <- fit.med.zqtl(effect = zqtl.data.eqtl$gwas.theta,
                          effect.se = zqtl.data.eqtl$gwas.se,
                          effect.m = joint.qtl.theta,
                          effect.m.se = joint.qtl.se,
                          X = zqtl.data.eqtl$X,
                          options = vb.opt)

log.msg('Finished fitting the joint model\n')

joint.null.stat <- as.vector(joint.out$bootstrap$stat.mat)

joint.effect.melt <- melt.effect(joint.out$param.mediated, joint.mediators, 'gwas')


## 2. M -> T (with much reduced number of bootstrap)
vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 5000, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2, nsample = 10, print.interv = 100,
               nboot = 3, bootstrap.method = 2, med.finemap = FALSE)

m2t.out <- fit.med.zqtl(effect = zqtl.data.eqtl$qtl.theta,
                        effect.se = zqtl.data.eqtl$qtl.se,
                        effect.m = zqtl.data.mqtl$qtl.theta,
                        effect.m.se = zqtl.data.mqtl$qtl.se,
                        X = zqtl.data.eqtl$X,
                        options = vb.opt)

log.msg('Finished fitting the M->T model\n')

m2t.null.stat <- as.vector(m2t.out$bootstrap$stat.mat)

m2t.effect.melt <- melt.effect(m2t.out$param.mediated, mqtl.mediators$med.id, eqtl.mediators$med.id)

################################################################
cat(joint.null.stat, file = gzfile(joint.null.file), sep = '\n')
cat(m2t.null.stat, file = gzfile(m2t.null.file), sep = '\n')

gene.out.tab <- joint.effect.melt %>% dplyr::select(-Var2) %>%
    filter(Var1 %in% eqtl.mediators$med.id) %>%
        dplyr::rename(med.id = Var1) %>%
            left_join(eqtl.mediators, by = 'med.id') %>%
                left_join(best.eqtl.tab, by = 'med.id') %>%
                    mutate(chr = ld.info[1,1], ld.lb = ld.info[1,2], ld.ub = ld.info[1,3]) %>%
                        dplyr::select(chr, ld.lb, ld.ub, med.id, theta, theta.var, lodds, hgnc,
                                      tss, tes, strand, best.gwas.z, best.gwas.rs, best.gwas.loc,
                                      best.qtl.z, best.qtl.rs, best.qtl.loc)

cpg.out.tab <- joint.effect.melt %>% dplyr::select(-Var2) %>%
    filter(Var1 %in% mqtl.mediators$med.id) %>%
        dplyr::rename(med.id = Var1) %>%
            left_join(mqtl.mediators, by = 'med.id') %>%
                left_join(best.mqtl.tab, by = 'med.id') %>%
                    mutate(chr = ld.info[1,1], ld.lb = ld.info[1,2], ld.ub = ld.info[1,3]) %>%
                        dplyr::select(chr, ld.lb, ld.ub, med.id, theta, theta.var, lodds,
                                      cg.loc, best.gwas.z, best.gwas.rs, best.gwas.loc,
                                      best.qtl.z, best.qtl.rs, best.qtl.loc)


write.tab(gene.out.tab, file = gzfile(joint.gene.file))
write.tab(cpg.out.tab, file = gzfile(joint.cpg.file))
write.tab(m2t.effect.melt, file = gzfile(m2t.out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished M -> T\n\n\n')
