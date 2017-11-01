#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) q()

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/5.ld.gz'
mqtl.stat.file <- argv[2]               # e.g., mqtl.stat.file = 'stat/IGAP/hic-merged/hs-lm/5.mqtl_bed.gz'
eqtl.stat.file <- argv[3]               # e.g., eqtl.stat.file = 'stat/IGAP/hic-merged/hs-lm/5.eqtl_bed.gz'
plink.hdr <- argv[4]                    # e.g., plink.hdr = 'geno/rosmap1709-chr5'
ld.idx <- as.integer(argv[5])           # e.g., ld.idx = 11

qtl.cutoff <- 0
if(length(argv) > 6) {
    qtl.cutoff <- as.numeric(argv[6])
    out.hdr <- argv[7]
} else {
    out.hdr <- argv[6]
}


## Fit two models:
## G -> M + T -> AD
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
pve.file <- out.hdr %&&% '.pve.gz'
joint.null.file <- out.hdr %&&% '.joint-null.gz'

files <- c(joint.gene.file, joint.cpg.file, joint.null.file, pve.file)

if(all(sapply(files, file.exists))) {
    log.msg('Files exists : %s\n\n\n', paste(files, collapse = ', '))
    q()
}

library(zqtl)
library(dplyr)
library(methods)
library(broom)
library(readr)
source('mediation.R')

n.snp.cutoff <- 10

## temporary directory
temp.dir <- system('mktemp -d ' %&&% out.hdr %&&% '-temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

ld.tab <- read_tsv(ld.file, col_names = c('chr', 'ld.lb', 'ld.ub', 'n.snp.ld', 'n.qtl.ld'))
ld.info <- ld.tab[ld.idx, ]

plink <- subset.plink(ld.info, temp.dir, plink.hdr)
x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))

eqtl.stat <- extract.sum.stat(ld.info, eqtl.stat.file, x.bim, temp.dir, is.eqtl = TRUE, qtl.cutoff)
eqtl.sum.stat <- eqtl.stat$sum.stat
eqtl.mediators <- eqtl.stat$mediators

## best QTL for each gene
best.eqtl.tab <- find.argmax.snps(eqtl.sum.stat)

## Only work on sufficiently large LD blocks
n.snps <- eqtl.sum.stat %>% select(snp.loc) %>% unique() %>% nrow()

if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n', n.snps, n.snp.cutoff)
    system('printf "" | gzip > ' %&&% joint.gene.file)
    system('printf "" | gzip > ' %&&% joint.cpg.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

mqtl.stat <- extract.sum.stat(ld.info, mqtl.stat.file, x.bim, temp.dir, is.eqtl = FALSE, qtl.cutoff)
mqtl.sum.stat <- mqtl.stat$sum.stat
mqtl.mediators <- mqtl.stat$mediators

n.snps <- mqtl.sum.stat %>% select(snp.loc) %>% unique() %>% nrow()

if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n', n.snps, n.snp.cutoff)
    system('printf "" | gzip > ' %&&% joint.gene.file)
    system('printf "" | gzip > ' %&&% joint.cpg.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

## best QTL for each CpG
best.mqtl.tab <- find.argmax.snps(mqtl.sum.stat)

zqtl.data.eqtl <- make.zqtl.data(plink, eqtl.sum.stat, eqtl.mediators)
zqtl.data.mqtl <- make.zqtl.data(plink, mqtl.sum.stat, mqtl.mediators)

.snps.eqtl <- zqtl.data.eqtl$snps %>% dplyr::select(rs) %>% mutate(eqtl.pos = 1:n())
.snps.mqtl <- zqtl.data.mqtl$snps %>% dplyr::select(rs) %>% mutate(mqtl.pos = 1:n())

.union <- .snps.eqtl %>% dplyr::full_join(.snps.mqtl, by = 'rs')

idx.e <- .union[, 'eqtl.pos'] %>% na.omit()
idx.m <- (.union %>% filter(is.na(eqtl.pos)))[, 'mqtl.pos'] %>% na.omit()

xx <- cbind(zqtl.data.eqtl$X %c% idx.e, zqtl.data.mqtl$X %c% idx.m) %>%
    scale()

n.snps <- ncol(xx)
n.genes <- nrow(eqtl.mediators)
n.cpgs <- nrow(mqtl.mediators)

joint.qtl.theta <- matrix(data = NA , nrow = n.snps, ncol = n.genes + n.cpgs)
joint.qtl.se <- matrix(data = NA , nrow = n.snps, ncol = n.genes + n.cpgs)
joint.mediators <- c(eqtl.mediators$med.id, mqtl.mediators$med.id)

joint.qtl.theta[!is.na(.union$eqtl.pos), 1:n.genes] <-
    zqtl.data.eqtl$qtl.theta %r% (.union$eqtl.pos %>% na.omit())

joint.qtl.se[!is.na(.union$eqtl.pos), 1:n.genes] <-
    zqtl.data.eqtl$qtl.se %r% (.union$eqtl.pos %>% na.omit())

joint.qtl.theta[!is.na(.union$mqtl.pos), -(1:n.genes)] <-
    zqtl.data.mqtl$qtl.theta %r% (.union$mqtl.pos %>% na.omit())

joint.qtl.se[!is.na(.union$mqtl.pos), -(1:n.genes)] <-
    zqtl.data.mqtl$qtl.se %r% (.union$mqtl.pos %>% na.omit())


joint.gwas <- rbind(zqtl.data.eqtl$gwas.theta %r% idx.e,
                    zqtl.data.mqtl$gwas.theta %r% idx.m)

joint.gwas.se <- rbind(zqtl.data.eqtl$gwas.se %r% idx.e,
                       zqtl.data.mqtl$gwas.se %r% idx.m)

log.msg('Have all data ready:\n%s mediators\n\n', length(joint.mediators))

if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n', n.snps, n.snp.cutoff)
    system('printf "" | gzip > ' %&&% joint.gene.file)
    system('printf "" | gzip > ' %&&% joint.cpg.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

## Joint model with short number of bootstrap steps (store them all)
vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 5000, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2, nsample = 10, print.interv = 100,
               nboot = 10, bootstrap.method = 2, med.finemap = FALSE,               
               out.residual = TRUE)

joint.out <- fit.med.zqtl(effect = joint.gwas,
                          effect.se = joint.gwas.se,
                          effect.m = joint.qtl.theta,
                          effect.m.se = joint.qtl.se,
                          X = xx,
                          options = vb.opt)

log.msg('Finished fitting the joint model\n')

joint.effect.melt <- melt.effect(joint.out$param.mediated, joint.mediators, 'gwas')

chr <- ld.info$chr[1]
ld.lb <- ld.info$ld.lb[1]
ld.ub <- ld.info$ld.ub[1]

gene.out.tab <- joint.effect.melt %>% dplyr::select(-Var2) %>%
    filter(Var1 %in% eqtl.mediators$med.id) %>%
        dplyr::rename(med.id = Var1) %>%
            left_join(eqtl.mediators, by = 'med.id') %>%
                left_join(best.eqtl.tab, by = 'med.id') %>%
                    mutate(chr, ld.lb, ld.ub) %>%
                        dplyr::select(chr, ld.lb, ld.ub, med.id, theta, theta.var, lodds, hgnc,
                                      tss, tes, strand, best.gwas.z, best.gwas.rs, best.gwas.loc,
                                      best.qtl.z, best.qtl.rs, best.qtl.loc, n.qtls)

cpg.out.tab <- joint.effect.melt %>% dplyr::select(-Var2) %>%
    filter(Var1 %in% mqtl.mediators$med.id) %>%
        dplyr::rename(med.id = Var1) %>%
            left_join(mqtl.mediators, by = 'med.id') %>%
                left_join(best.mqtl.tab, by = 'med.id') %>%
                    mutate(chr, ld.lb, ld.ub) %>%
                        dplyr::select(chr, ld.lb, ld.ub, med.id, theta, theta.var, lodds,
                                      cg.loc, best.gwas.z, best.gwas.rs, best.gwas.loc,
                                      best.qtl.z, best.qtl.rs, best.qtl.loc, n.qtls)


################################################################
## PVE calculation separate out genes and CpGs
##
## variance = theta' R theta
##

V.t <- joint.out$Vt
Y <- joint.out$Y
M <- joint.out$M
S.inv <- joint.out$S.inv.y
S <- 1/S.inv
S.inv.m <- joint.out$S.inv.m
D2 <- joint.out$D2
D <- sqrt(joint.out$D2)
Vd <- sweep(t(V.t), 2, D, `*`)
W <- sweep(t(V.t), 2, D, `/`)

## theta' R theta
eta.dir <- t(Vd) %*% joint.out$param.direct$theta
var.dir <- sum(eta.dir^2)

## estimated true qtl effect
## theta.hat        ~ S R inv(S) (aa * bb)
## inv(S) theta.hat ~ R inv(S) (aa * bb)
##                  ~ R inv(S) 

aa <- sweep(W %*% sweep(M, 1, D, `/`), 1, S, `*`)
bb <- joint.out$param.mediated$theta
ab <- aa %*% bb

eta.ab <- t(Vd) %*% ab
var.med <- sum(eta.ab^2)

## gene and Cpg separately
ab.e <- (aa %c% 1:n.genes) %*% (bb %r% 1:n.genes)
ab.m <- (aa %c% (-(1:n.genes))) %*% (bb %r% (-(1:n.genes)))

var.med.e <- sum((t(Vd) %*% ab.e)^2)
var.med.m <- sum((t(Vd) %*% ab.m)^2)

## residual variance
r.hat <- joint.out$resid$theta
rr <- sweep(W %*% (t(W) %*% r.hat), 1, S, `*`)
eta.r <- t(Vd) %*% rr
var.resid <- sum(eta.r^2)

## each mediation effect
var.each <- function(k) {
    ab.k <- (aa %c% k) %*% (bb %r% k)
    eta.ab.k <- t(Vd) %*% ab.k
    var.k <- sum(eta.ab.k^2)
    return(var.k)
}

var.med.vec <- sapply(1:(n.genes + n.cpgs), var.each)

pve.out.tab <- data.frame(chr,
                          ld.lb,
                          ld.ub,
                          med.id = c(eqtl.mediators$med.id, mqtl.mediators$med.id),
                          var.med.each = var.med.vec,
                          var.med.e,
                          var.med.m,
                          var.med,
                          var.dir,
                          var.resid)

################################################################
joint.null.stat <- as.vector(joint.out$bootstrap$stat.mat)

cat(joint.null.stat, file = gzfile(joint.null.file), sep = '\n')
write.tab(gene.out.tab, file = gzfile(joint.gene.file))
write.tab(cpg.out.tab, file = gzfile(joint.cpg.file))
write.tab(pve.out.tab, file = gzfile(pve.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished M -> T\n\n\n')
