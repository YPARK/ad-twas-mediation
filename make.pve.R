#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 8) q()

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/1.ld.gz'
sum.file <- argv[2]                     # e.g., sum.file = 'stat/IGAP/data/hs-lm/1.eqtl_bed.gz'
plink.hdr <- argv[3]                    # e.g., plink.hdr = 'geno/rosmap1709-chr1'
ld.idx <- as.integer(argv[4])           # e.g., ld.idx = 15
gwas.sample.size <- as.numeric(argv[5]) # e.g., gwas.sample.size = 74000
qtl.sample.size <- as.numeric(argv[6])  # e.g., qtl.sample.size = 356
is.eqtl <- as.logical(argv[7])          # e.g., is.eqtl = TRUE

################################################################
qtl.cutoff <- 0
if(length(argv) > 8) {
    qtl.cutoff <- as.numeric(argv[8])
    pve.out.file <- argv[9]
} else {
    pve.out.file <- argv[8]
}

if(file.exists(pve.out.file)) {
    log.msg('Fine %s exists\n\n', pve.out.file)
    q()
}

################################################################
dir.create(dirname(pve.out.file), recursive = TRUE)

source('util.R')
options(stringsAsFactors = FALSE)

library(zqtl)
library(dplyr)
library(methods)
library(broom)
source('mediation.R')

n.snp.cutoff <- 100

## temporary directory
temp.dir <- system('mktemp -d ' %&&% pve.out.file %&&% '-temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

ld.tab <- read.table(ld.file, col.names = c('chr', 'ld.lb', 'ld.ub', 'n.snp.ld', 'n.qtl.ld'))
ld.info <- ld.tab[ld.idx, ]

plink <- subset.plink(ld.info, temp.dir, plink.hdr)
x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))

sum.stat.out <- extract.sum.stat(ld.info, sum.file, x.bim, temp.dir, is.eqtl, qtl.cutoff)

sum.stat <- sum.stat.out$sum.stat
mediators <- sum.stat.out$mediators

## Only work on sufficiently large LD blocks
n.snps <- sum.stat %>% select(snp.loc) %>% unique() %>% nrow()
if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n\n\n', n.snps, n.snp.cutoff)
    system('printf "" | gzip > ' %&&% pve.out.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

zqtl.data <- make.zqtl.data(plink, sum.stat, mediators)

################################################################
## zQTL model for PVE calculation in mediation model
vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 7500, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2, nsample = 10, print.interv = 100,
               nboot = 0, med.finemap = FALSE, weight.y = TRUE, weight.m = TRUE,
               out.residual = TRUE)

med.out <- fit.med.zqtl(zqtl.data$gwas.theta, zqtl.data$gwas.se,
                        zqtl.data$qtl.theta, zqtl.data$qtl.se,
                        X = zqtl.data$X, n = gwas.sample.size,
                        n.med = qtl.sample.size, options = vb.opt)

V.t <- med.out$Vt
Y <- med.out$Y
M <- med.out$M
S.inv <- med.out$S.inv.y
S <- 1/S.inv
S.inv.m <- med.out$S.inv.m
D2 <- med.out$D2
D <- sqrt(med.out$D2)
Vd <- sweep(t(V.t), 2, D, `*`)
W <- sweep(t(V.t), 2, D, `/`)


## theta' R theta
eta.dir <- t(Vd) %*% med.out$param.direct$theta
var.dir <- sum(eta.dir^2)

## estimated true qtl effect
## theta.hat        ~ S R inv(S) (aa * bb)
## inv(S) theta.hat ~ R inv(S) (aa * bb)
##                  ~ R inv(S) 

aa <- sweep(W %*% sweep(M, 1, D, `/`), 1, S, `*`)
bb <- med.out$param.mediated$theta
ab <- aa %*% bb

eta.ab <- t(Vd) %*% ab
var.med <- sum(eta.ab^2)

## residual variance
r.hat <- med.out$resid$theta
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

var.med.vec <- sapply(1:ncol(zqtl.data$qtl.theta), var.each)

out.tab <- melt.effect(med.out$param.mediated, mediators$med.id, 1) %>%
    rename(med.id = Var1) %>%
        mutate(theta.se = sqrt(theta.var)) %>%
            mutate(chr = ld.info[1,1], ld.lb = ld.info[1,2], ld.ub = ld.info[1,3])

best.tab <- find.argmax.snps(sum.stat)

out.tab <- out.tab %>%
    left_join(mediators, by = 'med.id') %>%
        left_join(best.tab, by = 'med.id')

if(is.eqtl) {
    out.tab <- out.tab %>%
        dplyr::select(chr, ld.lb, ld.ub, med.id, theta, theta.se, lodds, hgnc,
                      tss, tes, strand, best.gwas.z, best.gwas.rs, best.gwas.loc,
                      best.qtl.z, best.qtl.rs, best.qtl.loc)
} else {
    out.tab <- out.tab %>%
        dplyr::select(chr, ld.lb, ld.ub, med.id, theta, theta.se, lodds,
                      cg.loc, best.gwas.z, best.gwas.rs, best.gwas.loc,
                      best.qtl.z, best.qtl.rs, best.qtl.loc)
}

out.tab <- data.frame(out.tab, v.med = var.med.vec, v.med.tot = var.med,
                      v.dir = var.dir, v.resid = var.resid)

write.tab(out.tab, file = gzfile(pve.out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished zQTL PVE\n\n\n')
