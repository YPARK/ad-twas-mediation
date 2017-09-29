#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 10) q()

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/19.ld.gz'
plink.hdr <- argv[2]                    # e.g., plink.hdr = 'geno/rosmap1709-chr19'
ld.idx <- as.integer(argv[3])           # e.g., ld.idx = 1

pve.qtl <- as.numeric(argv[4])           # e.g., pve.qtl <- 0.1
pve.med <- as.numeric(argv[5])           # e.g., pve.med <- 0.3
pve.dir <- as.numeric(argv[6])           # e.g., pve.dir <- 0.1

n.causal.med <- as.integer(argv[7])     # e.g., n.causal.med <- 1
n.causal.qtl <- as.integer(argv[8])     # e.g., n.causal.qtl <- 3
n.causal.direct <- as.integer(argv[9]) # e.g., n.causal.direct <- 3

out.file <- argv[10]

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

## temporary directory
temp.dir <- system('mktemp -d ' %&&% out.file %&&% '-temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

ld.tab <- read.table(ld.file, col.names = c('chr', 'ld.lb', 'ld.ub', 'n.snp.ld', 'n.qtl.ld'))
ld.info <- ld.tab[ld.idx, ]

plink <- subset.plink(ld.info, temp.dir, plink.hdr)
x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))
n.snps <- nrow(x.bim)
log.msg('num SNPs = %d\n\n', n.snps)

n.med <- ceiling(n.snps / 20)

## Only work on sufficiently large LD blocks
n.snp.cutoff <- 1000
if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n\n\n', n.snps, n.snp.cutoff)
    system('printf "" | gzip > ' %&&% out.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

################################################################
## generate synthetic data
X <- plink$BED %>% scale()
X[is.na(X)] <- 0
p <- ncol(X)
n <- nrow(X)

## make sure numbers
n.causal.med <- pmin(n.causal.med, n.med)
n.causal.qtl <- pmin(n.causal.qtl, p)
n.causal.direct <- pmin(n.causal.direct, p)

## (1) generate mediators
causal.qtls <- matrix(sample(p, n.causal.qtl * n.med, TRUE), nrow = n.causal.qtl)
theta.qtl <- .rnorm(n.causal.qtl, n.med) * sqrt(pve.qtl / n.causal.qtl)

sim.med <- function(j) {
    (X %c% causal.qtls[, j]) %*% (theta.qtl %c% j)
}

M <- sapply(1:ncol(causal.qtls), sim.med)
M <- M + sweep(.rnorm(n, n.med), 2, sqrt(1 - pve.qtl), `*`)

## genuine mediators
causal.med <- sample(n.med, n.causal.med, TRUE)

## generate pleiotropic direct effects
.temp <- unique(as.vector(causal.qtls %c% (-causal.med)))
pleiotropy.qtls <- unique(sample(.temp, n.causal.direct, TRUE))
pleiotropy.med <- apply(causal.qtls, 2, function(q) sum(q %in% pleiotropy.qtls))
pleiotropy.med <- which(pleiotropy.med > 0)
pleiotropy.med <- setdiff(pleiotropy.med, causal.med)

## (2) generate direct and mediation effects
n.causal.direct <- length(pleiotropy.qtls)
theta.dir <- .rnorm(n.causal.direct, 1) * sqrt(pve.dir / n.causal.direct)
theta.med <- .rnorm(n.causal.med, 1) * sqrt(pve.med / n.causal.med)

pve.error <- 1 - pve.dir - pve.med

y <- M %c% causal.med %*% theta.med +
    X %c% pleiotropy.qtls %*% theta.dir +
        sweep(.rnorm(n, 1), 2, sqrt(pve.error), `*`)

################################################################
## compute summary statistics
gwas.stat <- get.marginal.qtl(X, y, .melt = FALSE)
qtl.stat <- get.marginal.qtl(X, M, .melt = FALSE)

################################################################
## Just run the zQTL w/o bootstrapping w/o finemap
vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 5000, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = 0, nsample = 10, print.interv = 100,
               nboot = 0, med.finemap = FALSE, weight.y = TRUE, weight.m = TRUE)

z.out <- fit.med.zqtl(effect = gwas.stat$beta, effect.se = gwas.stat$beta.se,
                      effect.m = qtl.stat$beta, effect.m.se = qtl.stat$beta.se,
                      X = X, options = vb.opt)

z.out.tab <- melt.effect(z.out$param.mediated, 1:ncol(M), 'sim') %>%
    mutate(theta.se = sqrt(theta.var), gene = Var1) %>%
        dplyr::select(-Var2, -Var1, -theta.var) %>%
            dplyr::select(gene, theta, theta.se, lodds) %>%
                mutate(causal = ifelse(gene %in% causal.med, 1, 0)) %>%
                    mutate(pleiotropy = ifelse(gene %in% pleiotropy.med, -1, 0)) %>%
                        mutate(pve.med, pve.dir, pve.qtl, p, n.med, n.causal.qtl, n.causal.med, n.causal.direct)

log.msg('Finished zQTL\n\n')

################################################################
## run gene by gene TWAS by Mancuso et al. 2017

svd.out <- take.ld.svd(plink$BED, options = vb.opt)
W.t <- sweep(svd.out$V.t, 1, svd.out$D, `/`)
## V.t <- sweep(svd.out$V.t, 1, svd.out$D, `*`)
## LD.inv <- t(W.t) %*% W.t
## LD <- t(V.t) %*% V.t

gwas.z <- gwas.stat$beta / gwas.stat$beta.se

take.twas <- function(g) {

    qtl.z <- (qtl.stat$beta %c% g) / (qtl.stat$beta.se %c% g)

    ## This is too memory-intensive
    ## qtl.z.poly <- LD.inv %*% qtl.z
    ## num <- t(qtl.z.poly) %*% gwas.z
    ## denom <- t(qtl.z.poly) %*% LD %*% qtl.z.poly

    num <- t(W.t %*% gwas.z) %*% (W.t %*% qtl.z)
    denom <- t(W.t %*% qtl.z) %*% (W.t %*% qtl.z)
    log.msg('TWAS finished [%d / %d]\n', g, n.med)

    return(signif(as.numeric(num/sqrt(denom + 1e-16)), 4))
}

twas.tab <- data.frame(gene = as.character(1:n.med), twas = sapply(1:n.med, take.twas))

log.msg('TWAS Finished\n\n')

################################################################
## run MR egger for comparison

gwas.egger.stat <- get.marginal.qtl(X, y, .melt = TRUE) %>%    
    mutate(gwas.theta = beta, gwas.se = beta.z/beta) %>%
        dplyr::select(-gene, -beta, -beta.z)

colnames(M) <- as.character(1:n.med)

qtl.egger.stat <- get.marginal.qtl(X, M, .melt = TRUE) %>%
    mutate(qtl.theta = beta, qtl.se = beta.z/beta) %>%
        dplyr::select(-beta, -beta.z)

egger.tab <- qtl.egger.stat %>% left_join(gwas.egger.stat, by = 'snp') %>%
    group_by(gene) %>%
        do(egger = mr.egger(.)) %>%
            tidy(egger) %>%
                filter(.rownames != 'b.qtl') %>%
                    dplyr::select(-.rownames)

egger.tab <- egger.tab %>% as.data.frame() %>% mutate(gene = as.character(gene)) %>%
    rename(egger.theta=Estimate, egger.se=Std..Error, egger.t=t.value, egger.p.val=Pr...t..) %>%
        mutate(egger.p.val = pmin(-log10(egger.p.val), 100)) %>%
            mutate(egger.p.val = signif(egger.p.val, 4))

out.tab <- z.out.tab %>% left_join(egger.tab, by = 'gene') %>%
    left_join(twas.tab, by = 'gene') %>%
        dplyr::select(gene, causal, pleiotropy, pve.med, pve.dir, pve.qtl,
                      p, n.med, n.causal.qtl, n.causal.med, n.causal.med, n.causal.direct,
                      theta, theta.se, lodds,
                      egger.t, twas)                      

write.tab(out.tab, file = gzfile(out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished zQTL\n\n\n')
