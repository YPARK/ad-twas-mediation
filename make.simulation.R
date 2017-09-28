#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 12) q()

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/19.ld.gz'
sum.file <- argv[2]                     # e.g., sum.file = 'stat/IGAP/data/hs-lm/19.eqtl_bed.gz'
plink.hdr <- argv[3]                    # e.g., plink.hdr = 'geno/rosmap1709-chr19'
ld.idx <- as.integer(argv[4])           # e.g., ld.idx = 117
is.eqtl <- as.logical(argv[5])          # e.g., is.eqtl = TRUE

pve.qtl <- as.numeric(argv[6])           # e.g., pve.qtl <- 0.1
pve.med <- as.numeric(argv[7])           # e.g., pve.med <- 0.3
pve.dir <- as.numeric(argv[8])           # e.g., pve.dir <- 0.1

n.causal.med <- as.integer(argv[9])     # e.g., n.causal.med <- 1
n.causal.qtl <- as.integer(argv[10])    # e.g., n.causal.qtl <- 3
n.causal.direct <- as.integer(argv[11]) # e.g., n.causal.direct <- 3

out.file <- argv[12]

dir.create(dirname(out.file), recursive = TRUE)

source('util.R')
options(stringsAsFactors = FALSE)

if(file.exists(out.file)) {
    log.msg('File exists : %s\n\n\n', out.file)
    q()
}

stopifnot(pve.dir + pve.med <= 1)
stopifnot(pve.qtl <= 1)

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

sum.stat.out <- extract.sum.stat(ld.info, sum.file, x.bim, temp.dir, is.eqtl)
mediators <- sum.stat.out$mediators
sum.stat <- sum.stat.out$sum.stat
n.snps <- sum.stat %>% select(snp.loc) %>% unique() %>% nrow()

log.msg('num SNPs = %d\n\n', n.snps)

## Only work on sufficiently large LD blocks
n.snps <- sum.stat %>% select(snp.loc) %>% unique() %>% nrow()
n.snp.cutoff <- 500
if(n.snps < n.snp.cutoff) {
    log.msg('This LD block is too small : %d SNPs < %d\n\n\n', n.snps, n.snp.cutoff)
    system('printf "" | gzip > ' %&&% out.file)
    system('rm -r ' %&&% temp.dir)
    q()
}

################################################################
## generate synthetic data -- following Mancuso et al. 2017
X <- plink$BED %>% scale()
X[is.na(X)] <- 0
p <- ncol(X)
n <- nrow(X)
n.med <- nrow(mediators)

## make sure numbers
n.causal.med <- pmin(n.causal.med, n.med)
n.causal.qtl <- pmin(n.causal.qtl, p)
n.causal.direct <- pmin(n.causal.direct, p)

## (1) generate mediators : theta ~ N(0, pve/n.qtl.per.gene)
##                          error ~ N(0, (1-pve))
causal.qtls <- matrix(sample(p, n.causal.qtl * n.med), nrow = n.causal.qtl)
theta.qtl <- .rnorm(n.causal.qtl, n.med) * sqrt(pve.qtl / n.causal.qtl)

sim.med <- function(j) {
    (X %c% causal.qtls[, j]) %*% (theta.qtl %c% j)
}

M <- sapply(1:ncol(causal.qtls), sim.med)
M <- M + sweep(.rnorm(n, n.med), 2, sqrt(1 - pve.qtl), `*`)

## (2) generate direct and mediation effects
causal.direct <- sample(p, n.causal.direct)
theta.dir <- .rnorm(n.causal.direct, 1) * sqrt(pve.dir / n.causal.direct)

causal.med <- sample(n.med, n.causal.med)
theta.med <- .rnorm(n.causal.med, 1) * sqrt(pve.med / n.causal.med)

pve.error <- 1 - pve.dir - pve.med

y <- M %c% causal.med %*% theta.med +
    X %c% causal.direct %*% theta.dir +
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
                    mutate(pve.med, pve.dir, pve.qtl, p, n.causal.qtl, n.causal.med, n.causal.direct)

log.msg('Finished zQTL\n\n')

################################################################
## run gene by gene TWAS by Mancuso et al. 2017

svd.out <- take.ld.svd(plink$BED, options = vb.opt)
V.t.inv <- sweep(svd.out$V.t, 1, svd.out$D, `/`)
LD.inv <- t(V.t.inv) %*% V.t.inv
V.t <- sweep(svd.out$V.t, 1, svd.out$D, `*`)
LD <- t(V.t) %*% V.t

gwas.z <- gwas.stat$beta / gwas.stat$beta.se

take.twas <- function(g) {
    qtl.z <- (qtl.stat$beta %c% g) / (qtl.stat$beta.se %c% g)
    qtl.z.poly <- LD.inv %*% qtl.z
    num <- t(qtl.z.poly) %*% gwas.z
    denom <- t(qtl.z.poly) %*% LD %*% qtl.z.poly
    return(as.numeric(num/sqrt(denom + 1e-16)))
}

twas.tab <- data.frame(gene = as.character(1:n.med), twas = sapply(1:n.med, take.twas))

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
        dplyr::select(gene, causal, pve.med, pve.dir, pve.qtl,
                      p, n.causal.qtl, n.causal.med, n.causal.med, n.causal.direct,
                      theta, theta.se, lodds,
                      egger.t, twas)                      

write.tab(out.tab, file = gzfile(out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished zQTL\n\n\n')
