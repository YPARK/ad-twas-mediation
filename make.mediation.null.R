#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

ld.file <- argv[1]                      # e.g., 'stat/IGAP/ld/hs-fqtl/19.obs_ld.gz'
sum.file <- argv[2]                     # e.g., 'stat/IGAP/hs-fqtl/19.obs_bed.gz'
plink.hdr <- argv[3]                    # e.g., '1kg/chr19' # geno/rosmap1709-chr19'
ld.idx <- as.integer(argv[4])           # e.g., 117
gwas.sample.size <- as.numeric(argv[5]) # e.g., 74000
qtl.sample.size <- as.numeric(argv[6])  # e.g., 300
z.out.file <- argv[7]                   # e.g., 'temp'

dir.create(dirname(z.out.file), recursive = TRUE)

source('util.R')
options(stringsAsFactors = FALSE)

if(file.exists(z.out.file)) {
    log.msg('File exists : %s\n\n\n', z.out.file)
    q()
}

library(zqtl)
library(dplyr)
library(methods)
library(broom)
source('mediation.R')

################################################################
ld.tab <- read.table(ld.file, col.names = c('chr', 'ld.lb', 'ld.ub', 'n.snp.ld', 'n.qtl.ld'))
ld.info <- ld.tab[ld.idx, ]

################################################################
## temporary directory
temp.dir <- system('mktemp -d ' %&&% z.out.file %&&% 'temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

plink <- subset.plink(ld.info, temp.dir, plink.hdr)
x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))

sum.stat.out <- extract.sum.stat(ld.info, sum.file, x.bim, temp.dir)
sum.stat <- sum.stat.out$sum.stat
genes <- sum.stat.out$genes

zqtl.data <- make.zqtl.data(plink, sum.stat, genes)

################################################################
## Generate null data
vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e3,
               vbiter = 7500, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2,
               print.interv = 500, nboot = 0)

z.marg <- fit.zqtl(zqtl.data$gwas.theta, zqtl.data$gwas.se, X = zqtl.data$X,
                  n = gwas.sample.size, options = vb.opt)

svd.out <- take.ld.svd(zqtl.data$X)
V <- sweep(t(svd.out$V.t), 2, svd.out$D, `*`)
R <- V %*% t(V)

n <- nrow(z.marg$param$theta)
p <- ncol(z.marg$param$theta)
theta.rand <- z.marg$param$theta + .rnorm(n, p) * sqrt(z.marg$param$theta.var)
rand.idx <- sample(nrow(theta.rand))
theta.rand <- theta.rand %r% rand.idx
gwas.se <- zqtl.data$gwas.se %r% rand.idx
gwas.theta.z <- sweep(R, 2, gwas.se, `/`) %*% theta.rand
gwas.theta <- gwas.theta.z * gwas.se

################################################################
## Estimate mediation model on the generated data

z.out <- fit.med.zqtl(gwas.theta, gwas.se,
                      zqtl.data$qtl.theta, zqtl.data$qtl.se,
                      X = zqtl.data$X, n = gwas.sample.size,
                      n.med = qtl.sample.size, options = vb.opt)

out.tab <- melt.effect(z.out$param.mediated, genes$ensg, 'permuted') %>%
    dplyr::select(-Var2)

write.tab(out.tab, file = gzfile(z.out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished zQTL\n\n\n')
