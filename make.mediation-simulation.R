#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 11) q()

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/1.ld.gz'
sum.file <- argv[2]                     # e.g., sum.file = 'stat/IGAP/data/hs-lm/1.eqtl_bed.gz'
plink.hdr <- argv[3]                    # e.g., plink.hdr = 'geno/rosmap1709-chr1'
ld.idx <- as.integer(argv[4])           # e.g., ld.idx = 15
is.eqtl <- as.logical(argv[5])          # e.g., is.eqtl = TRUE

h2.qtl <- as.numeric(argv[6])           # h2.qtl <- 0.3
h2.med <- as.numeric(argv[7])           # h2.med <- 0.5
n.causal.med <- as.integer(argv[8])     # n.causal.med <- 2
n.causal.qtl <- as.integer(argv[9])     # n.causal.qtl <- 3
n.causal.direct <- as.integer(argv[10]) # e.g, n.causal.direct <- 1

out.file <- argv[11]

dir.create(dirname(out.hdr), recursive = TRUE)

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
temp.dir <- system('mktemp -d ' %&&% out.hdr %&&% '-temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

ld.tab <- read.table(ld.file, col.names = c('chr', 'ld.lb', 'ld.ub', 'n.snp.ld', 'n.qtl.ld'))
ld.info <- ld.tab[ld.idx, ]

plink <- subset.plink(ld.info, temp.dir, plink.hdr)
x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))

sum.stat.out <- extract.sum.stat(ld.info, sum.file, x.bim, temp.dir, is.eqtl)
sum.stat <- sum.stat.out$sum.stat
mediators <- sum.stat.out$mediators

################################################################
## generate synthetic data
X <- plink$BED %>% scale()
X[is.na(X)] <- 0
p <- ncol(X)
n <- nrow(X)
n.med <- nrow(mediators)

## make sure numbers
n.causal.med <- pmin(n.causal.med, n.med)
n.causal.qtl <- pmin(n.causal.qtl, p)
n.causal.direct <- pmin(n.causal.direct, p)
theta.pop = c(-1, 1, -0.5, 0.5)

## generate mediators
causal.qtls <- matrix(sample(p, n.causal.qtl * n.med), nrow = n.causal.qtl)
theta.qtl <- .sample(n.causal.qtl, n.med, pop = theta.pop)

med.mat <- do.call(cbind, lapply(1:ncol(causal.qtls),
                                 function(j) (X %c% causal.qtls[,j]) %*% (theta.qtl %c% j)))

v0 <- apply(med.mat, 2, var) * (1/h2.qtl - 1)
med.mat <- med.mat + sweep(.rnorm(n, n.med), 2, sqrt(v0), `*`)

## generate direct effects
causal.direct <- sample(p, n.causal.direct)
theta.dir <- .sample(n.causal.direct, 1, pop = theta.pop)

causal.med <- sample(n.med, n.causal.med)
theta.med <- .sample(n.causal.med, 1, theta.pop)

y <- (med.mat %c% causal.med) %*% theta.med +
    (X %c% causal.direct) %*% theta.dir

y <- y + sweep(.rnorm(n, 1), 2, sqrt(var(y) * (1/h2.med - 1)), `*`)

################################################################
## compute summary statistics
gwas.stat <- get.marginal.qtl(X, y, .melt = FALSE)
qtl.stat <- get.marginal.qtl(X, med.mat, .melt = FALSE)

## Just run the zQTL w/o bootstrapping w/o finemap
vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 5000, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2, nsample = 10, print.interv = 100,
               nboot = 0, med.finemap = FALSE,
               weight.y = TRUE, weight.m = TRUE)

z.out <- fit.med.zqtl(gwas.stat$beta, gwas.stat$beta.se,
                      qtl.stat$beta, qtl.stat$beta.se,
                      X = X, options = vb.opt)

z.out.tab <- melt.effect(z.out$param.mediated, 1:n.med, 'sim') %>%
    mutate(theta.se = sqrt(theta.var), gene = Var1) %>%
        dplyr::select(-Var2, -Var1, -theta.var) %>%
            dplyr::select(gene, theta, theta.se, lodds) %>%
                mutate(causal = ifelse(gene %in% causal.med, 1, 0)) %>%
                    mutate(h2.qtl, h2.med, p, n.causal.qtl, n.causal.med, n.causal.direct)

################################################################
## run MR egger for comparison
gwas.stat <- get.marginal.qtl(X, y, .melt = TRUE) %>%    
    mutate(gwas.theta = beta, gwas.se = beta.z/beta) %>%
        dplyr::select(-Var2, -beta, -beta.z)

colnames(med.mat) <- as.character(1:n.med)

qtl.stat <- get.marginal.qtl(X, med.mat, .melt = TRUE) %>%
    mutate(qtl.theta = beta, qtl.se = beta.z/beta) %>%
        dplyr::select(-beta, -beta.z)

egger.tab <- qtl.stat %>% left_join(gwas.stat, by = 'snp') %>%
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
    dplyr::select(gene, causal, h2.qtl, h2.med,
                  p, n.causal.qtl, n.causal.med, n.causal.med, n.causal.direct,
                  theta, theta.se, lodds,
                  egger.theta, egger.se, egger.t, egger.p.val)

write.tab(out.tab, file = gzfile(out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished zQTL\n\n\n')
