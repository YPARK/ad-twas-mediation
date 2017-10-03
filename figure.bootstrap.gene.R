#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

gene.idx <- as.integer(argv[1])

options(stringsAsFators = FALSE)
source('util.R')
source('mediation.R')
library(readr)
library(zqtl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(latex2exp)
library(Hmisc)
source('figure.util.R')

gene.tab.file <- 'tables/bootstrap_gene_significant.txt.gz'
gene.tab <- read.table(gene.tab.file, header = TRUE) %>% na.omit()

## Figure 2. gene by gene scatter plots
out.dir <- 'figures/genes/'
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

draw.gene <- function(gene.idx, gene.tab, temp.dir, out.hdr) {

    file.name <- out.hdr %&&% sprintf('gene_%d_%05d_%s.pdf', gene.tab[gene.idx, ]$chr,
                                      gene.idx, gene.tab[gene.idx, ]$hgnc)

    if(file.exists(file.name)) { return(NULL) }

    plink.hdr <- 'geno/rosmap1709-chr'
    sum.file.dir <- 'stat/IGAP/data/hs-lm/'

    ld.info <- as.vector(gene.tab[gene.idx, 1:3])
    chr <- ld.info[1]
    plink <- subset.plink(ld.info, temp.dir, plink.hdr %&&% chr)
    x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))

    sum.stat.obj <- extract.sum.stat(ld.info = ld.info,
                                     sum.file = sum.file.dir %&&% chr %&&% '.eqtl_bed.gz',
                                     x.bim = x.bim,
                                     temp.dir = temp.dir,
                                     is.eqtl = TRUE)

    sum.stat <- sum.stat.obj$sum.stat
    genes <- sum.stat.obj$mediators

    zqtl.data <- make.zqtl.data(plink, sum.stat, genes)

    y.pos <- which(genes$med.id == as.character(gene.tab[gene.idx, 'ensg']))

    param <- gene.tab[gene.idx, ] %>% select(theta, theta.se)

    ## re-estimate zQTL model
    vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 10000, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -0.1, nsample = 10, print.interv = 10,
               nboot = 0, med.finemap = FALSE, weight.y = TRUE, weight.m = TRUE)

    z.out <- fit.med.zqtl(zqtl.data$gwas.theta,
                          zqtl.data$gwas.se,
                          zqtl.data$qtl.theta %c% y.pos,
                          zqtl.data$qtl.se %c% y.pos,
                          X = zqtl.data$X,
                          options = vb.opt)

    z.x <- (zqtl.data$qtl.theta / zqtl.data$qtl.se) %c% y.pos
    z.y <- (zqtl.data$gwas.theta / zqtl.data$gwas.se)
    valid <- which(!is.na(z.x) & !is.na(z.y))

    svd.out <- take.ld.svd(zqtl.data$X, options = list(do.stdize = TRUE, eigen.tol = 1e-2))
    
    V.t <- svd.out$V.t

    zz.x <- (V.t %c% valid) %*% (z.x %r% valid)
    zz.y <- V.t %*% z.y

    ww <- 1/svd.out$D^2
    .wtd.std <- function(x) (x - wtd.mean(x, ww)) / sqrt(wtd.var(x, ww))

    zz.y.direct <- z.out$param.direct$theta / zqtl.data$gwas.se
    zz.y.direct <- sweep(V.t %*% zz.y.direct, 1, svd.out$D^2, `*`)

    ## rr.y <- matrix(lm(zz.y ~ zz.y.direct, weights = ww)$residuals, ncol = 1)
    rr.y <- zz.y - zz.y.direct

    p0 <- ggplot(data.frame(QTL = z.x[valid], GWAS = z.y[valid])) + theme_bw() +
        theme(legend.position = 'bottom', legend.background = element_blank(),
              title = element_text(size=8)) +
            geom_point(aes(x = QTL, y = GWAS), alpha = 0.5, size = 1) +
                ggtitle(genes[y.pos, ]$med.id)

    .df <- data.frame(QTL = .wtd.std(zz.x), GWAS = .wtd.std(zz.y), w = ww)

    .thm <- theme_bw() +
        theme(legend.position = 'none', axis.title.y=element_blank(), title = element_text(size=8))

    p1 <- ggplot(.df) + .thm +
        geom_point(aes(x = QTL, y = GWAS, size = w, alpha = w)) +
            scale_size_continuous(range = c(0, 1))

    p1 <- p1 + coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))
    
    .df <- data.frame(QTL = zz.x, GWAS = rr.y, w = ww)

    .wlm <- lm(GWAS ~ QTL + 1, data = .df, weights = w)
    .slope <- signif(coefficients(.wlm)[2], 2)
    .inter <- signif(coefficients(.wlm)[1], 2)

    .title <- sprintf('%.2f %s + %.2f', .slope, gene.tab[gene.idx, ]$hgnc, .inter)

    .df <- data.frame(QTL = .wtd.std(zz.x), GWAS = .wtd.std(rr.y), w = ww)

    .thm <- theme_bw() +
        theme(axis.title.y=element_blank(), title = element_text(size=8))

    p2 <- ggplot(.df) + .thm +
        geom_abline(slope = .slope, intercept = .inter, color = 'red', size = .5) +
            ggtitle(.title) +
                scale_size_continuous(name = TeX('$1/\\lambda$'), range = c(0, 1)) +
                    geom_point(aes(x = QTL, y = GWAS, size = w, alpha = w)) +
                        scale_alpha_continuous(guide = FALSE)
    
    p2 <- p2 + coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3))

    .df <- data.frame(llik = z.out$llik[, 1]) %>% mutate(iter = 1:n())

    p3 <- ggplot(.df) + theme_classic() + geom_line(aes(x=iter, y=llik)) +
        theme(axis.title.y=element_blank()) +
            xlab('iteration (x10)') + ggtitle('log-likelihood')

    pdf(file = file.name, width = 8, height = 2.5, useDingbats = FALSE)
    print(grid.hcat(list(p0, p1, p2, p3), widths = c(.9, .9, 1.3, .9)))
    dev.off()

}

## temporary directory
temp.dir <- system('mktemp -d ' %&&% out.dir %&&% '/temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

names(gene.tab)[2] <- 'ld.lb'
names(gene.tab)[3] <- 'ld.ub'

draw.gene(gene.idx, gene.tab, temp.dir, out.hdr = out.dir %&&% '/')

system('rm -r ' %&&% temp.dir)
log.msg('Finished figure drawing\n\n\n')
