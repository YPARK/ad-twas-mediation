#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

gene.idx <- as.integer(argv[1])

options(stringsAsFators = FALSE)
source('util.R')
source('mediation.R')
source('figure.util.R')
library(readr)
library(zqtl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(latex2exp)
library(reshape2)

ld.tab.file <- 'tables/bootstrap_ld_strict.txt.gz'
ld.tab <- read.table(ld.tab.file, header = TRUE)
ld.key <- ld.tab %>% select(chr, ld.lb, ld.ub) %>% unique()

## Draw all genes within the same LD block
ld.idx <- 1

temp.dir <- 'temp'
dir.create(temp.dir)

plink.hdr <- 'geno/rosmap1709-chr'
sum.file.dir <- 'stat/IGAP/data/hs-lm/'

ld.info <- ld.key[ld.idx, ]

plink <- subset.plink(ld.info, temp.dir, plink.hdr %&&% chr)
x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))

sum.stat.obj <- extract.sum.stat(ld.info = ld.info,
                                 sum.file = sum.file.dir %&&% chr %&&% '.eqtl_bed.gz',
                                 x.bim = x.bim,
                                 temp.dir = temp.dir,
                                 is.eqtl = TRUE)

sum.stat <- sum.stat.obj$sum.stat
genes <- sum.stat.obj$mediators

## load('temp.rdata')
## genes <- sum.stat %>% dplyr::select(med.id, hgnc, tss, tes, y.pos) %>% unique() %>% arrange(y.pos)

model.tab <- ld.tab %>%
    right_join(ld.info, by = c('chr', 'ld.lb', 'ld.ub'))

zqtl.data <- make.zqtl.data(plink, sum.stat, genes)

## just run on genes that exist on the list
# y.pos <- match(model.tab$ensg, genes$med.id)
# zqtl.data$qtl.theta <- zqtl.data$qtl.theta %c% y.pos
# zqtl.data$qtl.se <- zqtl.data$qtl.se %c% y.pos

## Re-estimate zQTL model for visualization
vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e4,
               vbiter = 5000, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2, nsample = 10, print.interv = 100,
               nboot = 0, med.finemap = FALSE, weight = TRUE, out.residual = FALSE)

med.out <- fit.med.zqtl(zqtl.data$gwas.theta,
                        zqtl.data$gwas.se,
                        zqtl.data$qtl.theta,
                        zqtl.data$qtl.se,
                        X = zqtl.data$X,
                        n = 74046,
                        n.med = 356,
                        options = vb.opt)

## predicted z-score trace
Vd <- sweep(t(med.out$Vt), 2, sqrt(med.out$D2), `*`)

z.dir <- Vd %*% (sweep(t(Vd), 2, 1/zqtl.data$gwas.se, `*`) %*% med.out$param.direct$theta)

qtl.theta.nz <- zqtl.data$qtl.theta
qtl.theta.nz[is.na(zqtl.data$qtl.theta)] <- 0
z.qtl.nz <- zqtl.data$qtl.theta / zqtl.data$qtl.se
z.qtl.nz[is.na(zqtl.data$qtl.theta)] <- 0

z.gwas <- zqtl.data$gwas.theta / zqtl.data$gwas.se
eta.dir <- t(Vd) %*% med.out$param.direct$theta
z.dir <- Vd %*% (t(Vd) %*% (med.out$param.direct$theta / zqtl.data$gwas.se))
var.dir <- sum(eta.dir^2)

get.pve <- function(k) {
    path.med <- (qtl.theta.nz %c% k) %*% (med.out$param.mediated$theta %r% k)
    var.med <- sum(path.med^2)
    pve <- var.med / (var.med + var.dir)
}

pve.vec <- sapply(1:ncol(z.qtl.nz), get.pve)

model.df <- cbind(model.tab, pve = pve.vec)



.scale.fill.count <- scale_fill_gradientn(colours = c('#7777FF', 'yellow'), trans = 'log',
                                         breaks = c(1, 10, 100))

## QTL information
gwas.z <- zqtl.data$gwas.theta[, 1] / zqtl.data$gwas.se[, 1]
gwas.df <- data.frame(snp.loc = zqtl.data$snps$snp.loc, gwas = gwas.z) 

qtl.z <- zqtl.data$qtl.theta / zqtl.data$qtl.se
colnames(qtl.z) <- model.tab$hgnc
rownames(qtl.z) <- zqtl.data$snps$snp.loc

qtl.melt <- melt(qtl.z, varnames = c('snp.loc', 'hgnc'))
qtl.melt$hgnc <- as.character(qtl.melt$hgnc)
qtl.melt <- qtl.melt %>% left_join(gwas.df, by = 'snp.loc')
qtl.melt$hgnc <- factor(qtl.melt$hgnc, model.tab$hgnc)

plt.qtl <- ggplot(qtl.melt %>% na.omit(), aes(x = value, y = gwas)) + theme_bw() +
    geom_hex() +
    .scale.fill.count +
    geom_smooth(method = 'lm', color = 'red', size = .5) +
    facet_grid(.~hgnc, scales = 'free') +
    xlab('QTL z-score') + ylab('GWAS z-score')


## After LD-adjustment
colnames(med.out$M) <- model.tab$hgnc
rownames(med.out$M) <- 1:nrow(ld.df)

Y <- data.frame(gwas = scale(med.out$Y / sqrt(med.out$D2))) %>%
    mutate(component = 1:n())

M <- sweep(med.out$M, 1, sqrt(med.out$D2), `/`) %>% scale()

m.df <- melt(M, varnames = c('component', 'hgnc'))
m.df$hgnc <- as.character(m.df$hgnc)
m.df <- m.df %>% left_join(Y, by = 'component')
m.df$hgnc <- factor(m.df$hgnc, model.tab$hgnc)

plt.qtl.adj <- ggplot(m.df, aes(x=value, y=gwas)) + theme_bw() +
    geom_hex() +
    .scale.fill.count +
    geom_smooth(method = 'lm', color = 'red', size = .5) +
    facet_grid(.~hgnc, scales = 'free') +
    xlab('adjusted QTL z') + ylab('adjusted GWAS z')



## TWAS vs mediation
blk.width <- 1e4

gwas.block.df <- gwas.df %>%
    mutate(blk.loc = floor(snp.loc/blk.width)*blk.width) %>%
    group_by(blk.loc) %>% slice(which.max(abs(gwas))) %>%
    mutate(gwas.p = l10.p.two(abs(gwas)))

plt.manhattan <- ggplot() + theme_bw() + ylab('-log10 P') + xlab('genomic location') +
    geom_hline(yintercept = -log10(5e-8), color = 'gray') + 
    geom_hline(yintercept = -log10(2.5e-6), lty = 2) + 
    geom_point(data=gwas.df, aes(x=snp.loc, y = l10.p.two(abs(gwas))), size = .1) +
    geom_rect(data=gwas.block.df, aes(xmin=blk.loc, xmax=blk.loc+blk.width, ymin=0, ymax=gwas.p),
              fill = 'gray', alpha = 0.5)

plt.manhattan +
    geom_rect(data = model.df,
              aes(xmin = tss, xmax = tes, ymin = 0, ymax = -log10(p.val.nwas)),
              fill = 'green', size = .5, color = 'green', alpha = .5) +
    geom_segment(data = model.df,
                 aes(x = (tes+tss)/2, xend = (tes+tss)/2, y = 0, yend = -log10(pmax(p.val.direct, p.val.marginal))), color = 'red', size = 1, arrow = arrow(length = unit(.05, 'inches')))







    ## LD information
    ld.pairs <- zqtl::take.ld.pairs(zqtl.data$X, cutoff = 0.1, stdize = TRUE) %>%
        as.data.frame()

    x.min <- min(zqtl.data$snps$snp.loc)
    x.max <- max(zqtl.data$snps$snp.loc)

    x.ld.min <- min(ld.pairs$x.pos + .5)
    x.ld.max <- max(ld.pairs$x.pos - .5)
    x.ld.len <- x.ld.max - x.ld.min
    x.len <- x.max - x.min

    ld.df <- ld.pairs %>% mutate(xx.pos = x.min + x.len * (x.pos - x.ld.min)/x.ld.len) %>%
        filter(x %in% qtl.df$x | y %in% qtl.df$x) %>%
            mutate(cov = pmin(pmax(cov, -0.5), 0.5))

    rm(ld.pairs)

    x.scale.ld <- scale_x_continuous(limits = range(ld.df$xx.pos) + c(-.5, .5),
                                      expand = c(0,0))

    ld.df <- ld.df %>% mutate(g = paste(x,y,sep='-'))




zz.y.direct <- med.out$param.direct$theta / zqtl.data$gwas.se
zz.y.direct <- sweep(svd.out$V.t %*% zz.y.direct, 1, svd.out$D^2, `*`)

p <- ncol(zqtl.data$X)
Vt.C <- svd.out$V.t %*% matrix(1/p, p, 1)
I.d <- matrix(1/p, nrow = nrow(Vt.C), ncol = 1)
zz.y.cov.1 <- sweep(Vt.C, 1, svd.out$D^2, `*`) %*% med.out$param.covariate.eta$theta
zz.y.cov.2 <- I.d %*% med.out$param.covariate.delta$theta

plot(zz.y.cov.1 + zz.y.cov.2)



gwas.z <- zqtl.data$gwas.theta / zqtl.data$gwas.se
qtl.z <- (zqtl.data$qtl.theta / zqtl.data$qtl.se)

mm <- lapply(1:ncol(qtl.z), function(j) {
    .q <- qtl.z %c% j
    .v <- which(!is.na(.q[, 1]))
    .x <- svd.out$V.t %c% .v
    .q <- .q %r% .v
    return(.x %*% .q)
})

M <- do.call(cbind, mm) %>% as.matrix()



y <- svd.out$V.t %*% gwas.z
y.m <- y - zz.y.direct

y.hat <- M %*% med.out$param.mediated$theta
r <- y.m - y.hat





lambda <- svd.out$D^2

.show <- which(1/lambda > 10)

plot(y.hat[.show], y[.show])

plot(y.hat[.show], y.m[.show])




med.theta <- model.genes %>% select(med.id) %>% rename(ensg = med.id) %>%
    left_join(model.tab, by = 'ensg') %>% select(theta) %>% as.matrix()




zz <- sweep(svd.out$V.t %*% qtl.z.nz, 1, zqtl.data$gwas.se, `*`)

y.hat <- zz %*% med.theta

plot(y.hat, y)



plot((qtl.z %c% 15) %*% (med.theta %r% 15), gwas.z)

cor.test((qtl.z %c% 15) %*% (med.theta %r% 15), gwas.z)

plot(qtl.z %c% 15, gwas.z)




## Figure : gene by gene LD block and statistics
out.dir <- 'figures/genes/'
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

draw.gene <- function(gene.idx, gene.tab, temp.dir, out.hdr) {

    file.name <- out.hdr %&&% sprintf('qtl_%d_%05d_%s.pdf', gene.tab[gene.idx, ]$chr,
                                      gene.idx, gene.tab[gene.idx, ]$hgnc)

    if(file.exists(file.name)) { return(NULL) }

    plink.hdr <- 'geno/rosmap1709-chr'
    sum.file.dir <- 'stat/IGAP/data/hs-lm/'

    ld.info <- as.vector(gene.tab[gene.idx, 1:3])
    chr <- ld.info[1]
    names(ld.info) <- c('chr', 'ld.lb', 'ld.ub')
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

    ensg <- as.character(gene.tab[gene.idx, ]$ensg)
    hgnc <- as.character(gene.tab[gene.idx, ]$hgnc)
    y.pos <- which(genes$med.id == ensg)

    ## QTL information
    qtl.z <- zqtl.data$qtl.theta[, y.pos] / zqtl.data$qtl.se[, y.pos]

    qtl.df <- data.frame(snp.loc = zqtl.data$snps$snp.loc, qtl.z) %>%
        mutate(x = 1:n()) %>%
            mutate(ln.p = -log10(2 * pnorm(abs(qtl.z), lower.tail = FALSE)))

    ## GWAS information
    gwas.z <- zqtl.data$gwas.theta[, 1] / zqtl.data$gwas.se[, 1]

    gwas.df <- data.frame(snp.loc = zqtl.data$snps$snp.loc, gwas.z) %>%
        mutate(x = 1:n()) %>%
            mutate(ln.p = -log10(2 * pnorm(abs(gwas.z), lower.tail = FALSE)))

    ## LD information
    ld.pairs <- zqtl::take.ld.pairs(zqtl.data$X, cutoff = 0.05, stdize = TRUE) %>%
        as.data.frame()

    x.min <- min(zqtl.data$snps$snp.loc)
    x.max <- max(zqtl.data$snps$snp.loc)

    x.ld.min <- min(ld.pairs$x.pos + .5)
    x.ld.max <- max(ld.pairs$x.pos - .5)
    x.ld.len <- x.ld.max - x.ld.min
    x.len <- x.max - x.min

    ld.df <- ld.pairs %>% mutate(xx.pos = x.min + x.len * (x.pos - x.ld.min)/x.ld.len) %>%
        filter(x %in% qtl.df$x | y %in% qtl.df$x) %>%
            mutate(cov = pmin(pmax(cov, -0.5), 0.5))

    rm(ld.pairs)

    x.scale.ld <- scale_x_continuous(limits = range(ld.df$xx.pos) + c(-.5, .5),
                                      expand = c(0,0))

    ld.df <- ld.df %>% mutate(g = paste(x,y,sep='-'))
    
    p1 <- ggplot(ld.df) + x.scale.ld +
        theme(axis.text.y = element_blank(), axis.title = element_blank()) +
            geom_polygon(aes(group = g, fill = cov, x = xx.pos, y = -y.pos)) +
                theme(legend.position = c(0, 1),
                      legend.justification = c(0, 1),
                      legend.title = element_text(size = 8),
                      legend.background = element_blank(),
                      axis.text.x = element_text(size = 8)) +
                          scale_fill_continuous(low = 'blue', high = 'yellow')
    
    p2 <- ggplot(qtl.df) + theme_classic() + ylab('eQTL') +
        geom_point(aes(x = snp.loc, y = ln.p), size = .5) +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_text(size = 8)) +
                x.scale.ld

    p3 <- ggplot(gwas.df) + theme_classic() + ylab('GWAS') +
        geom_point(aes(x = snp.loc, y = ln.p), size = .5) +
            xlab('Genomic location') +
                x.scale.ld +
                    theme(axis.text.x = element_text(size = 8))

    rr <- (max(abs(ld.df$y.pos)) + 3) / max(ld.df$x.pos)

    ld.file.name <- out.hdr %&&% sprintf('ld_%d_%05d_%s.png', gene.tab[gene.idx, ]$chr,
                                      gene.idx, gene.tab[gene.idx, ]$hgnc)

    png(file = ld.file.name,
        width = 100 * (1 + (x.len / 1e5) * .25),
        height = 200 * (.5 + 2 * rr))
    print(p1)
    dev.off()

    heights <- c(.5, .5)
    pdf(file = file.name, width = 1 + (x.len / 1e5) * .25, height = 2 * sum(heights) + 1,
        useDingbats = FALSE)
    print(grid.vcat(list(p2, p3), heights = heights))
    dev.off()
}

## temporary directory
temp.dir <- system('mktemp -d ' %&&% out.dir %&&% '/temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

draw.gene(gene.idx, gene.tab, temp.dir, out.hdr = out.dir %&&% '/')

system('rm -r ' %&&% temp.dir)
log.msg('Finished figure drawing\n\n\n')
