#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

ld.idx <- as.integer(argv[1])
out.dir <- argv[2]

                                        # ld.idx <- 1

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
library(scales)

ld.tab.file <- 'tables/bootstrap_ld_significant.txt.gz'
ld.tab <- read.table(ld.tab.file, header = TRUE)
ld.key <- ld.tab %>% select(chr, ld.lb, ld.ub) %>% unique()

if(ld.idx > nrow(ld.key)) q()

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

plink.hdr <- 'geno/rosmap1709-chr'
sum.file.dir <- 'stat/IGAP/data/hs-lm/'

ld.info <- ld.key[ld.idx, ]
chr <- as.integer(ld.info[1])

ld.out.file <- sprintf('%s/chr%d_ld%d_heatmap.png', out.dir, chr, ld.idx)

if(file.exists(ld.out.file)) q()

## Draw all genes within the same LD block
temp.dir <- system('mktemp -d ' %&&% out.dir %&&% '/temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

plink <- subset.plink(ld.info, temp.dir, plink.hdr %&&% chr)
x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))

sum.stat.obj <- extract.sum.stat(ld.info = ld.info,
                                 sum.file = sum.file.dir %&&% chr %&&% '.eqtl_bed.gz',
                                 x.bim = x.bim,
                                 temp.dir = temp.dir,
                                 is.eqtl = TRUE)

sum.stat <- sum.stat.obj$sum.stat
genes <- sum.stat.obj$mediators

model.tab <- ld.tab %>%
    right_join(ld.info, by = c('chr', 'ld.lb', 'ld.ub'))

zqtl.data <- make.zqtl.data(plink, sum.stat, genes)

## just run on genes that exist on the list
y.pos <- match(model.tab$ensg, genes$med.id)
zqtl.data$qtl.theta <- zqtl.data$qtl.theta %c% y.pos
zqtl.data$qtl.se <- zqtl.data$qtl.se %c% y.pos

################################################################
## Re-estimate zQTL model for PVE calculation
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

z.dir.rot <- sweep(t(Vd), 2, 1/zqtl.data$gwas.se, `*`) %*% med.out$param.direct$theta
z.dir <- Vd %*% z.dir.rot

qtl.theta.nz <- zqtl.data$qtl.theta
qtl.theta.nz[is.na(zqtl.data$qtl.theta)] <- 0
z.qtl.nz <- zqtl.data$qtl.theta / zqtl.data$qtl.se
z.qtl.nz[is.na(zqtl.data$qtl.theta)] <- 0

eta.dir <- t(Vd) %*% med.out$param.direct$theta
var.dir <- sum(eta.dir^2)

Wt <- sweep(med.out$Vt, 1, sqrt(med.out$D2), `/`)

get.med.var <- function(k) {
    path.med <- (qtl.theta.nz %c% k) %*% (med.out$param.mediated$theta %r% k)
    var.med <- t(Wt %*% path.med) %*% (Wt %*% path.med)
}

med.var.vec <- sapply(1:ncol(z.qtl.nz), get.med.var)
var.tot <- sum(med.var.vec) + var.dir
pve.vec <- med.var.vec / var.tot
pve.dir <- var.dir / var.tot

model.df <- cbind(model.tab, pve = pve.vec) %>%
    mutate(l10.lodds = pmin(-log10(pmax(p.val.direct, p.val.marginal)), 10))

################################################################
## scatter plots
.scale.fill.count <- scale_fill_gradientn(colours = c('#7777FF', 'yellow'), trans = 'log',
                                          breaks = c(1, 10, 100))

.scatter.w <- .5 + nrow(model.df) * .8
.scatter.h <- 1 + .5

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
            geom_smooth(method = 'lm', se = FALSE,  color = 'red', size = .5) +
                facet_grid(.~hgnc, scales = 'free') +
                    xlab('QTL z-score') + ylab('GWAS z-score')

plt.qtl.file <- sprintf('%s/chr%d_ld%d_qtl.pdf', out.dir, chr, ld.idx)

ggsave(plt.qtl, file = plt.qtl.file, width = .scatter.w, height = .scatter.h,
       useDingbats = FALSE, limitsize = FALSE)



## After LD-adjustment
colnames(med.out$M) <- model.tab$hgnc
rownames(med.out$M) <- 1:nrow(med.out$M)

## z.dir.rot = D^{-1} D^2 V' S^{-1} * theta
yy <- med.out$Y / sqrt(med.out$D2)
yy.1 <- yy - z.dir.rot

Y <- data.frame(gwas = scale(yy), gwas.1 = scale(yy.1)) %>% mutate(component = 1:n())

M <- sweep(med.out$M, 1, sqrt(med.out$D2), `/`) %>% scale()

m.df <- melt(M, varnames = c('component', 'hgnc'))
m.df$hgnc <- as.character(m.df$hgnc)
m.df <- m.df %>% left_join(Y, by = 'component')
m.df$hgnc <- factor(m.df$hgnc, model.tab$hgnc)

plt.qtl.adj <- ggplot(m.df, aes(x=value, y=gwas)) + theme_bw() +
    geom_hex() +
        .scale.fill.count +
            geom_smooth(method = 'lm', se = FALSE,  color = 'red', size = .5) +
                facet_grid(.~hgnc, scales = 'free') +
                    xlab('adjusted QTL z') + ylab('adjusted GWAS z')

plt.qtl.adj.file <- sprintf('%s/chr%d_ld%d_qtl_adj.pdf', out.dir, chr, ld.idx)

ggsave(plt.qtl.adj, file = plt.qtl.adj.file, width = .scatter.w, height = .scatter.h,
       useDingbats = FALSE, limitsize = FALSE)

plt.qtl.adj.1 <- ggplot(m.df, aes(x=value, y=gwas.1)) + theme_bw() +
    geom_hex() +
        .scale.fill.count +
            geom_smooth(method = 'lm', se = FALSE,  color = 'red', size = .5) +
                facet_grid(.~hgnc, scales = 'free') +
                    xlab('adjusted QTL z') + ylab('adjusted GWAS z')

plt.qtl.adj.file <- sprintf('%s/chr%d_ld%d_qtl_adj_nod.pdf', out.dir, chr, ld.idx)

ggsave(plt.qtl.adj.1, file = plt.qtl.adj.file, width = .scatter.w, height = .scatter.h,
       useDingbats = FALSE, limitsize = FALSE)



################################################################
## GWAS vs TWAS vs mediation
blk.width <- 1e4

gwas.block.df <- gwas.df %>%
    mutate(blk.loc = floor(snp.loc/blk.width)*blk.width) %>%
        group_by(blk.loc) %>% slice(which.max(abs(gwas))) %>%
            mutate(gwas.p = l10.p.two(abs(gwas)))

blk.range <- range(c(gwas.block.df$blk.loc, model.df$tss, model.df$tes, max(gwas.block.df$blk.loc) + blk.width))
blk.length <- blk.range[2] - blk.range[1]
.gwas.w <- blk.length / blk.width * 0.02 + 1

.gwas.kb <- function() {
    function(x) format(x/1e3, big.mark=',')
}

.gwas.x.scale <- scale_x_continuous(limits = blk.range + c(-1000, 1000), expand = c(0, 0),
                                    labels = .gwas.kb())

## -- GWAS
plt.manhattan <- ggplot() + theme_bw() + ylab('-log10 P') + xlab('genomic location') +
    .gwas.x.scale + 
        geom_hline(yintercept = -log10(5e-8), color = 'gray') + 
            geom_hline(yintercept = -log10(2.5e-6), lty = 2) + 
                geom_point(data=gwas.df, aes(x=snp.loc, y = l10.p.two(abs(gwas))), size = .1) +
                    geom_rect(data=gwas.block.df, aes(xmin=blk.loc, xmax=blk.loc+blk.width, ymin=0, ymax=gwas.p),
                              fill = 'gray', alpha = 0.5)

out.file <- sprintf('%s/chr%d_ld%d_gwas.pdf', out.dir, chr, ld.idx)

ggsave(plt.manhattan, file = out.file, width = .gwas.w, height = 3,
       useDingbats = FALSE, limitsize = FALSE)


## -- TWAS
plt.manhattan <- plt.manhattan +
    geom_rect(data = model.df,
              aes(xmin = tss, xmax = tes, ymin = 0, ymax = pmin(-log10(p.val.nwas), 10)),
              fill = 'green', size = .5, color = 'green', alpha = .3)

out.file <- sprintf('%s/chr%d_ld%d_gwas_twas.pdf', out.dir, chr, ld.idx)

ggsave(plt.manhattan, file = out.file, width = .gwas.w, height = 3,
       useDingbats = FALSE, limitsize = FALSE)


## -- mediation
plt.manhattan <- plt.manhattan +
    geom_segment(data = model.df,
                 aes(x = (tes+tss)/2, xend = (tes+tss)/2, y = 0, yend = l10.lodds),
                 color = 'red', size = .5, arrow = arrow(length = unit(.05, 'inches')))

if(model.df %>% filter(lodds > -2) %>% nrow() > 0) {

    plt.manhattan <- plt.manhattan +
        geom_text_repel(data = model.df %>% filter(lodds > -2),
                        aes(x = (tss + tes)/2, y = l10.lodds, label = hgnc), size = 3,
                        nudge_y = 2, segment.alpha = 0.5, segment.color = '#FFCCCC',
                        segment.size = 0.1, angle = 30, direction = 'y')

}

if(model.df %>% filter(lodds <= -2) %>% nrow() > 0) {

    plt.manhattan <- plt.manhattan +
        geom_text_repel(data = model.df %>% filter(lodds <= -2),
                        aes(x = (tss + tes)/2, y = 0, label = hgnc), size = 2,
                        segment.color = '#FFCCCC', segment.size = 0.1, angle = 30, direction = 'y')

}

out.file <- sprintf('%s/chr%d_ld%d_gwas_twas_med.pdf', out.dir, chr, ld.idx)

ggsave(plt.manhattan, file = out.file, width = .gwas.w, height = 3,
       useDingbats = FALSE, limitsize = FALSE)



## PVE plot
plt.pve <- ggplot(model.df) + theme_bw() + .gwas.x.scale +
    ylab('PVE') + theme(axis.title.x = element_blank()) +
        geom_hline(yintercept = pve.dir, lty = 2) +
            geom_rect(aes(xmin = tss, xmax = tes, ymin = 0, ymax = pve),
                      fill = 'orange', size = .5, color = 'orange', alpha = 0.5)


.pve.logit <- c(unique(round(log(pve.vec) - log(1 - pve.vec))), 0, 2, 4)
.pve.brk <- unique(signif(c(pve.dir, 1/(1+exp(-.pve.logit)))))

.logit.kb <- function() {
    function(x) format(x/1e3, big.mark=',')
}

.logit <- trans_new(".logit",
                    function(x) signif(log(x) - log(1-x), 2),
                    function(y) signif(1/(1+exp(-y)), 2))

.pve.scale <- scale_y_continuous(breaks = .pve.brk, trans = .logit)

plt.pve <- plt.pve + .pve.scale

if(model.df %>% filter(pve > .01) %>% nrow() > 0) {
    plt.pve <- plt.pve +
        geom_text_repel(data = model.df %>% filter(pve > .01),
                        aes(x = (tss + tes)/2, y = pve, label = hgnc),
                        size = 4, nudge_y = .1, segment.alpha = 0.5, angle = 30)
}

if(model.df %>% filter(pve <= .01) %>% nrow() > 0) {
    plt.pve <- plt.pve +
        geom_text_repel(data = model.df %>% filter(pve <= .01),
                        aes(x = (tss + tes)/2, y = 0, label = hgnc),
                        size = 2, nudge_y = -1, segment.size = 0.1, angle = 30)
}

out.file <- sprintf('%s/chr%d_ld%d_pve.pdf', out.dir, chr, ld.idx)

ggsave(plt.pve, file = out.file, width = .gwas.w, height = 3,
       useDingbats = FALSE, limitsize = FALSE)




################################################################
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
    mutate(cov = pmin(pmax(cov, -0.5), 0.5))

rm(ld.pairs)

ld.df <- ld.df %>% mutate(g = paste(x,y,sep='-'))

.thm <- theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_blank(),
              legend.position = c(0, 1),
              legend.justification = c(0, 1),
              legend.title = element_text(size = 8),
              legend.background = element_blank(),
              axis.text.x = element_text(size = 8))

plt.ld <- ggplot(ld.df) + .gwas.x.scale + .thm +
    theme(axis.text.y = element_blank(), axis.title = element_blank()) +
        geom_polygon(aes(group = g, fill = cov, x = xx.pos, y = -y.pos)) +
            scale_fill_continuous(low = 'blue', high = 'yellow')

ld.out.file <- sprintf('%s/chr%d_ld%d_heatmap.png', out.dir, chr, ld.idx)

rr <- (max(abs(ld.df$y.pos)) + 3) / max(ld.df$x.pos)

png(file = ld.out.file,
    width = 100 * .gwas.w,
    height = 100 * (.5 + 2 * rr),
    bg = 'transparent')
print(plt.ld)
dev.off()

system('rm -r ' %&&% temp.dir)
log.msg('Finished figure drawing\n\n\n')
