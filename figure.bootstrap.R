#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)
source('util.R')
source('mediation.R')
library(readr)
library(zqtl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(latex2exp)
source('figure.util.R')

gene.tab.file <- 'tables/bootstrap_gene.txt.gz'
gene.tab <- read.table(gene.tab.file, header = TRUE) %>% na.omit()
gene.tab.sig.file <- 'tables/bootstrap_gene_significant.txt.gz'
gene.tab.sig <- read.table(gene.tab.sig.file, header = TRUE) %>% na.omit()

## Figure 1. global
draw.chr <- function(.df) {

    .df.sig <- .df %>% filter(p.val.marginal < 1e-6, p.val.direct < 1e-6, lodds > 2)

    plt <- ggplot() + theme_bw() +
        geom_linerange(data = .df %>% filter(lodds > 0),
                       aes(x = (tss+tes)/2, ymin = theta - 2 * theta.se, ymax = theta + 2* theta.se),
                       color = 'gray40', alpha = 0.5)

    plt <- plt +
        geom_point(data = .df, aes(x = (tss+tes)/2, y = theta, size = 1/(1+exp(-lodds))),
                   alpha = 0.8, color = 'gray40', show.legend = FALSE)

    plt <- plt +
        geom_point(data = .df.sig, aes(x = (tss+tes)/2, y = theta),
                   color = 'red', pch = 4, size = 1) +
                       facet_grid(.~chr, scales = 'free', space = 'free') +
                           scale_size_continuous(limits = c(0, 1), range = c(0, 1))
    
    plt <- plt + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    plt <- plt +
        geom_text_repel(data = .df.sig %>% filter(theta > 0), aes(x = (tss+tes)/2, y = theta+2*theta.se, label = hgnc), nudge_y = .1, size = 2, segment.color = 'green', segment.alpha = 0.5)

    plt <- plt +
        geom_text_repel(data = .df.sig %>% filter(theta < 0), aes(x = (tss+tes)/2, y = theta-2*theta.se, label = hgnc), nudge_y = -.1, size = 2, segment.color = 'green', segment.alpha = 0.5) +
            ylab('Mediation effect') + xlab('Genomic location')
    
    return(plt)
}

p.list <- list(draw.chr(gene.tab %>% filter(chr <= 3)),
               draw.chr(gene.tab %>% filter(chr >= 4, chr <= 7)),
               draw.chr(gene.tab %>% filter(chr >= 7, chr <= 11)),
               draw.chr(gene.tab %>% filter(chr >= 12, chr <= 16)),
               draw.chr(gene.tab %>% filter(chr >= 17)))

out.dir <- 'figures/genes/'
out.files <- out.dir %&&% 'global_' %&&% 1:length(p.list) %&&% '.pdf'

sapply(1:length(p.list), function(j) ggsave(plot = p.list[[j]], filename = out.files[j],
                                            useDingbats = FALSE, width = 10, height = 3))

## Figure 2. gene by gene scatter plots
out.dir <- 'figures/genes/'
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

draw.gene <- function(gene.idx, gene.tab, temp.dir, out.hdr) {

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

    sum.stat.tab <- sum.stat.obj$sum.stat
    genes <- sum.stat.obj$mediators

    zqtl.data <- make.zqtl.data(plink, sum.stat.tab, genes)

    y.pos <- which(genes$med.id == as.character(gene.tab[gene.idx, 'ensg']))

    param <- gene.tab[gene.idx, ] %>% select(theta, theta.se)

    z.x <- (zqtl.data$qtl.theta / zqtl.data$qtl.se) %c% y.pos
    z.y <- (zqtl.data$gwas.theta / zqtl.data$gwas.se)

    valid <- which(!is.na(z.x) & !is.na(z.y))

    svd.out <- take.ld.svd(zqtl.data$X, options = list(do.stdize = TRUE, eigen.tol = 1e-2))
    
    V.t <- svd.out$V.t

    zz.x <- V.t %c% valid %*% (z.x %r% valid)
    zz.y <- V.t %c% valid %*% (z.y %r% valid)

    p0 <- ggplot(data.frame(QTL = z.x[valid], GWAS = z.y[valid])) + theme_bw() +
        theme(legend.position = 'bottom', legend.background = element_blank()) +
            geom_point(aes(x = QTL, y = GWAS), alpha = 0.75, size = 1) +
                ggtitle(genes[y.pos, ]$med.id)

    .df <- data.frame(QTL = zz.x, GWAS = zz.y, w = 1/svd.out$D^2)
    viz.cutoff <- median(.df$w) * 0.20

    .title <- sprintf('[%s] %.4f +/- %.2e', gene.tab[gene.idx, ]$hgnc, param[1], param[2])

    p1 <- ggplot(.df %>% filter(w > viz.cutoff)) + theme_bw() +
        geom_smooth(aes(x = QTL, y = GWAS, weight = w), method = 'lm', size = .5) +
            ggtitle(.title) +
                theme(legend.background = element_blank()) +
                    scale_size_continuous(name = TeX('$1/\\lambda$'), range = c(0, 2)) +
                        geom_point(aes(x = QTL, y = GWAS, size=w), alpha = 0.75)
    
    file.name <- out.hdr %&&% sprintf('%d_%05d_%s.pdf', gene.tab[gene.idx, ]$chr,
                                      gene.idx, gene.tab[gene.idx, ]$hgnc)

    pdf(file = file.name, width = 7, height = 4, useDingbats = FALSE)
    print(grid.hcat(list(p0, p1), widths = c(1, 1.55)))
    dev.off()

}

## temporary directory
temp.dir <- system('mktemp -d ' %&&% out.dir %&&% '/temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

names(gene.tab.sig)[2] <- 'ld.lb'
names(gene.tab.sig)[3] <- 'ld.ub'

for(g in 1:nrow(gene.tab.sig)){
    draw.gene(g, gene.tab.sig, temp.dir, out.hdr = out.dir %&&% '/')
    log.msg('Finished : %d / %d\n', g, nrow(gene.tab.sig))
}

system('rm -r ' %&&% temp.dir)
log.msg('Finished figure drawing\n\n\n')
