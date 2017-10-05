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
source('figure.util.R')

gene.tab.file <- 'tables/bootstrap_gene_significant.txt.gz'
gene.tab <- read.table(gene.tab.file, header = TRUE) %>% na.omit()

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
        width = 300 * (1 + (x.len / 1e5) * .25),
        height = 300 * (.5 + 2 * rr))
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
