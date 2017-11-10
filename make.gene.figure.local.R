#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 5) q()

ld.tab.file <- argv[1]        # e.g., ld.tab.file <- 'tables/genes/full-5/significant_LD.txt.gz'
chr <- as.integer(argv[2])
ld.lb <- as.integer(argv[3])
ld.ub <- as.integer(argv[4])
out.file <- argv[5]

if(file.exists(out.file)) q()

options(stringsAsFators = FALSE)
source('util.R')
source('mediation.R')
source('figure.util.R')
library(readr)
library(zqtl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(scales)

ld.tab <- read.table(ld.tab.file, header = TRUE)

med.cutoff <- ld.tab %>% filter(qval < 1e-2) %>%
    slice(which.max(pval))

med.cutoff <- med.cutoff$pval

dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

plink.hdr <- 'geno/rosmap1709-chr'
sum.file.dir <- 'stat/IGAP/data/hs-lm/'

ld.info <- data.frame(chr = chr, ld.lb = ld.lb, ld.ub = ld.ub)
gwas.file <- 'IGAP/chr' %&&% chr %&&% '.txt.gz'

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

################################################################
## should include additional GWAS SNPs not matched with QTL
missing.gwas <- find.missing.gwas(gwas.file, ld.info, sum.stat)

sum.stat.combined <-
    rbind(sum.stat %>% dplyr::select(x.pos, y.pos, gwas.theta, gwas.se, qtl.theta, qtl.se),
          missing.gwas %>%
              dplyr::select(x.pos, y.pos, gwas.theta, gwas.se, qtl.theta, qtl.se))


################################################################
model.tab <- ld.tab %>%
    right_join(ld.info, by = c('chr', 'ld.lb', 'ld.ub'))

model.tab <- model.tab %>%
    mutate(med.p = -log10(pval))

gwas.df <- read.gwas.tab(gwas.file, ld.info) %>%
    mutate(gwas.p = l10.p.two(abs(gwas.z)))

################################################################
## Figure out the range
blk.width <- 1e4

gwas.block.df <- gwas.df %>%
    mutate(blk.loc = floor(snp.loc/blk.width)*blk.width) %>%
        group_by(blk.loc) %>% slice(which.max(abs(gwas.z)))

gwas.range <- range(gwas.block.df$snp.loc)
gene.range <- range(c(model.tab$tss, model.tab$tes))
blk.range <- range(c(gwas.range, gene.range))
blk.range <- range(c(min(blk.range - blk.width), blk.range, max(gwas.block.df$blk.loc) + blk.width))
blk.length <- blk.range[2] - blk.range[1]
.gwas.w <- blk.length / blk.width * 0.02 + 1

.gwas.kb <- function() {
    function(x) format(x/1e3, big.mark=',')
}

.gwas.x.scale <- scale_x_continuous(limits = blk.range + c(-1000, 1000), expand = c(0, 0),
                                    labels = .gwas.kb())

.gwas.x.scale.top <- scale_x_continuous(limits = blk.range + c(-1000, 1000), expand = c(0, 0),
                                    labels = .gwas.kb(), position = 'top')

plt.gwas <-
    gg.plot() + ylab('-log10 GWAS p-value') +
    theme(axis.title.x = element_blank()) +
    geom_rect(data=gwas.block.df, aes(xmin=blk.loc, xmax=blk.loc+blk.width, ymin=0, ymax=gwas.p),
              fill = 'gray80') +
    geom_point(data=gwas.df, aes(x=snp.loc, y = gwas.p), size = .1)

if(max(gwas.df$gwas.p) > -log10(5e-8)) {
    plt.gwas <- plt.gwas + geom_hline(yintercept = -log10(5e-8), color = 'red')
}

plt.gwas <- plt.gwas + scale_y_reverse() + .gwas.x.scale.top

################################################################
## [1] PVE of genes (uniformly distributed)
## [2] grid - gene location
## [3] genes - SNPs
## [4] SNPs - linked SNPs
## [5] GWAS upside down

model.tab.sorted <- model.tab %>% arrange((tss+tes)/2) %>%
    mutate(.gene.loc = seq(blk.range[1] + blk.width/2, blk.range[2] - blk.width/2,
               by = blk.length / n()))

empty.theme <- theme(axis.text = element_blank(), 
                     line = element_blank(), axis.line = element_blank(),
                     axis.ticks = element_blank(), panel.border = element_blank(),
                     legend.position = c(1,1), legend.justification = c(1, 1))

p1 <- gg.plot(model.tab.sorted) + ylab('% variance explained') + xlab('') + .gwas.x.scale +
    geom_segment(aes(x = .gene.loc, xend = .gene.loc, y = 0, yend = 100 * PVE, color = theta/theta.se), size = 3) +
    geom_text(aes(x = .gene.loc, y = 100 * PVE, label = 100 * signif(PVE, 2)), size = 3,
              hjust=0, vjust=0) +
    scale_color_gradient2(low = 'blue', high = 'red', mid = 'gray40', guide = FALSE) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

top.genes <- model.tab %>% arrange(desc(PVE)) %>% head(1)
strong.genes <- model.tab %>% filter(PVE > .01 | qval < 5e-2)
show.genes <- unique(c(top.genes$med.id, strong.genes$med.id))

p2 <- gg.plot(model.tab.sorted) +
    ylab('') + xlab('') + .gwas.x.scale +
    geom_segment(aes(x = tss, xend = tes, y=0, yend=0), color = '#005500', size = 3) +
    geom_segment(aes(x = (tss+tes)/2, xend = .gene.loc, y = 0, yend =1), color = 'gray') +
    geom_text(aes(x = .gene.loc, y = 1, label = format(pval, digits=2, scientific=TRUE)),
              vjust = 1, hjust = 0, size = 3) +
    geom_text(aes(x = .gene.loc, y = 1, label = hgnc),
              vjust = 0, hjust = 0, size = 3) +
    empty.theme + ylim(c(0, 2))

.qtl.aes <- aes(x = snp.loc, xend = (tss+tes)/2, y = 0, yend = 1,
                color = pmin(pmax(qtl.z, -5), 5), size = abs(qtl.z))

sum.stat.significant <- sum.stat %>% filter(med.id %in% show.genes)
sum.stat.significant$rs <- as.character(sum.stat.significant$rs)

p3 <- gg.plot(sum.stat.significant) + ylab('') + xlab('') + .gwas.x.scale +
    geom_segment(.qtl.aes) +
    scale_size_continuous(range = c(0, 1), guide = FALSE) +
    scale_color_gradientn('QTL', colors = c('blue', 'white', 'red')) +
    empty.theme +
    theme(legend.position = c(1,1), legend.justification = c(1, 1))

## Who are significant LD partners?
x.bim.matched <- 
    sum.stat.significant %>% dplyr::select(rs) %>% unique() %>%
    left_join(x.bim)

svd.out <- take.ld.svd(plink$BED, options = list(do.stdize = TRUE))

Vd <- sweep(t(svd.out$V.t), 2, svd.out$D, `*`)
Vd.sub <- Vd %r% (x.bim.matched$x.pos)
LD <- Vd.sub %*% t(Vd)

take.ld.bound <- function(r) {
    .temp <- which(abs(LD[r, ]) > .1)
    data.frame(r = r, lb = min(.temp), ub = max(.temp))
}

ld.bound <- do.call(rbind, lapply(1:nrow(LD), take.ld.bound))

.ld.bound.tab <- data.frame(rs = as.character(x.bim.matched$rs),
                           snp.loc = x.bim.matched$snp.loc,
                           ld.lb = plink$BIM[ld.bound[, 2], 4],
                           ld.ub = plink$BIM[ld.bound[, 3], 4])

ld.bound.tab <- sum.stat.significant %>% left_join(.ld.bound.tab) %>% na.omit() %>%
    group_by(med.id) %>%
    summarize(snp.lb = min(snp.loc), snp.ub = max(snp.loc),
              ld.lb = min(ld.lb), ld.ub = max(ld.ub))

ld.poly.tab <- ld.bound.tab %>%
    gather(key = 'y.loc', value = 'x', -med.id) %>%
    mutate(y.loc = gsub(y.loc, pattern = '.lb', replacement = '')) %>%
    mutate(y.loc = gsub(y.loc, pattern = '.ub', replacement = '')) %>%
    mutate(y = ifelse(y.loc == 'snp', 1, 0)) %>%
    arrange(x)

p4 <-
    gg.plot() + empty.theme + .gwas.x.scale +
    geom_polygon(data = ld.poly.tab,
                 aes(x = x, y = y, group = med.id),
                 alpha = 0.25, fill = 'green', show.legend = FALSE) +
    geom_segment(data = ld.bound.tab,
                 aes(x = snp.lb, xend = ld.lb, y = 1, yend = 0, color = med.id),
                 show.legend = FALSE) +
    geom_segment(data = ld.bound.tab,
                 aes(x = snp.ub, xend = ld.ub, y = 1, yend = 0, color = med.id),
                 show.legend = FALSE) +
    theme(axis.title = element_blank())

gg <- grid.vcat(list(p1, p2, p3, p4, plt.gwas), heights = c(2, .5, .5, .5, 2))

ggsave(filename = out.file, plot = gg, width = 8, height = 8, units = 'in', limitsize = FALSE)

system('rm -r ' %&&% temp.dir)
log.msg('Finished figure drawing\n\n\n')
