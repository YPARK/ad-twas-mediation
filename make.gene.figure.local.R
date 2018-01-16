#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 5) q()

ld.tab.file <- argv[1]        # e.g., ld.tab.file <- 'tables/genes/full-0/significant_LD.txt.gz'
chr <- as.integer(argv[2])    # e.g., chr <- 1 
ld.lb <- as.integer(argv[3])  # e.g., ld.lb <- 64105056
ld.ub <- as.integer(argv[4])  # e.g., ld.ub <- 65696189
out.file <- argv[5]           # e.g., out.file <- 'temp.pdf'

if(file.exists(out.file)) q()

options(stringsAsFators = FALSE)
source('util.R')
source('mediation.R')
source('figure.util.R')
library(readr)
library(zqtl)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(reshape2)
library(scales)
library(ggrepel)

ld.tab <- read.table(ld.tab.file, header = TRUE)

med.cutoff <- ld.tab %>% filter(qval < 1e-4) %>%
    slice(which.min(lodds))

med.cutoff <- med.cutoff$lodds

dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)

plink.hdr <- 'geno/rosmap1709-chr'
sum.file.dir <- 'stat/IGAP/data/hs-lm/'
nwas.file <- 'nwas/IGAP_rosmap_eqtl_hs-lm_' %&&% chr %&&% '.nwas.gz'

ld.info <- data.frame(chr = chr, ld.lb = ld.lb, ld.ub = ld.ub)
gwas.file <- 'IGAP/chr' %&&% chr %&&% '.txt.gz'

## Draw all genes within the same LD block
out.dir <- dirname(out.file)
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
## Read TWAS results for comparison
nwas.tab <- read_tsv(nwas.file,
                     col_names = c('med.id', 'twas', 'chr', 'ld.lb', 'ld.ub', 'n.qtl.ld'))

twas.p <- 2 * pnorm(abs(nwas.tab$twas), lower.tail = FALSE)
twas.q <- p.adjust(twas.p, method = 'fdr')
twas.cutoff <- max(twas.p[twas.q < 1e-4])

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

################################################################
## Combine mediation and TWAS

model.tab.sorted <- model.tab %>%
    arrange((tss+tes)/2) %>%
        mutate(.gene.loc = seq(blk.range[1] + blk.width/2, blk.range[2] - blk.width/2, by = blk.length / n()))

model.tab.sorted <- model.tab.sorted %>%
    left_join(nwas.tab) %>%
        mutate(twas.p = 2 * pnorm(abs(twas), lower.tail = FALSE))
        
top.genes <- model.tab.sorted %>% arrange(desc(PVE.glob)) %>% head(1)
top.twas.gene <- model.tab.sorted %>% arrange(twas.p) %>% head(1)
strong.genes <- model.tab.sorted %>% filter(lodds > 0 | twas.p < twas.cutoff)
show.genes <- unique(c(top.genes$med.id, top.twas.gene$med.id, strong.genes$med.id))

## LD boundary for each gene
x.bim.matched <- 
    sum.stat %>% dplyr::select(rs) %>% unique() %>%
    left_join(x.bim)

svd.out <- take.ld.svd(plink$BED, options = list(do.stdize = TRUE))
Vd <- sweep(t(svd.out$V.t), 2, svd.out$D, `*`)
Vd.sub <- Vd %r% (x.bim.matched$x.pos)
LD <- Vd.sub %*% t(Vd)

take.ld.bound <- function(r) {
    .temp <- which(abs(LD[r, ]) >= .1)
    data.frame(r = r, lb = min(.temp), ub = max(.temp))
}

ld.bound <- do.call(rbind, lapply(1:nrow(LD), take.ld.bound))

.ld.bound.tab <- data.frame(rs = as.character(x.bim.matched$rs),
                           snp.loc = x.bim.matched$snp.loc,
                           ld.lb = plink$BIM[ld.bound[, 2], 4],
                           ld.ub = plink$BIM[ld.bound[, 3], 4])

ld.bound.tab <- sum.stat %>% left_join(.ld.bound.tab) %>% na.omit() %>%
    group_by(med.id) %>%
        summarize(snp.lb = min(snp.loc), snp.ub = max(snp.loc),
                  lb = min(ld.lb), ub = max(ld.ub))

## Gene frequency in 10kb bin
func.footprint <- function(tab, res = 1e4) {
    ret <- seq(floor(min(tab$lb)/res), ceiling(max(tab$ub)/res) + 1) * res
}

take.footprint.freq <- function(tab) {
    tab %>% group_by(.id) %>%
        do(footprint = func.footprint(.)) %>%
            tidy(footprint) %>%
                dplyr::rename(.bin = x) %>%
                    group_by(.bin) %>%
                        summarize(.freq = n())
}

gene.freq.tab <-
    model.tab.sorted %>% rename(.id = med.id) %>%
        select(.id, tss, tes) %>% unique() %>% dplyr::rename(lb = tss, ub = tes) %>%
            take.footprint.freq()

ld.freq.tab <-
    ld.bound.tab %>% rename(.id = med.id) %>%
        take.footprint.freq()

snp.freq.tab <- 
    sum.stat %>% select(rs, snp.loc) %>% unique() %>%
        mutate(.id = rs, lb = snp.loc, ub = snp.loc) %>%
            take.footprint.freq()


################################################################

empty.theme <- theme(axis.text = element_blank(), 
                     line = element_blank(), axis.line = element_blank(),
                     axis.ticks = element_blank(), panel.border = element_blank(),
                     legend.position = c(1,1), legend.justification = c(1, 1))

empty.x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())

color.grad2 <- scale_color_gradient2(low = 'blue', high = 'red', mid = 'gray40', guide = FALSE)

## Fig a. gene frequency
.freq.aes <- aes(xmin = .bin, xmax = .bin + 1e4, ymin = 0, ymax = .freq)

.tab <- ld.bound.tab %>% filter(med.id %in% show.genes) %>%
    left_join(model.tab.sorted) %>%
        as.data.frame() %>%
            arrange(tss) %>%
                mutate(y.loc = 1:n()) %>%
                    mutate(y.loc = -y.loc)

p1.a <-
    gg.plot() + ylab('LD blks') + .gwas.x.scale +
    geom_rect(data = ld.freq.tab, .freq.aes, fill = '#99FF66') +
    geom_segment(data = .tab, aes(x = lb, xend = ub, y = y.loc, yend = y.loc), size = 2, color = '#005500') +
    geom_text(data = .tab, aes(x = ub, y = y.loc, label = hgnc), hjust = 0, size = 3) +
    empty.x + scale_y_continuous(breaks = seq(1, max(ld.freq.tab$.freq), 3))

p1.b <-
    gg.plot() + .gwas.x.scale + ylab('genes') +
    geom_rect(data = gene.freq.tab, .freq.aes, fill = 'gray50', alpha = .7) +
    geom_segment(data = .tab, aes(x = snp.lb, xend = snp.ub, y = y.loc, yend = y.loc), size = 1, color = 'gray20') +
    geom_text(data = .tab, aes(x = snp.ub, y = y.loc, label = hgnc), hjust = 0, size = 3) +
    empty.x + scale_y_continuous(breaks = seq(1, max(gene.freq.tab$.freq), 2))

## Fig b. PVE 
.pve.aes <- aes(x = (tss+tes)/2, xend = (tss+tes)/2, y = 0, yend = PVE.glob, color = theta/theta.se)
.pve.aes.lab <- aes(x = (tss+tes)/2, y = PVE.glob, label = signif(PVE.glob, 2)) 

p2 <-
    gg.plot() + ylab('var. exp.') + xlab('') + .gwas.x.scale +
    geom_segment(data = model.tab.sorted, .pve.aes, size = 2) +
    geom_text_repel(data = model.tab.sorted %>% filter(med.id %in% show.genes), .pve.aes.lab, size = 3) +
    color.grad2 +
    empty.x +
    scale_y_continuous(breaks = c(1e-6, 1e-5, 1e-4, 1e-3), trans = 'sqrt')

## Fig b. mediation posterior probability
.med.aes <- aes(x = (tss+tes)/2, xend = (tss+tes)/2, y = 0, yend = lodds, color = theta/theta.se)

p3 <- 
    gg.plot(model.tab.sorted) + ylab('logit mediation prob') + xlab('') + .gwas.x.scale +
    geom_hline(yintercept = med.cutoff, lty = 2) +
    geom_hline(yintercept = 0, lty = 1)+
    geom_segment(.med.aes, size = 2) +
    geom_text_repel(data = model.tab.sorted %>% filter(med.id %in% show.genes), aes(x = (tss+tes)/2, y = lodds, label = hgnc), size = 3) +
    color.grad2 +
    empty.x

## Fig c. TWAS
.twas.aes <- aes(x = (tss+tes)/2, xend = (tss+tes)/2, y = 0,
                 yend = pmin(-log10(2 * pnorm(abs(twas), lower.tail = FALSE)), 20),
                 color = twas, label = hgnc)

.twas.aes.lab <- aes(x = (tss+tes)/2,
                     y = pmin(-log10(2 * pnorm(abs(twas), lower.tail = FALSE)), 20),
                     color = twas, label = hgnc)

p4 <-
    gg.plot(model.tab.sorted) + ylab('-log10 TWAS P') + xlab('') + .gwas.x.scale +
    geom_hline(yintercept = -log10(twas.cutoff), lty = 2) +
    geom_segment(.twas.aes, size = 2) +
    geom_text_repel(data = model.tab.sorted %>% filter(med.id %in% show.genes), .twas.aes.lab, size = 3) +
    color.grad2 +
    empty.x

## Fig d. snp frequency
p5 <- 
    gg.plot(snp.freq.tab) + ylab('QTLs') +
    geom_rect(.freq.aes, color = 'gray50') +
    empty.x  + scale_y_reverse() + .gwas.x.scale.top


## Fig e. GWAS
plt.gwas <-
    gg.plot() + ylab('-log10 GWAS P') +
    theme(axis.title.x = element_blank()) +
    geom_rect(data=gwas.block.df, aes(xmin=blk.loc, xmax=blk.loc+blk.width, ymin=0, ymax=gwas.p),
              fill = 'gray80') +
    geom_point(data=gwas.df, aes(x=snp.loc, y = gwas.p), size = .1)

if(max(gwas.df$gwas.p) > -log10(5e-8)) {
    plt.gwas <- plt.gwas + geom_hline(yintercept = -log10(5e-8), color = 'red')
}

plt.gwas <- plt.gwas + scale_y_reverse() + .gwas.x.scale.top

ld.hh <- (length(show.genes) + max(ld.freq.tab$.freq)) * .05
gene.hh <- pmin((length(show.genes) + max(gene.freq.tab$.freq)) * .2, .5)

hh.vec <- c(ld.hh, gene.hh, 1, 1.5, 1, 1.5, .5)

hh <- sum(hh.vec)  +  1
ww <- pmax(pmin(ceiling((gwas.range[2] - gwas.range[1]) / 1e6 * 5), 15), 5)

gg <- grid.vcat(list(p1.a, p1.b, p2, p3, p4, plt.gwas, p5), heights = hh.vec)

ggsave(filename = out.file, plot = gg, width = ww, height = hh, units = 'in', limitsize = FALSE)

system('rm -r ' %&&% temp.dir)
log.msg('Finished figure drawing\n\n\n')
