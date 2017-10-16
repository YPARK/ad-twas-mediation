#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

ld.idx <- as.integer(argv[1]) # e.g., ld.idx <- 1
ld.tab.file <- argv[2] # e.g., ld.tab.file <- 'tables/genes_ld_significant.txt.gz'
out.dir <- argv[3]     # e.g., out.dir <- 'temp' 

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

ld.tab <- read.table(ld.tab.file, header = TRUE)
ld.key <- ld.tab %>% select(chr, ld.lb, ld.ub) %>% unique()

if(ld.idx > nrow(ld.key)) q()

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

plink.hdr <- 'geno/rosmap1709-chr'
sum.file.dir <- 'stat/IGAP/data/hs-lm/'

ld.info <- ld.key[ld.idx, ]
chr <- as.integer(ld.info[1])

ld.out.file <- sprintf('%s/chr%d_ld%d_gwas_med_twas.pdf', out.dir, chr, ld.idx)

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

model.tab <- model.tab %>%
    mutate(med.p = -log10(pmax(p.val.direct, p.val.marginal)))

model.tab <- model.tab %>%
    mutate(pve = v.med / (v.med.tot + v.resid + v.dir))

zqtl.data <- make.zqtl.data(plink, sum.stat, genes)

## just run on genes that exist on the list
y.pos <- match(model.tab$ensg, genes$med.id)
zqtl.data$qtl.theta <- zqtl.data$qtl.theta %c% y.pos
zqtl.data$qtl.se <- zqtl.data$qtl.se %c% y.pos

## QTL information
gwas.z <- zqtl.data$gwas.theta[, 1] / zqtl.data$gwas.se[, 1]
gwas.df <- data.frame(snp.loc = zqtl.data$snps$snp.loc, gwas = gwas.z) %>%
    mutate(gwas.p = l10.p.two(abs(gwas)))

qtl.z <- zqtl.data$qtl.theta / zqtl.data$qtl.se
colnames(qtl.z) <- model.tab$hgnc
rownames(qtl.z) <- zqtl.data$snps$snp.loc

qtl.melt <- melt(qtl.z, varnames = c('snp.loc', 'hgnc'))
qtl.melt$hgnc <- as.character(qtl.melt$hgnc)
qtl.melt <- qtl.melt %>% left_join(gwas.df, by = 'snp.loc')
qtl.melt <- qtl.melt %>% na.omit() %>%
    mutate(qtl.p = l10.p.two(abs(value)))
qtl.melt$hgnc <- as.character(qtl.melt$hgnc)

blk.width <- 1e4

gwas.block.df <- gwas.df %>%
    mutate(blk.loc = floor(snp.loc/blk.width)*blk.width) %>%
        group_by(blk.loc) %>% slice(which.max(abs(gwas))) %>%
            mutate(gwas.p = l10.p.two(abs(gwas)))

blk.range <- range(c(gwas.block.df$blk.loc, model.tab$tss, model.tab$tes, max(gwas.block.df$blk.loc) + blk.width))

blk.length <- blk.range[2] - blk.range[1]

.gwas.w <- blk.length / blk.width * 0.02 + 1

.gwas.kb <- function() {
    function(x) format(x/1e3, big.mark=',')
}

.gwas.x.scale <- scale_x_continuous(limits = blk.range + c(-1000, 1000), expand = c(0, 0),
                                    labels = .gwas.kb())

gg.plot <- function(...) {
    ggplot(...) + theme_bw() + theme(plot.background = element_blank(),
                                     panel.background = element_blank(),
                                     strip.background = element_blank(),
                                     legend.background = element_blank())    
}

################################################################
## -- GWAS

plt.gwas <- gg.plot() + ylab('-log10 GWAS p-value') +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    .gwas.x.scale + 
    geom_point(data=gwas.df, aes(x=snp.loc, y = gwas.p), size = .1) +
    geom_rect(data=gwas.block.df, aes(xmin=blk.loc, xmax=blk.loc+blk.width, ymin=0, ymax=gwas.p),
              fill = 'gray', alpha = 0.5)

if(max(gwas.df$gwas.p) > -log10(5e-8)) {
    plt.gwas <- plt.gwas + geom_hline(yintercept = -log10(5e-8), color = 'red')
}

out.file <- sprintf('%s/chr%d_ld%d_gwas.pdf', out.dir, chr, ld.idx)
ggsave(file = out.file, plot = plt.gwas, width = .gwas.w, height = 3,
       useDingbats = FALSE, limitsize = FALSE)

################################################################
## -- QTLs
lb <- min(gwas.block.df$blk.loc)
ub <- max(gwas.block.df$blk.loc)

plt.qtl <- gg.plot(qtl.melt) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    geom_point(aes(x = snp.loc, y = qtl.p), size = .1, color = 'gray') +
    ylab('-log10 QTL p') +
    .gwas.x.scale


plt.edge <- ggplot() + theme_void() + .gwas.x.scale +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())

.temp <- model.tab %>% filter(med.p > - log10(2.5e-6)) %>%
    dplyr::select(hgnc, tss, tes)

if(nrow(.temp) > 0) {
    qtl.edges <- qtl.melt %>% right_join(.temp, by = 'hgnc') %>% na.omit()
    plt.edge <- plt.edge +
        geom_segment(data = qtl.edges, aes(x = snp.loc, xend = (tss+tes)/2, y = 1, yend = 0),
                     size = .1, color = 'gray', alpha = 0.5)
}


plt.edge <- plt.edge +
    geom_segment(data = model.tab, aes(x = (tss+tes)/2, xend = (tss+tes)/2, y = 0, yend = pve, size = pve),
                 color = 'orange') +
    geom_text_repel(data = model.tab, aes(x = (tss+tes)/2, label = hgnc, y = pve, size = pve),
                    angle = 45, segment.color = 'green', alpha = 0.5) +
    scale_size_continuous(guide = FALSE, range=c(2, 4))


v.tot <- model.tab %>% mutate(v.tot = v.med.tot + v.resid + v.dir) %>%
    slice(which.max(v.tot)) %>% dplyr::select(v.tot, v.med.tot, v.resid, v.dir) %>%
        melt()

v.tot.val <- v.tot %>% filter(variable == 'v.tot') %>% dplyr::select(value) %>% as.numeric()

v.tot <- v.tot %>% filter(variable != 'v.tot') %>% mutate(pve = value / v.tot.val)

get.plt <- function() {
    ret <- gg.plot() +
        theme(axis.text.x = element_blank(), axis.title = element_blank())
    return(ret)
}

p1 <- get.plt() +
    geom_bar(data = v.tot, aes(x='0', y = pve, fill = variable), stat = 'identity') +
        ylim(c(0, 1))

p2 <- get.plt() +
    geom_bar(data = model.tab %>% arrange(desc(pve)) %>% head(7),
             aes(x='1', y = pve, fill = hgnc), stat = 'identity') +
                 ylim(c(0, 1))

out.file <- sprintf('%s/chr%d_ld%d_pve.pdf', out.dir, chr, ld.idx)
pdf(file = out.file, width = 4, height = 3, useDingbats = FALSE)
grid.hcat(list(p1, p2))
dev.off()

plt.med <- gg.plot(model.tab) +
    geom_hline(yintercept = -log10(2.5e-6), lty = 2, color = 'gray') + 
    geom_segment(aes(x = (tss + tes)/2, xend = (tss + tes)/2, y = 0, yend = med.p, color = theta),
                 arrow = arrow(length = unit(.05, 'inches')),
                 size = 1) +
                     scale_y_reverse() + scale_color_gradient2(low = 'blue', high = 'red', mid = 'gray40', guide = FALSE) +
                         .gwas.x.scale + ylab('-log10 mediation p-value') +
    xlab('genomic location (kb)')


plt.med.twas <- gg.plot(model.tab) +
    geom_hline(yintercept = -log10(2.5e-6), lty = 2, color = 'gray') + 
    geom_segment(aes(x = (tss + tes)/2, xend = (tss + tes)/2, y = 0, yend = pmin(-log10(p.val.nwas), 10)),
                 color = '#009900', size = 2, alpha = 0.5) +
    geom_segment(aes(x = (tss + tes)/2, xend = (tss + tes)/2, y = 0, yend = med.p, color = theta),
                 arrow = arrow(length = unit(.05, 'inches')),
                 size = 1) +
                     scale_y_reverse() + scale_color_gradient2(low = 'blue', high = 'red', mid = 'gray40', guide = FALSE) +
                         .gwas.x.scale + ylab('-log10 mediation p-value') +
    xlab('genomic location (kb)')


out.file <- sprintf('%s/chr%d_ld%d_gwas_med_twas.pdf', out.dir, chr, ld.idx)
pdf(file = out.file, width = 8, height = 6, useDingbats = FALSE)
grid.vcat(list(plt.gwas, plt.edge, plt.med.twas), heights = c(.25, .25, .5))
dev.off()          


out.file <- sprintf('%s/chr%d_ld%d_gwas_med.pdf', out.dir, chr, ld.idx)
pdf(file = out.file, width = 8, height = 6, useDingbats = FALSE)
grid.vcat(list(plt.gwas, plt.edge, plt.med), heights = c(.25, .25, .5))
dev.off()          

system('rm -r ' %&&% temp.dir)
log.msg('Finished figure drawing\n\n\n')
