#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
library(dplyr)
library(qvalue)
library(readr)

dir.create('tables', recursive = TRUE, showWarnings = FALSE)

read.chr <- function(in.hdr, in.cols, file.name) {
    .files <- in.hdr %&&% 1:22 %&&% '.' %&&% file.name %&&% '.gz'
    .dat <- lapply(.files, read_tsv, col_names = in.cols)
    return(do.call(rbind, .dat))
}

herit.cpgs.tab <- read_tsv('heritable.cpgs.qvalue.txt.gz', col_names = TRUE)
h.cpgs <- herit.cpgs.tab %>% filter(q.val < 0.01)
h.cpgs <- h.cpgs$cg

herit.genes.tab <- read_tsv('heritable.genes.qvalue.txt.gz', col_names = TRUE)
h.genes <- herit.genes.tab %>% filter(q.val < 0.01) %>% dplyr::select(ENSG)
h.genes <- h.genes$ENSG

################################################################
qtl.cutoff <- abs(qnorm(1e-4 / 2)) ## 10^-4

################################################################
meth.cols <- c('chr', 'ld.lb', 'ld.ub', 'cg', 'meth.theta',
               'meth.theta.var', 'meth.lodds', 'cg.loc',
               'meth.best.gwas.z', 'meth.best.gwas.rs', 'meth.best.gwas.loc',
               'meth.best.qtl.z', 'meth.best.qtl.rs', 'meth.best.qtl.loc')

meth.tab <- read.chr('joint/IGAP_rosmap_hs-lm_', meth.cols, 'cpg') %>%
    filter(cg %in% h.cpgs, abs(meth.best.qtl.z) > qtl.cutoff) %>%
    mutate(meth.theta.var = sqrt(meth.theta.var)) %>%
    rename(meth.theta.se = meth.theta.var)

joint.null.tab <- read.chr('joint/IGAP_rosmap_hs-lm_', 'lodds.null', 'null')

gene.cols <- c('chr', 'ld.lb', 'ld.ub', 'ensg', 'gene.theta',
               'gene.theta.var', 'gene.lodds', 'hgnc', 'tss', 'tes', 'strand',
               'gene.best.gwas.z', 'gene.best.gwas.rs', 'gene.best.gwas.loc',
               'gene.best.qtl.z', 'gene.best.qtl.rs', 'gene.best.qtl.loc')

gene.tab <- read.chr('joint/IGAP_rosmap_hs-lm_', gene.cols, 'genes') %>%
    filter(ensg %in% h.genes, abs(gene.best.qtl.z) > qtl.cutoff) %>%
    mutate(gene.theta.var = sqrt(gene.theta.var)) %>%
    rename(gene.theta.se = gene.theta.var)


################################################################
p.cutoff <- 0.05 / (nrow(meth.tab) + nrow(gene.tab))

lodds.stat <- c(meth.tab$meth.lodds, gene.tab$gene.lodds)
p.val <- empPvals(lodds.stat, joint.null.tab$lodds.null)
q.val <- qvalue(p.val)

meth.idx <- 1:nrow(meth.tab)
gene.idx <- -meth.idx

meth.tab <- cbind(meth.tab, p.val = p.val[meth.idx], q.val = q.val$qvalue[meth.idx])
gene.tab <- cbind(gene.tab, p.val = p.val[gene.idx], q.val = q.val$qvalue[gene.idx])
                  

################################################################
## Identify CpGs linked to Transcript
read.m2t <- function(in.hdr = 'bootstrap/rosmap_hs-lm_',
                     heritable.cpgs = h.cpgs,
                     heritable.genes = h.genes) {
    in.cols <- c('cg', 'ensg', 'm2t.theta', 'm2t.theta.var', 'm2t.lodds')    
    ret <- read.chr(in.hdr, in.cols, 'm2t') %>%
        filter(cg %in% heritable.cpgs, ensg %in% heritable.genes)

    ret <- ret %>%
        mutate(m2t.theta.se = sqrt(m2t.theta.var))

    return(ret)
}

read.m2t.null <- function(in.hdr = 'bootstrap/rosmap_hs-lm_') {
    in.cols <- 'm2t.null'
    read.chr(in.hdr, in.cols, 'm2t-null')
}

m2t.tab <- read.m2t(heritable.cpgs = unique(meth.tab$cg),
                    heritable.genes = unique(gene.tab$ensg))

m2t.null.tab <- read.m2t.null()

pval.m2t <- empPvals(m2t.tab$m2t.lodds, m2t.null.tab$m2t.null)

qval.m2t.obj <- qvalue(pval.m2t)

m2t.tab <- cbind(m2t.tab, m2t.pval = qval.m2t.obj$pvalues, m2t.qval = qval.m2t.obj$qvalues)

## Take links at FDR 10%
.temp.gene <- gene.tab %>%
    dplyr::select(chr, ld.lb, ld.ub, ensg, hgnc, tss, tes,
                  gene.best.qtl.z, gene.best.qtl.rs, gene.best.qtl.loc)

.temp.meth <- meth.tab %>%
    dplyr::select(cg, cg.loc, meth.best.qtl.z, meth.best.qtl.rs, meth.best.qtl.loc)

m2t.tab.fdr10 <- m2t.tab %>% filter(m2t.qval < 0.1) %>%
    right_join(.temp.gene, by = 'ensg') %>% na.omit() %>%
    right_join(.temp.meth, by = 'cg') %>% na.omit() %>%
    group_by(chr, ld.lb, ld.ub, cg, ensg) %>%
    slice(which.max(m2t.lodds))

################################################################
write.tab.named(meth.tab, file = gzfile('tables/joint_meth.txt.gz'))
write.tab.named(gene.tab, file = gzfile('tables/joint_gene.txt.gz'))
write.tab.named(m2t.tab.fdr10, file = gzfile('tables/joint_m2t_fdr10.txt.gz'))

