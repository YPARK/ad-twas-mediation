#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

ld.file <- argv[1]                      # e.g., ld.file = 'stat/IGAP/ld/hs-fqtl/19.obs_ld.gz'
sum.file <- argv[2]                     # e.g., sum.file = 'stat/IGAP/hs-fqtl/19.obs_bed.gz'
plink.hdr <- argv[3]                    # e.g., plink.hdr = '1kg/chr19' # geno/rosmap1709-chr19'
ld.idx <- as.integer(argv[4])           # e.g., ld.idx = 117
gwas.sample.size <- as.numeric(argv[5]) # e.g., gwas.sample.size = 74000
qtl.sample.size <- as.numeric(argv[6])  # e.g., qtl.sample.size = 300
out.hdr <- argv[7]                      # e.g., 'temp'

dir.create(dirname(out.hdr), recursive = TRUE)

zqtl.num.boot <- 111

source('util.R')
options(stringsAsFactors = FALSE)

z.out.file <- out.hdr %&&% '.zqtl.gz'
egger.out.file <- out.hdr %&&% '.egger.gz'
rot.egger.out.file <- out.hdr %&&% '.rot-egger.txt.gz'

if(file.exists(z.out.file)) {
    log.msg('File exists : %s\n\n\n', z.out.file)
    q()
}

library(zqtl)
library(dplyr)
library(methods)
library(broom)
source('mediation.R')

################################################################
## temporary directory
temp.dir <- system('mktemp -d ' %&&% out.hdr %&&% '-temp.XXXX',
                   intern = TRUE, ignore.stderr = TRUE)

ld.tab <- read.table(ld.file, col.names = c('chr', 'ld.lb', 'ld.ub', 'n.snp.ld', 'n.qtl.ld'))
ld.info <- ld.tab[ld.idx, ]

plink <- subset.plink(ld.info, temp.dir)
x.bim <- data.frame(plink$BIM, x.pos = 1:nrow(plink$BIM))

sum.stat.out <- extract.sum.stat(ld.info, sum.file, x.bim, temp.dir)
sum.stat <- sum.stat.out$sum.stat
genes <- sum.stat.out$genes

zqtl.data <- make.zqtl.data(plink, sum.stat, genes)

################################################################
## 1. run gene-by-gene Egger regression
egger.tab <- sum.stat %>% group_by(ensg) %>% do(egger = mr.egger(.)) %>%
    tidy(egger) %>% filter(.rownames != 'b.qtl') %>% dplyr::select(-.rownames) %>%
        rename(egger.theta=Estimate, egger.se=Std..Error, egger.t=t.value, egger.p.val=Pr...t..)

egger.tab <- egger.tab %>% left_join(genes, by = 'ensg')

log.msg('Finished EGGER\n\n\n')

################################################################
## 2. run rotated Egger regression
rot.egger.tab <- take.rot.egger(zqtl.data, genes)

log.msg('Finished Rotated EGGER\n\n\n')

################################################################
## 3. zqtl
vb.opt <- list(pi = -2, tau = -4, do.hyper = FALSE, tol = 1e-8, gammax = 1e3,
               vbiter = 7500, do.stdize = TRUE, eigen.tol = 1e-2,
               rate = 1e-2, decay = -1e-2,
               print.interv = 1500, nboot = zqtl.num.boot)

z.out <- fit.med.zqtl(zqtl.data$gwas.theta, zqtl.data$gwas.se,
                      zqtl.data$qtl.theta, zqtl.data$qtl.se,
                      X = zqtl.data$X, n = gwas.sample.size,
                      n.med = qtl.sample.size, options = vb.opt)

boot.tab <- melt.effect(list(lodds.boot.mean = z.out$boot.mediated.lodds.mean,
                             lodds.boot.se = sqrt(z.out$boot.mediated.lodds.var)),
                        genes$ensg, 'boot') %>%
                            rename(ensg = Var1) %>% dplyr::select(-Var2)

out.tab <- melt.effect(z.out$param.mediated, genes$ensg, 'obs') %>%
    rename(ensg = Var1) %>% dplyr::select(-Var2) %>%
        left_join(genes, by = 'ensg') %>%
            left_join(boot.tab, by = 'ensg') %>%
                left_join(data.frame(ensg = genes$ensg, emp.pval = z.out$mediation.pvalue),
                          by = 'ensg')

out.tab <- out.tab %>%
    mutate(z.pval = pnorm((lodds - lodds.boot.mean) / (1e-2 + lodds.boot.se), lower.tail = FALSE))

write.tab(out.tab, file = gzfile(z.out.file))
write.tab(egger.tab, file = gzfile(egger.out.file))
write.tab(rot.egger.tab, file = gzfile(rot.egger.out.file))

system('rm -r ' %&&% temp.dir)
log.msg('Finished zQTL\n\n\n')
