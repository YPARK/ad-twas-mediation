#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
library(dplyr)
library(tidyr)
library(qvalue)
library(readr)
library(pander)
library(goseq)
library(xlsx)
source('figure.util.R')
source('mediation.R')

read.chr <- function(in.hdr, in.cols, file.name) {
    .files <- in.hdr %&&% 1:22 %&&% '.' %&&% file.name %&&% '.gz'
    .dat <- lapply(.files, read_tsv, col_names = in.cols)
    return(do.call(rbind, .dat))
}


read.mediation <- function(in.hdr) {

    in.cols <-
        c('chr', 'ld.lb', 'ld.ub', 'med.id',
          'theta', 'theta.se', 'lodds', 'hgnc',
          'tss', 'tes', 'strand',
          'best.gwas.z', 'best.gwas.rs', 'best.gwas.loc',
          'best.qtl.z', 'best.qtl.rs', 'best.qtl.loc', 'n.qtls',
          'v.med', 'v.med.tot', 'v.dir', 'v.resid')

    ## It is unstable to use mediation based on less than 50 QTLs
    ret <- read.chr(in.hdr, in.cols, 'mediation')

    return(ret)
}

read.null <- function(in.hdr) {
    in.cols <- 'null.lodds'
    return(read.chr(in.hdr, in.cols, 'null'))
}

make.pandoc.tab <- function(med.tab) {

    .out <- med.tab %>%
        mutate(theta = signif(theta, 2),
               theta.se = signif(theta.se, 2),
               pip = signif(1/(1+exp(-lodds)), 1),
               ld.gwas.z = signif(ld.gwas.z, 2),
               best.gwas.z = signif(best.gwas.z, 2),
               best.qtl.z = signif(best.qtl.z, 2),
               Gene = hgnc,
               pval = signif(pval, 1),
               qval = signif(qval, 1))

    .out <- .out %>%
        mutate(TSS = ifelse(strand == '+',
                   format(round(tss/1e3), big.mark=','),
                   format(round(tes/1e3), big.mark=',')),
               TES = ifelse(strand == '-',
                   format(round(tss/1e3), big.mark=','),
                   format(round(tes/1e3), big.mark=',')),
               ld.gwas.loc = format(round(ld.gwas.loc/1e3), big.mark=','),
               best.qtl.loc = format(round(best.qtl.loc/1e3), big.mark=','),
               PVE.loc = signif(PVE, 2),
               PVE.glob = signif(PVE.glob, 2))

    .out <- .out %>%
        mutate(Mediation = theta %&&% ' (' %&&% theta.se %&&% ', ' %&&% pip %&&% ')')
    .out <- .out %>%
        mutate(GWAS = ld.gwas.rs %&&% ' (' %&&% ld.gwas.z %&&% ', ' %&&% ld.gwas.loc %&&% ')') %>%
            mutate(QTL = best.qtl.rs %&&% ' (' %&&% best.qtl.z %&&% ', ' %&&% best.qtl.loc %&&% ')')

    .out <- .out %>%
        dplyr::select(Gene, chr, TSS, Mediation, GWAS, PVE.loc, PVE.glob)

    ret <- pandoc.table.return(.out, row.names = FALSE, style = 'simple',
                               split.tables = 200, digits = 2)
    return(ret)
}


take.goseq <- function(gene.test.tab,
                       qval.cutoff = 1e-2,
                       pve.cutoff = 1e-2,
                       go.pval.cutoff = 1e-2) {

    de.genes <- gene.test.tab %>% dplyr::filter(qval < qval.cutoff, PVE > pve.cutoff)
    genes.tot <- gene.test.tab$med.id %>% unique()
    gene.vector <- as.integer(genes.tot %in% de.genes$med.id)
    print(mean(gene.vector))

    names(gene.vector) <- genes.tot
    pwf <- nullp(gene.vector, 'hg19', 'ensGene', plot.fit = FALSE)
    go.out <- goseq(pwf, 'hg19', 'ensGene')

    go.analysis.tab <- go.out %>%
        dplyr::rename(go.pval = over_represented_pvalue) %>%
            arrange(go.pval, ontology) %>% as.data.frame()

    take.ensg.cat <- function(ensg, go.tab) {
        go.cat <- go.tab$category
        ensg2cat <- getgo(ensg, 'hg19', 'ensGene')

        .parse <- function(idx, cat.list, ensg.name) {
            if(is.na(ensg.name[idx])) return(NULL)
            .cat <- intersect(cat.list[[idx]], go.cat)
            if(length(.cat) < 1) return(NULL)

            return(data.frame(ensg = ensg.name[idx],  category = .cat))
        }

        ensg2cat.list <- lapply(seq_along(ensg2cat),
                                .parse,
                                cat.list = ensg2cat,
                                ensg.name = names(ensg2cat))

        .merge <- function(a, b) {
            return(rbind(a, b))
        }

        ret <- do.call(rbind, Filter(Negate(is.null), ensg2cat.list))
        return(ret)
    }

    go.tab.sig <- go.analysis.tab %>% dplyr::filter(go.pval < go.pval.cutoff)
    ensg2cat.tab <- take.ensg.cat(de.genes$med.id, go.tab.sig)

    return(list(go.tab = go.analysis.tab, ensg2cat = ensg2cat.tab))
}

draw.go.tab <- function(go.summary, goseq.result, go.pdf.file) {

    go.orders <- order.pair(go.summary %>% dplyr::select(hgnc, term) %>%
                            dplyr::rename(row = hgnc, col = term) %>%
                            dplyr::mutate(weight = 1))

    term.short <- sapply(go.orders$cols, function(s) substr(s, start = 1, stop = 90))

    plt.df <- go.summary
    plt.df$hgnc <- factor(plt.df$hgnc, go.orders$rows)
    plt.df$term <- factor(plt.df$term, go.orders$cols, term.short)

    tab.df <- goseq.result$go.tab %>% filter(term %in% go.summary$term)
    tab.df$term <- factor(tab.df$term, go.orders$cols, term.short)

    p1 <- gg.plot(plt.df, aes(x = hgnc, y = term)) +
        geom_tile(fill = 'black') +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.y = element_blank()) + xlab('Genes')

    p2 <-
        gg.plot(tab.df, aes(y = term, yend = term, x = 0, xend = -log10(go.pval))) +
            geom_segment(color = 'orange', size = 2) + xlab('-log10 P') +
                theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank())

    gg <- grid.hcat(list(p1, p2), widths = c(8, 2))
    hh <- length(go.orders$cols) * .1 + 1
    ww <- max(sapply(go.orders$cols, nchar)) * .05 + hh
    ggsave(filename = go.pdf.file, plot = gg, height = hh, width = ww, units = 'in', limitsize=FALSE)
}

write.tables <- function(qtl.data, qval.cutoff = 1e-2, pve.cutoff = 1e-2) {

    gwas.summary.stat <- read_tsv('stat/IGAP/ld.summary.txt.gz', col_names = TRUE) %>%
        dplyr::rename(ld.gwas.z = gwas.z, ld.gwas.rs = rs, ld.gwas.p = gwas.p, ld.gwas.loc = snp.loc)

    med.hdr <- 'result/mediation/' %&&% qtl.data %&&% '_IGAP_eqtl_hs-lm_'
    null.hdr <- 'result/null/' %&&% qtl.data %&&% '_IGAP_eqtl_hs-lm_'

    out.dir <- 'tables/genes/' %&&% qtl.data
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    out.significant.file <- out.dir %&&% '/significant_genes.txt.gz'
    out.ld.file <- out.dir %&&% '/significant_LD.txt.gz'
    out.priority.pandoc <- out.dir %&&% '/priority_genes.md'
    out.priority.xlsx <- out.dir %&&% '/priority_genes.xlsx'
    out.tot.file <- out.dir %&&% '/total_genes.txt.gz'

    out.subthreshold.pandoc <- out.dir %&&% '/subthreshold_genes.md'
    out.subthreshold.xlsx <- out.dir %&&% '/subthreshold_genes.xlsx'

    out.figure.pval.file <- out.dir %&&% '/pvalues.pdf'
    out.figure.qval.file <- out.dir %&&% '/qvalues.pdf'

    med.tab <- read.mediation(med.hdr)
    null.tab <- read.null(null.hdr)

    pval <- empPvals(stat = med.tab$lodds, stat0 = null.tab$null.lodds)
    q.obj <- qvalue(pval)

    ## calculate global total variance
    ld.tab <- med.tab %>% dplyr::group_by(chr, ld.lb, ld.ub) %>%
        slice(which.max(v.med)) %>%
            arrange(chr, ld.lb) %>%
                dplyr::select(chr, ld.lb, ld.ub, v.med, v.med.tot, v.dir, v.resid) %>%
                    arrange(desc(v.med.tot)) %>%
                        as.data.frame()
    
    sum.tab <- ld.tab %>%
        summarize(v.med = sum(v.med.tot), v.dir = sum(v.dir), v.resid = sum(v.resid))

    v.tot <- sum(sum.tab)

    med.tab <- cbind(med.tab, pval = pval, qval = q.obj$qvalues) %>%
        mutate(gwas.p = pmin(l10.p.two(abs(best.gwas.z)),20)) %>%
            mutate(PVE = signif(v.med / (max(v.med, v.med.tot) + v.dir + v.resid), 2)) %>%
                mutate(PVE.glob = signif(v.med / v.tot, 4)) %>%
                    left_join(gwas.summary.stat)

    pdf(file = out.figure.pval.file)
    hist(pval, 30)
    dev.off()

    pdf(file = out.figure.qval.file)
    plot(q.obj)
    dev.off()

    sig.tab <- med.tab %>% dplyr::filter(qval < qval.cutoff) %>%
        arrange(chr, tss)

    if(nrow(sig.tab) < 1) {
        return(NULL)
    }

    for(go.cutoff in c(1e-2, 1e-4)) {

        goseq.result <- take.goseq(med.tab,
                                   qval.cutoff = go.cutoff,
                                   pve.cutoff = 0,
                                   go.pval.cutoff = 1e-2)
        
        if(!is.null(goseq.result$ensg2cat)) {
            go.summary <-
                sig.tab %>%
                    dplyr::select(med.id, chr, hgnc) %>%
                        dplyr::rename(ensg = med.id) %>%
                            left_join(goseq.result$ensg2cat) %>% na.omit() %>%
                                left_join(goseq.result$go.tab) %>%
                                    dplyr::select(chr, ensg, hgnc, category, term, ontology) %>%
                                        unique() %>%
                                            arrange(chr)
            
            go.tab.file <- out.dir %&&% '/GO_' %&&% as.character(go.cutoff) %&&% '.txt.gz'
            write_tsv(goseq.result$go.tab %>% dplyr::filter(numInCat <= 100), path = go.tab.file)
            
            ## show top 50 GO tables
            .tab <- goseq.result$go.tab %>% dplyr::filter(numInCat <= 100) %>% head(50) %>%
                dplyr::mutate(category = gsub(pattern = 'GO:', replacement = '', category)) %>%
                    dplyr::rename(pval=go.pval) %>%
                        dplyr::mutate(pval=signif(pval,2)) %>%
                            dplyr::mutate(numDEInCat = numDEInCat %&&% "/" %&&% numInCat) %>%
                                dplyr::rename(num=numDEInCat, "-"=ontology) %>%
                                    dplyr::select(-under_represented_pvalue,-numInCat)
            
            
            .ret <- pandoc.table.return(.tab, row.names = FALSE, style = 'simple',
                                        justify = 'right',
                                        split.tables = 200, digits = 2)
            
            cat(.ret, file = out.dir %&&% '/GO_top50_' %&&% as.character(go.cutoff) %&&% '.md')
            
            if(nrow(go.summary) > 0) {
                go.pdf.file <- out.dir %&&% '/GO' %&&% as.character(go.cutoff) %&&% '.pdf'
                draw.go.tab(go.summary, goseq.result, go.pdf.file)
            }
        }
    }

    write_tsv(sig.tab, path = out.significant.file)

    ## priority
    .temp.tab <- med.tab %>%
        dplyr::filter(qval < qval.cutoff, PVE.glob >= pve.cutoff) %>%
            arrange(chr, tss)

    .temp <- make.pandoc.tab(.temp.tab)

    cat(.temp, file = out.priority.pandoc)

    if(nrow(.temp.tab) > 0) {
        write.xlsx(.temp.tab %>% as.data.frame(), out.priority.xlsx,
                   sheetName = sprintf('Significant q < 1e-4 and PVEg > %e', pve.cutoff))
    }

    ## overall subthreshold
    .temp.tab <- med.tab %>%
        dplyr::filter(qval < 1e-4, ld.gwas.p < 1e-4) %>%
            arrange(chr, tss)

    .temp <- make.pandoc.tab(.temp.tab)

    cat(.temp, file = out.subthreshold.pandoc)

    if(nrow(.temp.tab) > 0) {
        write.xlsx(.temp.tab %>% as.data.frame(), out.subthreshold.xlsx,
                   sheetName = 'Significant q < 1e-4 and GWAS p < 1e-4')
    }

    .temp <- med.tab %>% dplyr::filter(qval < qval.cutoff) %>%
        dplyr::select(chr, ld.lb, ld.ub) %>% unique()

    write_tsv(med.tab %>% right_join(.temp),
              path = out.ld.file)

    write_tsv(med.tab, path = out.tot.file)
}

write.tables('full-0', qval.cutoff = 1e-4, pve.cutoff = 5e-4)

##qtl.data.vec <- c('full-' %&&% c(5, 3, 0), 'hic-' %&&% c(5, 3, 0))
##for(qtl.data in qtl.data.vec) {
##    write.tables(qtl.data, qval.cutoff = 5e-2, pve.cutoff = 5e-2)
##}
