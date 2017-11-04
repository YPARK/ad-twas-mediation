#!/usr/bin/env Rscript

options(stringsAsFators = FALSE)

source('util.R')
library(dplyr)
library(tidyr)
library(qvalue)
library(readr)
library(pander)
library(goseq)
source('figure.util.R')

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

    return(read.chr(in.hdr, in.cols, 'mediation'))
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
               best.gwas.z = signif(best.gwas.z, 1),
               best.qtl.z = signif(best.gwas.z, 1),
               PVE = signif(v.med / (v.med.tot + v.dir + v.resid), 2),
               Gene = hgnc)

    .out <- .out %>%
        mutate(TSS = ifelse(strand == '+',
                   format(tss, big.mark=','),
                   format(tes, big.mark=',')),
               TES = ifelse(strand == '-',
                   format(tss, big.mark=','),
                   format(tes, big.mark=',')))

    .out <- .out %>%
        mutate(Mediation = theta %&&% ' (' %&&% theta.se %&&% ', ' %&&% pip %&&% ')')
    .out <- .out %>%
        mutate(GWAS = best.gwas.rs %&&% ' (' %&&% best.gwas.z %&&% ')')

    .out <- .out %>%
        dplyr::select(Gene, chr, TSS, TES, Mediation, GWAS, PVE, pval, qval)

    ret <- pandoc.table.return(.out, row.names = FALSE, style = 'simple',
                               split.tables = 200, digits = 2)
    return(ret)
}


take.goseq <- function(gene.test.tab, qval.cutoff = 1e-2,
                       go.pval.cutoff = 1e-2) {

    de.genes <- gene.test.tab %>% dplyr::filter(qval < qval.cutoff)
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

write.tables <- function(qtl.data, qval.cutoff = 1e-2) {

    med.hdr <- 'result/mediation/' %&&% qtl.data %&&% '_IGAP_eqtl_hs-lm_'
    null.hdr <- 'result/null/' %&&% qtl.data %&&% '_IGAP_eqtl_hs-lm_'

    out.dir <- 'tables/genes/' %&&% qtl.data
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    out.significant.file <- out.dir %&&% '/significant_genes.txt.gz'
    out.ld.file <- out.dir %&&% '/significant_LD.txt.gz'
    out.significant.pandoc <- out.dir %&&% '/significant_genes.md'

    out.figure.pval.file <- out.dir %&&% '/pvalues.pdf'
    out.figure.qval.file <- out.dir %&&% '/qvalues.pdf'

    med.tab <- read.mediation(med.hdr)
    null.tab <- read.null(null.hdr)

    pval <- empPvals(stat = med.tab$lodds, stat0 = null.tab$null.lodds)
    q.obj <- qvalue(pval)
    med.tab <- cbind(med.tab, pval = pval, qval = q.obj$qvalues)

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

    goseq.result <- take.goseq(med.tab, qval.cutoff = qval.cutoff,
                               go.pval.cutoff = 1e-2)

    go.summary <-
        sig.tab %>%
            dplyr::select(med.id, chr, hgnc) %>%
                dplyr::rename(ensg = med.id) %>%
                    left_join(goseq.result$ensg2cat) %>% na.omit() %>%
                        left_join(goseq.result$go.tab) %>%
                            dplyr::select(chr, ensg, hgnc, category, term, ontology) %>%
                                unique() %>%
                                    arrange(chr)

    go.tab.file <- out.dir %&&% '/GO.txt.gz'
    write_tsv(goseq.result$go.tab, path = go.tab.file)

    ## show top 20 GO tables
    .tab <- goseq.result$go.tab %>% head(20) %>%
        dplyr::select(-under_represented_pvalue)
    .ret <- pandoc.table.return(.tab, row.names = FALSE, style = 'simple',
                                split.tables = 200, digits = 2)

    cat(.ret, file = out.dir %&&% '/GO_top20.md')

    if(nrow(go.summary) > 0) {
        go.pdf.file <- out.dir %&&% '/GO.pdf'
        draw.go.tab(go.summary, goseq.result, go.pdf.file)
    }

    write_tsv(sig.tab, path = out.significant.file)

    .temp.tab <- med.tab %>% dplyr::filter(qval < qval.cutoff) %>%
        arrange(chr, tss)
    .temp <- make.pandoc.tab(.temp.tab)
    cat(.temp, file = out.significant.pandoc)

    .temp <- med.tab %>% dplyr::filter(qval < qval.cutoff) %>%
        dplyr::select(chr, ld.lb, ld.ub) %>% unique()

    write_tsv(med.tab %>% right_join(.temp),
              path = out.ld.file)

}

qtl.data.vec <- c('full-' %&&% seq(2, 8, 2), 'hic-' %&&% seq(2, 8, 2))

for(qtl.data in qtl.data.vec) {
    write.tables(qtl.data)
}
