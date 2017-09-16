## subset summary stat data
extract.sum.stat <- function(ld.info, sum.file, x.bim, temp.dir, qtl.cutoff = 2) {

    bedtools.result <- temp.dir %&&% '/stat.ucsc_bed.gz'
    bedtools.hdr <- 'source /broad/software/scripts/useuse > /dev/null; reuse -q BEDTools'
    bedtools.cmd <- paste('bedtools', 'intersect', '-a', sum.file, '-b', 'stdin -wa | gzip > ',
                          bedtools.result)
    sum.cmd <- sprintf('%s; printf "%s" | %s', bedtools.hdr, paste(ld.info, collapse = '\t'),
                       bedtools.cmd)
    
    system(sum.cmd, intern = TRUE, ignore.stderr = TRUE)
    
    sum.stat.cols <- c('chr', 'snp.loc.1', 'snp.loc', 'rs',
                       'qtl.a1', 'qtl.a2', 'qtl.theta', 'qtl.z',
                       'gwas.theta', 'gwas.se', 'gwas.z',
                       'ld', 'ensg', 'hgnc', 'tss', 'tes', 'strand')
    
    sum.stat <- read.table(bedtools.result, col.names = sum.stat.cols, sep = '\t') %>%
        mutate(qtl.se = qtl.theta / qtl.z) %>%
            left_join(x.bim, by = c('chr', 'rs', 'snp.loc'))

    sum.stat <- sum.stat %>%
        filter(((plink.a1 == qtl.a1) & (plink.a2 == qtl.a2)) | ((plink.a1 == qtl.a2) & (plink.a2 == qtl.a1))) %>%
            mutate(gwas.z.flip = if_else(qtl.a1 != plink.a1, -gwas.z, gwas.z)) %>%
                mutate(gwas.theta.flip = if_else(qtl.a1 != plink.a1, -gwas.theta, gwas.theta)) %>%
                    mutate(qtl.z.flip = if_else(qtl.a1 != plink.a1, -qtl.z, qtl.z)) %>%
                        mutate(qtl.theta.flip = if_else(qtl.a1 != plink.a1, -qtl.theta, qtl.theta))
    
    sum.stat <- sum.stat %>%
        rename(.gwas.z = gwas.z, .gwas.theta = gwas.theta, .qtl.z = qtl.z, .qtl.theta = qtl.theta) %>%
            rename(gwas.z = gwas.z.flip, gwas.theta = gwas.theta.flip,
                   qtl.z = qtl.z.flip, qtl.theta = qtl.theta.flip)

    sum.stat <- sum.stat %>% filter(abs(qtl.z) >= qtl.cutoff)

    genes <-
        sum.stat %>% group_by(ensg) %>%
            slice(which.max(abs(gwas.z))) %>%
                dplyr::select(ensg, hgnc, tss, tes, strand, gwas.theta, gwas.z) %>%
                    arrange(tss)
    
    sum.stat <- sum.stat %>%
        mutate(y.pos = match(ensg, genes$ensg)) %>%
            na.omit()

    return(list(sum.stat = sum.stat, genes = genes))
}

subset.plink <- function(ld.info, temp.dir) {

    ## read plink
    plink.lb <- ld.info$ld.lb
    plink.ub <- ld.info$ld.ub
    chr <- ld.info$chr

    plink.cmd <- sprintf('./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --chr %d --from-bp %d --to-bp %d --out %s', plink.hdr, chr, plink.lb, plink.ub, temp.dir %&&% '/plink')
    system(plink.cmd)

    plink <- read.plink(temp.dir %&&% '/plink')
    colnames(plink$BIM) <- c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')
    return(plink)
}

make.zqtl.data <- function(plink, sum.stat, genes) {

    X <- plink$BED
    qtl.theta <- matrix(NA, nrow = ncol(X), ncol = nrow(genes))
    qtl.se <- matrix(NA, nrow = ncol(X), ncol = nrow(genes))

    qtl.idx <- sum.stat %>% dplyr::select(x.pos, y.pos) %>% as.matrix()
    qtl.theta[qtl.idx] <- sum.stat$qtl.theta
    qtl.se[qtl.idx] <- sum.stat$qtl.se

    gwas.stat <- sum.stat %>% group_by(x.pos) %>%
        summarize(gwas.theta = mean(gwas.theta), gwas.se = mean(gwas.se))

    gwas.theta <- matrix(0, nrow = ncol(X), ncol = 1)
    gwas.se <- matrix(0, nrow = ncol(X), ncol = 1)
    gwas.theta[as.integer(gwas.stat$x.pos), 1] <- gwas.stat$gwas.theta
    gwas.se[as.integer(gwas.stat$x.pos), 1] <- gwas.stat$gwas.se

    valid.x.loc <- sort(unique(sum.stat$x.pos))

    qtl.theta <- qtl.theta %r% valid.x.loc
    qtl.se <- qtl.se %r% valid.x.loc
    gwas.theta <- gwas.theta %r% valid.x.loc
    gwas.se <- gwas.se %r% valid.x.loc

    X <- X %c% valid.x.loc

    return(list(X = X, qtl.theta = qtl.theta, qtl.se = qtl.se,
                gwas.theta = gwas.theta, gwas.se = gwas.se))
}

melt.effect <- function(effect.obj, .rnames, .cnames) {
    require(reshape2)
    melt.list <- lapply(seq_along(effect.obj),
                        function(mat.list, val.names, i) {
                            mat <- signif(mat.list[[i]], 4)
                            val <- val.names[[i]]
                            rownames(mat) <- .rnames
                            colnames(mat) <- .cnames
                            melt(mat, value.name = val) },
                        mat.list = effect.obj, val.names = names(effect.obj))

    ret <- Reduce(function(...) left_join(..., by = c('Var1', 'Var2')), melt.list) %>%
        mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
    return(ret)
}

l10.p.two <- function(...) pmin(- log10(2 * pnorm(..., lower.tail = FALSE)), 100)

l10.p.one <- function(...) pmin(- log10(pnorm(..., lower.tail = FALSE)), 100)

mr.egger <- function(tab, do.perm = FALSE) {

    qtl.theta <- tab$qtl.theta
    qtl.se <- tab$qtl.se
    gwas.theta <- tab$gwas.theta
    gwas.se <- tab$gwas.se
    
    ## Do the same thing as Bowden, Smith, Burgess, IJE (2015)

    b.gwas <- gwas.theta * sign(qtl.theta)
    b.qtl <- abs(qtl.theta)
    if(do.perm){
        b.qtl <- b.qtl[sample(length(b.qtl))]
    }

    coefficients(summary(lm(b.gwas ~ b.qtl, weights = 1/gwas.se^2)))
}

take.rot.egger <- function(zqtl.data, genes) {

    svd.out <- take.ld.svd(zqtl.data$X)
    V.t <- svd.out$V.t
    weight <- svd.out$D^2

    gwas.rot <- V.t %*% (zqtl.data$gwas.theta / zqtl.data$gwas.se)
    qtl.rot <- (zqtl.data$qtl.theta / zqtl.data$qtl.se)
    qtl.rot[is.na(zqtl.data$qtl.theta)] <- 0
    qtl.rot <- V.t %*% qtl.rot

    wlm.out <- apply(qtl.rot, 2, function(qrot, grot, d) {
        ret <- coefficients(summary(lm(grot ~ qrot, weights = 1/d^2)))
        ret[2, ]
    }, grot = gwas.rot[, 1], d = svd.out$D)

    data.frame(genes, t(wlm.out)) %>%
        rename(egger.theta=Estimate, egger.se=Std..Error, egger.t=t.value, egger.p.val=Pr...t..)        
}

