options(stringsAsFactors = FALSE)

log.msg <- function(...) {
    cat('[', date() ,']', sprintf(...), file = stderr())
}

`%&&%` <- function(a, b) paste(a, b, sep = '')

`%c%` <- function(a, b) a[, b, drop = FALSE]

`%r%` <- function(a, b) a[b, , drop = FALSE]

.NA <- function(nrow, ncol) {
    matrix(NA, nrow, ncol)
}

.rnorm <- function(nrow, ncol) {
    matrix(rnorm(nrow * ncol), nrow, ncol)
}

.sample <- function(nrow, ncol, pop = c(1,-1)) {
    matrix(sample(pop, nrow * ncol, TRUE), nrow, ncol)
}

fast.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / n.obs
    return(ret)
}

fast.z.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / sqrt(n.obs)
    return(ret)
}

get.marginal.qtl <- function(xx, yy, .melt = TRUE, TOL = 1e-8) {

    .xx <- scale(xx)
    .yy <- scale(yy)

    beta <- fast.cov(.xx, .yy)
    beta.z <- fast.z.cov(.xx, .yy)

    if(.melt) {
        require(reshape2)
        require(dplyr)

    beta <- beta %>% as.matrix() %>% melt()
    beta.z <- beta.z %>% as.matrix() %>% melt()

    colnames(beta) <- c('snp', 'gene', 'beta')
    colnames(beta.z) <- c('snp', 'gene', 'beta.z')
    out <- beta %>% left_join(beta.z, by = c('snp', 'gene')) %>%
        as.data.frame()
    return(out)
    }
    return(list(beta = beta, beta.se = beta / (beta.z +  TOL)))
}

fast.cor <- function(x, y) {
    x.sd <- apply(x, 2, sd, na.rm = TRUE)
    y.sd <- apply(y, 2, sd, na.rm = TRUE)
    ret <- fast.cov(scale(x, scale = FALSE), scale(y, scale = FALSE))
    ret <- sweep(sweep(ret, 1, x.sd, `/`), 2, y.sd, `/`)    
    return(ret)
}

write.tab <- function(x, ...) {
    write.table(x, sep = '\t', quote = FALSE,
                row.names = FALSE, col.names = FALSE, ...)
}

write.tab.named <- function(x, ...) {
    write.table(x, sep = '\t', quote = FALSE,
                row.names = FALSE, col.names = TRUE, ...)
}

write.mat <- function(x, digits = 4, ...) {
    write.tab(round(x, digits), ...)
}

take.theta <- function(qtl.out) {
    if('mean.left' %in% names(qtl.out)){
        theta <- qtl.out$mean.left$theta %*% t(qtl.out$mean.right$theta)
    } else {
        theta <- qtl.out$mean$theta
    }
    return(theta)
}

take.pve <- function(qtl.out, xx, cc) {

    theta <- take.theta(qtl.out)
    theta.cov <- qtl.out$mean.cov$theta
    resid <- qtl.out$resid$theta

    xx[is.na(xx)] <- 0
    theta[is.na(theta)] <- 0
    theta.cov[is.na(theta.cov)] <- 0
    cc[is.na(cc)] <- 0
    eta.1 <- xx %*% theta
    eta.2 <- cc %*% theta.cov

    data.frame(v1 = apply(eta.1, 2, var, na.rm = TRUE),
               v2 = apply(eta.2, 2, var, na.rm = TRUE),
               vr = apply(resid, 2, var, na.rm = TRUE))
}
