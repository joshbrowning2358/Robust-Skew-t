library(sn)
library(ggplot2)

setwd("~/Professional Files/Mines/Research/Robust Estimators/")

sn2 <- makeSECdistr(dp=list(xi=c(2,1), Omega=diag(2), alpha=c(3,-5))
        ,family="SN", name="Mandy", compNames=c("m1","m2"))

show(sn2)
summary(sn2)
plot(sn2)

#------------
# http://tesi.cab.unipd.it/7115/

st2 <- makeSECdistr(dp=list(xi=c(2,1), Omega=diag(2), alpha=c(3,-5),df=6)
        ,family="ST", name="Skew-T", compNames=c("m1","m2"))

plot(st2)


makeSECdistr <- function (dp, family, name, compNames) 
{
  #Data quality checks:
    if (!(toupper(family) %in% c("SN", "ESN", "SC", "ST"))) 
        stop("unknown family")
    family <- toupper(family)
    ndp <- if (family %in% c("SN", "SC")) 
        3
    else 4
    if (length(dp) != ndp) 
        stop("wrong number of dp components")
    if (family == "ST") {
        nu <- as.numeric(dp[4])
        if (nu <= 0) 
            stop("'nu' for ST family must be positive")
        if (nu == Inf) {
            warning("ST family with 'nu==Inf' is changed to SN family")
            family <- "SN"
            dp <- dp[-4]
        }
    }

    if (is.numeric(dp)) {
        if (dp[2] <= 0) 
            stop("omega parameter must be positive")
        fourth <- switch(family, SN = NULL, ESN = "tau", SC = NULL, 
            ST = "nu")
        names(dp) <- c("xi", "omega", "alpha", fourth)
        name <- if (!missing(name)) 
            as.character(name)[1]
        else paste("Unnamed-", toupper(family), sep = "")
        obj <- new("SECdistrUv", dp = dp, family = family, name = name)
    }
    else {
        if (is.list(dp)) {
            names(dp) <- rep(NULL, ndp)
            d <- length(dp[[3]])
            if (any(abs(dp[[3]]) == Inf)) 
                stop("Inf in alpha not allowed")
            if (length(dp[[1]]) != d) 
                stop("mismatch of parameters size")
            Omega <- matrix(dp[[2]], d, d)
            if (any(Omega != t(Omega))) 
                stop("Omega matrix must be symmetric")
            if (min(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values) <= 
                0) 
                stop("Omega matrix must be positive definite")
            dp0 <- list(xi = as.vector(dp[[1]]), Omega = Omega, 
                alpha = dp[[3]])
            name <- if (!missing(name)) 
                as.character(name)[1]
            else paste("Unnamed-", toupper(family), "[d=", as.character(d), 
                "]", sep = "")
            if (family == "ST") 
                dp0$nu <- nu
            if (family == "ESN") 
                dp0$tau <- dp[[4]]
            if (d == 1) 
                warning(paste("A multivariate distribution with dimension=1 is a near-oxymoron.", 
                  "\nConsider using a 'dp' vector to define a univariate distribution.", 
                  "\nHowever, I still build a multivariate distribution for you."))
            if (missing(compNames)) {
                compNames <- if (length(colnames(dp[[1]])) == 
                  d) 
                  colnames(dp[[1]])
                else as.vector(outer("V", as.character(1:d), 
                  paste, sep = ""))
            }
            else {
                if (length(compNames) != d) 
                  stop("Wrong length of 'compNames'")
                compNames <- as.character(as.vector(compNames))
            }
            obj <- new("SECdistrMv", dp = dp0, family = family, 
                name = name, compNames = compNames)
        }
        else stop("'dp' must be either a numeric vector or a list")
    }
    obj
}

edit( sn:::plot.SECdistrMv )
st2 <- makeSECdistr(dp=list(xi=c(2,1), Omega=diag(2), alpha=c(3,-5),df=6)
        ,family="ST", name="MandyST", compNames=c("m1","m2"))
x = st2
landmarks="auto"
data=NULL; data.par=NULL; gap=0.5

plot.SECdistrMv <- function (x, range, probs, npt, landmarks = "auto", main, comp, 
    compLabs, data = NULL, data.par = NULL, gap = 0.5, ...) 
{
    obj <- x
    if (slot(obj, "class") != "SECdistrMv") 
        stop("object of wrong class")
    compNames <- slot(obj, "compNames")
    d <- length(compNames)
#    if (missing(probs)) 
        probs <- c(0.25, 0.5, 0.75, 0.95)
    if (any(probs <= 0 || probs >= 1)) 
        stop("probs must be within (0,1)")
    if (sum(probs > 0 && probs < 1) == 0) 
        stop("invalid probs")
#    if (missing(npt)) 
        npt <- rep(101, d)
#    if (missing(main)) {
        main <- if (d == 2) 
            paste("Density function of", slot(obj, "name"))
        else paste("Bivariate densities of", slot(obj, "name"))
#    }
#    if (missing(comp)) 
        comp <- seq(1, d)
#    if (missing(compLabs)) 
        compLabs <- compNames
    if (length(compLabs) != d) 
        stop("wrong length of 'compLabs' or 'comp' vector")
    family <- toupper(obj@family)
    lc.family <- tolower(family)
    if (lc.family == "esn") 
        lc.family <- "sn"
    dp <- slot(obj, "dp")
#    if (missing(range)) {
        range <- matrix(NA, 2, d)
        q.fn <- get(paste("q", lc.family, sep = ""), inherits = TRUE)
        for (j in 1:d) {
            marg <- marginalSECdistr(obj, comp = j, drop = TRUE)
            q <- q.fn(c(0.05, 0.25, 0.75, 0.95), dp = marg@dp)
            dq <- diff(q)
            range[, j] <- c(q[1] - 1.5 * dq[1], q[length(q)] + 
                1.5 * dq[length(dq)])
            if (!is.null(data)) {
                range[1, j] <- min(range[1, j], min(data[, j]))
                range[2, j] <- max(range[2, j], max(data[, j]))
            }
        }
#    }
    dots <- list() #dots <- list(...)
    nmdots <- names(dots)
    if (d == 1) {
        message("Since dimension=1, plot as a univariate distribution")
        objUv <- marginalSECdistr(obj, comp = 1, drop = TRUE)
        out <- plot(objUv, data = data, ...)
    }
    if (d == 2) 
        out <- plot.SECdistrBv(x, range, probs, npt, compNames, 
            compLabs, landmarks, data, data.par, main)
    if (d > 2) {
        textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
            y, txt, cex = cex, font = font)
        localAxis <- function(side, x, y, xpd, bg, main, oma, 
            ...) {
            if (side%%2 == 1) 
                Axis(x, side = side, xpd = NA, ...)
            else Axis(y, side = side, xpd = NA, ...)
        }
        localPlot <- function(..., oma, font.main, cex.main) plot.SECdistrBv(...)
        text.diag.panel <- compLabs
        oma <- if ("oma" %in% nmdots) 
            dots$oma
        else NULL
        if (is.null(oma)) {
            oma <- c(4, 4, 4, 4)
            if (!is.null(main)) 
                oma[3L] <- 6
        }
        opar <- par(mfrow = c(length(comp), length(comp)), mar = rep(c(gap, 
            gap/2), each = 2), oma = oma)
        on.exit(par(opar))
        out <- list()
        count <- 0
        for (i in comp) for (j in comp) {
            count <- count + 1
            if (i == j) {
                plot(1, type = "n", xlab = "", ylab = "", axes = FALSE)
                text(1, 1, text.diag.panel[i], cex = 2)
                box()
                out[[count]] <- paste("diagonal component", 
                  compNames[i])
            }
            else {
                ji <- c(j, i)
                marg <- marginalSECdistr(obj, comp = ji)
                out[[count]] <- localPlot(x = marg, range = range[, 
                  ji], probs = probs, npt = npt[ji], compNames = compNames[ji], 
                  compLabs = compLabs[ji], landmarks = landmarks, 
                  data = data[, ji], data.par = data.par, main = "", 
                  yaxt = "n", xaxt = "n", ...)
                if (i == comp[1]) 
                  axis(3)
                if (j == length(comp)) 
                  axis(4)
                if (j == comp[1]) 
                  axis(2)
                if (i == length(comp)) 
                  axis(1)
                box()
            }
        }
        par(new = FALSE)
        if (!is.null(main)) {
            font.main <- if ("font.main" %in% nmdots) 
                dots$font.main
            else par("font.main")
            cex.main <- if ("cex.main" %in% nmdots) 
                dots$cex.main
            else par("cex.main")
            mtext(main, side = 3, TRUE, line = 5, outer = TRUE, 
                at = NA, cex = cex.main, font = font.main, adj = 0.5)
        }
    }
    invisible(out)
}


edit( sn:::plot.SECdistrBv )
x <- makeSECdistr(dp=list(xi=c(2,1), Omega=diag(2), alpha=c(3,-5),df=6)
        ,family="ST", name="MandyST", compNames=c("m1","m2"))
npt = rep(101,2); compNames=c("a","b"); compLabs=c("a","b")
landmarks="auto"
range = matrix(c(-1,3.3,6.2,2.6), nrow=2, byrow=T)
probs = c(1:3/4,.95)
data = NULL
data.par= NULL
main = "Title"
plot.SECdistrBv <- function (x, range, probs, npt = rep(101, 2), compNames, compLabs, 
    landmarks, data = NULL, data.par, main, ...) 
{
    obj <- x
    dp <- slot(obj, "dp")
    family <- slot(obj, "family")
    lc.family <- tolower(family)
    if (lc.family == "esn") 
        lc.family <- "sn"
    d.fn <- get(paste("dm", lc.family, sep = ""), inherits = TRUE)
    n1 <- npt[1]
    n2 <- npt[2]
    x1 <- seq(min(range[, 1]), max(range[, 1]), length = n1)
    x2 <- seq(min(range[, 2]), max(range[, 2]), length = n2)
    x1.x2 <- cbind(rep(x1, n2), as.vector(matrix(x2, n1, n2, 
        byrow = TRUE)))
    X <- matrix(x1.x2, n1 * n2, 2, byrow = FALSE)
    pdf <- matrix(d.fn(X, dp = dp), n1, n2)
    Omega <- dp[[2]]
    Omega.bar <- cov2cor(Omega)
    alpha <- dp[[3]]
    alpha.star <- sqrt(sum(alpha * as.vector(Omega.bar %*% alpha)))
    omega <- sqrt(diag(Omega))
    if (lc.family == "sn") {
        k.tau <- if (length(dp) == 4) 
            (zeta(2, dp[[4]]) * pi)^2/4
        else 1
        log.levels <- (log(1 - probs) - log(2 * pi) - 0.5 * 
            log(1 - Omega.bar[1, 2]^2) + k.tau * log(1 + exp(-1.544/alpha.star))) - 
            sum(log(omega))
    }
    if (lc.family == "st" | lc.family == "sc") {
        nu <- if (lc.family == "st") 
            obj@dp[[4]]
        else 1
        l.nu <- (-1.3/nu - 4.93)
        h <- 100 * log(exp(((1.005 * alpha.star - 0.045) * l.nu - 
            1.5)/alpha.star) + 1)
        K <- h * (1.005 * alpha.star - 0.1) * (1 + nu)/(alpha.star * 
            nu)
        qF <- qf(probs, 2, nu)
        log.levels <- (lgamma(nu/2 + 1) - lgamma(nu/2) - log(pi * 
            nu) - 0.5 * log(1 - Omega.bar[1, 2]^2) - (nu/2 + 
            1) * log(2 * qF/nu + 1) + K - sum(log(omega)))
    }
    oo <- options()
    options(warn = -1)
    contour(x1, x2, pdf, levels = exp(log.levels), labels = paste("p=", 
        as.character(probs), sep = ""), main = main, xlab = compLabs[1], 
        ylab = compLabs[2])
    if (!is.null(data)) {
        col <- if (!is.null(data.par$col)) 
            data.par$col
        else par()$col
        pch <- if (!is.null(data.par$pch)) 
            data.par$pch
        else par()$pch
        cex <- if (!is.null(data.par$cex)) 
            data.par$cex
        else par()$cex
        points(data, col = col, pch = pch, cex = cex)
    }
    if (landmarks != "") {
        if (landmarks == "auto") {
            mean.type <- "proper"
            if (lc.family == "sc") 
                mean.type <- "pseudo"
            if (lc.family == "st") {
                if (dp[[4]] <= 1) 
                  mean.type <- "pseudo"
            }
        }
        else mean.type <- landmarks
        landmarks.label <- c("origin", "mode", if (mean.type == 
            "proper") "mean" else "mean~")
        cp <- dp2cpMv(dp, family, cp.type = mean.type, upto = 1)
        mode <- modeSECdistrMv(dp, family)
        x.pts <- c(dp$xi[1], mode[1], cp[[1]][1])
        y.pts <- c(dp$xi[2], mode[2], cp[[1]][2])
        points(x.pts, y.pts, ...)
        text(x.pts, y.pts, landmarks.label, pos = 2, offset = 0.3, 
            ...)
        lines(x.pts, y.pts, lty = 2)
    }
    options(oo)
    return(list(x = x1, y = x2, names = compNames, density = pdf))
}

getPval = function(data, dp, probs=100:0/100){
  probs = sort(probs, decreasing=T)

  nu <- dp$nu
  Omega <- dp[[2]]
  Omega.bar <- cov2cor(Omega)
  alpha <- dp[[3]]
  alpha.star <- sqrt(sum(alpha * as.vector(Omega.bar %*% alpha)))
  omega <- sqrt(diag(Omega))
  l.nu <- (-1.3/nu - 4.93)
  h <- 100 * log(exp(((1.005 * alpha.star - 0.045) * l.nu - 
      1.5)/alpha.star) + 1)
  K <- h * (1.005 * alpha.star - 0.1) * (1 + nu)/(alpha.star * 
      nu)
  qF <- qf(probs, 2, nu)
  log.levels <- (lgamma(nu/2 + 1) - lgamma(nu/2) - log(pi * 
      nu) - 0.5 * log(1 - Omega.bar[1, 2]^2) - (nu/2 + 
      1) * log(2 * qF/nu + 1) + K - sum(log(omega)))

  dens = dmst(data, dp=dp)
  pval = rep(0, NROW(data))
  for(i in 1:length(probs)){
    filt = dens>exp(log.levels[i])
    pval[filt] = 1-probs[i]
  }
  return(pval)
}

getContour = function(dp, probs=c(0.99, 0.98, 0.95, 0.9, 0.8)){
  nu <- dp$nu
  Omega <- dp[[2]]
  Omega.bar <- cov2cor(Omega)
  alpha <- dp[[3]]
  alpha.star <- sqrt(sum(alpha * as.vector(Omega.bar %*% alpha)))
  omega <- sqrt(diag(Omega))
  l.nu <- (-1.3/nu - 4.93)
  h <- 100 * log(exp(((1.005 * alpha.star - 0.045) * l.nu - 
      1.5)/alpha.star) + 1)
  K <- h * (1.005 * alpha.star - 0.1) * (1 + nu)/(alpha.star * 
      nu)
  qF <- qf(probs, 2, nu)
  log.levels <- (lgamma(nu/2 + 1) - lgamma(nu/2) - log(pi * 
      nu) - 0.5 * log(1 - Omega.bar[1, 2]^2) - (nu/2 + 
      1) * log(2 * qF/nu + 1) + K - sum(log(omega)))
  return(exp(log.levels))
}

out = rmst( 10000, dp=st2@dp )
pvalues = getPval( out, dp=st2@dp, probs=0:100/100 )
qplot( pvalues, binwidth=.01 )
sum( pvalues==1 )

for(p in c(.5, .75, .9, .95, .99)){
  png(paste0("Results/Skew_CI/alpha=",p,".png"), width=4000, height=4000)
  plot(st2, probs=p, landmarks="")
  points(out[,1], out[,2], col=as.numeric(pvalues<1-p)+1, pch=16, cex=.1)
  dev.off()
}

levels = c(.5, .75, .9, .95, .99)
png("Results/Skew_CI/alpha=multiple_zoomed.png", width=4000, height=4000)
plot(st2, probs=levels, landmarks="")
points(out[,1], out[,2], col=sapply(pvalues, function(x) sum(x<1-levels))+1, pch="+", cex=.6)
dev.off()

levels = c(.5, .75, .9, .95, .99)
png("Results/Skew_CI/alpha=multiple.png", width=400, height=400)
plot(st2, probs=levels, landmarks="")
points(out[,1], out[,2], col=sapply(pvalues, function(x) sum(x<1-levels))+1, pch="+", cex=.6)
dev.off()