gen.lms <-
function (dat, p = 0.5, bet, epsilon = 0, k.u = 0) 
{
    n <- nrow(dat)
    if (k.u == 0) {
        k.u <- kl.ku(n, p, bet, epsilon)[2]
    }
    else {
        if (k.u < 0) {
            stop("k.u must be a positive integer!\n")
        }
    }
    dat.b <- dat[rowSums(abs(dat)) < Inf, ]
    n.b <- nrow(dat.b)
    dat.b0 <- dat[rowSums(abs(dat[, c(3, 4)])) < Inf, ]
    n.b0 <- nrow(dat.b0)
    if (n.b0 < k.u) {
        stop("Too many unbounded data!\n")
    }
    else {
        if (n.b < k.u) {
            b.pot <- 0
        }
        else {
            dat.unique <- unique(dat.b, MARGIN = 1)
            m <- nrow(dat.unique)
            ind <- 1:m
            ind.vi <- 0
            for (i in 1:(m - 1)) {
                ind.vi <- c(ind.vi, rep(i, (m - i)))
            }
            ind.vj <- 0
            for (j in 2:m) {
                ind.vj <- c(ind.vj, ind[j:m])
            }
            ind.dat <- cbind(ind.vi, ind.vj)[-1, ]
            b.pot <- matrix(0, choose(m, 2) + 1, 4)
            for (i in 1:choose(m, 2)) {
                b.pot[i, 1] <- (dat.unique[ind.dat[i, 1], 3] - 
                  dat.unique[ind.dat[i, 2], 3])/(dat.unique[ind.dat[i, 
                  1], 1] - dat.unique[ind.dat[i, 2], 1])
                b.pot[i, 2] <- (dat.unique[ind.dat[i, 1], 3] - 
                  dat.unique[ind.dat[i, 2], 3])/(dat.unique[ind.dat[i, 
                  1], 2] - dat.unique[ind.dat[i, 2], 2])
                b.pot[i, 3] <- (dat.unique[ind.dat[i, 1], 4] - 
                  dat.unique[ind.dat[i, 2], 4])/(dat.unique[ind.dat[i, 
                  1], 1] - dat.unique[ind.dat[i, 2], 1])
                b.pot[i, 4] <- (dat.unique[ind.dat[i, 1], 4] - 
                  dat.unique[ind.dat[i, 2], 4])/(dat.unique[ind.dat[i, 
                  1], 2] - dat.unique[ind.dat[i, 2], 2])
            }
            b.pot <- unique(as.vector(round(b.pot, 10)))
            b.pot <- na.omit(b.pot)
            if (max(b.pot) == Inf) {
                b.pot <- b.pot[-which(b.pot == Inf)]
            }
            if (min(b.pot) == -Inf) {
                b.pot <- b.pot[-which(b.pot == -Inf)]
            }
            max.b <- max(abs(b.pot))
        }
    }
    a.pot <- rep(0, length(b.pot))
    q.pot <- rep(0, length(b.pot))
    for (i in 1:length(b.pot)) {
        if (b.pot[i] == 0) {
            d <- dat.b0[order(dat.b0[, 4]), c(3, 4)]
            d.u <- d[, 2]
            d.l <- d[, 1]
            d.l <- sort(d.l)
            diff <- matrix(0, (n.b0 - k.u + 1), 3)
            for (k in 1:nrow(diff)) {
                d.u.k <- d[d[, 1] >= d.l[k], 2]
                diff[k, 1] <- d.u.k[k.u] - d.l[k]
                diff[k, 2] <- (d.u.k[k.u] + d.l[k])/2
                diff[k, 3] <- (d.u.k[k.u] - d.l[k])/2
            }
            diff <- unique(diff, MARGIN = 1)
            k.opt <- which(diff[, 1] == min(diff[, 1]))
            if (length(k.opt) > 1) {
                for (l in 2:length(k.opt)) {
                  b.pot <- c(b.pot, b.pot[i])
                  a.pot <- c(a.pot, diff[k.opt[l], 2])
                  q.pot <- c(q.pot, diff[k.opt[l], 3])
                }
            }
            a.pot[i] <- diff[k.opt[1], 2]
            q.pot[i] <- diff[k.opt[1], 3]
        }
        else {
            d <- matrix(0, n.b, 2)
            if (b.pot[i] > 0) {
                d[, 1] <- dat.b[, 3] - b.pot[i] * dat.b[, 2]
                d[, 2] <- dat.b[, 4] - b.pot[i] * dat.b[, 1]
            }
            else {
                d[, 1] <- dat.b[, 3] - b.pot[i] * dat.b[, 1]
                d[, 2] <- dat.b[, 4] - b.pot[i] * dat.b[, 2]
            }
            d <- d[order(d[, 2]), ]
            d.u <- d[, 2]
            d.l <- d[, 1]
            d.l <- sort(d.l)
            diff <- matrix(0, (n.b - k.u + 1), 3)
            for (k in 1:nrow(diff)) {
                d.u.k <- d[d[, 1] >= d.l[k], 2]
                diff[k, 1] <- d.u.k[k.u] - d.l[k]
                diff[k, 2] <- (d.u.k[k.u] + d.l[k])/2
                diff[k, 3] <- (d.u.k[k.u] - d.l[k])/2
            }
            diff <- unique(diff, MARGIN = 1)
            k.opt <- which(diff[, 1] == min(diff[, 1]))
            if (length(k.opt) > 1) {
                for (l in 2:length(k.opt)) {
                  b.pot <- c(b.pot, b.pot[i])
                  a.pot <- c(a.pot, diff[k.opt[l], 2])
                  q.pot <- c(q.pot, diff[k.opt[l], 3])
                }
            }
            a.pot[i] <- diff[k.opt[1], 2]
            q.pot[i] <- diff[k.opt[1], 3]
        }
    }
    preresult <- cbind(a.pot, b.pot, q.pot)
    preresult <- unique(preresult)
    attr(preresult, "dimnames")[[2]] <- c("a.lrm", "b.lrm", "q.lrm")
    list(lrm = preresult[preresult[, 3] == min(preresult[, 3]), 
        ], max.b)
}
idf.create <-
function (dat, var.labels = NULL) 
{
    if (!(is.data.frame(dat))) 
        stop("The data must be provided as a data frame! \n")
    m <- ncol(dat)/2
    if (ceiling(m) != m) 
        stop("The data frame must contain two neighboring columns (with lower and upper endpoints, respectively) for each variable! \n")
    idf <- vector("list", m)
    var.lab <- NULL
    for (i in 1:m) {
        idf[[i]] <- dat[, c((2 * i - 1), 2 * i)]
    }
    if (!is.null(var.labels)) {
        for (i in 1:m) {
            attr(idf[[i]], "names") <- c(paste(var.labels[i], 
                "l", sep = "."), paste(var.labels[i], "u", sep = "."))
        }
        attr(idf, "names") <- var.labels
    }
    else {
        for (i in 1:m) {
            attr(idf[[i]], "names") <- c(paste(paste("var", i, 
                sep = ""), "l", sep = "."), paste(paste("var", 
                i, sep = ""), "u", sep = "."))
            var.lab <- c(var.lab, paste("var", i, sep = ""))
        }
        attr(idf, "names") <- var.lab
    }
    idf$n <- nrow(dat)
    class(idf) <- "idf"
    idf
}
kl.ku <-
function (n, p = 0.5, bet, epsilon = 0) 
{
    if ((max(p, 1 - p) + epsilon)^n > bet) {
        stop("k.l and k.u are not defined ! \n")
    }
    lambda <- function(s, t) {
        (s/t)^(-s) * ((1 - s)/(1 - t))^(s - 1)
    }
    k.l <- Find(function(k) {
        lambda(k/n, p - epsilon) <= bet^(1/n)
    }, 0:floor((p - epsilon) * n), right = T)
    k.u <- Find(function(k) {
        lambda(k/n, p + epsilon) <= bet^(1/n)
    }, (floor((p + epsilon) * n) + 1):n)
    c(k.l, k.u)
}
plot.idf <-
function (x, y = NULL, ..., var = NULL, k.x = 1, k.y = 1, p.cex = 1, 
    x.adj = 0.5, x.padj = 3, y.las = 0, y.adj = 1, y.padj = 0, 
    x.lim = c(0, 0), y.lim = c(0, 0), x.lab = "X", y.lab = "Y") 
{
    dat.idf <- x
    if (!is.null(var)) {
        dat <- cbind(dat.idf[[which(names(dat.idf) == var[1])]], 
            dat.idf[[which(names(dat.idf) == var[2])]])
    }
    else {
        dat <- cbind(dat.idf[[1]], dat.idf[[2]])
    }
    if (x.lim[1] == 0 & x.lim[2] == 0) {
        x.min <- floor(min(dat[dat[, 1] != -Inf, 1]))
        x.max <- ceiling(max(dat[dat[, 2] != Inf, 2]))
    }
    else {
        x.min <- x.lim[1]
        x.max <- x.lim[2]
    }
    if (y.lim[1] == 0 & y.lim[2] == 0) {
        y.min <- floor(min(dat[dat[, 3] != -Inf, 3]))
        y.max <- ceiling(max(dat[dat[, 4] != Inf, 4]))
    }
    else {
        y.min <- y.lim[1]
        y.max <- y.lim[2]
    }
    dat[dat[, 1] == -Inf, 1] <- min(dat[dat[, 1] != -Inf, 1]) - 
        (x.max - x.min)/10
    dat[dat[, 2] == Inf, 2] <- max(dat[dat[, 2] != Inf, 2]) + 
        (x.max - x.min)/10
    dat[dat[, 3] == -Inf, 3] <- min(dat[dat[, 3] != -Inf, 3]) - 
        (y.max - y.min)/10
    dat[dat[, 4] == Inf, 4] <- max(dat[dat[, 4] != Inf, 4]) + 
        (y.max - y.min)/10
    dat[, 1] <- floor(dat[, 1] * k.x)/k.x
    dat[, 2] <- ceiling(dat[, 2] * k.x)/k.x
    dat[, 3] <- floor(dat[, 3] * k.y)/k.y
    dat[, 4] <- ceiling(dat[, 4] * k.y)/k.y
    n <- dat.idf$n
    dat <- dat[order(dat[, 1]), ]
    x.steps <- rep(0, n)
    y.steps <- rep(0, n)
    for (i in 1:nrow(dat)) {
        x.steps[i] <- round((dat[i, 2] - dat[i, 1]) * k.x, 10)
        y.steps[i] <- round((dat[i, 4] - dat[i, 3]) * k.y, 10)
    }
    plot(mean(dat[, 1]), mean(dat[, 3]), type = "p", pch = 15, 
        col = 0, xlim = c(x.min, x.max), ylim = c(y.min, y.max), 
        xlab = " ", ylab = " ", las = y.las)
    mtext(x.lab, side = 1, las = 1, adj = x.adj, padj = x.padj)
    mtext(y.lab, side = 2, las = y.las, adj = y.adj, padj = y.padj)
    for (i in 1:nrow(dat)) {
        if (x.steps[i] > 0 & y.steps[i] > 0) {
            Z.i <- matrix(0, nrow = (x.steps[i] * y.steps[i]), 
                2)
            if (x.steps[i] == 1) {
                Z.i[, 1] <- rep(dat[i, 1] + 1/(2 * k.x), times = y.steps[i])
            }
            else {
                Z.i[, 1] <- rep(seq(dat[i, 1] + 1/(2 * k.x), 
                  dat[i, 2] - 1/(2 * k.x), 1/k.x), times = y.steps[i])
            }
            if (y.steps[i] == 1) {
                Z.i[, 2] <- rep(dat[i, 3] + 1/(2 * k.y), each = x.steps[i])
            }
            else {
                Z.i[, 2] <- rep(seq(dat[i, 3] + 1/(2 * k.y), 
                  dat[i, 4] - 1/(2 * k.y), 1/k.y), each = x.steps[i])
            }
            points(Z.i[, 1], Z.i[, 2], pch = 15, col = "lightgray", 
                cex = p.cex)
        }
        points(dat[i, 2], dat[i, 4], pch = 20, cex = 0.25)
        points(dat[i, 2], dat[i, 3], pch = 20, cex = 0.25)
        points(dat[i, 1], dat[i, 4], pch = 20, cex = 0.25)
        points(dat[i, 1], dat[i, 3], pch = 20, cex = 0.25)
        segments(dat[i, 1], dat[i, 4], dat[i, 2], dat[i, 4], 
            lty = 1, lwd = 2)
        segments(dat[i, 1], dat[i, 3], dat[i, 2], dat[i, 3], 
            lty = 1, lwd = 2)
        segments(dat[i, 1], dat[i, 3], dat[i, 1], dat[i, 4], 
            lty = 1, lwd = 2)
        segments(dat[i, 2], dat[i, 3], dat[i, 2], dat[i, 4], 
            lty = 1, lwd = 2)
    }
}
plot.s.linlir <-
function (x, y = NULL, ..., typ, para.typ = "polygon", b.range = c(-1e-05, 
    1e-05), b.grid = 1000, nb.func = 1000, seed.func = NULL, 
    pl.lrm = TRUE, pl.band = FALSE, pl.dat = FALSE, k.x = 1, 
    k.y = 1, p.cex = 1, x.adj = 0.5, x.padj = 3, y.las = 0, y.adj = 1, 
    y.padj = 0, x.lim = c(0, 0), y.lim = c(0, 0), x.lab = " ", 
    y.lab = " ") 
{
    x.s.linlir <- x
    if (typ == "para") {
        if (x.lim[1] == 0 & x.lim[2] == 0) {
            x.min <- floor(x.s.linlir$b.undom[1])
            x.max <- ceiling(x.s.linlir$b.undom[2])
        }
        else {
            x.min <- x.lim[1]
            x.max <- x.lim[2]
        }
        if (y.lim[1] == 0 & y.lim[2] == 0) {
            y.min <- floor(x.s.linlir$a.undom[1])
            y.max <- ceiling(x.s.linlir$a.undom[2])
        }
        else {
            y.min <- y.lim[1]
            y.max <- y.lim[2]
        }
        if (para.typ == "polygon") {
            if (b.range[1] == -1e-05 & b.range[2] == 1e-05) {
                b.range <- c(x.min, x.max)
            }
            b.pot <- seq(b.range[1], b.range[2], by = (b.range[2] - 
                b.range[1])/b.grid)
            n <- x.s.linlir$n
            k.l <- x.s.linlir$config$k.l
            a.l.plot <- matrix(NA, nrow = (n - k.l), ncol = length(b.pot))
            a.u.plot <- matrix(NA, nrow = (n - k.l), ncol = length(b.pot))
            for (i in 1:length(b.pot)) {
                a.int <- undom.a(x.s.linlir$dat, b = b.pot[i], 
                  x.s.linlir$q.lrm, x.s.linlir$config$p, x.s.linlir$config$bet, 
                  x.s.linlir$config$epsilon)
                a.l.plot[, i] <- a.int[[1]][, 1]
                a.u.plot[, i] <- a.int[[1]][, 2]
            }
            if (x.lab == " ") {
                x.lab <- "b"
            }
            if (y.lab == " ") {
                y.lab <- "a"
            }
            plot(x.min, y.max, type = "n", xlim = c(x.min, x.max), 
                ylim = c(y.min, y.max), xlab = " ", ylab = " ", 
                las = y.las)
            mtext(x.lab, side = 1, las = 1, adj = x.adj, padj = x.padj)
            mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                padj = y.padj)
            for (j in 1:(n - k.l)) {
                polygon(c(b.pot[a.l.plot[j, ] <= a.u.plot[j, 
                  ]], rev(b.pot[a.l.plot[j, ] <= a.u.plot[j, 
                  ]])), c(a.l.plot[j, a.l.plot[j, ] <= a.u.plot[j, 
                  ]], rev(a.u.plot[j, a.l.plot[j, ] <= a.u.plot[j, 
                  ]])), col = "darkgrey", lty = 0)
            }
        }
        else {
            undom.para <- x.s.linlir$undom.para
            plot(undom.para[, 2], undom.para[, 1], pch = 19, 
                col = "darkgrey", xlim = c(x.min, x.max), ylim = c(y.min, 
                  y.max), xlab = " ", ylab = " ", las = y.las)
            mtext(x.lab, side = 1, las = 1, adj = x.adj, padj = x.padj)
            mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                padj = y.padj)
        }
        if (pl.lrm == TRUE) {
            points(x.s.linlir$f.lrm[2], x.s.linlir$f.lrm[1], 
                pch = 19, col = 4, cex = 1)
        }
    }
    else {
        dat <- x.s.linlir$dat
        if (x.lim[1] == 0 & x.lim[2] == 0) {
            x.min <- floor(min(dat[dat[, 1] != -Inf, 1]))
            x.max <- ceiling(max(dat[dat[, 2] != Inf, 2]))
        }
        else {
            x.min <- x.lim[1]
            x.max <- x.lim[2]
        }
        if (y.lim[1] == 0 & y.lim[2] == 0) {
            y.min <- floor(min(dat[dat[, 3] != -Inf, 3]))
            y.max <- ceiling(max(dat[dat[, 4] != Inf, 4]))
        }
        else {
            y.min <- y.lim[1]
            y.max <- y.lim[2]
        }
        x.d <- (x.max - x.min)/10
        if (x.lab == " ") {
            x.lab <- "X"
        }
        if (y.lab == " ") {
            y.lab <- "Y"
        }
        if (typ == "lrm") {
            if (pl.dat == TRUE) {
                dat.idf <- idf.create(dat)
                plot(dat.idf, k.x = k.x, k.y = k.y, p.cex = p.cex, 
                  x.adj = x.adj, x.padj = x.padj, y.las = y.las, 
                  y.adj = y.adj, y.padj = y.padj, x.lim = c(x.min, 
                    x.max), y.lim = c(y.min, y.max), x.lab = x.lab, 
                  y.lab = y.lab)
                curve(x.s.linlir$f.lrm[1] + x.s.linlir$f.lrm[2] * 
                  x, x.min - x.d, x.max + x.d, add = T, lty = 1, 
                  col = 4, lwd = 2)
            }
            else {
                plot(min(dat[, 1]), max(dat[, 4]), type = "n", 
                  xlim = c(x.min, x.max), ylim = c(y.min, y.max), 
                  xlab = " ", ylab = " ", las = y.las)
                mtext(x.lab, side = 1, las = 1, adj = x.adj, 
                  padj = x.padj)
                mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                  padj = y.padj)
                curve(x.s.linlir$f.lrm[1] + x.s.linlir$f.lrm[2] * 
                  x, x.min - x.d, x.max + x.d, add = T, lty = 1, 
                  col = 4, lwd = 2)
            }
        }
        else {
            if (!is.null(seed.func)) {
                set.seed(seed.func)
            }
            undom <- x.s.linlir$undom.para[sample(1:nrow(x.s.linlir$undom.para), 
                nb.func), ]
            if (pl.dat == TRUE) {
                dat.idf <- idf.create(dat)
                plot(dat.idf, k.x = k.x, k.y = k.y, p.cex = p.cex, 
                  x.adj = x.adj, x.padj = x.padj, y.las = y.las, 
                  y.adj = y.adj, y.padj = y.padj, x.lim = c(x.min, 
                    x.max), y.lim = c(y.min, y.max), x.lab = x.lab, 
                  y.lab = y.lab)
                for (i in 1:nrow(undom)) {
                  curve(undom[i, 1] + undom[i, 2] * x, x.min - 
                    x.d, x.max + x.d, add = T, lty = 1, col = "darkgrey")
                }
            }
            else {
                plot(min(dat[, 1]), max(dat[, 4]), type = "n", 
                  xlim = c(x.min, x.max), ylim = c(y.min, y.max), 
                  xlab = " ", ylab = " ", las = y.las)
                mtext(x.lab, side = 1, las = 1, adj = x.adj, 
                  padj = x.padj)
                mtext(y.lab, side = 2, las = y.las, adj = y.adj, 
                  padj = y.padj)
                for (i in 1:nrow(undom)) {
                  curve(undom[i, 1] + undom[i, 2] * x, x.min - 
                    x.d, x.max + x.d, add = T, lty = 1, col = "darkgrey")
                }
            }
            if (pl.lrm == TRUE) {
                curve(x.s.linlir$f.lrm[1] + x.s.linlir$f.lrm[2] * 
                  x, x.min - x.d, x.max + x.d, add = T, lty = 1, 
                  col = 4, lwd = 2)
            }
        }
        if (pl.band == TRUE) {
            curve(x.s.linlir$f.lrm[1] + x.s.linlir$q.lrm + x.s.linlir$f.lrm[2] * 
                x, x.min - x.d, x.max + x.d, add = T, lty = 2, 
                col = 4, lwd = 2)
            curve(x.s.linlir$f.lrm[1] - x.s.linlir$q.lrm + x.s.linlir$f.lrm[2] * 
                x, x.min - x.d, x.max + x.d, add = T, lty = 2, 
                col = 4, lwd = 2)
        }
    }
}
print.s.linlir <-
function (x, ...) 
{
    x.s.linlir <- x
    cat("\nSimple linear LIR analysis\n")
    cat("\nCall:\n")
    print(x.s.linlir$call)
    cat("\nLIR settings:\n")
    cat(paste("p:", eval(x.s.linlir$config$p), "   beta:", eval(x.s.linlir$config$bet), 
        "   epsilon:", eval(x.s.linlir$config$epsilon), "   k.l:", 
        eval(x.s.linlir$config$k.l), "   k.u:", eval(x.s.linlir$config$k.u), 
        " \n"))
    cat(paste("confidence level of each confidence interval:", 
        eval(round(x.s.linlir$config$conf.lev.ci * 100, 2)), 
        "% \n"))
}
s.linlir <-
function (dat.idf, var = NULL, p = 0.5, bet, epsilon = 0, b.grid = 1000) 
{
    if (class(dat.idf) != "idf") {
        stop("The data must be provided as an *idf* object !\n")
    }
    if (!is.null(var)) {
        dat <- cbind(dat.idf[[which(names(dat.idf) == var[1])]], 
            dat.idf[[which(names(dat.idf) == var[2])]])
    }
    else {
        dat <- cbind(dat.idf[[1]], dat.idf[[2]])
    }
    lrm <- gen.lms(dat, p, bet, epsilon, k.u = 0)
    s.linlir <- vector("list", 2)
    if (!is.vector(lrm[[1]])) {
        s.linlir[[1]] <- lrm[[1]][, 1:2]
        attr(s.linlir, "names")[1] <- "f.lrm"
        s.linlir[[2]] <- lrm[[1]][1, 3]
        attr(s.linlir, "names")[2] <- "q.lrm"
    }
    else {
        s.linlir[[1]] <- lrm[[1]][1:2]
        attr(s.linlir, "names")[1] <- "f.lrm"
        s.linlir[[2]] <- lrm[[1]][3]
        attr(s.linlir, "names")[2] <- "q.lrm"
    }
    b.search <- lrm[[2]][1]
    if (!is.vector(s.linlir$f.lrm)) {
        for (j in 1:nrow(s.linlir$f.lrm)) {
            i <- 5
            prepara <- matrix(NA, 1, 2)
            while (i <= 15) {
                b.range <- c(s.linlir$f.lrm[j, 2] - i * b.search, 
                  s.linlir$f.lrm[j, 2] + i * b.search)
                par <- undom.para(dat, b.range, b.extra = 0, 
                  b.grid, q.lrm = s.linlir$q.lrm, p, bet, epsilon)
                if (par$b.undom[1] > b.range[1] & par$b.undom[2] < 
                  b.range[2]) {
                  i <- 99
                }
                else {
                  i <- i + 1
                }
            }
            if (i == 16) {
                b.extra <- runif(b.grid/2, -1e+09, 1e+09)
                para <- undom.para(dat, b.range, b.extra, b.grid, 
                  q.lrm = s.linlir$q.lrm, p, bet, epsilon)
            }
            else {
                para <- undom.para(dat, c(floor(par$b.undom[1]) - 
                  0.1, ceiling(par$b.undom[2]) + 0.1), b.extra = 0, 
                  b.grid, q.lrm = s.linlir$q.lrm, p, bet, epsilon)
            }
            prepara <- rbind(prepara, para[[3]])
            prepara <- unique(prepara)
        }
        preresult <- prepara[is.na(prepara[, 1]) == F, ]
        result1 <- c(min(preresult[, 1]), max(preresult[, 1]))
        attr(result1, "names") <- c("a.min", "a.max")
        result2 <- c(min(preresult[, 2]), max(preresult[, 2]))
        attr(result2, "names") <- c("b.min", "b.max")
        if (round(result2[1], 10) == b.range[1] | round(result2[2], 
            10) == b.range[2]) {
            print("b.range too small or possibly unbounded !")
        }
        result3 <- preresult
        para <- list(a.undom = result1, b.undom = result2, undom.para = result3)
    }
    else {
        i <- 5
        while (i <= 15) {
            b.range <- c(s.linlir$f.lrm[2] - i * b.search, s.linlir$f.lrm[2] + 
                i * b.search)
            par <- undom.para(dat, b.range, b.extra = 0, b.grid, 
                q.lrm = s.linlir$q.lrm, p, bet, epsilon)
            if (par$b.undom[1] > b.range[1] & par$b.undom[2] < 
                b.range[2]) {
                i <- 99
            }
            else {
                i <- i + 1
            }
        }
        if (i == 16) {
            b.extra <- runif(b.grid/2, -1e+09, 1e+09)
            para <- undom.para(dat, b.range, b.extra, b.grid, 
                q.lrm = s.linlir$q.lrm, p, bet, epsilon)
        }
        else {
            para <- undom.para(dat, c(floor(par$b.undom[1]) - 
                0.1, ceiling(par$b.undom[2]) + 0.1), b.extra = 0, 
                b.grid, q.lrm = s.linlir$q.lrm, p, bet, epsilon)
        }
    }
    s.linlir$a.undom <- para[[1]]
    s.linlir$b.undom <- para[[2]]
    s.linlir$undom.para <- para[[3]]
    klku <- kl.ku(dat.idf$n, p, bet, epsilon)
    if (epsilon <= 0) {
        conf.l <- sum(dbinom((klku[1] + 1):klku[2], dat.idf$n, 
            p))
    }
    else {
        if (p <= 0.5) {
            conf.l <- sum(dbinom((klku[1] + 1):klku[2], dat.idf$n, 
                (p + epsilon)))
        }
        else {
            conf.l <- sum(dbinom((klku[1] + 1):klku[2], dat.idf$n, 
                (p - epsilon)))
        }
    }
    as.conf <- pchisq(q = (-2 * log(bet)), df = 1)
    s.linlir$config <- list(p = p, beta = bet, epsilon = epsilon, 
        k.l = klku[1], k.u = klku[2], conf.lev.ci = conf.l, as.conf.lev.ci = as.conf)
    s.linlir$dat <- dat
    s.linlir$n <- dat.idf$n
    s.linlir$call <- match.call()
    class(s.linlir) <- "s.linlir"
    s.linlir
}
summary.idf <-
function (object, ...) 
{
    dat.idf <- object
    cat("\nSummary of interval data frame\n\n")
    cat(dat.idf$n, "observations\n")
    cat((length(dat.idf) - 1), "interval-valued variables \n")
    for (i in 1:(length(dat.idf) - 1)) {
        cat(paste("\nVariable", names(dat.idf)[i], ": \n\n", 
            sep = " "))
        print(summary(dat.idf[[i]]))
    }
}
summary.s.linlir <-
function (object, ...) 
{
    x.s.linlir <- object
    cat("\nSimple linear LIR analysis results \n")
    cat("\nCall:\n")
    print(x.s.linlir$call)
    cat("\nEstimated parameters of the function f.lrm:\n")
    if (!is.vector(x.s.linlir$f.lrm)) {
        cat("intercepts of f.lrm: ", eval(x.s.linlir$f.lrm[, 
            1]), "\n")
        cat("slopes of f.lrm: ", eval(x.s.linlir$f.lrm[, 2]), 
            "\n")
    }
    else {
        cat("intercept of f.lrm: ", eval(x.s.linlir$f.lrm[1]), 
            "\n")
        cat("slope of f.lrm: ", eval(x.s.linlir$f.lrm[2]), "\n")
    }
    cat("\nRanges of parameter values of the undominated functions:\n")
    cat("intercept of f in [", eval(x.s.linlir$a.undom[1]), ",", 
        eval(x.s.linlir$a.undom[2]), "]\n", sep = "")
    cat("slope of f in [", eval(x.s.linlir$b.undom[1]), ",", 
        eval(x.s.linlir$b.undom[2]), "]\n", sep = "")
    cat("\nNumber of observations:", x.s.linlir$n, "\n")
    cat("\nLIR settings:\n")
    cat(paste("p:", eval(x.s.linlir$config$p), "   beta:", eval(x.s.linlir$config$bet), 
        "   epsilon:", eval(x.s.linlir$config$epsilon), "   k.l:", 
        eval(x.s.linlir$config$k.l), "   k.u:", eval(x.s.linlir$config$k.u), 
        " \n"))
    cat(paste("confidence level of each confidence interval:", 
        eval(round(x.s.linlir$config$conf.lev.ci * 100, 2)), 
        "% \n"))
}
undom.a <-
function (dat, b, q.lrm, p = 0.5, bet, epsilon = 0) 
{
    n <- nrow(dat)
    k.l <- kl.ku(n, p, bet, epsilon)[1]
    dat.b <- dat[rowSums(abs(dat)) < Inf, ]
    n.b <- nrow(dat.b)
    dat.b0 <- dat[rowSums(abs(dat[, c(3, 4)])) < Inf, ]
    n.b0 <- nrow(dat.b0)
    d <- matrix(0, n, 2)
    if (b > 0) {
        d[, 1] <- dat[, 3] - b * dat[, 2]
        d[, 2] <- dat[, 4] - b * dat[, 1]
    }
    else {
        d[, 1] <- dat[, 3] - b * dat[, 1]
        d[, 2] <- dat[, 4] - b * dat[, 2]
    }
    d.l <- sort(d[, 1])
    d.u <- sort(d[, 2])
    a.undom <- matrix(0, (n - k.l), 2)
    for (j in 1:nrow(a.undom)) {
        a.undom[j, ] <- c(d.l[k.l + j] - q.lrm, d.u[j] + q.lrm)
    }
    result1 <- a.undom[order(a.undom[, 1]), ]
    attr(result1, "dimnames")[[2]] <- c("a.l", "a.u")
    if (nrow(matrix(a.undom[a.undom[, 1] <= a.undom[, 2], ], 
        ncol = 2)) < 1) {
        result2 <- "There is no undominated line with the chosen slope b!"
    }
    else {
        a.undom <- unique(a.undom, MARGIN = 1)
        preresult2 <- a.undom[a.undom[, 1] <= a.undom[, 2], ]
        if (is.vector(preresult2) == T) {
            result2 <- matrix(preresult2, 1, 2)
            attr(result2, "dimnames")[[2]] <- c("a.l", "a.u")
        }
        else {
            result2 <- preresult2
            i <- 1
            while (i < nrow(result2)) {
                for (k in 1:(nrow(result2) - i)) {
                  if (result2[i, 2] >= result2[i + k, 1]) {
                    result2[i, 2] <- max(result2[i, 2], result2[i + 
                      k, 2])
                  }
                }
                if (result2[i, 2] == max(result2[, 2])) {
                  i <- nrow(result2)
                }
                else {
                  i <- i + 1
                }
            }
            for (i in 2:nrow(result2)) {
                if (result2[i, 2] <= result2[i - 1, 2]) {
                  result2[i, 1] <- result2[i - 1, 1]
                  result2[i, 2] <- result2[i - 1, 2]
                }
            }
            result2 <- unique(result2)
            attr(result2, "dimnames")[[2]] <- c("a.l", "a.u")
        }
    }
    list(result1, result2)
}
undom.para <-
function (dat, b.range, b.extra = 0, b.grid = 1000, q.lrm, p = 0.5, 
    bet, epsilon = 0) 
{
    n <- nrow(dat)
    if (length(b.range) != 2) {
        stop("b.range must be given as a vector containing the lower and upper limit of the considered range !\n")
    }
    if (b.range[2] <= b.range[1]) {
        stop("b.range must be given as b.range[1] < b.range[2] !\n")
    }
    b.pot <- c(seq(b.range[1], b.range[2], by = (b.range[2] - 
        b.range[1])/b.grid), b.extra)
    para.undom <- matrix(NA, 1, 2)
    attr(para.undom, "dimnames")[[2]] <- c("a", "b")
    for (i in 1:length(b.pot)) {
        a.undom <- undom.a(dat, b = b.pot[i], q.lrm, p, bet, 
            epsilon)[[2]]
        if (is.matrix(a.undom) == T) {
            for (k in 1:nrow(a.undom)) {
                if (round(a.undom[[k, 2]], 10) == round(a.undom[[k, 
                  1]], 10)) {
                  para.undom <- rbind(para.undom, matrix(c(round(a.undom[k, 
                    1], 10), b.pot[i]), nrow = 1, ncol = 2))
                }
                else {
                  k.grid <- b.grid/2
                  para.undom <- rbind(para.undom, matrix(c(seq(a.undom[[k, 
                    1]], a.undom[[k, 2]], by = (a.undom[[k, 2]] - 
                    a.undom[[k, 1]])/k.grid), rep(b.pot[i], (k.grid + 
                    1))), nrow = (k.grid + 1), ncol = 2))
                }
            }
        }
    }
    if (which(is.na(para.undom[, 1])) > 1) {
        stop("Too many NA's !\n")
    }
    preresult <- para.undom[is.na(para.undom[, 1]) == F, ]
    result1 <- c(min(preresult[, 1]), max(preresult[, 1]))
    attr(result1, "names") <- c("a.min", "a.max")
    result2 <- c(min(preresult[, 2]), max(preresult[, 2]))
    attr(result2, "names") <- c("b.min", "b.max")
    if (round(result2[1], 10) == b.range[1] | round(result2[2], 
        10) == b.range[2]) {
        print("b.range too small or possibly unbounded !")
    }
    result3 <- preresult
    list(a.undom = result1, b.undom = result2, undom.para = result3)
}
