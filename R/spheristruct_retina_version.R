##' Optimise the mapping from the flat outline to the sphere
##'
##' @title Optimise mapping
##' @param r reconstructedOutline object
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cutoff coefficient
##' @param nu Power to which to raise area
##' @param method Method to pass to \code{optim}
##' @param plot.3d If \code{TRUE} make a 3D plot in an RGL window
##' @param dev.flat Device handle for plotting grid to
##' @param dev.polar Device handle for plotting polar plot to
##' @param ... Extra arguments to pass to \code{\link{fire}}
##' @return reconstructedOutline object
##' @author David Sterratt
##' @export
solveMappingCart <- function(r, alpha=4, x0=0.5, nu=1, method="BFGS",
                               plot.3d=FALSE, dev.flat=NA, dev.polar=NA, ...) {
  phi <- r$phi
  lambda <- r$lambda
  R <- r$R
  phi0 <- r$phi0
  lambda0 <- r$lambda0
  Tt <- r$Tt
  A <- r$A
  Cut <- r$Cut
  Ct <- r$Ct
  Pt <- r$Pt
  Lt <- r$Lt
  Bt <- r$Bt
  Rsett <- r$Rsett
  i0t <- r$i0t
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)
  
  ## Optimisation and plotting 
  opt <- list()
  opt$x <- sphere.spherical.to.sphere.cart(phi, lambda, R)
  opt$conv <- 1

  ## Compute "mass" for each node
  minL <- rep(Inf, nrow(Pt))
  for (i in 1:nrow(Cut)) {
    minL[Cut[i,1]] <- min(minL[Cut[i,1]], Lt[i])
    minL[Cut[i,2]] <- min(minL[Cut[i,2]], Lt[i])
  }
  m <- 1/minL
  m <- m/mean(m)
  count <- 50
  
  while (opt$conv && count) {
    ## Optimise
    opt <- fire(opt$x,
                force=function(x) {Fcart(x, Ct, Lt, Tt, A, R, alpha, x0, nu)},
                restraint=function(x) {Rcart(x, R, Rsett, i0t, phi0, lambda0)},
                dt=1,
                nstep=200,
                m=m, verbose=TRUE, ...) 
    count <- count - 1
    ## Report
    E.tot <- Ecart(opt$x, Cu=Cut, L=Lt, R=R, T=Tt, A=A,
                   alpha=alpha, x0=x0, nu=nu)
    E.l <- Ecart(opt$x, Cu=Cut, L=Lt, R=R, T=Tt, A=A,
                 alpha=0, x0=x0, nu=0)

    s <- sphere.cart.to.sphere.spherical(opt$x, R)
    phi <-    s[,"phi"]
    lambda <- s[,"lambda"]
    ft <- flipped.triangles(phi, lambda, Tt, R)
    nflip <- sum(ft$flipped)

    ## Plot
    if (plot.3d) {
      sphericalplot(list(phi=phi, lambda=lambda, R=R,
                          Tt=Tt, Rsett=Rsett, gb=r$gb, ht=r$ht),
                     datapoints=FALSE)
    }

    if (!is.na(dev.flat)) {
      dev.set(dev.flat)
      flatplot(r, grid=TRUE, strain=TRUE,
                datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE)
    }

    if (!is.na(dev.polar)) {
      dev.set(dev.polar)
      r$phi <- phi
      r$lambda <- lambda
      projection(r)
    }
  }

  o <- merge(list(phi=phi, lambda=lambda, opt=opt, nflip=sum(ft$flipped),
                  E.tot=E.tot, E.l=E.l),
             r)
  o$mean.strain    <- mean(abs(getStrains(o)$spherical$strain))
  o$mean.logstrain <- mean(abs(getStrains(o)$spherical$logstrain))
  class(o) <- class(r)
  return(o)
}




##' Optimise the mapping from the flat outline to the sphere
##'
##' @title Optimise mapping
##' @param r reconstructedOutline object
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cut-off coefficient
##' @param nu Power to which to raise area
##' @param method Method to pass to \code{optim}
##' @param plot.3d If \code{TRUE} make a 3D plot in an RGL window
##' @param dev.flat Device handle for plotting flatplot updates to. If
##' \code{NA} don't make any flat plots
##' @param dev.polar Device handle for plotting polar plot updates
##' to. If \code{NA} don't make any polar plots.
##' @param control Control argument to pass to \code{optim}
##' @return reconstructedOutline object
##' @author David Sterratt
##' @export
optimiseMapping <- function(r, alpha=4, x0=0.5, nu=1, method="BFGS",
                             plot.3d=FALSE, dev.flat=NA, dev.polar=NA,
                             control=list()) {
  phi <- r$phi
  lambda <- r$lambda
  R <- r$R
  phi0 <- r$phi0
  lambda0 <- r$lambda0
  Tt <- r$Tt
  A <- r$A
  Cut <- r$Cut
  Ct <- r$Ct
  Pt <- r$Pt
  Lt <- r$Lt
  Bt <- r$Bt
  Rsett <- r$Rsett
  i0t <- r$i0t
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)
  
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi[-Rsett], lambda[-i0t])
  opt$conv <- 1
  count <- 0
  while (opt$conv) {
    ## Optimise
    opt <- stats::optim(opt$p, E, gr=dE,
                        method=method,
                        T=Tt, A=A, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R,
                        alpha=alpha,  N=Nt, x0=x0, nu=nu,
                        Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi,
                        verbose=FALSE, control=control)
    
    ## Report
    E.tot <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               alpha=alpha,  N=Nt, x0=x0, nu=nu,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
    E.l <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               alpha=0,  N=Nt, x0=x0, nu=nu,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)

    ft <- flipped.triangles(phi, lambda, Tt, R)
    nflip <- sum(ft$flipped)
    ## Decode p vector
    phi          <- rep(phi0, Nt)
    phi[-Rsett]  <- opt$p[1:Nphi]
    lambda       <- rep(lambda0, Nt)
    lambda[-i0t] <- opt$p[Nphi+1:(Nt-1)]

    ## Plot
    if (plot.3d) {
      sphericalplot(list(phi=phi, lambda=lambda, R=R,
                          Tt=Tt, Rsett=Rsett, gb=r$gb, ht=r$ht),
                     datapoints=FALSE)
    }

    if (!is.na(dev.flat)) {
      dev.set(dev.flat)
      flatplot(r, grid=TRUE, strain=TRUE,
               datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE)
    }

    if (!is.na(dev.polar)) {
      dev.set(dev.polar)
      r$phi <- phi
      r$lambda <- lambda
      projection(r)
    }
  }

  o <- merge(list(phi=phi, lambda=lambda, opt=opt, nflip=sum(ft$flipped),
                  E.tot=E.tot, E.l=E.l),
             r)
  o$mean.strain    <- mean(abs(getStrains(o)$spherical$strain))
  o$mean.logstrain <- mean(abs(getStrains(o)$spherical$logstrain))
  class(o) <- class(r)
  return(o)
}
