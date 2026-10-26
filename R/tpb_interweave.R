# A stepping-out slice update on an unconstrained scalar.
tpb_slice <- function(x, log_density, width = 1) {
  current <- log_density(x)
  if (!is.finite(current)) stop("Non-finite starting slice density")
  height <- current - rexp(1)
  left <- x - runif(1) * width
  right <- left + width
  while (log_density(left) > height) left <- left - width
  while (log_density(right) > height) right <- right + width
  repeat {
    proposal <- runif(1, left, right)
    if (log_density(proposal) >= height) return(proposal)
    if (proposal < x) left <- proposal else right <- proposal
  }
}

# Gaussian likelihood with beta integrated out: covariance sigmaSq*I + X D X'.
tpb_log_marginal_y <- function(X,y,log_var,sigmaSq) {
  n <- nrow(X)
  log_sigma <- log(sigmaSq)
  shift <- max(log_sigma,log_var)
  A <- sweep(X,2,exp((log_var-shift)/2),"*")
  gram <- tcrossprod(A)
  noise <- exp(log_sigma-shift)
  loses_noise <- any(diag(gram)+noise==diag(gram))
  R <- if (loses_noise) NULL else
    tryCatch(chol(gram+diag(noise,n)),error=function(e)NULL)
  if (!is.null(R)) {
    residual <- forwardsolve(t(R),y*exp(-shift/2))
    return(-.5*(n*shift+2*sum(log(diag(R)))+sum(residual^2)))
  }
  decomp <- svd(A,nv=0)
  rotated <- drop(crossprod(decomp$u,y))
  residual <- y-drop(decomp$u %*% rotated)
  log_eigenvalue <- shift+2*log(decomp$d)
  log_variance <- log_sigma+tpb_softplus(log_eigenvalue-log_sigma)
  -.5*((n-length(rotated))*log_sigma+sum(residual^2)/sigmaSq+
        sum(log_variance+rotated^2*exp(-log_variance)))
}

tpb_collapsed_target <- function(X,y,log_var,a,b,phi,sigmaSq,hp) {
  x <- log_var-log(phi)
  negative <- x<0
  z <- x*ifelse(negative,a,b)
  transform <- function(lp) {
    param <- exp(lp)
    x <- z/ifelse(negative,param[1],param[2])
    list(param=param,x=x,log_var=lp[3]+x)
  }
  target <- function(lp) {
    value <- transform(lp)
    a <- value$param[1]; b <- value$param[2]
    if (any(!is.finite(value$log_var)) || any(!is.finite(value$param)) || any(value$param<=0))
      return(-Inf)
    log_z <- sum(-abs(z)-(a+b)*tpb_softplus(-abs(value$x))-
                   lbeta(a,b)-log(ifelse(negative,a,b)))
    tpb_log_marginal_y(X,y,value$log_var,sigmaSq)+log_z+
      tpb_log_hyperprior(lp,hp)+sum(lp)
  }
  list(target=target,transform=transform)
}

tpb_interweave <- function(X,y,log_var,a,b,phi,sigmaSq,hp,
                           woodbury=TRUE,diagX=FALSE,slice_width=rep(1,3)) {
  representation <- tpb_collapsed_target(X,y,log_var,a,b,phi,sigmaSq,hp)
  lp <- log(c(a,b,phi))
  for (j in 1:3) {
    target <- function(value) {
      proposal <- lp
      proposal[j] <- value
      representation$target(proposal)
    }
    lp[j] <- tpb_slice(lp[j],target,slice_width[j])
  }
  value <- representation$transform(lp)
  # beta was discarded in this block. Redraw it immediately under the new D.
  gaussian <- Gibbs_beta(X,y,sigmaSq=sigmaSq,log_local_var=value$log_var,
                         woodbury=woodbury,diagX=diagX,return_log=TRUE)
  list(a=value$param[1],b=value$param[2],phi=value$param[3],log_var=value$log_var,
       beta=gaussian$beta,log_abs_beta=gaussian$log_abs_beta,sign_beta=gaussian$sign_beta)
}
