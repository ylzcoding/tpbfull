Gibbs_phi = function(p, nu, lambda, w, beta){
  return(1/rgamma(1, shape = (p+1)/2, rate = 1/w+0.5*sum(beta^2 * lambda / nu)))
}

# w|phi_hyper ~ InvGamma(1.5, 0.5)
Gibbs_w = function(phi, phi_hyper = 1){return(1/rgamma(1, 1, 1/phi+1/phi_hyper))}

Gibbs_sigmaSq = function(n, X, y, beta){
  return(1/rgamma(1, shape = n/2+3/2, rate = 1/2 + 0.5*(t(y - X%*%beta)%*%(y - X%*%beta))))
}
