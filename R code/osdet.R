##### phi(theta) #####
gauss_phi <- function(theta, mu, sigma2){
  phi = mu*theta + 0.5*sigma2*theta^2
  return(phi)
}

##### osdet algorithm #####
gauss_osdet <- function(n, b, l, s = 0, mu, sigma2){
  ## to compute estimates a_n = P(S_n > n*b)
  ### with S_n = X1 +...+ Xn
  
  k = 0
  L = 1
  w = b
  mu1 = mu
  sigma21 = sigma2
  
  ## perform state-dependent OET while event of interest is rare
  repeat{
    
    # update parameters by exp tilt
    theta <- max((w - mu1) / sigma21, 0)
    mu1 = mu1 + theta*sigma21
    
    # sample from distn
    X = rnorm(1)
    
    # update IS measure
    L = exp(theta*X - gauss_phi(theta, mu = mu1, sigma2 = sigma21))*L
    
    # update RW
    s <- s + X*L
    k = k + 1
    w <- max((n*b - s)/(n - k), 0)
    
    if (n == k | w <= 1/sqrt(n-k) | w > l) {
      break
    }
  }
  
  ### 'while' loop ended implies event of interest is no longer rare
  #### so finish algorithm with standard OET
  barS = 0
  sbar = 0
  if (k < n){
    size = length((k+1):n)
    theta = max((w - mu1) / sigma21,0)
    mu1 = mu1 + theta*sigma21
    
    barS = sum(rnorm(size))
    
    L = exp(theta*barS - (n-k)*gauss_phi(theta, mu = mu1, sigma2 = sigma21))*L
    sbar = barS*L
  } 
  
  Ind <- 0+(s + sbar > n*b)
  Y = L*Ind
  return(Y)
}
