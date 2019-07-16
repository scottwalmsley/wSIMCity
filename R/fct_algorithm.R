


#' Mixture model: Laplace and Uniform
#'
#' Performs expectation maximization algorithm on the dM values determined from the data.
#'
#' @param x vector of numbers.
#' @param tol numeric number indicating tolerance in ppm.
#' @param maxiter maximum number of iterations.
#' @param boost A weight value applied to Laplace distribution.
#' @param instrument.tol numeric value assigned to
#' @param x 
#' @param dM 
#' @param tol 
#' @param maxiter 
#' @param boost 
#' @param instrument_tol 
#' 
#' @return list of data model
#' @export
laplace_unif_EM = function(x,dM,tol = 1e-5, maxiter = 1000,boost = 2,instrument_tol = .01) {
  
  
  
  weights= NULL
  
  ox = x
  wo = seq(1,length(x)) # index the original data
  dens = density(x, bw = 0.001, n =2048)
  m = dens$x[dens$y == max(dens$y, na.rm=T)]
  correction_ppm = m
  
  dens_dm = density(dM, bw = 0.001, n =2048)
  m.d = dens_dm$x[dens_dm$y == max(dens_dm$y, na.rm=T)]
  correction_dM = m.d
  ## adjust the data to the mode
  dM = dM - m.d
  x = x - m
  dens = density(x, bw = 0.001, n =2048)
  m = dens$x[dens$y == max(dens$y, na.rm=T)]
  
  
  
  #if(boost > 0){
  #  rn = rnorm(length(o.x)*boost, mean=m, sd = 1)
  #  x = c(o.x, rn)
  
  #}
  
  n = length(x)
  
  b = 1  ## initial value
  z = rep(0.15, n)
  pi_match = 0.5
  
  par_diff = 1e5
  iter = 0
  
  ### mode m is fixed, just update the variance of laplace 2b^2
  i = 0
  while ((iter < maxiter) & (abs(par_diff) > tol)) {
    new_v = sum(z * (x - m) ^ 2) / sum(z)
    new_b = sqrt(new_v / 2)
    
    dens_x = dlaplace(x, m, sqrt(new_b))
    if(boost > 0){
      #weights = weight_laplace(x,m,instrument_tol,boost)
      weights = weight_gauss(dM,m,instrument_tol,boost)
      dens_x = dens_x * weights
    }
    
    p1 = pi_match * dens_x #dlaplace(x, m, sqrt(new_b))
    
    p2 = (1 - pi_match) * dunif(x, min(x) - 0.01, max(x) + 0.01)
    z = p1 / (p1 + p2)
    new_pi_match = mean(z)
    
    par_diff = max(abs(c(new_pi_match - pi_match, new_b - b)))
    iter = iter + 1
    if (iter %% 10 == 0)
      cat(paste(iter,"\n"))
    
    b = new_b
    pi_match = new_pi_match
    i=i+1
  }
  cat("\n")
  o = order(ox)
  ox = ox[o]
  nz = z
  z = z[wo]
  z = z[o]
  
  
  res = list("cr.dM" =correction_dM, "cr.ppm" = correction_ppm,
    mu = m, b = b,o.x= ox,x= x, z=z ,weights = weights,n.z = nz,dens_x = dens_x, xmin = min(x) - 0.01, xmax = max(x) + 0.01, pos = pi_match,unif = mean(p2,na.rm =TRUE)
  )
  res
}



#' Score similarity using Gaussian function
#'
#' Computes the similarity score using the gaussian assumption.  Use for Da data or RT data.
#'
#' @param v0 numeric the reference value
#' @param v1 numeric the observed value
#' @param tol numeric the tolerance window in generic units
#'
#' @return numeric the score
#' @export
#'
#' @examples
#' val1 <- 0.2
#' val2 <- 0.21
#' similarityScore_gauss(val2-val1,tol= 0.1)
#' 0.9950125
#'
#' val1 <- 0.35
#' similarityScore_gauss(val2-val1,tol= 0.1)
#' 0.3246525
#'
#'
similarityScore_gauss = function(deltaVal,tol){
  exp(-0.5 * (abs(deltaVal)/tol)^2)
}


#' Score similarity using the Laplace function
#'
#' Computes the similarity score usingg the laplace assumption. Used for data measured with ppm errors.
#'
#'
#' @param v0 numeric the reference value
#' @param v1 numeric the observed value
#' @param tol numeric the tolerance window in generic units
#'
#' @return numeric the score
#' @export
#'
#' @examples
#' #' val2 <- 0.21
#' similarityScore_laplace(val2-val1,tol= 0.1)
#' 0.9048374
#'
#' val1 <- 0.35
#' similarityScore_gauss(val2-val1,tol= 0.1)
#' 0.2231302
similarityScore_laplace <- function(deltaVal,tol){
  exp(-1* (abs(deltaVal)/tol))
}


#' weight_laplace
#'
#' Assigned a weight to the distribution of X values given expected instrument performance.
#' For example, dM is the input X values, tol is a tolerance in ppm, and boost is the multiplier.
#'
#' @param dM vector or numeric the vector of measurement errors
#' @param tol numeric the ppm tolerance of the instrument
#' @param boost numeric the boost value (0-2)
#' @param mu numeric the mean of the laplace distribution
#'
#' @return numeric or vector the weighted value(S).
#' @export
#'
#' @examples 
#' weight_laplace(dM,mu=0.1, tol = 5, boost = 2)
weight_laplace <- function(dM,mu,tol,boost){
  boost*exp(-1* (abs(mu-dM)/tol))
}



#' Apply weights based on mass decimal
#'
#' @param dM numeric delta mass
#' @param mu numeric mean
#' @param tol numeric tolerance
#' @param boost numeric integer value to boost signal
#'
#' @return vector of weights
#' @export
#'
# @examples
weight_gauss <- function(dM,mu,tol,boost){
  boost*exp(-0.5 * (abs(mu-dM)/tol)^2)
}



#' Compute Join Aligner Score
#'
#' Performs join aligner score to produce the feature score
#'
#' @param m0 numeric reference mass
#' @param m1 numeric query mass
#' @param aM numeric the mass score contribution (a+b must add up to one)
#' @param r0 numeric reference rt
#' @param r1 numeric query rt
#' @param aR numeric the RT score contribition (see aM)
#'
#' @return numeric the feature score
#' @export
#'
#' @examples
#' score_feature(dM = 1.2,aM = 0.5,tolM = 5, drt = 0.2, aR = 0.5, tolR = 0.3)
#' 0.7936826
#'
#' score_feature(dM = 0.2,aM = 0.5,tolM = 5, drt = 0.1, aR = 0.5, tolR = 0.3)
#' 0.9533745
#'
score_feature  =  function(dM,aM,tolM,drt,aR,tolR){
  aM * similarityScore_laplace(dM,tolM) + aR * similarityScore_gauss(drt,tolR)
}




#' Get a mass tolerance window
#'
#' @param m numeric mass
#' @param ppm numeric ppm tolerance
#'
#' @return numeric vector containing minimum and maximum data values
#' @export
#' 
# @examples
getMassTolRange <- function(m,ppm){
  
  dM <- ppm*m/1e6
  return(c(m-dM,m+dM))

}
