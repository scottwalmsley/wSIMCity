#' ppmErr
#'
#' @param m the mass
#' @param m0 the reference mass
#'
#' @return the ppm error
#' @export
#'
ppmErr= function(m,m0){

   (m-m0)/m0*1e6

}

#' getMassTolRange
#'
#' @param m the mass / mz
#' @param ppm the ppm tolerance
#'
#' @return vector of minimum and maximum range
#' @export
#'
getMassTolRange <- function(m,ppm){

   dM <- ppm*m/1e6

   return(c(m-dM,m+dM))

}

#' Score similarity using Gaussian function
#'
#' Computes the similarity score using the gaussian assumption.  Use for Da data or RT data.
#' @param deltaVal the delta decimal value to score
#' @param tol the max decimal tolerance \(eg 0.01 for Da or 3 min for RT\)
#' @return decimal score 0-1
#' @export
similarityScore_gauss = function(deltaVal,tol){

   exp(-0.5 * (abs(deltaVal)/tol)^2)

}


#' Score similarity using Gaussian function
#'
#' @param deltaVal the delta decimal value in ppm to score
#' @param tol the max ppm tolerance
#'
#' @return decimal score 0-1
#' @export
similarityScore_laplace <- function(deltaVal,tol){
   exp(-1 * (abs(deltaVal)/tol))
}

