#' @title A method of figure out precision matrix
#' @description A method of figure out precision matrix
#' @param sigmahat the matrix to be inverse
#' @param rho0 the radius
#' @return a list of precision matrix,parameter gamma, and eigens of precision matrix
#' @importFrom stats uniroot
#' @examples
#' \dontrun{
#' #a <- matrix(c(2,0,0,3),nrow=2)
#' b <- robpre(a,0.002)
#' s <- svd(b)
#' s$u
#' s$d
#' }
#' @export
robpre <- function(sigmahat, rho0) {
  S <- svd(sigmahat)
  U <- S$u
  lambda <- S$d
  f1 <- function(x,gam){
    sqrt(x^2*gam^2+4*x*gam)
  }
  target <- function(gamma,lamb,rhor){
    p <- length(lamb)
    gamma*(rhor^2-sum(lamb)/2)-p+sum(sapply(lamb,f1,gam=gamma))/2
  }
  minum <- 0
  maxum <- 100
  while(target(gamma=maxum,lamb=lambda,rhor=rho0)<0){
    minum <- maxum
    maxum <- 2*maxum
  }
  gamma0 <- uniroot(target,interval=c(minum,maxum),lamb=lambda,rhor=rho0)$root
  f2 <- function(x,gam){
    gam*(1-(sqrt(x^2*gam^2+4*x*gam)-x*gam)/2)
  }
  x <- sapply(lambda,f2,gam=gamma0)
  list(premat=U%*%diag(x)%*%t(U),gamma=gamma0,est=x)
}
