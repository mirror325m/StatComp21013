#' @title A method of PCA
#' @description A method of PCA
#' @param hatsigma the matrix
#' @param r the number of factors to be calculated
#' @param T max times bound of the iteration
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' a <- matrix(1:9,nrow=3)
#' Heteropca(a,2,100)
#' }
#' @export
Heteropca <- function(hatsigma, r, T) {
  N <- nrow(hatsigma)
  e <- 0.01
  matN <- hatsigma
  matN0 <- hatsigma
  for(i in 1:N)
    matN[i,i] <- 0;
  for(t in 1:T){
    Nsvd <- svd(matN);
    U <-  Nsvd$u[,1:r]
    D <- diag(Nsvd$d[1:r]);
    V <- Nsvd$v[,1:r];
    matN0 <- U %*% D %*% t(V);
    
    if(max(abs(matN-matN0))<e) break;
    for(i in 1:N)
      matN[i,i] <- matN0[i,i];
  }
  U
}
