#' @title A Gibbs sampler using R
#' @description A Gibbs sampler using R
#' @param N the number of samples
#' @param thin the number of between-sample random numbers
#' @return a random sample of size \code{n}
#' @importFrom stats rbeta rbinom
#' @examples
#' \dontrun{
#' rnR <- gibbsR(100,10)
#' par(mfrow=c(2,1));
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' }
#' @export
gibbsR <- function(N, thin) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- 0
  y <- 0.5
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rbinom(1, 3, y)
      y <- rbeta(1, x + 1, 3 - x + 2)
    }
    mat[i, ] <- c(x, y)
  }
  mat
}
