## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library("StatComp21013")
#source("robpre.R")
#a <- matrix(c(2,0,0,3),nrow=2)
#b <- robpre(a,0.002)
#s <- svd(b)
#s$u
#s$d

lambda <- sapply(1:5,function(x) 10^(x-3))
hatsigma <- diag(lambda)
robpre(hatsigma,0.001)
q <- 10^(seq(-2,-0.8,0.01))
gammastar <- sapply(q,function(x,sigma) robpre(hatsigma,x)$gamma,sigma=hatsigma)
x <- sapply(q,function(x,sigma) robpre(hatsigma,x)$est,sigma=hatsigma)
plot(gammastar,type="l",log="y")
plot(x[5,]/x[1,],type="l",log="y")
plot(x[5,],type="l")

## -----------------------------------------------------------------------------
library("StatComp21013")
#source("HeteroPCA.R")
library(MASS)

set.seed(123)
p <- 5
n <- 400
r <- 5
w <- runif(p)
sigma2 <- runif(p)
U0 <- matrix(runif(p*r),nrow=p)
U <- qr.Q(qr(diag(w)%*%U0))
sigma0 <- U%*%diag(1:r)%*%t(U)
x <- mvrnorm(n,mu=rep(0,r),Sigma = sigma0)
e <- mvrnorm(n,mu=rep(0,r),Sigma = diag(sigma2*sigma2))
y <- x+e
for(i in 1:r){
  y[,i]=y[,i]-mean(y[,i])
}
hat <- t(y)%*%y/(n-1)
U0 <- Heteropca(hat,2,500)
1-svd(U[-(1:2),]%*%U0)$d[1]

