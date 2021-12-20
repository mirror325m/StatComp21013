## -----------------------------------------------------------------------------
set.seed(1234)
f <- function(x,n=100){##The original is larger.
u <- runif(n,0,x)
est <- x*mean(dbeta(u,3,3))
est
}
x <- 1:9/10
esti <- sapply(x,f)
theo <-pbeta(x,3,3)
cbind(esti,theo)

## -----------------------------------------------------------------------------
n <- 100##The original is larger.
sigma <- 1
u1 <- runif(n)
x1 <- sqrt(-2*sigma^2*log(1-u1))
x2 <- sqrt(-2*sigma^2*log(u1))
u2 <- runif(n)
x3 <- sqrt(-2*sigma^2*log(1-u2))
(sd((x1+x3)/2)^2-sd((x1+x2)/2)^2)/sd((x1+x3)/2)^2

## -----------------------------------------------------------------------------
c1 <- 1/(1-pnorm(1))
c2 <- 1/(1-pexp(1))
n <- 100##The original is larger.
u1 <- rnorm(n,0,1)
x1 <- u1[which(u1>1)]
y1 <- 1/c1*x1^2
est1 <- mean(y1)
u2 <- rexp(n,1)
x2 <- u2[which(u2>1)]
y2 <- 1/c2*x2^2/sqrt(2*pi)*exp(-x2^2/2+x2)
est2 <- mean(y2)
cbind(est1,est2)
cbind(sd(x1)^2,sd(x2)^2)

## -----------------------------------------------------------------------------
set.seed(1234)
m <- 100##The original is larger.
n <- 20
alpha <- 0.05
y1 <- numeric(m)
y2 <- numeric(m)
for(j in 1:m){
x <- rchisq(n,2)
a <- mean(x)-sd(x)/sqrt(n)*qt(1-alpha/2,df=n-1)
b <- mean(x)+sd(x)/sqrt(n)*qt(1-alpha/2,df=n-1)
y1[j] <- as.numeric(2>a & 2<b)
UCL <- (n-1) * var(x) / qchisq(alpha, df=n-1)
y2[j] <- as.numeric(4<UCL)
}
cbind(mean(y1),mean(y2))

## -----------------------------------------------------------------------------
m <- 100##The original is larger.
n <- 20
alpha <- 0.05
y1 <- numeric(m)
y2 <- numeric(m)
y3 <- numeric(m)
for(j in 1:m){
x1 <- rchisq(n,1)
t1 <- sqrt(n)*(mean(x1)-1)/sd(x1)
y1[j] <- as.numeric(abs(t1)>qt(1-alpha/2,df=n-1))
x2 <- runif(n,0,2)
t2 <- sqrt(n)*(mean(x2)-1)/sd(x2)
y2[j] <- as.numeric(abs(t2)>qt(1-alpha/2,df=n-1))
x3 <- rexp(n,1)
t3 <- sqrt(n)*(mean(x3)-1)/sd(x3)
y3[j] <- as.numeric(abs(t3)>qt(1-alpha/2,df=n-1))
}
cbind(mean(y1),mean(y2),mean(y3))

## -----------------------------------------------------------------------------
set.seed(123)
library(MASS)
alpha <- 0.05#6.8
d <- 2
mu <- c(1,1)
Sigma <- matrix(c(2,1,1,3),2,2)
n <- c(10,20,30,50,100)
cv <- qchisq(0.95,d*(d+1)*(d+2)/6)*6/n
stat <- function(x) {
  xm <- cbind(x[,1]-mean(x[,1]),x[,2]-mean(x[,2]))
  n <- nrow(x)
  sigma <- t(xm)%*%xm/n
  ins <- solve(sigma)
  s <- xm%*%ins%*%t(xm)
  y <- sum(s^3)/n^2
  y
}
p <- numeric(length(n))
m <- 50##The original is larger.
for (i in 1:length(n)) {
  sktests <- numeric(m)
  for(j in 1:m){
    x <- mvrnorm(n=n[i],mu,Sigma)
    y <- stat(x)
    sktests[j] <- as.numeric(stat(x)>=cv[i]) 
  }
  p[i] <- mean(sktests)
}
p
#6.10
alpha <- 0.1
n <- 30
m <- 50##The original is larger.
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qchisq(0.95,d*(d+1)*(d+2)/6)*6/n
for (j in 1:N) { #for each epsilon
e <- epsilon[j]
sktests <- numeric(m)
gener <- function(c){
  mu <- c(1,2)
  Sigma1 <- matrix(c(10,0,0,10),2,2)
  Sigma2 <- matrix(c(1,0,0,1),2,2)
  y <- c*mvrnorm(n=1,mu=c(1,2),Sigma=Sigma1)+(1-c)*mvrnorm(n=1,mu=c(1,2),Sigma=Sigma2)
  y
}
for (i in 1:m) { #for each replicate
r <- sample(c(0, 1), replace = TRUE, size = n, prob = c(1-e, e))
x <- t(sapply(r,gener))
sktests[i] <- as.integer(abs(stat(x)) >= cv)
}
pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",
xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
set.seed(123)
library('bootstrap')
n <- nrow(scor)
T <- cov(scor)
T1 <- svd(T)$d
theta.hat <- T1[1]/(T1[1]+T1[2]+T1[3]+T1[4]+T1[5])
B <- 200##The original is larger.
n <- nrow(scor)
theta.b <- numeric(B)
for (b in 1:B){
    i <- sample(1:n, size = n, replace = TRUE)
    scor1 <- scor[i,]
    T <- cov(scor1)
    T1 <- svd(T)$d
    theta.b[b] <- T1[1]/(T1[1]+T1[2]+T1[3]+T1[4]+T1[5])
}
bias <- mean(theta.b) - theta.hat
print(bias)
print(se.theta.b <- sd(theta.b))

## -----------------------------------------------------------------------------
theta.jack <- numeric(n)
for(i in 1:n){
    T <- cov(scor[-i,])
    T1 <- svd(T)$d
    theta.jack[i] <- T1[1]/(T1[1]+T1[2]+T1[3]+T1[4]+T1[5])}
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
print(bias)

se <- sqrt((n-1) *
mean((theta.jack - mean(theta.jack))^2))
print(se)

## -----------------------------------------------------------------------------
alpha <- c(.025, .975)
print(quantile(theta.b, alpha))
conf <- 0.95
zalpha <- qnorm(alpha)
z0 <- qnorm(sum(theta.b < theta.hat) / B)
L <- mean(theta.jack) - theta.jack
a <- sum(L^3)/(6 * sum(L^2)^1.5)
adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
limits <- quantile(theta.b, adj.alpha)
print(limits)

## -----------------------------------------------------------------------------
B <- 200##The original is larger.
n <- 40##The original is larger.
sk <- function(x) {
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}

m <- 50##The original is larger.

rx1sktests <- numeric(m)
lx1sktests <- numeric(m)
rx2sktests <- numeric(m)
lx2sktests <- numeric(m)
for (j in 1:m) {
x1 <- rnorm(n)
x2 <- rchisq(n,5)
theta.b <- numeric(B)
for (b in 1:B){
    i <- sample(1:n, size = n, replace = TRUE)
    x11 <- x1[i]
    x21 <- x2[i]
}
theta.x1 <- mean(x11)
cv <- qnorm(.975, 0, sqrt(6/B))
rx1sktests[j] <- as.integer(sk(x11) >= cv )
lx1sktests[j] <- as.integer(sk(x11) <= -cv ) 
rx2sktests[j] <- as.integer(sk(x21) >= cv )
lx2sktests[j] <- as.integer(sk(x21) <= -cv )
}
r1p.reject <- mean(rx1sktests)
l1p.reject <- mean(lx1sktests)
r2p.reject <- mean(rx2sktests)
l2p.reject <- mean(lx2sktests)
cbind(l1p.reject,r1p.reject)
cbind(l2p.reject,r2p.reject)

## -----------------------------------------------------------------------------
library(boot)
set.seed(123)
n <- 50##The original is larger.
x <- rnorm(n,0,1)
y <- rexp(n,1)
z <- cbind(x,y)
cor2 <- function(a,ix,dims){
p <- dims[1]
q <- dims[2]
d <- p + q
x <- z[ , 1:p]
y <- z[ix, -(1:p)]
return(cor(x,y,method="spearman"))
}
boot.obj <- boot(data = z, statistic = cor2, R = 999, sim = "permutation",dims = c(1, 1))
tb <- c(boot.obj$t0, boot.obj$t)
pcor <- mean(tb>=tb[1])
prefer <- cor.test(x,y,method="spearman")$p.value
c(pcor,prefer)

## -----------------------------------------------------------------------------
n <- 50##The original is larger.
x1 <- rnorm(n,2,1)
y1 <- rexp(n,2)
z1 <- c(x1,y1)
x2 <- rnorm(n,1,1)
y2 <- rexp(n,2)
z2 <- c(x2,y2)
x3 <- rt(n,1)
y3 <- 0.5*rnorm(n,0,2)+0.5*rnorm(n,2,1)
z3 <- c(x3,y3)

library(RANN)
library(energy)
library(Ball)
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1)
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
N <- c(n, n)

boot.obj <- boot(data = z1, statistic = Tn, R = 999,
sim = "permutation", sizes = N,k=3)
ts <- c(boot.obj$t0,boot.obj$t)
N1p.value <- mean(ts>=ts[1])

boot.obs <- eqdist.etest(z1, sizes=N, R=999)
E1p.value <- boot.obs$p.value

B1p.value = bd.test(x = x1, y = y1, num.permutations=999)$p.value

c(N1p.value,E1p.value,B1p.value)

boot.obj <- boot(data = z2, statistic = Tn, R = 999,
sim = "permutation", sizes = N,k=3)
ts <- c(boot.obj$t0,boot.obj$t)
N2p.value <- mean(ts>=ts[1])

boot.obs <- eqdist.etest(z2, sizes=N, R=999)
E2p.value <- boot.obs$p.value

B2p.value = bd.test(x = x2, y = y2, num.permutations=999)$p.value

c(N2p.value,E2p.value,B2p.value)

boot.obj <- boot(data = z3, statistic = Tn, R = 999,
sim = "permutation", sizes = N,k=3)
ts <- c(boot.obj$t0,boot.obj$t)
N3p.value <- mean(ts>=ts[1])

boot.obs <- eqdist.etest(z3, sizes=N, R=999)
E3p.value <- boot.obs$p.value

B3p.value = bd.test(x = x3, y = y3, num.permutations=999)$p.value

c(N3p.value,E3p.value,B3p.value)

## -----------------------------------------------------------------------------
set.seed(123)
f <- function(x, location = 0, scale = 1) {
stopifnot(scale > 0)
return(1/(scale*pi*(1+((x-location)/scale)^2)))
}
m <- 100##The original is larger.
x <- numeric(m)
x[1] <- rnorm(1, mean=1)
k <- 0
u <- runif(m)
for (i in 2:m) {
xt <- x[i-1]
y <- rnorm(1, mean = xt)
num <- f(y) * dnorm(xt, mean = y)
den <- f(xt) * dnorm(y, mean = xt)
if (u[i] <= num/den) x[i] <- y else {
x[i] <- xt
k <- k+1
}
}
bx <- x[-(1:10)]##The original is larger.
r <- 1:9/10
qx <- quantile(bx,r)
q <- qcauchy(r)
t(cbind(qx,q))

## -----------------------------------------------------------------------------
N <- 500 #length of chain    ##The original is larger.
burn <- 100 #burn-in length  ##The original is larger.
X <- matrix(0, N, 2) #the chain, a bivariate sample
rho <- -.75 #correlation
a <- 1
b <- 2
n <- 3
X[1, ] <- c(0,0.5) #initialize
for (i in 2:N) {
y1 <- X[i-1, 2]
X[i, 1] <- rbinom(n = 1, size = n, prob = y1)
x1 <- X[i, 1]
c1 <- x1+a
c2 <- n-x+b
X[i, 2] <- rbeta(1, c1, c2)
}
b <- burn + 1
x <- X[b:N, ]

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

## -----------------------------------------------------------------------------
cauchy.chain <- function(N,X1){
f <- function(x, location = 0, scale = 1) {
stopifnot(scale > 0)
return(1/(scale*pi*(1+((x-location)/scale)^2)))
}
m <- N
x <- numeric(m)
x[1] <- X1
k <- 0
u <- runif(m)
for (i in 2:m) {
xt <- x[i-1]
y <- dnorm(1, mean = xt)
num <- f(y) * dnorm(xt, mean = y)
den <- f(xt) * dnorm(y, mean = xt)
if (u[i] <= num/den) x[i] <- y else {
x[i] <- xt
k <- k+1
}
}
return(x)
}
k <- 4
n <- 1500##The original is larger.
b <- 100##The original is larger.
x0 <- c(-0.75,-0.25,0.25,0.75)
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- cauchy.chain(n, x0[i])
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
#plot psi for the four chains
#par(mfrow=c(2,2))
for (i in 1:k)
plot(psi[i, (b+1):n], type="l",
xlab=i, ylab=bquote(psi))
#par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
f <- function(k,a){
  d <- length(a)
  l <- (k+1)*log(sum(a*a))-lgamma(k+1)-k*log(2)-log(2*k+1)-log(2*k+2)+lgamma((d+1)/2)+lgamma(k+1.5)-lgamma(k+d/2+1)
  exp(l)*(-1)^k
}

a <- c(1:20)
k <- 100
f(k,a)

S <- function(a,n=100){#the amount of terms calculated  ##The original is larger.
  sum(sapply(n:0,function(k) f(k,a)))
}

a <- c(1,2)
S(a)

## -----------------------------------------------------------------------------
k <- c(4,25)
m <- length(k)
asolve <- numeric(m)

f <- function(k,a){
  lgamma(k/2)-lgamma((k-1)/2)-log(k-1)/2+log(integrate(function(u) (1+u^2/(k-1))^(-k/2),lower=0,upper=sqrt(a^2*(k-1)/(k-a^2)))$value)-lgamma((k+1)/2)+lgamma((k)/2)+log(k)/2-log(integrate(function(u) (1+u^2/(k))^(-(k+1)/2),lower=0,upper=sqrt(a^2*(k)/(k+1-a^2)))$value)
}

for(i in 1:m){
  asolve[i] <- uniroot(function(a) f(k[i],a),c(0.1,k[i]^0.5-0.01))$root
}

cbind(k,asolve)

## -----------------------------------------------------------------------------
y <- c(0.54,0.48,0.33,0.43,1.00,1.00,0.91,1.00,0.21,0.85)
print(sum(y)/(length(y)-sum(y[which(y==1)])))

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
data <- mtcars
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
modl <- lapply(formulas,lm,data = data)
r2list <- lapply(modl,rsq)
r2list

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
modl <- lapply(bootstraps,lm,formula = mpg ~ disp)
r2list_lapply <- lapply(modl,rsq)
r2list_lapply

n <- length(bootstraps)
modl <- vector("list", n)
for (i in 1:n) {
modl[[i]] <- lm(formula = mpg ~ disp,data = bootstraps[[i]])
}
r2list_for <- lapply(modl,rsq)
r2list_for

## -----------------------------------------------------------------------------
data <- mtcars
vapply(1:ncol(data),function(i) sd(data[,i]),FUN.VALUE=0)
data <- InsectSprays
vapply(1:ncol(data),function(i) if(typeof(data[,i])=="double")sd(data[,i]) else 0,FUN.VALUE=0)

## -----------------------------------------------------------------------------
library("parallel")
mcsapply <- function(mc.cores,X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE){
  cl <- makeCluster(getOption("cl.cores", mc.cores))
  y <- parSapply(cl, X, FUN, ..., simplify = TRUE,
          USE.NAMES = TRUE, chunk.size = NULL)
  stopCluster(cl)
  y
}

mcvapply <- function(mc.cores,X, FUN, FUN.VALUE, ..., USE.NAMES = TRUE){
  cl <- makeCluster(getOption("cl.cores", mc.cores))
  y <- parSapply(cl, X, FUN, ..., simplify = TRUE,
          USE.NAMES = TRUE, chunk.size = NULL)
  a <- which(typeof(y)!=typeof(FUN.VALUE))
  y[a] <- rep(NA,length(a))
  stopCluster(cl)
  y
}

## -----------------------------------------------------------------------------
library("StatComp21013")
suppressWarnings(library('Rcpp'))
suppressWarnings(library('microbenchmark'))
#source("gibbsR.R")
#sourceCpp("gibbsC.cpp")
gibbR <- gibbsR(100,10)
gibbC <- gibbsC(100,10)
xR <- gibbR[,1]
yR <- gibbR[,2]
xC <- gibbC[,1]
yC <- gibbC[,2]
qqplot(xR,xC)
qqplot(yR,yC)
ts <- microbenchmark(gibbR=gibbsR(100,10),gibbC=gibbsC(100,10))
summary(ts)[,c(1,3,5,6)]

