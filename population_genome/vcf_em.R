library(mixtools)
library(depmixS4)
state2col <- function(x) {
  col <- rep("black", length(x))
  for (i in 1:length(x)) {
    if (x[i]==1) {
      col[i] <- c("fill_color=yellow")
    } else {
      col[i] <- c("fill_color=purple")
    }
  }
  return(col)
}

# fit 2 speed model
## with EM
vcf = read.table("/Users/alexwang/0data/0mango/genome/2speed/intraspecific/hn47_density10k.txt")
vcfem = normalmixEM(vcf$V4/10,k=2)
#### model assessment
mu <- mean(vcf$V4/10)
sigma <- sd(vcf$V4/10)
# D'Test and KS test check normality
ks.test(vcf$V4/10, "pnorm",mu, sigma)
library(fBasics)
dagoTest(x = vcf$V4/10)
## log likehood function of Gaussian
# methods 1
normal <- function(theta,x){
  mu <- theta[1]
  sigma <- theta[2]
  n <- length(x)
  logL <- -0.5*n*log(2*pi)-0.5*n*log(sigma^2)-(1/(2*sigma^2))*sum((x-mu)**2)
  return (-logL)
}
res <- optim(c(mu,sigma),normal,x=vcf$V4/10)  # optimization
# method 2
LL <- function(mu, sigma) {
  R = dnorm(vcf$V4/10, mu, sigma)
  -sum(log(R))
}
library(stats4)
fit <- mle(LL, start = list(mu = mu, sigma = sigma))
# summary(fit)
# logLik(fit)
aic.gaussian <- 2*2 - (-2*res$value)
aic.mixGaussian <- 2*5 - 2*vcfem$loglik
# Likelihood Rate Test judge whether normal distribution can substitude mixGaussian model
LR = 2*((vcfem$loglik)- (-res$value))
dchisq(x = LR, df = 2)
# expand y lim for plot.mixEM
plot.em <- function(x,
           whichplots = 1,
           loglik = 1 %in% whichplots,
           density = 2 %in%
             whichplots,
           xlab1 = "Iteration",
           ylab1 = "Log-Likelihood",
           main1 = "Observed Data Log-Likelihood",
           col1 = 1,
           lwd1 = 2,
           xlab2 = NULL,
           ylab2 = NULL,
           main2 = NULL,
           col2 = NULL,
           lwd2 = 2,
           alpha = 0.05,
           marginal = FALSE,
           ...) {
    mix.object <- x
    k <- ncol(mix.object$posterior)
    x <- sort(mix.object$x)
    a <- hist(x, plot = FALSE)
    maxy <-
      max(max(a$density), 0.3989 * mix.object$lambda / mix.object$sigma) + 0.01
    if (is.null(main2)) {
      main2 <- "Density Curves"
    }
    if (is.null(xlab2)) {
      xlab2 <- "Data"
    }
    if (is.null(col2)) {
      col2 <- 2:(k + 1)
    }
    hist(
      x,
      prob = TRUE,
      main = main2,
      xlab = xlab2,
      ylim = c(0, maxy),
      ...
    )
    if (length(mix.object$mu) == 1) {
      arbvar <- TRUE
      mix.object$sigma <- mix.object$scale * mix.object$sigma
      arbmean <- FALSE
    }
    if (length(mix.object$mu) == k && length(mix.object$sigma) ==
        1) {
      arbmean <- TRUE
      arbvar <- FALSE
    }
    if (length(mix.object$sigma) == k && length(mix.object$mu) == k) {
      arbmean <- TRUE
      arbvar <- TRUE
    }
    for (i in 1:k) {
      lines(
        x,
        mix.object$lambda[i] * dnorm(x, mean = mix.object$mu[i *arbmean + (1 - arbmean)], sd = mix.object$sigma[i *arbvar + (1 - arbvar)]),
        col = col2[i],
        lwd = lwd2
      )
    }
  }
par(mar=c(4.5,4.5,1,1))
plot.em(vcfem,breaks=50,lwd2 = 2, w=1.1, xlab2="Variants/kb",main2="",
        col2=c("#1B9E77","#D95F02"))
# str(vcfem)
library(shape)
Arrows(5.5, 0.05, 4.5, 0.057, arr.length = 0.2, segment = T, code = 1, arr.adj = 0.5, col="purple")
Arrows(18, 0.055, 19, 0.062, arr.length = 0.2, segment = T, code = 1, arr.adj = 0.5, col="orange" )
text(4.4, 0.060, expression(paste("slow: ", mu, " = 7.78,  ", sigma," = 3.68")))
text(21, 0.066, expression(paste("fast: ", mu, " = 14.78,  ", sigma," = 2.59")))

# mtext("A", adj=0.01, line=-2, outer=T)

## decode with veterbi
msp <- depmix(V4/4~1,nstates=2,data=vcf)
set.seed(1)
fmsp <- fit(msp)
#plot(ts(posterior(fmsp)))  
state = posterior(fmsp)$state

# check accuracy
ss <- vcfem$posterior
aa <- ifelse(ss[,1]>ss[,2],1,2)
xx <- aa == state
length(xx[which(xx==TRUE)])/length(xx)
# write.table(paste(vcf$V1, vcf$V2, vcf$V3, state),file="/Users/alexwang/Downloads/density/state.tsv", sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
# write.table(paste(vcf$V1, vcf$V2, vcf$V3, state2col(state)),file="/Users/alexwang/Downloads/density/state.txt", sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

