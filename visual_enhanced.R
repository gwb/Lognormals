require(ggplot2)
require(reshape2)
require(grid)
require(PerformanceAnalytics)
source('/Users/gwb/Hacks/Projects/myRUtils/graphics.R')



get.convo.emp.cdf <- function(m1, sd1, m2, sd2, n.points, n.sample){

  v1 <- log(exp(sd1)^2)
  v2 <- log(exp(sd2)^2)
  
  s1 <- sqrt(log(1 + exp(v1-2*m1)))
  u1 <- m1 - (s1^2)/2

  s2 <- sqrt(log(1 + exp(v2-2*m2)))
  u2 <- m2 - (s2^2)/2
  
  samples.cdf = list()
  full.data = NULL
  # Defines samples
  for (it in seq(1, n.sample)){
    y1 <- rnorm(n.points, u1, s1)
    y2 <- rnorm(n.points, u2, s2)
    simul.sample <- y1 + log(1 + exp(y2 - y1))

    #browser()
    full.data <- c(full.data, simul.sample)
    samples.cdf <- cbind(samples.cdf, ecdf(simul.sample))
  }

  #browser()
  qtile <- quantile(x=full.data, probs=c(0.01, 0.99))
  data.x <- seq(qtile[1], qtile[2], length.out=4000)

  u.bound <- NULL
  m <- NULL
  l.bound <- NULL
  for (i in data.x){
    #browser()
    samples.results <- do.call("cbind", lapply(samples.cdf, function(x) x(i)))
    m <- c(m, rowMeans(samples.results))
    samples.qtile <- quantile(x=samples.results, probs=c(0.05, 0.95))
    u.bound <- c(u.bound, samples.qtile[2])
    l.bound <- c(l.bound, samples.qtile[1])
  }

  return(rbind(data.x, l.bound, m, u.bound))
  
}


test.analytical.closeness <- function(m1, sd1, m2, sd2, n.sim, n.points=100, n.samples=200){

  v1 <- log(exp(sd1)^2)
  v2 <- log(exp(sd2)^2)
  
  s1 <- sqrt(log(1 + exp(v1-2*m1)))
  u1 <- m1 - (s1^2)/2

  s2 <- sqrt(log(1 + exp(v2-2*m2)))
  u2 <- m2 - (s2^2)/2
  
  dt <- get.convo.emp.cdf(m1, sd1, m2, sd2, n.points, n.samples)
  dimnames(dt) <- NULL
  dt <- t(dt)
  dt <- data.frame(dt)
  names(dt) <- c("x", "l", "m", "u")
  
  y1 <- rnorm(n.sim, u1, s1)
  y2 <- rnorm(n.sim, u2, s2)
  simul.sample <- y1 + log(1 + exp(y2-y1))
  
  m <- mean(simul.sample)
  s <- sqrt(var(simul.sample))

  MC.cdf <- function(x) pnorm(x, m, s)

  data.x <- dt$x
  dt$approx <- MC.cdf(data.x)

  return(dt)
}


plot.analytical.closeness <- function(m1, sd1, m2, sd2, n.sim, n.points=100, n.samples=200){

  dt <- test.analytical.closeness(m1, sd1, m2, sd2, n.sim)
  dt <- melt(foo, measure.vars=c("l", "m", "u", "approx"))
  g <- ggplot(data=dt, aes(x=x)) + geom_line(aes(y=value, color=variable))

  return(g)
}

m1 <- 3
sd1 <- 1
m2 <- 4
sd2 <- 2

foo <- test.analytical.closeness(m1, sd1, m2, sd2, 2000000)
mfoo <- melt(foo, measure.vars=c("l", "m", "u", "approx"))
ggplot(data=mfoo, aes(x=x)) + geom_line(aes(y=value, color=variable))


m1 <- 1
sd1 <- 3
m2 <- 9
sd2 <- 3

foo <- test.analytical.closeness(m1, sd1, m2, sd2, 2000000)
mfoo <- melt(foo, measure.vars=c("l", "m", "u", "approx"))
ggplot(data=mfoo, aes(x=x)) + geom_line(aes(y=value, color=variable))


m1 <- 1
sd1 <- 1
m2 <- 9
sd2 <- 1

foo <- test.analytical.closeness(m1, sd1, m2, sd2, 2000000)
mfoo <- melt(foo, measure.vars=c("l", "m", "u", "approx"))
ggplot(data=mfoo, aes(x=x)) + geom_line(aes(y=value, color=variable))





