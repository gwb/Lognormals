
args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]
row.index <- as.integer(args[2])
filepath_out <- args[3]

params <- read.csv(filename)


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
  data.x <- seq(qtile[1], qtile[2], length.out=2000)

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


get.cdf.data <- function(m1, sd1, m2, sd2, n.sim, n.points=100, n.samples=200){

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


test.cdf.closeness <- function(m1, sd1, m2, sd2, n.sim, n.points=100, n.samples=200){
  dt <- get.cdf.data(m1, sd1, m2, sd2, n.sim)
  res <- with(dt, (approx < l) | (u < approx))
  return(sum(res)/length(res))
}


m1 <- params$m1[row.index]
m2 <- params$m2[row.index]
sd1 <- params$sd1[row.index]
sd2 <- params$sd2[row.index]

res <- test.cdf.closeness(m1, sd1, m2, sd2, 2000000, n.points=200)

dt.out <- data.frame(m1=m1, m2=m2, sd1=sd1, sd2=sd2, distance=res)

write.csv(dt.out, filepath_out)

