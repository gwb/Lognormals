# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# The script is an ad-hoc method for measuring the distance,
# or rather the acceptability of the distance between two cdf.
# More specifically, the algorithm computes the empirical cdf
# using a specified number of points, and repeats the simulation
# many times in order to obtain confidence bounds around the cdf.
# Them, the cdf of the approximation is computed, and we count the
# proportion of points for which the approximation is NOT within
# the bounds.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Parameters
param.file <- "params-reversed-3.csv"
base.filepath.out <- "res/out"
row.number <- 10

# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
row.index <- as.integer(args[1])
filepath_out <- paste(base.filepath.out, row.index, ".csv", sep="")




get.convo.emp.cdf <- function(m1, sd1, m2, sd2, n.points, n.sample){

  v1 <- log(exp(sd1)^2)
  v2 <- log(exp(sd2)^2)
  
  s1 <- sqrt(log1p(exp(v1-2*m1)))
  u1 <- m1 - (s1^2)/2

  s2 <- sqrt(log1p(exp(v2-2*m2)))
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
  
  s1 <- sqrt(log1p(exp(v1-2*m1)))
  u1 <- m1 - (s1^2)/2

  s2 <- sqrt(log1p(exp(v2-2*m2)))
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


# Note that the following could be written much more compactly using a data.table,
# but given the recent segfault issues I had using data.tables, I don't feel like
# spending hours debugging it on the cluster, so I stick to a good ol' loop. Plus
# with ten iterations, it's not gonna be critically slow..

params <- read.csv(param.file)

m1.list <- NULL
m2.list <- NULL
sd1.list <- NULL
sd2.list <- NULL
distance.list <- NULL


# takes care of the edge cases
if ( (row.index - 1) * row.number + row.number > dim(params)[1] ){
  row.count <- row.number - (row.index - 1) * row.number - row.number + dim(params)[1]
}else{
  row.count <- row.number
}


for (i in seq(1, row.count)){
  print(i)
  
  ind <- (row.index - 1) * row.number + i
  
  m1 <- params$m1[ind]
  m2 <- params$m2[ind]
  sd1 <- params$sd1[ind]
  sd2 <- params$sd2[ind]
  distance <- test.cdf.closeness(m1, sd1, m2, sd2, 2000000, n.points=200)

  m1.list <- c(m1.list, m1)
  m2.list <- c(m2.list, m2)
  sd1.list <- c(sd1.list, sd1)
  sd2.list <- c(sd2.list, sd2)

  distance.list <- c(distance.list, distance)
}
  
dt.out <- data.frame(m1=m1.list, m2=m2.list, sd1=sd1.list, sd2=sd2.list, distance=distance.list)

write.csv(dt.out, filepath_out)

