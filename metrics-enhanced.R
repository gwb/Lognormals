require(numDeriv)
require(fastGHQuad)


# Parameters
param.file <- "params.csv"
base.filepath.out <- "/n/airoldifs2/lab/gbasse/lognormals/res/out"
number_of_simulations = 4
row.number <- 10
GH.n.nodes <- 20#100
numiter = 200#2000


# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
row.index <- as.integer(args[1])
filepath_out <- paste(base.filepath.out, row.index, ".csv", sep="")


get.log.integrand <- function(t, u1, u2, s1, s2){

  fn <- function(x){
    return(dlnorm(t-exp(x), u1, s1) * dlnorm(exp(x), u2, s2) * exp(x))
  }
  return(fn)
}

get.log.first.deriv <- function(t, u1, u2, s1, s2){
  fn <- function(x){
    a <- 1 / ( (t - exp(x)) * sqrt(2*pi*s1^2) )
    b <- exp( -( log(t-exp(x)) - u1 )^2 / (2*s1^2) )
    c <- 1 / sqrt(2*pi*s2^2) * exp( -(x-u2)^2 / (2*s2^2) )
    
    a.p <-  exp(x) / ( (t-exp(x))^2 * sqrt(2*pi*s1^2) )
    b.p <- (log(t - exp(x)) - u1) / ( (t*exp(-x) - 1) * s1^2 ) * b
    c.p <- - (x-u2) / s2^2 * c

    res <- a.p * b * c + a * b.p * c + a * b * c.p
    
    return(res)
  }
  return(fn)
}

get.log.first.deriv.num <- function(t, u1, u2, s1, s2){

  intg <- get.log.integrand(t, u1, u2, s1, s2)
  fn <- function(x){
    return(grad(intg, x))
  }
  return(fn)
}

get.log.second.deriv <- function(t, u1, u2, s1, s2){

  fn <- function(x){
    a <- 1 / ( (t - exp(x)) * sqrt(2*pi*s1^2) )
    b <- exp( -( log(t-exp(x)) - u1 )^2 / (2*s1^2) )
    c <- 1 / sqrt(2*pi*s2^2) * exp( -(x-u2)^2 / (2*s2^2) )
    
    a.p <-  exp(x) / ( (t-exp(x))^2 * sqrt(2*pi*s1^2) )
    b.p <- (log(t - exp(x)) - u1) / ( (t*exp(-x) - 1) * s1^2 ) * b
    c.p <- - (x-u2) / s2^2 * c

    a.pp <- (  exp(x) * (t-exp(x)) + 2*exp(2*x) ) / ( (t-exp(x))^3 * sqrt(2*pi*s1^2) )
    b.pp <- (  -exp(x) / (t - exp(x)) * 1 / ( (t*exp(-x) - 1) * s1^2 )  + (log(t-exp(x)) - u1) * t*exp(-x) / ( (t*exp(-x) - 1)^2 * s1^2 ) + ( (log(t-exp(x)) - u1)/( (t*exp(-x) -1) * s1^2 ) )^2  ) * b
    c.pp <- (  -1 / ( s2^2 * sqrt(2*pi*s2^2) )  + 1 / sqrt(2*pi*s2^2) * ( - (x-u2)/s2^2 )^2  ) * exp( -(x-u2)^2 / (2*s2^2) )
    
    res <- a.pp* b * c + a.p * b.p * c + a.p * b * c.p +
      a.p * b.p * c + a * b.pp * c + a * b.p * c.p +
        a.p * b * c.p + a * b.p * c.p + a * b * c.pp
    return(res)
  }
  return(fn)
}

get.log.second.deriv.num <- function(t, u1, u2, s1, s2){

  intg <- get.log.integrand(t, u1, u2, s1, s2)
  fn <- function(x){
    return(hessian(intg, x))
  }
  return(fn)
}


get.mu.hat<- function(t, u1, u2, s1, s2, eps = 0.0000001){

  fn <- function(x){
    return( exp(x)/( t-exp(x) )  + exp(x) / (t - exp(x)) * ( log(t-exp(x)) - u1 ) / s1^2 - (x - u2) / s2^2 - 1)
  }

  a0 = b0 = ifelse(log(t) >0, log(t) / 2, log(t) * 3/2) #min(max(u1,u2), abs(log(t))/2) # prevents a0 and b0 from being bigger than t

  if (is.na(fn(a0))){
    print(a0)
    print(t)
    print(u1)
    print(u2)
    print(s1)
    print(s2)
  }

  #browser()
  while(fn(a0) <= 0){
    #a0 <- max(a0 - 10, a0 / 2)
    tmp <- a0
    a0 <- a0 - abs(log(t)) / 2
    count <- 1
    while(is.na(fn(a0))){
      a0 <- tmp - abs(log(t)) / i
      i = i+1
    }
    #if (is.na(fn(a0))){print(a0)}
  }
  #browser()
  while(fn(b0) >= 0){
    b0 <- min(b0 + 10, b0 + (log(t)-b0) / 2)
  }

  
  # use the bisection method
  #while( abs( fn( (a0+b0)/2 ) ) > eps ){
  #  u = (a0 + b0) / 2
  #  if (fn(u) < 0){
  #    b0 = u
  #  }else{
  #    a0 = u
  #  }
  #}

  intg <- get.log.integrand(t, u1, u2, s1, s2)
  return(optimize(intg, c(a0, b0), maximum=T)$maximum)
  #return( (a0 + b0) / 2 )
    
  
}

# get.mu.hat.alt(0.1915875, -3.656631, -2.063464, 1.145976, 0.3562696) => -0.191
# which corresponds to -t, the inferior bound. So it's clearly wrong

#get.mu.hat.alt <- function(t, u1, u2, s1, s2){
#  intg <- get.log.integrand(t, u1, u2, s1, s2)
#  return(optimize(intg, c(-t,log(t)), maximum=T)$maximum)
#}

get.convol.density <- function(u1, u2, s1, s2){
  fn <- function(t){
    intg <- get.log.integrand(t, u1, u2, s1, s2)
    u.hat <- get.mu.hat(t, u1, u2, s1, s2)
    s.hat <- 1/sqrt(-get.log.second.deriv(t, u1, u2, s1, s2)(u.hat))
    #s.hat <- 1 / sqrt(- hessian(intg, u.hat))
    rule <- gaussHermiteData(GH.n.nodes)
    res <- aghQuad(intg, u.hat, s.hat,rule)

    # note: as far as I can see, this only happens
    # when res << 1
    if(is.na(res)){
      res = 0
    }
    return(res)
  }
  return(fn)
}



run.mcmc <- function(numiter, x0, draw.proposal, eval.proposal, eval.target){

  rejections <- 0
  xs <- NULL
  
  x.tm1 <- x0
  
  it <- 1
  while( it < numiter ){
    x.star <- draw.proposal()
    log.ratio <- log(eval.target(x.star)) - log(eval.proposal(x.star)) - ( log(eval.target(x.tm1)) - log(eval.proposal(x.tm1)))

    if (is.na(log.ratio)){
      #it = it+1
      #next
      return(list("error", NA))
    }
    
    if(log(runif(1)) < min(log.ratio,0)){
      x.tm1 <- x.star
    }else{
      rejections <- rejections + 1
    }
    xs <- c(xs, x.tm1)
    it <- it + 1
  }
  return(list(xs, rejections/numiter))
}



test.mcmc <- function(m1, m2, sd1, sd2, x0=NULL, numiter=100){

  v1 <- log(exp(sd1)^2)
  v2 <- log(exp(sd2)^2)
  
  s1 <- sqrt(log(1 + exp(v1-2*m1)))
  u1 <- m1 - (s1^2)/2

  s2 <- sqrt(log(1 + exp(v2-2*m2)))
  u2 <- m2 - (s2^2)/2
  
  # Get MOM parameters
  sim <- log(rlnorm(100000, u1, s1) + rlnorm(100000,u2, s2))
  n.u <- mean(sim)
  n.s <- sqrt(var(sim))

  # get eval and proposal functions
  e.target <- get.convol.density(u1, u2, s1, s2)
  e.proposal <- function(x) dlnorm(x, n.u, n.s)
  d.proposal <- function() rlnorm(1, n.u, n.s)

  if (length(x0) == 0){
    x0 <- d.proposal()
  }
  
  res <- run.mcmc(numiter, x0, d.proposal, e.proposal, e.target)

  return(res)
}



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

  # transforming params
  v1 <- log(exp(sd1)^2)
  v2 <- log(exp(sd2)^2)
  
  s1 <- sqrt(log(1 + exp(v1-2*m1)))
  u1 <- m1 - (s1^2)/2

  s2 <- sqrt(log(1 + exp(v2-2*m2)))
  u2 <- m2 - (s2^2)/2
  
  # simulates from convolution
  sim <- log(rlnorm(100000, u1, s1) + rlnorm(100000,u2, s2))

  # extract quantiles 
  qtile = quantile(sim, probs=seq(0.01, 0.99, length.out=number_of_simulations))
  
  tmp.distance.vect = NULL
  for (i in seq(1, number_of_simulations)){
    distance <- test.mcmc(m1, m2, sd1, sd2, x0=exp(qtile[i]), numiter=numiter)[[2]]
    tmp.distance.vect <- c(tmp.distance.vect, distance)
  }

  distance = mean(tmp.distance.vect, na.rm=T)

  if (is.na(distance)){ next }
  
  m1.list <- c(m1.list, m1)
  m2.list <- c(m2.list, m2)
  sd1.list <- c(sd1.list, sd1)
  sd2.list <- c(sd2.list, sd2)

  distance.list <- c(distance.list, distance)
}
  
dt.out <- data.frame(m1=m1.list, m2=m2.list, sd1=sd1.list, sd2=sd2.list, distance=distance.list)

#print(dt.out)

write.csv(dt.out, filepath_out)



#res <- test.mcmc(-1, 1, 0, -2)
#print(res[[2]])
