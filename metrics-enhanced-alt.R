


# Parameters
param.file <- "params.csv"
base.filepath.out <- "/n/airoldifs2/lab/gbasse/lognormals/res/out"
#base.filepath.out <- "res/out"
number_of_simulations = 4
row.number <- 10
GH.n.nodes <- 20#100
numiter = 200#2000


# Command line arguments
args <- commandArgs(trailingOnly=TRUE)
row.index <- as.integer(args[1])
filepath_out <- paste(base.filepath.out, row.index, ".csv", sep="")



get.convol.density <- function(u1, u2, s1, s2, probs=c(0.0001, 0.9999), n.sim=5000000, n.points=2^13){

  sim <- rlnorm(n.sim, u1, s1) + rlnorm(n.sim, u2, s2)
  lsim <- log(sim)

  qtile <- quantile(lsim, probs=probs)
  
  dens <- density(lsim, from=qtile[1], to=qtile[2], n=n.points)
  dens.fn <- approxfun(dens$x, dens$y, yleft=0, yright=0)

  transf.dens.fn <- function(y) { return(dens.fn(log(y)) * 1 / y) }
  transf.dens.fn <- Vectorize(transf.dens.fn)

  return(transf.dens.fn)

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
  sim <- log(rlnorm(1000000, u1, s1) + rlnorm(1000000,u2, s2))
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
  sim <- log(rlnorm(1000000, u1, s1) + rlnorm(1000000,u2, s2))

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

