require(Rmpfr)
require(ggplot2)
require(reshape2)
require(fastGHQuad)


get.integrand <- function(t, u1, u2, s1, s2, log.opts= FALSE){
  fn <- function(y2){
    if(y2 >= t){
      return(0)
    }else if(y2 <= 0){
      return(0)
    }else{
      if(log.opts){
        return(dlnorm(t-y2, u1, s1, log= TRUE) + dlnorm(y2, u2, s2, log= TRUE))
      }else{
        return(exp(dlnorm(t-y2, u1, s1, log= TRUE) + dlnorm(y2, u2, s2, log= TRUE)))
        }
      #return(1/( (t-y2) * sqrt(2*pi*s1^2) * y2 * sqrt(2*pi*s2^2) ) * exp( -(log(t-y2) - u1)^2 / (2*s1^2) - (log(y2) - u2)^2 / (2*s2^2)))
    }
  }
    return(Vectorize(fn))
}

get.first.deriv <- function(t, u1, u2, s1, s2){
    fn <- function(x){
    a <- 1 / ( (t-x) * sqrt(2*pi*s1^2) )
    b <- exp( -(log(t-x) - u1)^2/(2*s1^2) )
    c <- 1 / ( x * sqrt(2*pi*s2^2) )
    d <- exp( -(log(x) - u2)^2 / (2*s2^2) )

    a.p <- 1 / ( (t-x)^2 * sqrt(2*pi*s1^2) )
    b.p <- 1 / (t - x)  *  (log(t-x) - u1) / s1^2 * exp( -(log(t-x) - u1)^2 / (2*s1^2) )
    c.p <- - 1 / ( x^2 * sqrt(2*pi*s2^2) )
    d.p <- - 1 / x * (log(x) - u2) / s2^2 * exp( -(log(x) - u2)^2 / (2*s2^2) )
    return( a.p * b * c * d + a * b.p * c * d + a * b * c.p * d + a * b * c * d.p)
  }
    return(fn)
}


get.mu.hat<- function(t, u1, u2, s1, s2, eps = 0.0000000001){

  fn <- function(x){
    return(1/(t-x) * (s1^2 + log(t-x) - u1) / s1^2 - 1/x * (s2^2 + log(x) - u2) /s2^2)
  }

  a0 = b0 = min(exp(max(u1,u2)), t/2) # prevents a0 and b0 from being bigger than t

  while(fn(a0) <= 0){
    a0 <- max(a0 - 10, a0 / 2) # prevents a0 from being assigned negative values
  }
  #browser()
  while(fn(b0) >= 0){
    b0 <- min(b0 + 10, b0 + (t-b0) / 2)
  }

  
  # use the bisection method
  while( abs( fn( (a0+b0)/2 ) ) > eps ){
    u = (a0 + b0) / 2
    if (fn(u) < 0){
      b0 = u
    }else{
      a0 = u
    }
  }

  return( (a0 + b0) / 2 )
    
  
}



get.mu.hat.alt <- function(t, u1, u2, s1, s2){
  intg <- get.integrand(t, u1, u2, s1, s2, log= TRUE)
  return(optimize(intg, c(0,t), maximum=T)$maximum)
}


get.sigma.hat <- function(x, t, u1, u2, s1, s2){
  a <- 1 / ( (t-x) * sqrt(2*pi*s1^2) )
  b <- exp( -(log(t-x) - u1)^2/(2*s1^2) )
  c <- 1 / ( x * sqrt(2*pi*s2^2) )
  d <- exp( -(log(x) - u2)^2 / (2*s2^2) )

  a.p <- 1 / ( (t-x)^2 * sqrt(2*pi*s1^2) )
  b.p <- 1 / (t - x)  *  (log(t-x) - u1) / s1^2 * exp( -(log(t-x) - u1)^2 / (2*s1^2) )
  c.p <- - 1 / ( x^2 * sqrt(2*pi*s2^2) )
  d.p <- - 1 / x * (log(x) - u2) / s2^2 * exp( -(log(x) - u2)^2 / (2*s2^2) )

  a.p2 <- 2 / ( (t-x)^3 * sqrt(2*pi*s1^2) )
  b.p2 <-  1 / (t-x)^2 * (log(t-x) - u1) / s1^2 * exp( -(log(t-x) - u1)^2 / (2*s1^2) ) +
          - 1 / (t-x)   * 1 / ((t - x) * s1^2)      * exp( -(log(t-x) - u1)^2 / (2*s1^2) ) +
           (1 / (t-x)   *  (log(t-x) - u1) / s1^2)^2 * exp( -(log(t-x) - u1)^2 / (2*s1^2) )
  c.p2 <- 2 / ( x^3 * sqrt(2*pi*s2^2) )
  d.p2 <- - 1 / x^2 * -(log(x) - u2) / s2^2 * exp( -(log(x) - u2)^2 / (2*s2^2) ) +
            1 / x   * - 1 / (x * s2^2)      * exp( -(log(x) - u2)^2 / (2*s2^2) ) +
           (1 / x   *  -(log(x) - u2) / s2^2)^2 * exp( -(log(x) - u2)^2 / (2*s2^2) )

  res <- a.p2* b  * c  * d  + a.p* b.p * c  * d  + a.p* b  * c.p * d  + a.p* b  * c  * d.p +
         a.p * b.p* c  * d  + a  * b.p2* c  * d  + a  * b.p* c.p * d  + a  * b.p* c  * d.p +
         a.p * b  * c.p* d  + a  * b.p * c.p* d  + a  * b  * c.p2* d  + a  * b  * c.p* d.p +
         a.p * b  * c  * d.p+ a  * b.p * c  * d.p+ a  * b  * c.p * d.p+ a  * b  * c  * d.p2

  return( 1/ sqrt(-res))
}


get.convol.density <- function(u1, u2, s1, s2){
  fn <- function(t){
    intg <- get.integrand(t, u1, u2, s1, s2)
    u.hat <- get.mu.hat.alt(t, u1, u2, s1, s2)
    s.hat <- get.sigma.hat(u.hat,t, u1, u2, s1, s2)
    print(s.hat)
    #s.hat <- 1 / sqrt(- hessian(intg, u.hat))
    rule <- gaussHermiteData(50)
    res <- aghQuad(intg, u.hat, s.hat,rule)
    return(res)
  }
  return(fn)
}


# ---- Testing the density
# crappy results with u1=0.1, u2=5, s1=0.05, s2=1)

sim <- rlnorm(3000000, 0.1, 0.05) + rlnorm(3000000, 3, 1)
qtile <- quantile(x=sim, probs=c(0.01, 0.99))
foo <- density(sim, n=2000,from=qtile[1], to=qtile[2])

#bar <- get.convol.density(0.1, 3, 0.05, 1)
bli <- get.convol.density(0.1, 3, 0.05, 1)

#bar.y <- sapply(foo$x, bar)
bli.y <- sapply(foo$x, bli)

dt <- data.frame(convolGH=bli.y, mc=foo$y)
dt <- melt(dt, measure.vars=c("convolGH", "mc"))
dt$x <- c(foo$x, foo$x)

ggplot(data=dt, aes(x=x)) + geom_line(aes(y=value, color=variable))
ggplot(data=data.frame(y=bar.y, x=foo$x), aes(x=x,y=y)) + geom_line()


# ---- Testing the density 2.0
# check, for this set of parameters:
# - the mode
# - the second deriv


sim <- rlnorm(3000000, 3, 1) + rlnorm(3000000, 4, 2)
qtile <- quantile(x=sim, probs=c(0.01, 0.99))
foo <- density(sim, n=2000,from=qtile[1], to=qtile[2])

#bar <- get.convol.density(0.1, 3, 0.05, 1)
bli <- get.convol.density(3, 4, 1, 2)

#bar.y <- sapply(foo$x, bar)
bli.y <- sapply(foo$x, bli)

dt <- data.frame(convolGH=bli.y, mc=foo$y)
dt <- melt(dt, measure.vars=c("convolGH", "mc"))
dt$x <- c(foo$x, foo$x)

g1 <- ggplot(data=dt, aes(x=x)) + geom_line(aes(y=value, color=variable))
ggsave(plot=g1, filename="old_method.pdf", width=17, height=10)


sapply(foo$x, get.mu.hat, u1, u2, s1 ,s2)
sapply(foo$x, get.mu.hat.alt, u1, u2, s1 ,s2)


run.mcmc <- function(numiter, x0, draw.proposal, eval.proposal, eval.target){

  rejections <- 0
  xs <- NULL
  
  x.tm1 <- x0
  
  it <- 1
  while( it < numiter ){

    x.star <- draw.proposal()
    log.ratio <- log(eval.target(x.star)) - log(eval.proposal(x.star)) - ( log(eval.target(x.tm1)) - log(eval.proposal(x.tm1)))

    if(log(runif(1)) < min(log.ratio,0)){
      x.tm1 <- x.star
    }else{
      rejections <- rejections + 1
    }
    xs <- c(xs, x.tm1)
    it <- it + 1
  }
  return(list(xs, rejections))
}



u1 <- 0.1
s1 <- 0.05
u2 <- 4
s2 <- 1

d.prop <- function() rlnorm(1,4,1)
e.prop <- function(x) dlnorm(x, 4, 1)
e.target <- function(x) dlnorm(x,5,1)

res <- run.mcmc(100, 67, d.prop, e.prop, e.target)


test.mcmc <- function(u1, u2, s1, s2, x0=NULL, numiter=100){

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


u1 <- 0.1
s1 <- 0.05
u2 <- 4
s2 <- 1

foo <- test.mcmc(u1, u2, s1, s2, numiter=300)
foo <- test.mcmc(3, 3, 1, 1, numiter=100)
