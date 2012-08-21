require(rootSolve)

get.FW <- function(u1, u2, s1, s2){
  sq1 <- s1^2
  sq2 <- s2^2
  us <- c(u1,u2)
  sqs <- c(sq1,sq2)
  n.u1 <- sum(exp(us + sqs/2))
  t1 <- sum(exp(2*us + 2*sqs))
  t2 <- 2 * exp(u1 + u2) * exp((sq1 + sq2) / 2)
  n.u2 <- t1 + t2

  print(n.u1)
  print(n.u2)
  
  return( c(2 * log(n.u1) - 0.5 * log(n.u2), log(n.u2) - 2 * log(n.u1)) )
  
}



get.SY <- function(u1, u2, s1, s2){

  # AG's checked with python versions: same results
  AG1 = c(-0.1153,
       -0.5667,
       0.1151,
       -0.71624,
       0.13129,
       -0.16113,
       0.20842,
       -0.44997,
       0.27562,
       -0.51097,
       0.13451,
       -0.26701,
       0.59191,
       -0.36953,
       0.69127,
       0.80570,
       0.14297,
       -0.32259,
       0.20556,
       -0.38886,
       -0.31453,
       -0.27300,
       0.62442,
       -0.40482,
       0.77467)

AG1P = c(0,1,2,1,1,0,2,2,2,1,0,2,2,2,1,-1,2,2,2,1,-1,1,1,1,0)

AG2 = c(0.40128,
       -0.44832,
       0.73917,
       -0.37721,
       0.52622,
       -0.15791,
       0.16372,
       -0.29949,
       0.13378,
       -0.17474,
       0.21745,
       -0.21141,
       0.39767,
       -0.17768,
       0.21953,
       -0.75253,
       0.11471,
       -0.21895,
       0.10067,
       -0.12249,
       0.84479,
       -0.22228,
       0.42895,
       -0.20357,
       0.24995)

AG2P = c(-1,1,1,1,0,1,2,2,2,1,1,2,2,2,1,0,2,2,2,1,-1,1,1,1,0)

AG3 = c(-0.39586,
       -0.54549,
       0.11392,
       -0.71163,
       0.13152,
       0.68399,
       0.19625,
       -0.43140,
       0.26477,
       -0.49405,
       -0.47172,
       -0.24868,
       0.56185,
       -0.35018,
       0.65609,
       0.18193,
       0.13235,
       -0.30500,
       0.19377,
       -0.36583,
       -0.28174,
       -0.25182,
       0.58933,
       -0.38094,
       0.72601)

AG3P = c(1,1,2,1,1,1,2,2,2,1,1,2,2,2,1,1,2,2,2,1,0,1,1,1,0)
  
  # All the Gis have been checked against the python versions, same results
  G1.A <- AG1 * 10^AG1P
  G2.A <- AG2 * 10^AG2P
  G3.A <- AG3 * 10^AG3P

  G1 <- function(m, ssq){
    is <- rep(seq(0,4), each=5)
    js <- rep(seq(0,4), 5)
    return( 10^sum((G1.A * sqrt(ssq) ^ (is/2) * abs(m) ^ (js/2))) )
  }

  G2 <- function(m, ssq){
    is <- rep(seq(0,4), each=5)
    js <- rep(seq(0,4), 5)
    return( 10^sum((G2.A * sqrt(ssq) ^ (is/2) * abs(m) ^ (js/2))) )
  }

  G3 <- function(m, ssq){
    is <- rep(seq(0,4), each=5)
    js <- rep(seq(0,4), 5)
    return( 10^sum((G3.A * sqrt(ssq) ^ (is/2) * abs(m) ^ (js/2))) )
  }

  if ( u2 > u1 ){
    tmp = u1
    u1 = u2
    u2 = tmp

    tmp = s1
    s1 = s2
    s2 = tmp
  }

  ssq1 <- s1^2
  ssq2 <- s2^2
  
  mw <- u2 - u1
  ssqw <- ssq1 + ssq2
  
  VZ = ssq1 - G1(mw, ssqw) ^ 2 - 2 * ssq1 / ssqw * G3(mw, ssqw) + G2(mw, ssqw)
  EZ = u1 + G1(mw, ssqw)

  return(c(EZ, sqrt(VZ)))
}


# Tested against Python implementation, works
get.WU <- function(u1, u2, s1, s2, t1, t2, initv=c(0,1), maxiter=1000){
  get.npsi <- function(s){
    xs = c(0.24534,
          0.73747,
          1.23407,
          1.73853,
          2.25497,
          2.78880,
          3.34785,
          3.94476,
          4.60369,
          5.38748)
    ws = c(4.62243 * 10^(-1),
          2.86675 * 10^(-1),
          1.09017 * 10^(-1),
          2.48105 * 10^(-2),
          3.24377 * 10^(-3),
          2.28338 * 10^(-4),
          7.80255 * 10^(-6),
          1.08606 * 10^(-7),
          4.39934 * 10^(-10),
          2.22939 * 10^(-13))
    xs <- c(xs, -xs)
    ws <- c(ws, ws)
    fn <- function(mu, ssq){
      return(sum( ws / sqrt(pi) * exp( -s * exp( sqrt( 2*ssq ) * xs + mu ) ) ))
    }
    return(fn)
  }
  
  get.nfit <- function(nu1, nu2, ns1, ns2, nt1, nt2){
    mgfi1 <- get.npsi(nt1)
    mgfi2 <- get.npsi(nt2)

    mgfi1s <- mgfi1(nu1, ns1^2) * mgfi1(nu2, ns2^2)
    mgfi2s <- mgfi2(nu1, ns1^2) * mgfi2(nu2, ns2^2)

    fn <- function(x){
      mu <- x[1]
      s <- x[2]
      ssq <- s^2
      return(c(mgfi1(mu, ssq) - mgfi1s, mgfi2(mu, ssq) - mgfi2s))
    }
    return(fn)
  }

  f <- get.nfit(u1, u2, s1, s2, t1, t2)
  return( multiroot(f, initv, maxiter=maxiter) )
  
}


test.FW <- function(u1, u2, s1, s2){
  res <- get.FW(u1, u2, s1, s2)
  n.u <- res[1]
  n.s <- sqrt(res[2])

  sim <- rlnorm(3000000, u1, s1) + rlnorm(3000000, u2, s2)
  m <- log(mean(sim))
  v <- log(var(sim))

  s <- sqrt(log(1 + exp(v-2*m)))
  u <- m - (s^2)/2

  return(c(n.u, u, n.s, s))
}

                        
