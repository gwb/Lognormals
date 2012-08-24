require(rootSolve)

get.WU <- function(u1, u2, s1, s2, t1=0.01, t2=0.05, maxiter=1000){

  initv <- c(max(u1,u2), max(s1,s2))
  
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
  return( multiroot(f, initv, maxiter=maxiter)$root )
  
}
