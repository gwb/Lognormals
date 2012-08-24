

get.FW <- function(u1, u2, s1, s2){
  sq1 <- s1^2
  sq2 <- s2^2
  us <- c(u1,u2)
  sqs <- c(sq1,sq2)
  n.u1 <- sum(exp(us + sqs/2))
  t1 <- sum(exp(2*us + 2*sqs))
  t2 <- 2 * exp(u1 + u2) * exp((sq1 + sq2) / 2)
  n.u2 <- t1 + t2

  return( c(2 * log(n.u1) - 0.5 * log(n.u2), log(n.u2) - 2 * log(n.u1)) )
  
}
