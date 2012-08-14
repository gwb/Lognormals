require(ggplot2)
require(reshape2)
require(grid)
require(PerformanceAnalytics)
source('/Users/gwb/Hacks/Projects/myRUtils/graphics.R')



# ----- The convo matching only the first and second moments
get.convol.emp.CDFs <- function(u1,s1, u2, s2, nl){
  
# Defines sample
y1 <- rnorm(nl, u1, s1)
y2 <- rnorm(nl, u2, s2)
simul.sample = y1 + log(1 + exp(y2 - y1))

# Obtains the quantiles for the samples. This allows us to make the density
# estimation more robust to sampling variations
qtile <- quantile(x=simul.sample, probs=c(0.01, 0.99))

# Creates MC lognormal approx
m <- mean(simul.sample)
s <- sqrt(var(simul.sample))

MC.cdf <- function(x) pnorm(x, m, s)

# Creates CDFs data
data.x <- seq(qtile[1], qtile[2], length.out=2000)
th.emp.cdf <- ecdf(simul.sample)(data.x)
approx.cdf <- MC.cdf(data.x)


# Plots the results
dt <- data.frame(th.emp.cdf = th.emp.cdf, approx.cdf = approx.cdf)
m.dt <- melt(dt, measure.vars=c("th.emp.cdf", "approx.cdf"))
m.dt$x <- c(data.x, data.x)
g <- ggplot(data=m.dt, aes(x=x)) + geom_line(aes(y=value, color=variable)) +
  custom_opts() +
    xlab("x") +
      ylab("F(x)") +
        scale_color_hue("CDFs", labels=c("Truth (empirical)   ", "Approximation   ")) +
            opts(title=paste("Parameters: u1 = ", u1, ", s1 = ", s1, ", u2 = ", u2, ", s2 = ", s2, sep=""))
         

return(g)
}


## equal low sds, varying means
#g1 <- get.convol.emp.CDFs(2, 1, 5, 1, 1000000)
#g2 <- get.convol.emp.CDFs(5, 1, 5, 1, 1000000)
#g3 <- get.convol.emp.CDFs(5, 1, 10, 1, 1000000)
#g4 <- get.convol.emp.CDFs(10, 1, 10, 1, 1000000)
#
#stack_plots(list(g1,g2,g3,g4), 2,2, add_opts=function(x) opts(legend.position="bottom"))
#
## equal means, varying sds
#g1 <- get.convol.emp.CDFs(5, 1, 5, 1, 1000000)
#g2 <- get.convol.emp.CDFs(5, 3, 5, 1, 1000000)
#g3 <- get.convol.emp.CDFs(5, 5, 5, 1, 1000000)
#g4 <- get.convol.emp.CDFs(5, 2, 5, 2, 1000000)
#g5 <- get.convol.emp.CDFs(5, 3, 5, 3, 1000000)
#g6 <- get.convol.emp.CDFs(5, 4, 5, 4, 1000000)
#
#stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3, add_opts=function(x) opts(legend.position="bottom"))
#
## one mean low, one high, varying sds
#g1 <- get.convol.emp.CDFs(5, 1, 10, 1, 1000000)
#g2 <- get.convol.emp.CDFs(5, 3, 10, 1, 1000000)
#g3 <- get.convol.emp.CDFs(5, 5, 10, 1, 1000000)
#g4 <- get.convol.emp.CDFs(5, 1, 10, 2, 1000000)
#g5 <- get.convol.emp.CDFs(5, 1, 10, 3, 1000000)
#g6 <- get.convol.emp.CDFs(5, 1, 10, 6, 1000000)
#g7 <- get.convol.emp.CDFs(5, 2, 10, 2, 1000000)
#g8 <- get.convol.emp.CDFs(5, 3, 10, 3, 1000000)
#g9 <- get.convol.emp.CDFs(5, 5, 10, 5, 1000000)
#
#stack_plots(list(g1,g2,g3,g4,g5,g6,g7,g8,g9), 3,3, add_opts=function(x) opts(legend.position="bottom"))



# ----- The convo matching the 3rd and 4th moments


double_factorial <- function(n){
  a <- 2*seq(1,n) -1
  return(prod(a[a<=n]))
}


get.convol.emp.all.CDFs <- function(u1,s1, u2, s2, nl){

# Functions to convert u and s in log EY and Var EY
get.SDY <- function (u,s) return( log( sqrt( (exp(s^2) - 1) * exp(2 * u + (s^2)) ) ) )
get.EY <- function (u,s) return(u + s^2/2)

m1 <- get.EY(u1,s1)
sd1 <- get.SDY(u1,s1)
m2 <- get.EY(u2,s2)
sd2 <- get.SDY(u2,s2)


# Defines sample
y1 <- rnorm(nl, u1, s1)
y2 <- rnorm(nl, u2, s2)
simul.sample = y1 + log(1 + exp(y2 - y1))

# Obtains the quantiles for the samples. This allows us to make the density
# estimation more robust to sampling variations
qtile <- quantile(x=simul.sample, probs=c(0.01, 0.99))


# Computes the moments used for matching
# -- the classical mean and std
m <- mean(simul.sample)
s <- sqrt(var(simul.sample))

# -- the std obtained from the 4th central moment
#    for the normal approximation
s.4 <- (mean((simul.sample-m)^4) / double_factorial(3))^(1/4)

# -- the moments necessary for the t - approximation
V <- var(simul.sample)
mu <- mean(simul.sample)
K <- kurtosis(simul.sample-mu, method='excess')
sq <- (6 + 2*K)/(6+4*K) * V
nu <- 6/K + 4


# Generates the approximations
MC.reg.normal.cdf <- function(x) pnorm(x, m, s)
MC.tail.normal.cdf <- function(x) pnorm(x, m, s.4)
MC.t.cdf <- function(x) pt( (x-mu) / sqrt(sq) , nu)

# Creates CDFs data
data.x <- seq(qtile[1], qtile[2], length.out=2000)
th.emp.cdf <- ecdf(simul.sample)(data.x)
approx.reg.normal.cdf <- MC.reg.normal.cdf(data.x)
approx.tail.normal.cdf <- MC.tail.normal.cdf(data.x)
approx.t.cdf <- MC.t.cdf(data.x)

# Plots the results
dt <- data.frame(th.emp.cdf = th.emp.cdf, approx.reg.normal.cdf = approx.reg.normal.cdf, approx.tail.normal.cdf = approx.tail.normal.cdf, approx.t.cdf = approx.t.cdf)
m.dt <- melt(dt, measure.vars=c("th.emp.cdf", "approx.reg.normal.cdf", "approx.tail.normal.cdf", "approx.t.cdf"))
m.dt$x <- c(data.x, data.x, data.x, data.x)
g <- ggplot(data=m.dt, aes(x=x)) + geom_line(aes(y=value, color=variable)) +
  custom_opts() +
    xlab("x") +
      ylab("F(x)") +
        scale_color_hue("CDFs", labels=c("Truth (empirical)   ", "Normal (head)   ", "Normal (tail)   ", "t   ")) +
            opts(title=paste("Parameters: u1 = ", u1, ", s1 = ", s1, ", u2 = ", u2, ", s2 = ", s2,
                           "\n            m1 = ", round(m1,2), ", sd1 = ",round(sd1,2), ", m2 = ", round(m2,2), ", sd2 = ", round(sd2,2),
                   sep=""))
         

return(g)
}


get.convol.unlog.CDFs <- function(m1,sd1, m2, sd2, nl){

# Here, m and v are given as log(E[Y]) and log(Var[Y]) instead
# of u and s directly, so we need to explicitely convert them:

v1 <- log(exp(sd1)^2)
v2 <- log(exp(sd2)^2)
  
s1 <- sqrt(log(1 + exp(v1-2*m1)))
u1 <- m1 - (s1^2)/2

s2 <- sqrt(log(1 + exp(v2-2*m2)))
u2 <- m2 - (s2^2)/2

  
# Defines sample
y1 <- rnorm(nl, u1, s1)
y2 <- rnorm(nl, u2, s2)
simul.sample = y1 + log(1 + exp(y2 - y1))

# Obtains the quantiles for the samples. This allows us to make the density
# estimation more robust to sampling variations
qtile <- quantile(x=simul.sample, probs=c(0.01, 0.99))


# Computes the moments used for matching
# -- the classical mean and std
m <- mean(simul.sample)
s <- sqrt(var(simul.sample))

# -- the std obtained from the 4th central moment
#    for the normal approximation
s.4 <- (mean((simul.sample-m)^4) / double_factorial(3))^(1/4)

# -- the moments necessary for the t - approximation
V <- var(simul.sample)
mu <- mean(simul.sample)
K <- kurtosis(simul.sample-mu, method='excess')
sq <- (6 + 2*K)/(6+4*K) * V
nu <- 6/K + 4


# Generates the approximations
MC.reg.normal.cdf <- function(x) pnorm(x, m, s)
MC.tail.normal.cdf <- function(x) pnorm(x, m, s.4)
MC.t.cdf <- function(x) pt( (x-mu) / sqrt(sq) , nu)

# Creates CDFs data
data.x <- seq(qtile[1], qtile[2], length.out=2000)
th.emp.cdf <- ecdf(simul.sample)(data.x)
approx.reg.normal.cdf <- MC.reg.normal.cdf(data.x)
approx.tail.normal.cdf <- MC.tail.normal.cdf(data.x)
approx.t.cdf <- MC.t.cdf(data.x)

# Plots the results
dt <- data.frame(th.emp.cdf = th.emp.cdf, approx.reg.normal.cdf = approx.reg.normal.cdf, approx.tail.normal.cdf = approx.tail.normal.cdf, approx.t.cdf = approx.t.cdf)
m.dt <- melt(dt, measure.vars=c("th.emp.cdf", "approx.reg.normal.cdf", "approx.tail.normal.cdf", "approx.t.cdf"))
m.dt$x <- c(data.x, data.x, data.x, data.x)
g <- ggplot(data=m.dt, aes(x=x)) + geom_line(aes(y=value, color=variable)) +
  custom_opts() +
    xlab("x") +
      ylab("F(x)") +
        scale_color_hue("CDFs", labels=c("Truth (empirical)   ", "Normal (head)   ", "Normal (tail)   ", "t   ")) +
            opts(title=paste("Parameters: m1 = ", m1, ",    sd1 = ", sd1, ",    m2 = ", m2, ",    sd2 = ", sd2,
                         "\n              u1 = ", round(u1,2), ", s1 = ", round(s1,2), ", u2 = ", round(u2,2), ", s2 = ", round(s2,2),
                   sep=""))
         

return(g)
}

get.EY.hm <- function(min.u, max.u, step.u, min.s, max.s, step.s){

get.EY <- function (u,s) return(u + s^2/2)

param.space <- expand.grid(u=seq(min.u, max.u, step.u), s=seq(min.s, max.s, step.s))
EY <- apply(param.space, 1, function(x) get.EY(x[1], x[2]))

param.space <- data.frame(param.space)
param.space$EY <- EY

g <- ggplot(data=param.space, aes(x=u, y=s, fill=EY)) + geom_tile() + opts(title = "Log E[Y]")  + guides(fill=guide_legend(title="Log E[Y]"))

return(g)  
}

get.SDY.hm <- function(min.u, max.u, step.u, min.s, max.s, step.s){

get.SDY <- function (u,s) return( log( sqrt( (exp(s^2) - 1) * exp(2 * u + (s^2)) ) ) )

param.space <- expand.grid(u=seq(min.u, max.u, step.u), s=seq(min.s, max.s, step.s))
SDY <- apply(param.space, 1, function(x) get.SDY(x[1], x[2]))

param.space <- data.frame(param.space)
param.space$SDY <- SDY

g <- ggplot(data=param.space, aes(x=u, y=s, fill=SDY)) + geom_tile() + opts(title = "Log SD[Y]") + guides(fill=guide_legend(title="Log SD[Y]"))

return(g)  
}
