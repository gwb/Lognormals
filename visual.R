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


# equal low sds, varying means
g1 <- get.convol.emp.CDFs(2, 1, 5, 1, 1000000)
g2 <- get.convol.emp.CDFs(5, 1, 5, 1, 1000000)
g3 <- get.convol.emp.CDFs(5, 1, 10, 1, 1000000)
g4 <- get.convol.emp.CDFs(10, 1, 10, 1, 1000000)

stack_plots(list(g1,g2,g3,g4), 2,2, add_opts=function(x) opts(legend.position="bottom"))

# equal means, varying sds
g1 <- get.convol.emp.CDFs(5, 1, 5, 1, 1000000)
g2 <- get.convol.emp.CDFs(5, 3, 5, 1, 1000000)
g3 <- get.convol.emp.CDFs(5, 5, 5, 1, 1000000)
g4 <- get.convol.emp.CDFs(5, 2, 5, 2, 1000000)
g5 <- get.convol.emp.CDFs(5, 3, 5, 3, 1000000)
g6 <- get.convol.emp.CDFs(5, 4, 5, 4, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3, add_opts=function(x) opts(legend.position="bottom"))

# one mean low, one high, varying sds
g1 <- get.convol.emp.CDFs(5, 1, 10, 1, 1000000)
g2 <- get.convol.emp.CDFs(5, 3, 10, 1, 1000000)
g3 <- get.convol.emp.CDFs(5, 5, 10, 1, 1000000)
g4 <- get.convol.emp.CDFs(5, 1, 10, 2, 1000000)
g5 <- get.convol.emp.CDFs(5, 1, 10, 3, 1000000)
g6 <- get.convol.emp.CDFs(5, 1, 10, 6, 1000000)
g7 <- get.convol.emp.CDFs(5, 2, 10, 2, 1000000)
g8 <- get.convol.emp.CDFs(5, 3, 10, 3, 1000000)
g9 <- get.convol.emp.CDFs(5, 5, 10, 5, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6,g7,g8,g9), 3,3, add_opts=function(x) opts(legend.position="bottom"))



# ----- The convo matching the 3rd and 4th moments


double_factorial <- function(n){
  a <- 2*seq(1,n) -1
  return(prod(a[a<=n]))
}


get.convol.emp.all.CDFs <- function(u1,s1, u2, s2, nl){

# Temporary, just to debug
#u1 = 5
#s1 = 1
#u2 = 10
#s2 = 2
#nl = 100000
  
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
            opts(title=paste("Parameters: u1 = ", u1, ", s1 = ", s1, ", u2 = ", u2, ", s2 = ", s2, sep=""))
         

return(g)
}


# equal low sds, varying means
g1 <- get.convol.emp.all.CDFs(2, 1, 5, 1, 1000000)
g2 <- get.convol.emp.all.CDFs(5, 1, 5, 1, 1000000)
g3 <- get.convol.emp.all.CDFs(5, 1, 10, 1, 1000000)
g4 <- get.convol.emp.all.CDFs(10, 1, 10, 1, 1000000)

stack_plots(list(g1,g2,g3,g4), 2,2, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="equal_low_std.pdf")

# equal means, varying sds
g1 <- get.convol.emp.all.CDFs(5, 1, 5, 1, 1000000)
g2 <- get.convol.emp.all.CDFs(5, 3, 5, 1, 1000000)
g3 <- get.convol.emp.all.CDFs(5, 5, 5, 1, 1000000)
g4 <- get.convol.emp.all.CDFs(5, 2, 5, 2, 1000000)
g5 <- get.convol.emp.all.CDFs(5, 3, 5, 3, 1000000)
g6 <- get.convol.emp.all.CDFs(5, 4, 5, 4, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="equal_mean_varying_std.pdf")

# one mean low, one high, varying sds
g1 <- get.convol.emp.all.CDFs(5, 1, 10, 1, 1000000)
g2 <- get.convol.emp.all.CDFs(5, 3, 10, 1, 1000000)
g3 <- get.convol.emp.all.CDFs(5, 5, 10, 1, 1000000)
g4 <- get.convol.emp.all.CDFs(5, 1, 10, 2, 1000000)
g5 <- get.convol.emp.all.CDFs(5, 1, 10, 3, 1000000)
g6 <- get.convol.emp.all.CDFs(5, 1, 10, 6, 1000000)
g7 <- get.convol.emp.all.CDFs(5, 2, 10, 2, 1000000)
g8 <- get.convol.emp.all.CDFs(5, 3, 10, 3, 1000000)
g9 <- get.convol.emp.all.CDFs(5, 5, 10, 5, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6,g7,g8,g9), 3,3, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="mean_large_medium_varying_std.pdf")

# Adhoc tests

# low mean, varying variance

g1 <- get.convol.emp.all.CDFs(1, 0.1, 1, 0.1, 1000000)
g2 <- get.convol.emp.all.CDFs(1, 0.5, 1, 0.1, 1000000)
g3 <- get.convol.emp.all.CDFs(1, 1, 1, 0.1, 1000000)
g4 <- get.convol.emp.all.CDFs(1, 0.2, 1, 0.2, 1000000)
g5 <- get.convol.emp.all.CDFs(1, 0.6, 1, 0.6, 1000000)
g6 <- get.convol.emp.all.CDFs(1, 1, 1, 1, 1000000)
g7 <- get.convol.emp.all.CDFs(1, 0.7, 1, 0.6, 1000000)
g8 <- get.convol.emp.all.CDFs(1, 0.7, 1, 0.9, 1000000)
g9 <- get.convol.emp.all.CDFs(1, 0.6, 1, 0.9, 1000000)

#g7 <- get.convol.emp.all.CDFs(1, 0.1, 2, 0.1, 1000000)
#g8 <- get.convol.emp.all.CDFs(1, 0.6, 2, 0.1, 1000000)
#g9 <- get.convol.emp.all.CDFs(1, 0.1, 2, 0.6, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6,g7,g8,g9), 3,3, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="low_mean_varying_var.pdf")


# -- extension (2)
g1 <- get.convol.emp.all.CDFs(1, 0.2, 1, 0.1, 1000000)
g2 <- get.convol.emp.all.CDFs(1, 0.3, 1, 0.1, 1000000)
g3 <- get.convol.emp.all.CDFs(1, 0.4, 1, 0.1, 1000000)
g4 <- get.convol.emp.all.CDFs(1, 0.4, 1, 0.2, 1000000)
g5 <- get.convol.emp.all.CDFs(1, 0.6, 1, 0.2, 1000000)
g6 <- get.convol.emp.all.CDFs(1, 0.8, 1, 0.2, 1000000)
g7 <- get.convol.emp.all.CDFs(1, 0.6, 1, 0.4, 1000000)
g8 <- get.convol.emp.all.CDFs(1, 0.7, 1, 0.4, 1000000)
g9 <- get.convol.emp.all.CDFs(1, 0.9, 1, 0.4, 1000000)

#g7 <- get.convol.emp.all.CDFs(1, 0.1, 2, 0.1, 1000000)
#g8 <- get.convol.emp.all.CDFs(1, 0.6, 2, 0.1, 1000000)
#g9 <- get.convol.emp.all.CDFs(1, 0.1, 2, 0.6, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6,g7,g8,g9), 3,3, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="low_mean_varying_var_2.pdf")

# -- extension (3)
g1 <- get.convol.emp.all.CDFs(1, 0.5, 1, 0.2, 1000000)
g2 <- get.convol.emp.all.CDFs(1, 0.6, 1, 0.3, 1000000)
g3 <- get.convol.emp.all.CDFs(1, 0.7, 1, 0.3, 1000000)
g4 <- get.convol.emp.all.CDFs(1, 0.8, 1, 0.4, 1000000)
g5 <- get.convol.emp.all.CDFs(1, 0.8, 1, 0.5, 1000000)
g6 <- get.convol.emp.all.CDFs(1, 0.9, 1, 0.5, 1000000)
g7 <- get.convol.emp.all.CDFs(1, 1, 1, 0.6, 1000000)
g8 <- get.convol.emp.all.CDFs(1, 1, 1, 0.7, 1000000)
g9 <- get.convol.emp.all.CDFs(1, 1, 1, 0.8, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6,g7,g8,g9), 3,3, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="low_mean_varying_var_3.pdf")


# Medium mean, varying variance
g1 <- get.convol.emp.all.CDFs(5, 0.1, 5, 0.1, 1000000)
g2 <- get.convol.emp.all.CDFs(5, 0.5, 5, 0.1, 1000000)
g3 <- get.convol.emp.all.CDFs(5, 1, 5, 0.1, 1000000)
g4 <- get.convol.emp.all.CDFs(5, 2, 5, 0.1, 1000000)
g5 <- get.convol.emp.all.CDFs(5, 3, 5, 0.1, 1000000)
g6 <- get.convol.emp.all.CDFs(5, 4, 5, 0.1, 1000000)
g7 <- get.convol.emp.all.CDFs(5, 1, 5, 1, 1000000)
g8 <- get.convol.emp.all.CDFs(5, 2, 5, 2, 1000000)
g9 <- get.convol.emp.all.CDFs(5, 5, 5, 5, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6,g7,g8,g9), 3,3, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="medium_mean_varying_var.pdf")

# -- extension (2)

g1 <- get.convol.emp.all.CDFs(5, 0.4, 5, 0.1, 1000000)
g2 <- get.convol.emp.all.CDFs(5, 0.5, 5, 0.1, 1000000)
g3 <- get.convol.emp.all.CDFs(5, 0.5, 5, 0.2, 1000000)
g4 <- get.convol.emp.all.CDFs(5, 0.6, 5, 0.2, 1000000)
g5 <- get.convol.emp.all.CDFs(5, 0.6, 5, 0.3, 1000000)
g6 <- get.convol.emp.all.CDFs(5, 0.7, 5, 0.3, 1000000)
g7 <- get.convol.emp.all.CDFs(5, 0.7, 5, 0.4, 1000000)
g8 <- get.convol.emp.all.CDFs(5, 0.8, 5, 0.4, 1000000)
g9 <- get.convol.emp.all.CDFs(5, 1, 5, 0.4, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6,g7,g8,g9), 3,3, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="medium_mean_varying_var_2.pdf")



# -- extension (3)

g1 <- get.convol.emp.all.CDFs(5, 1.5, 5, 1, 1000000)
g2 <- get.convol.emp.all.CDFs(5, 2, 5, 1.5, 1000000)
g3 <- get.convol.emp.all.CDFs(5, 2.5, 5, 2, 1000000)
g4 <- get.convol.emp.all.CDFs(5, 3, 5, 2.5, 1000000)
g5 <- get.convol.emp.all.CDFs(5, 3.5, 5, 3, 1000000)
g6 <- get.convol.emp.all.CDFs(5, 4, 5, 3.5, 1000000)
g7 <- get.convol.emp.all.CDFs(5, 4.5, 5, 4, 1000000)
g8 <- get.convol.emp.all.CDFs(5, 5, 5, 4.5, 1000000)
g9 <- get.convol.emp.all.CDFs(5, 4, 5, 3, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6,g7,g8,g9), 3,3, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="medium_mean_varying_var_3.pdf")


# -- extension (4)
g1 <- get.convol.emp.all.CDFs(5, 1.5, 5, 0.5, 1000000)
g2 <- get.convol.emp.all.CDFs(5, 2.5, 5, 1.5, 1000000)
g3 <- get.convol.emp.all.CDFs(5, 3, 5, 2, 1000000)
g4 <- get.convol.emp.all.CDFs(5, 3.5, 5, 2.5, 1000000)
g5 <- get.convol.emp.all.CDFs(5, 4, 5, 3, 1000000)
g6 <- get.convol.emp.all.CDFs(5, 5, 5, 3.5, 1000000)

stack_plots(list(g1,g2,g3,g4,g5,g6), 3,2, title.common="CDF Comparison", legend.common=T, axis.common=c("x", "F(x)"), fout="medium_mean_varying_var_4.pdf")
