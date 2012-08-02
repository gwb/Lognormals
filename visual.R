require(ggplot2)
require(reshape2)
require(grid)
source('/Users/gwb/Hacks/Projects/myRUtils/graphics.R')


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
