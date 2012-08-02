require(ggplot2)
require(reshape2)
require(grid)
source('/Users/gwb/Hacks/Projects/myRUtils/graphics.R')


# --------------  CDF Experiments
"""
In this section, we try to understand when the different bandwidth parameters
should be use. For that, we plot the sample CDF and the theoretical CDF for different
sets of parameters, and see when/if it breaks

"""

#u = 3
#s = 0.1
#nl=10000
#ns=5000

get.CDFs <- function(u,s,nl, ns, n.nodes, quantile.low=0.01, quantile.high=0.99){
# For near perfect cdf matches, use: u=3, s=0.1, nl=10000, ns=5000
  
# Defines a large and a small sample
large.sample = exp(rnorm(nl, u, s))
small.sample = large.sample[seq(1,ns)]

# Obtains the quantiles for the samples. This allows us to make the density
# estimation more robust to sampling variations
qt.small <- quantile(x=small.sample, probs=c(quantile.low, quantile.high))
qt.large <- quantile(x=large.sample, probs=c(quantile.low, quantile.high))

# Computes the densities and an adjusted density
small.sample.density <- density(small.sample, n=n.nodes, bw="ucv", from=qt.small[1], to=qt.small[2])
large.sample.density <- density(large.sample, n=n.nodes, bw="nrd", from=qt.large[1], to=qt.large[2])

# Converts the densities into CDFs
small.sample.cdf <- cumsum(small.sample.density$y/sum(small.sample.density$y))
large.sample.cdf <- cumsum(large.sample.density$y/sum(large.sample.density$y))

# Creates the theoretica CDFs
# Note that we purposefully avoid using the plnorm method, which would
# give more accurate CDFs, but render the comparison useless, since ultimately
# we are not interested in computing the cdf but the kl divergence.

th.small.density <- dlnorm(small.sample.density$x, u, s)
th.large.density <- dlnorm(large.sample.density$x, u, s)

th.small.cdf <- cumsum(th.small.density/sum(th.small.density))
th.large.cdf <- cumsum(th.large.density/sum(th.large.density))

e.large.cdf <- ecdf(large.sample)(large.sample.density$x)

# Plots the results
dt <- data.frame(sample.small= small.sample.cdf, th.small = th.small.cdf, sample.large = large.sample.cdf, th.large = th.large.cdf, empirical.large = e.large.cdf)
m.dt <- melt(dt, measure.vars=c("sample.small", "th.small", "sample.large", "th.large", "empirical.large"))
m.dt$x <- c(small.sample.density$x, small.sample.density$x, large.sample.density$x, large.sample.density$x, large.sample.density$x)
g <- ggplot(data=m.dt, aes(x=x)) + geom_line(aes(y=value, color=variable)) + scale_x_log10() +
  custom_opts() +
    xlab("x") +
      ylab("F(x)") +
        scale_color_hue("", labels=c("ucv   ", "theoretical ucv   ", "nrd   ", "theoretical nrd   ", "empirical   ")) +
            opts(title=paste("Parameters: u = ", u, ", s = ", s, ", nodes = ", n.nodes, "\nlarge = ", formatC(nl), ", small = ", formatC(ns), sep=""))
         

return(g)
}


# Plot results - sample size modest for ucv  <= this is our prefered setup
g1 <- get.CDFs(1, 1, 1000000, 15000, 2^10)
g2 <- get.CDFs(1, 1.5, 1000000, 15000, 2^10)
g3 <- get.CDFs(1, 2, 1000000, 15000, 2^10)
g4 <- get.CDFs(2, 1, 1000000, 15000, 2^10)
g5 <- get.CDFs(2, 1.5, 1000000, 15000, 2^10)
g6 <- get.CDFs(2, 2, 1000000, 15000, 2^10)

stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3,fout="u1-2-n1024.pdf", function(x) opts(legend.position="bottom"))

# Plot results - sample size modest for ucv, nodes maxed
g1 <- get.CDFs(1, 1, 1000000, 15000, 2^11)
g2 <- get.CDFs(1, 1.5, 1000000, 15000, 2^11)
g3 <- get.CDFs(1, 2, 1000000, 15000, 2^11)
g4 <- get.CDFs(2, 1, 1000000, 15000, 2^11)
g5 <- get.CDFs(2, 1.5, 1000000, 15000, 2^11)
g6 <- get.CDFs(2, 2, 1000000, 15000, 2^11)

stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3,fout="u1-2-n2048.pdf", function(x) opts(legend.position="bottom"))

# Plot results - sample size maxed
g1 <- get.CDFs(1, 1, 1000000, 45000, 2^10)
g2 <- get.CDFs(1, 1.5, 1000000, 45000, 2^10)
g3 <- get.CDFs(1, 2, 1000000, 45000, 2^10)
g4 <- get.CDFs(2, 1, 1000000, 45000, 2^10)
g5 <- get.CDFs(2, 1.5, 1000000, 45000, 2^10)
g6 <- get.CDFs(2, 2, 1000000, 45000, 2^10)

stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3,fout="u1-2-n1024-ssmaxed.pdf", function(x) opts(legend.position="bottom"))

# Plot results - sample size maxed, nodes maxed 
g1 <- get.CDFs(1, 1, 1000000, 45000, 2^11)
g2 <- get.CDFs(1, 1.5, 1000000, 45000, 2^11)
g3 <- get.CDFs(1, 2, 1000000, 45000, 2^11)
g4 <- get.CDFs(2, 1, 1000000, 45000, 2^11)
g5 <- get.CDFs(2, 1.5, 1000000, 45000, 2^11)
g6 <- get.CDFs(2, 2, 1000000, 45000, 2^11)

stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3,fout="u1-2-n2048-ssmaxed.pdf", function(x) opts(legend.position="bottom"))

# Plot results - equal sample size
g1 <- get.CDFs(1, 1, 15000, 15000, 2^10)
g2 <- get.CDFs(1, 1.5, 15000, 15000, 2^10)
g3 <- get.CDFs(1, 2, 15000, 15000, 2^10)
g4 <- get.CDFs(2, 1, 15000, 15000, 2^10)
g5 <- get.CDFs(2, 1.5, 15000, 15000, 2^10)
g6 <- get.CDFs(2, 2, 15000, 15000, 2^10)

stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3,fout="u1-2-ssequal.pdf", function(x) opts(legend.position="bottom"))

# Plot results - tails exploration
g1 <- get.CDFs(1, 1, 1000000, 15000, 2^10, 0.001, 0.999)
g2 <- get.CDFs(1, 1.5, 1000000, 15000, 2^10, 0.001, 0.999)
g3 <- get.CDFs(1, 2, 1000000, 15000, 2^10, 0.001, 0.999)
g4 <- get.CDFs(2, 1, 1000000, 15000, 2^10, 0.001, 0.999)
g5 <- get.CDFs(2, 1.5, 1000000, 15000, 2^10, 0.001, 0.999)
g6 <- get.CDFs(2, 2, 1000000, 15000, 2^10, 0.001, 0.999)

stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3,fout="u1-2-tail.pdf", function(x) opts(legend.position="bottom"))

g1 <- get.CDFs(3, 1, 100000, 10000)
g2 <- get.CDFs(3, 1.5, 100000, 10000)
g3 <- get.CDFs(3, 2, 100000, 10000)
g4 <- get.CDFs(4, 1, 100000, 10000)
g5 <- get.CDFs(4, 1.5, 100000, 10000)
g6 <- get.CDFs(4, 2, 100000, 10000)

stack_plots(list(g1,g2,g3,g4,g5,g6), 2,3,fout="u3-4.pdf", function(x) opts(legend.position="bottom"))




# --------------  CDFs for convol of lognormals


get.convol.CDFs <- function(u1,s1, u2, s2,nl, ns){
# For near perfect cdf matches, use: u=3, s=0.1, nl=10000, ns=5000
  
# Defines a large and a small sample
large.sample = exp(rnorm(nl, u1, s1)) + exp(rnorm(nl, u2, s2))
small.sample = large.sample[seq(1,ns)]

# Obtains the quantiles for the samples. This allows us to make the density
# estimation more robust to sampling variations
qt.small <- quantile(x=small.sample, probs=c(0.01, 0.99))
qt.large <- quantile(x=large.sample, probs=c(0.01, 0.99))

# Computes the densities and an adjusted density
small.sample.density <- density(small.sample, bw="ucv", from=qt.small[1], to=qt.small[2])
large.sample.density <- density(large.sample, bw="nrd", from=qt.large[1], to=qt.large[2])

# Converts the densities into CDFs
small.sample.cdf <- cumsum(small.sample.density$y/sum(small.sample.density$y))
large.sample.cdf <- cumsum(large.sample.density$y/sum(large.sample.density$y))


# Plots the results
dt <- data.frame(sample.small= small.sample.cdf, sample.large = large.sample.cdf)
m.dt <- melt(dt, measure.vars=c("sample.small", "sample.large"))
m.dt$x <- c(small.sample.density$x, large.sample.density$x)
g <- ggplot(data=m.dt, aes(x=x)) + geom_line(aes(y=value, color=variable)) + scale_x_log10() +
  custom_opts() +
    xlab("x") +
      ylab("F(x)") +
        scale_color_hue("", labels=c("ucv   ", "nrd   ")) +
            opts(title=paste("Parameters: u1 = ", u1, ", s1 = ", s1, ", u2 = ", u2, ", s2 = ", s2, sep=""))
         

return(g)
}




# --------------  KL Divergence Experiments

get.kldivs <- function (u,s, nl, ns){

# Defines a large and a small sample
large.sample = exp(rnorm(nl, u, s))
small.sample = large.sample[seq(1,ns)]

# Obtains the quantiles for the samples. This allows us to make the density
# estimation more robust to sampling variations
qt.small <- quantile(x=small.sample, probs=c(0.05, 0.95))
qt.large <- quantile(x=large.sample, probs=c(0.05, 0.95))

# Computes the densities and an adjusted density
small.sample.density <- density(small.sample, bw="ucv", from=qt.small[1], to=qt.small[2])
large.sample.density <- density(large.sample, bw="nrd", from=qt.large[1], to=qt.large[2])
large.sample.adjusted.density <- density(large.sample, bw="nrd", from=qt.small[1], to=qt.small[2])

# Check that the evaluation points are the same
all(large.sample.density$x == small.sample.density$x)

# Computes the kl between each sample and the theoretical density.
# Note that two theoretical densities are computed to match the evaluation
# points of the different samples
th.small.density <- dlnorm(small.sample.density$x, u,s)
th.large.density <- dlnorm(large.sample.density$x, u,s)

# A control density is computed for each sample. It serves as a reference
# for what a *bad* kl is.
control.small.density <- dlnorm(small.sample.density$x, u+1, s+1)
control.large.density <- dlnorm(large.sample.density$x, u+1, s+1)

# Computes the different kl divergences
kl.small.th <- sum(th.small.density * (log(th.small.density) - log(small.sample.density$y)))
kl.large.th <- sum(th.small.density * (log(th.large.density) - log(large.sample.density$y)))
kl.large.adjusted.th <- sum(th.small.density * (log(th.small.density) - log(large.sample.adjusted.density$y)))
kl.control.small <- sum(th.small.density * (log(th.small.density) - log(control.small.density)))
kl.control.large <- sum(th.large.density * (log(th.large.density) - log(control.large.density)))

return(c(kl.small.th, kl.control.small, kl.large.th, kl.control.large, kl.large.adjusted.th))
}


# parameters
u <- seq(1,8)
s <- c(0.1, 0.5, 1, 5, 10)
parameters <- expand.grid(u=u,s=s)


# Runs the simulations
res <- apply(parameters, 1, function(x) get.kldivs(x[1],x[2], 100000, 10000))
res <- cbind(parameters, t(res))
dt <- data.frame(res)
colnames(dt) <- c("u","s","ucv","bli","ndr","foo","bar")

# Plots the results
ggplot(data=dt, aes(x=u, y=s)) + geom_tile(aes(fill=ucv-ndr))
