


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
th.small.density <- dlnorm(small.sample.density$x)
th.large.density <- dlnorm(large.sample.density$x)

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
