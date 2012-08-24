require(data.table)

# First Simulation parameters
#params <- expand.grid(m1 = seq(-3, 6, 1), m2 = seq(-3, 6, 1), sd1 = seq(-3, 8, .5), sd2 = seq(-3, 8, .5))
#params <- data.table(params)
#reduced.params <- params[m2>m1, ]
#write.csv(reduced.params, "params.csv")

# Second simulation
params <- expand.grid(m1 = c(-1,1,3,6) , m2 = c(-1, 1, 3, 6), sd1 = seq(-3, 8, .2), sd2 = seq(-3, 8, .2))
params <- data.table(params)
reduced.params <- params[m2>m1, ]
write.csv(reduced.params, "params-2.csv")


# Third simulation
params <- expand.grid(m1 = seq(-1, 30, 5) , m2 = seq(-1, 30, 5), sd1 = seq(-3, 30, 3), sd2 = seq(-3, 30, 3))
params <- data.table(params)
reduced.params <- params[m2>m1, ]
write.csv(reduced.params, "params-3.csv")


# First reversed simulation
params <- expand.grid(m1= seq(-1, 30, 1) , m2 = seq(-1, 30, 1), sd1 = c(1 , 2, 3), sd2 = c(1, 2, 3))
params <- data.table(params)
reduced.params <- params[sd2>sd1, ]

write.csv(reduced.params, "params-reversed-1.csv")


# Second reversed simulation
params <- expand.grid(m1= seq(-1, 30, 1) , m2 = seq(-1, 30, 1), sd1 = c(0.5 , 5, 10), sd2 = c(0.5, 5, 10))
params <- data.table(params)
reduced.params <- params[sd2>sd1, ]

write.csv(reduced.params, "params-reversed-2.csv")


# Second reversed simulation
params <- expand.grid(m1= seq(0, 15, 0.25) , m2 = seq(0, 15, 0.25), sd1 = c(5), sd2 = c(10))
params <- data.table(params)
reduced.params <- params[sd2>sd1, ]

write.csv(reduced.params, "params-reversed-2.csv")

