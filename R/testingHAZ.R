library(fdapace)

source('CreateHazard.R')

# Exponential data
expDat <- rexp(1000, rate = 15)
expCensorIndicator <- rbinom(expDat, size = 1, prob = 0.1)
expHazard <- CreateHazard(y = expDat, c = expCensorIndicator)
plot(expHazard$x, expHazard$haz)

# Normal data
normDat <- rnorm(500)
normCensorIndicator <- rbinom(normDat, size = 1, prob = 0.2)
normHazard <- CreateHazard(y = normDat, c = normCensorIndicator)
plot(normHazard$x, normHazard$haz)

# Uniform data (add custom output grid)
uniDat <- runif(2000)
uniCensorIndicator <- rbinom(uniDat, size = 1, prob = 0.15)
uniHazard <- CreateHazard(y = uniDat, c = uniCensorIndicator,
                          optns = list(outputGrid = seq(0, 1, length.out = 500)))
plot(uniHazard$x, uniHazard$haz)
