# Temporary testing for dens2lhz.R
# Matt Dawson
# 2017

library(fdapace)
source('dens2lhz.R')
source('lhz2dens.R')
source('CreateDensity.R')

#===========================================================#
# Exponential Data
#===========================================================#

# Create data for testing function
tmpData <- rexp(10000, rate = 40)

# Estimate density function from data
tmpDensity <- CreateDensity(tmpData, optns = list(nRegGrid = 1000))

# Log Hazard transformation (specify delta)
tmpLHZ <- dens2lhz(dens = tmpDensity$y,
                   dSup = tmpDensity$x,
                   lhSup = seq(0, 1, length.out = 1000),
                   delta = .2)

# Check to see if the log hazard function makes sense
# (note that for an exponential distribution the hazard function is
# constant at lambda (rate))
plot(tmpLHZ)
log(40)

# Inverse transformation, log-hazard to density
plot(lhz2dens(lhz = tmpLHZ$y,
              lhSup = tmpLHZ$x,
              delta = tmpLHZ$delta,
              dSup = tmpDensity$x),
     type = 'l', lwd = 2)

lines(tmpDensity, col = 'red')

#===========================================================#
# Beta Data
#===========================================================#

# Create data for testing function
tmpData2 <- rbeta(1000, shape1 = 2, shape2 = 3)

# Estimate density function from data
tmpDensity2 <- CreateDensity(tmpData2, optns = list(nRegGrid = 1000))

# Log Hazard transformation (no delta specified)
tmpLHZ2 <- dens2lhz(dens = tmpDensity2$y,
                    dSup = tmpDensity2$x,
                    lhSup = seq(0, 1, length.out = 1000),
                    delta = .1)

# Plot the log hazard function
plot(tmpLHZ2)

# Inverse transformation, log-hazard to density
plot(lhz2dens(lhz = tmpLHZ2$y,
              lhSup = tmpLHZ2$x,
              delta = tmpLHZ2$delta,
              dSup = tmpDensity2$x),
     type = 'l', lwd = 2)

lines(tmpDensity2, col = 'red')

#===========================================================#
# Normal Data
#===========================================================#

# Create data for testing function
tmpData3 <- rnorm(1000, mean = 3, sd = 5)

# Estimate density function from data
tmpDensity3 <- CreateDensity(tmpData3, optns = list(nRegGrid = 2000))

# Log Hazard transformation (no delta specified)
tmpLHZ3 <- dens2lhz(dens = tmpDensity3$y,
                    dSup = tmpDensity3$x,
                    lhSup = seq(0, 1, length.out = 1000),
                    delta = .1)

# Plot the log hazard function
plot(tmpLHZ3)

# Inverse transformation, log-hazard to density
plot(lhz2dens(lhz = tmpLHZ3$y,
              lhSup = tmpLHZ3$x,
              delta = tmpLHZ3$delta,
              dSup = tmpDensity3$x),
     type = 'l', lwd = 2)

lines(tmpDensity3, col = 'red')

