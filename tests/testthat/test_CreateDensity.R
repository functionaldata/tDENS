
makePNG = FALSE
if(makePNG) png( "HadesEffortNorm.png" , width = 16.92, height = 9, units = 'in', res=213 )
par(mfrow=c(2,3))
makeComparisonPlotN <- function(N, mySeed = 123, twoGaussians = FALSE){
  
  set.seed(mySeed)
  N = N
  if(! twoGaussians){
    asdf2  = (rnorm(N))
    plot(density(asdf2, bw = "SJ"), main= "Gaussian N(0,1)", xlab = paste0(collapse = '', c( "N = ", as.character(N))))
    Unormal = CreateDensity(y = asdf2)
    
    lines(col='magenta' , Unormal$x, dnorm(Unormal$x))
  } else {
    asdf2  = c( (rnorm(N)) , rnorm(N, mean = 6) ) 
    plot(density(asdf2, bw = "SJ"), main= "Bimodal Gaussian N(0,1) + N(0,6)", xlab = paste0(collapse = '', c( "N = ", as.character(2*N))))
    Unormal = CreateDensity(y = asdf2)
    lines(col='magenta' , Unormal$x, dnorm(Unormal$x) * 0.5 + dnorm(Unormal$x, mean = 6) * 0.5)
  }
  lines(density(asdf2), col='red') 
  lines(col='green', x = Unormal$x, y = Unormal$y)
  abline(v = min(asdf2))
  abline(v = max(asdf2))
  legend(legend = c("Marron", "R-default", "HADES-like", "True PDF"), lwd= 2, col=c("black", "red", "green", "magenta"), bty = 'n', 'topright')
}
makeComparisonPlotN(100)
makeComparisonPlotN(2000)
makeComparisonPlotN(40000)
makeComparisonPlotN(100, twoGaussians = TRUE)
makeComparisonPlotN(2000, twoGaussians = TRUE)
makeComparisonPlotN(40000, twoGaussians = TRUE)
if(makePNG) dev.off()



if(makePNG) png( "HadesEffortExp.png" , width = 16.92, height = 4.5, units = 'in', res=213 )
par(mfrow=c(1,3))
makeComparisonPlotE <- function(N, mySeed = 123){
  
  set.seed(mySeed)
  N = N 
  asdf2  = (rexp(N, rate = 1.5))
  plot(density(asdf2, bw = "SJ"), main= "Exponential (rate=1.5)", xlab = paste0(collapse = '', c( "N = ", as.character(N))) , ylim = c(0, 1.45))
  lines(density(asdf2), col='red')
  Unormal = CreateDensity(y = asdf2, list( infSupport=FALSE))
  lines(col='green', x = Unormal$x, y = Unormal$y)
  lines(col='magenta' , Unormal$x, dexp(Unormal$x, rate = 1.5))
  abline(v = min(asdf2))
  abline(v = max(asdf2))
  legend(legend = c("Marron", "R-default", "HADES-like", "True PDF"), lwd= 2, col=c("black", "red", "green", "magenta"), bty = 'n', 'topright')
}
makeComparisonPlotE(100)
makeComparisonPlotE(2000)
makeComparisonPlotE(40000) 
if(makePNG) dev.off()



if(makePNG) png( "HadesEffortUnif.png" , width = 16.92, height = 4.5, units = 'in', res=213 )
par(mfrow=c(1,3))
makeComparisonPlotE <- function(N, mySeed = 123){
  
  set.seed(mySeed)
  N = N 
  asdf2  = (runif(N ))
  plot(density(asdf2, bw = "SJ"), main= "Uniform (0,1)", xlab = paste0(collapse = '', c( "N = ", as.character(N))) , ylim = c(0, 1.45))
  lines(density(asdf2), col='red')
  Unormal = CreateDensity(y = asdf2,  list( infSupport=FALSE) )
  lines(col='green', x = Unormal$x, y = Unormal$y)
  lines(col='magenta' , Unormal$x, dunif(Unormal$x))
  abline(v = min(asdf2))
  abline(v = max(asdf2))
  legend(legend = c("Marron", "R-default", "HADES-like", "True PDF"), lwd= 2, col=c("black", "red", "green", "magenta"), bty = 'n', 'topright')
}
makeComparisonPlotE(100)
makeComparisonPlotE(2000)
makeComparisonPlotE(40000) 
if(makePNG) dev.off()
