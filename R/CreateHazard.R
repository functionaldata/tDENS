CreateHazard = function(y, c, optns = list()) {
  
  if(is.null(optns$kernel)) {
    kernel = 'gauss'
  } else {
    kernel =  optns$kernel
  }
  
  if(is.null(optns$delta)) {
    delta = max(c(diff(range(y))/1000, min(diff(sort(unique(y))))))
  } else {
    delta = optns$delta
  }
  
  if(is.null(optns$outputGrid)){
    outputGrid = NULL
  } else {
    outputGrid = optns$outputGrid
  }
  
  if(is.null(outputGrid)) {
    histgrid = seq(min(y) - delta*0.5,
                   max(y) + delta*0.5,
                   by = delta)
  } else {
    histgrid = seq(min(outputGrid) - delta*0.5,
                   max(outputGrid) + delta*0.5,
                   by = delta)
  }
  
  # Number of observations
  N <- length(y)
  d <- y[c == 0]
  c <- y[c == 1]
  
  # Grid for binning
  histgrid = seq(min(y) - delta*0.5,
                 max(y) + delta*0.5,
                 by = delta)
  if((max(y) + delta*0.5) > histgrid[length(histgrid)]) {
    histgrid[length(histgrid)] =  max(y) + delta*0.5;
  }
  
  # Number of bins
  m <- length(histgrid)
  
  # Histogram counts
  histObjD =  hist(d, breaks = histgrid, plot = FALSE)
  histObjC = hist(c, breaks = histgrid, plot = FALSE)
  x = histObjD$mids
  
  # Hazard estimation
  dHaz = histObjD$counts
  cHaz = histObjC$counts
  events = histObjD$counts + histObjC$counts
  
  n = N - cumsum(events)
  n_1 = head(c(N, n), -1)
  
  yHaz = 2*dHaz / (delta * (n + n_1))
  bwHaz = fdapace:::CVLwls1D(y = yHaz,
                             t = x,
                             kernel = kernel,
                             npoly = 1,
                             nder = 0,
                             dataType = 'Dense',
                             kFolds = 5)
  
  if(is.null(outputGrid)) {
    outputGrid = x
  }
  
  haz1 = fdapace:::Lwls1D(bw = bwHaz,
                          kernel_type = 'gauss',
                          win = rep(1, m-1),
                          xin = x,
                          yin = yHaz,
                          xout = outputGrid)
  
  # Final estimator and normalization
  haz = delta^(-1) * log((2 + delta * haz1) / (2 - delta * haz1))
  haz[haz < 0] = 0
  
  # Return list of output grid and hazard estimator
  out = list(x = outputGrid, haz = haz)
  return(out)
}
