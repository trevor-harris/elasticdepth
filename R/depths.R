#' Computes the amplitude and phase depths of the R1 trajectories in \code{fmat}
#' @param fmat Matrix of trajectories in R1. Each column is a function.
#'
#' @return amplitude and phase depths in a list
#'
#' @export
depth.R1 = function(fmat) {
  # save the dimensions for convenience
  obs = nrow(fmat)
  fns = ncol(fmat)
  
  amp_dist = matrix(0, fns, fns)
  phs_dist = matrix(0, fns, fns)
  
  time = seq(0, 1, length.out = nrow(fmat))
  
  for (f1 in 1:(fns-1)) {
    
    dist = future_sapply(f1:ncol(fmat), function(y) {
      unlist(elastic.distance(fmat[,f1], fmat[,y], time))
    })
    
    phs_dist[f1, f1:fns] = dist[2,]
    amp_dist[f1, f1:fns] = dist[1,]
  }
  
  amp_dist = amp_dist + t(amp_dist)
  phs_dist = phs_dist + t(phs_dist)
  
  amp = 1 / (1 + apply(amp_dist, 1, median))
  phase = 1 / (1 + apply(phs_dist, 1, median))
  phase = ((2+pi)/pi) * (phase - 2/(2+pi))
  return(list(amplitude = amp, phase = phase))
}


#' Computes the amplitude and phase depths of the R2 trajectories in \code{tmat}
#' @param tmat Array of trajectories in R2. First axis is the X, Y coordinates. Second axis
#' is the trajectory value, third axis indexes the trajetories.
#'
#' @return amplitude and phase depths in a list
#'
#' @export
depth.R2 = function(tmat) {
  
  fns = dim(tmat)[3]
  
  amp_dist = matrix(0, fns, fns)
  phs_dist = matrix(0, fns, fns)
  
  for (f in 1:(fns-1)) {
    
    dist = future_sapply(f:fns, function(y) {
      unlist(elastic.distance.r2(tmat[,,f], tmat[,,y])) 
    })
    
    amp_dist[f, f:fns] = dist[2,]
    phs_dist[f, f:fns] = dist[1,]
  }
  
  amp_dist = amp_dist + t(amp_dist)
  phs_dist = phs_dist + t(phs_dist)
  
  amp = 1 / (1 + apply(amp_dist, 1, median))
  phase = 1 / (1 + apply(phs_dist, 1, median))
  phase = ((2+pi)/pi) * (phase - 2/(2+pi))
  return(list(amplitude = amp, phase = phase))
}

#' Computes the amplitude and phase depths of the S2 trajectories in \code{tmat}
#' @param tmat Array of trajectories in S2. First axis is the X, Y, Z coordinates. Second axis
#' is the trajectory value, third axis indexes the trajetories.
#'
#' @note 3D euclidean coordinates are automatically converted to spherical coordinates.
#'
#' @return amplitude and phase depths in a list
#'
#' @export
depth.S2 = function(tmat) {
  
  fns = dim(tmat)[3]
  
  amp_dist = matrix(0, fns, fns)
  phs_dist = matrix(0, fns, fns)
  
  for (f in 1:(fns-1)) {
    dist = future_sapply(f:fns, function(y) {
      unlist(elastic.distance.s2(tmat[,,f], tmat[,,y])) 
    })
    
    phs_dist[f, f:fns] = dist[1,]
    amp_dist[f, f:fns] = dist[2,]
  }
  
  amp_dist = amp_dist + t(amp_dist)
  phs_dist = phs_dist + t(phs_dist)
  
  amp = 1 / (1 + apply(amp_dist, 1, median))
  phase = 1 / (1 + apply(phs_dist, 1, median))
  phase = ((2+pi)/pi) * (phase - 2/(2+pi))
  return(list(amplitude = amp, phase = phase))
}




