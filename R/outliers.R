#' Computes the amplitude and phase depths of the S2 trajectories in \code{tmat}
#' @param elastic_depths elastic depths produced by one of the depth.R1, depth.R2, or 
#' depth.S2 functions. Both amplitude and phase depths need to be included.
#' @param ka Controls the outlier bound / whisker depth. Higher numbers make outlier
#' bound further away from the median. Defaults to 1.8
#' @param thresh Controls the amount of depth thresholding. Set from 0 to 1, where
#' 0 means only the top 0% of functions can be outliers and 1 means the top 100% of 
#' functions can be outliers. Defaults to 0.5.
#'
#' @return amplitude and phase depths in a list
#'
#' @export
elastic_outliers = function(elastic_depths, ka = 1.8, thresh = 0.5) {
  amp = elastic_depths$amplitude
  phs = elastic_depths$phase
  
  amp.100 = max(amp)
  phs.100 = max(phs)
  
  amp.50 = as.numeric(quantile(amp, 0.5))
  phs.50 = as.numeric(quantile(phs, 0.5))
  
  amp.iqr = amp.100 - amp.50
  phs.iqr = phs.100 - phs.50
  
  amp.lim = max(amp.50 - ka*amp.iqr, 0)
  phs.lim = max(phs.50 - ka*phs.iqr, 0)
  
  amp.thre = as.numeric(quantile(amp, thresh))
  phs.thre = as.numeric(quantile(phs, thresh))
  
  amp.out = (amp < amp.lim)*(amp < amp.thre)
  phs.out = (phs < phs.lim)*(phs < phs.thre)
  
  return(list(amp = amp.out, phs = phs.out))
}