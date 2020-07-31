##### outliers
#' @export
elastic_outliers = function(elastic_depths, ka = 1.8, thresh = 1) {
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