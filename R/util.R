##### Shapes / curves in R2
srsf = function(f) {
  f_to_srvf(f, seq(0, 1, length.out = length(f)))
}

elastic.distance.r2 = function(f1, f2) {
  align = curve_pair_align(f1, f2)
  
  # amplitude distance
  # amp = sqrt(mean(colSums((align$q1 - align$q2n)^2)))
  amp = calc_shape_dist(f1, f2)
  
  # phase distance
  phs = acos(min(1, mean(rep(1, length(align$gam)) * srsf(align$gam))))
  
  return(list(Dx = phs, Dy = amp))
}



##### trajectories on S
rescale = function(x) (x - min(x)) / (max(x) - min(x))


inner_s2 = function(f1, f2) {
  sum(diag(f1 %*% t(f2)))
}

norm_s2 = function(f) {
  sqrt(inner_s2(f, f))
}

inv_exp = function(f1, f2) {
  theta = acos(inner_s2(f1, f2))
  
  if ((is.nan(theta))|| (theta < 1e-10)){
    exp_inv = rep(0,length(f1))
  } else {
    exp_inv = theta / sin(theta) * (f2 - cos(theta)*f1)
  }
  return(exp_inv)
}


parallel_transport = function(h, p1, p2) {
  h - 2 * inner_s2(h, p2) * (p1 + p2) / norm_s2(p1 + p2)^2
}

backward_parallel_translation = function(p, v) {
  qp = v
  
  for(i in rev(2:nrow(p))) {
    qp = parallel_transport(qp, p[i, ], p[i-1, ])
  }
  
  qp
}

tsrvf = function(f) {
  t = nrow(f)
  n = ncol(f)
  
  adot = t(sapply(1:(t-1), function(x) (t-1) * inv_exp(f[x,], f[x+1,])))
  q = rbind(adot / sqrt(norm_s2(adot)), 0)
  
  q[t,] = parallel_transport(q[(t-1), ], f[t-1, ], f[t, ])
  
  q[2:t,] = sapply(2:t, function(x) backward_parallel_translation(f[1:x, ], q[x, ]))
  
  return(q)
}

warp_p_gam = function(p, gam) {
  pn = p
  t = seq(0, 1, length.out = nrow(p))
  for(i in 1:ncol(p)) {
    pn[,i] = approx(t, p[,i], gam)$y
  }
  return(pn)
}

pair_align_path = function(p1, p2) {
  h1 = tsrvf(p1)
  h2 = tsrvf(p2)
  
  n = nrow(h1)
  
  gam = .Call('DPQ', PACKAGE = 'fdasrvf', h1, h2, 3, n, 0, 0, rep(0,n))
  gam = rescale(invertGamma(gam))
  
  p2.warp = warp_p_gam(p2, gam)
  return(list(gam = gam, p2.warp = p2.warp))
}

elastic.distance.s2 = function(f, g) {
  
  # convert from the Rn format to the S2 format
  # really should have made these consistent...
  f = t(f)
  g = t(g)
  align = pair_align_path(f, g)
  
  # amplitude dist
  amp = sqrt(mean((tsrvf(f) - tsrvf(align$p2.warp))^2))
  
  # phase dist
  ident = rep(1, nrow(f))
  phs = acos(min(mean(ident * srsf(align$gam)), 1))
  
  return(list(Dx = phs, Dy = amp))
}
