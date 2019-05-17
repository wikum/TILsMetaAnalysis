
global_p_adj = function(tblist, method="BH"){
  P = lapply(tblist, function(x) x$P.Value)
  pv = Reduce(c, P)
  pv.adj = p.adjust(pv, method=method)
  j = 0
  for(i in 1:length(tblist)){
    a = j + 1
    b = length(P[[i]]) + j
    tblist[[i]]$adj.P.Val = pv.adj[a:b]
    j = b
  }
  tblist
}
