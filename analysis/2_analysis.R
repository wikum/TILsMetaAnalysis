

sink("log/2.log.txt", split=TRUE)

tryCatch({

  set.seed(1)
  
  library(knitr)
  library(ggplot2)
  library(gridExtra)
  library(ggrepel)
  library(rutils)
  library(RColorBrewer)
  library(gplots)
  library(limma)
  library(edgeR)
  library(WriteXLS)
  
  pal.default = c("#F8766D", "#619CFF", "#00BA38", "#C77CFF")
  log_fn = function(x) log2(x+1) 
  
  l1 = load("../data/data.rda")
  l2 = load("../data/pheno.rda")
  
  colnames(mat.list$PRAD) = pheno.list$PRAD[colnames(mat.list$PRAD), "sample"]

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
  
  cat("Running CD8 differential expression analysis..\n")
  
  table.list = lapply(names(pheno.list), function(tissue){
    
    x = mat.list[[tissue]]
    y = pheno.list[[tissue]]
    
    y$class = factor(y$class)
    
    sel_exp = which(y$class == "TUMOR_CD8_EXP")
    sel_act = which(y$class == "PBMC_CD8_ACTIVE")
    
    x = x[, c(sel_exp, sel_act)]
    y = y[c(sel_exp, sel_act), ]
    
    y$class = factor(as.character(y$class))
    
    mod = model.matrix(~ 0 + patient + class, y)
    
    con = makeContrasts(group=classTUMOR_CD8_EXP, levels=mod)
    
    dge = DGEList(counts=x, group=y$class, remove.zeros=TRUE, genes=rownames(x))
    keep = rowSums(cpm(dge) > 1) >= 1
    dge = dge[keep, ]
    dge$samples$lib.size = colSums(dge$counts)
    
    dge_n = calcNormFactors(dge)
    
    V = voom(dge_n, mod, plot=FALSE)
    lmfit_voom = lmFit(V, mod)
    lmfit_con = contrasts.fit(lmfit_voom, contrasts=con)
    
    top_table = topTable(eBayes(lmfit_con), coef=1, number=Inf, adjust.method="bon")
    
    top_table = top_table[order(top_table$adj.P.Val, decreasing=FALSE), ]
    
    dt = decideTests(eBayes(lmfit_con), adjust.method="bon", method="global",lfc=1, p.value=0.05)
    
    list(table=top_table, dt=dt)
    
  })
  names(table.list) = names(pheno.list)
  
  table.list.global = lapply(table.list, function(x) x$table)
  
  table.list.global = global_p_adj(table.list.global, method="BH")
  table.list.global = lapply(table.list.global, function(x) x[order(x$adj.P.Val), ])
  
  venn.list = lapply(table.list.global, function(x) rownames(x)[which(x$adj.P.Val <= 0.00001)])
  
  int = Reduce(intersect, venn.list)
  
  mat.sign = sapply(table.list.global, function(x) sign(x[int, "logFC"]) )
  
  plot.list = lapply(table.list, function(x){
    
    T = x$table
    dt = x$dt
    
    df = data.frame(T, differential_expression=factor(dt[rownames(T), 1], levels=c(-1, 0, 1)),
                    label=rep(FALSE, nrow(T)))
    
    #df$label[1:min(20, sum(df$differential_expression != 0))] = TRUE
    df$label[which(df$gene %in% int)] = TRUE
    
    df$col = rep("black", nrow(df))
    df$col[df$label == TRUE & df$differential_expression == 1] = "darkred"
    df$col[df$label == TRUE & df$differential_expression == -1] = "darkblue"
    
    ggplot(df, aes(x=logFC, y=B, col=differential_expression, label=label))+
      geom_point(size=3, alpha=.8)+
      scale_color_manual(values=c("steelblue", "darkgray", "salmon"), 
                         breaks=c(-1, 1),
                         labels=c(sprintf("logFC < -1, P < 0.05\n         (n=%d)\n", sum(df$differential_expression == -1)), 
                                  sprintf("logFC > 1, P < 0.05\n         (n=%d)", sum(df$differential_expression == 1)))
      )+
      geom_text_repel(data=subset(df, df$label == TRUE), 
                      aes(x=logFC, y=B, label=genes), 
                      nudge_x=.1, nudge_y=0, color=df$col[df$label],
                      show.legend=FALSE, fontface = 'bold')+
      ylim(min(df$B), max(df$B)+8)+
      xlab("logFC")+ylab("log Odds")+
      gtheme.GENERIC()+
      theme(legend.title=element_blank(),
            legend.position=c(.5, .9),
            legend.text = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
            legend.box.background=element_rect(fill=alpha("white", .9), size=10, color="white"))
    
  })
  
  # ==============================================================================
  
  png("../fig/vol_GBM.png", width=800, height=1200)
  tryCatch({
    print(plot.list$GBM)
  }, error = function(e){print(e)})
  dev.off()  
  
  png("../fig/vol_PRAD.png", width=800, height=1200)
  tryCatch({
    print(plot.list$PRAD)
  }, error = function(e){print(e)})
  dev.off()  

  png("../fig/vol_RCC.png", width=800, height=1200)
  tryCatch({
    print(plot.list$RCC)
  }, error = function(e){print(e)})
  dev.off()  

  png("../fig/vol_BLCA.png", width=800, height=1200)
  tryCatch({
    print(plot.list$BLCA)
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  
}, error = function(e){
  print(e)
})

print(sessionInfo())

rm(list=ls())
gc()

sink()


