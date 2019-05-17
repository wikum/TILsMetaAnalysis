

sink("log/3.log.txt", split=TRUE)

tryCatch({


  set.seed(1)
  
  library(knitr)
  library(ggplot2)
  library(gridExtra)
  library(ggrepel)
  library(rutils)
  library(mclust)
  library(RColorBrewer)
  library(WriteXLS)
  
  pal.default = c("#F8766D", "#619CFF", "#00BA38", "#C77CFF")
  log_fn = function(x) log2(x+1) 
  
  l1 = load("../data/data.rda")
  
  cat("Running immune checkpoint marker analysis..\n")
  
  # plot genes of interest in boxplot form
  gt0 = theme(strip.text.x=element_text(size=20))
  
  bx.genes = read.csv("../data/gene_list.txt", stringsAsFactors=FALSE, header=FALSE)[, 1]
  bx.genes.2 = bx.genes
  bx.genes.2[which(bx.genes == "PDCD1")] = "PD1"
  bx.genes.2[which(bx.genes == "HAVCR2")] = "TIM3"
  
  bx.list = utils.lapply_c(names(mat.list), function(tissue){
    M = log_fn(mat.list[[tissue]])
    P = pheno.list[[tissue]]
    
    
    M = M[bx.genes, ]
    rownames(M) = bx.genes.2
    t1 = utils.meltMat(M)
    t2 = data.frame(t1, class=sapply(t1$colKey, function(x) P[x, "class"]))
    t2$gene = t2
    
    ggplot(t2, aes(x=class, y=value, fill=class))+
      geom_boxplot()+ylab("expression")+
      gtheme.GENERIC()+gtheme.NO_X_LABS()+
      facet_wrap(~rowKey, ncol=3)+
      gt0
    
  })
  
  # ==============================================================================

  png("../fig/bx_GBM.png", width=800, height=600)
  tryCatch({
    print(bx.list$GBM)
  }, error = function(e){print(e)})
  dev.off()  

  png("../fig/bx_PRAD.png", width=800, height=600)
  tryCatch({
    print(bx.list$PRAD)
  }, error = function(e){print(e)})
  dev.off()  
  
  png("../fig/bx_RCC.png", width=800, height=600)
  tryCatch({
    print(bx.list$RCC)
  }, error = function(e){print(e)})
  dev.off()
  
  png("../fig/bx_BLCA.png", width=800, height=600)
  tryCatch({
    print(bx.list$BLCA)
  }, error = function(e){print(e)})
  dev.off()
  
  # ==============================================================================
  
  palette.density = c("violetred2", "maroon",
                      "lightblue", "steelblue",
                      "lightgreen", "turquoise4",
                      "yellow2", "orange",
                      "magenta", "purple")
  
  palette.density.2 = c("maroon", "steelblue", "turquoise4", "orange")
  
  
  n = length(mat.list)
  
  df.exp = Reduce(rbind, utils.lapply_i(mat.list, function(x, i, y){
    Reduce(rbind, utils.lapply_c(c("PDCD1", "HAVCR2", "LAG3"), function(g){
      t0 = log_fn(x[g, ])
      t1 = data.frame(tissue=y, expression=t0, gene=g)
      t1
    }))
  }))
  
  df.exp$label = factor(df.exp$gene, levels=c("PDCD1", "HAVCR2", "LAG3"), labels=c("PD1", "TIM3", "LAG3"))
  df.exp$sample = rownames(df.exp)
  
  # ==============================================================================
  
  png("../fig/densities.png", width=1200, height=800)
  tryCatch({
    
    print(
      
    ggplot(df.exp, aes(x=expression, fill=tissue))+
      geom_density()+
      scale_fill_manual(name="", values=palette.density.2)+
      gtheme.GENERIC()+
      xlab("expression")+
      facet_wrap(~tissue+label, ncol=3, scales="free")+
      theme(strip.text.x=element_text(size=20))
    
    )
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  
  
  gt = theme(plot.title = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="bold"))
  gt2 = theme(legend.justification="center", legend.position="bottom", legend.direction="vertical")
  gt3 = guides(shape=FALSE, col=FALSE)
  
  df.GBM = data.frame(
    PDCD1=log_fn(mat.list$GBM["PDCD1", ]),
    TIM3=log_fn(mat.list$GBM["HAVCR2", ]),
    LAG3=log_fn(mat.list$GBM["LAG3", ]),
    pheno.list$GBM
  )
  
  df.PRAD = data.frame(
    PDCD1=log_fn(mat.list$PRAD["PDCD1", ]),
    TIM3=log_fn(mat.list$PRAD["HAVCR2", ]),
    LAG3=log_fn(mat.list$PRAD["LAG3", ]),
    pheno.list$PRAD
  )
  
  df.RCC = data.frame(
    PDCD1=log_fn(mat.list$RCC["PDCD1", ]),
    TIM3=log_fn(mat.list$RCC["HAVCR2", ]),
    LAG3=log_fn(mat.list$RCC["LAG3", ]),
    pheno.list$RCC
  )
  
  df.BLCA = data.frame(
    PDCD1=log_fn(mat.list$BLCA["PDCD1", ]),
    TIM3=log_fn(mat.list$BLCA["HAVCR2", ]),
    LAG3=log_fn(mat.list$BLCA["LAG3", ]),
    pheno.list$BLCA
  )
  
  
  # GBM
  g1_1 =  ggplot(df.GBM[df.GBM$type == "CD8", ], aes(x=PDCD1, y=TIM3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("PD1 vs. TIM3")+
    xlab("PD1 expression")+ylab("TIM3 expression")+
    geom_text_repel(aes(x=PDCD1, y=TIM3, label=patient), size=5)
  
  g1_2 = ggplot(df.GBM[df.GBM$type == "CD8", ], aes(x=PDCD1, y=LAG3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("PD1 vs. LAG3")+
    xlab("PD1 expression")+ylab("LAG3 expression")+
    geom_text_repel(aes(x=PDCD1, y=LAG3, label=patient), size=5)
  
  g1_3 = ggplot(df.GBM[df.GBM$type == "CD8", ], aes(x=TIM3, y=LAG3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("TIM3 vs. LAG3")+
    xlab("TIM3 expression")+ylab("LAG3 expression")+
    geom_text_repel(aes(x=TIM3, y=LAG3, label=patient), size=5)
  
  # PRAD
  g2_1 = ggplot(df.PRAD[df.PRAD$type == "CD8", ], aes(x=PDCD1, y=TIM3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("PD1 vs. TIM3")+
    xlab("PD1 expression")+ylab("TIM3 expression")+
    geom_text_repel(aes(x=PDCD1, y=TIM3, label=patient), size=5)
  
  g2_2 = ggplot(df.PRAD[df.PRAD$type == "CD8", ], aes(x=PDCD1, y=LAG3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("PD1 vs. LAG3")+
    xlab("PD1 expression")+ylab("LAG3 expression")+
    geom_text_repel(aes(x=PDCD1, y=LAG3, label=patient), size=5)
  
  g2_3 = ggplot(df.PRAD[df.PRAD$type == "CD8", ], aes(x=TIM3, y=LAG3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("TIM3 vs. LAG3")+
    xlab("TIM3 expression")+ylab("LAG3 expression")+
    geom_text_repel(aes(x=TIM3, y=LAG3, label=patient), size=5)
  
  # RCC
  g3_1 = ggplot(df.RCC[df.RCC$type == "CD8", ], aes(x=PDCD1, y=TIM3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("PD1 vs. TIM3")+
    xlab("PD1 expression")+ylab("TIM3 expression")+
    geom_text_repel(aes(x=PDCD1, y=TIM3, label=patient), size=5)
  
  g3_2 =  ggplot(df.RCC[df.RCC$type == "CD8", ], aes(x=PDCD1, y=LAG3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("PD1 vs. LAG3")+
    xlab("PD1 expression")+ylab("LAG3 expression")+
    geom_text_repel(aes(x=PDCD1, y=LAG3, label=patient), size=5)
  
  g3_3 = ggplot(df.RCC[df.RCC$type == "CD8", ], aes(x=TIM3, y=LAG3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("TIM3 vs. LAG3")+
    xlab("TIM3 expression")+ylab("LAG3 expression")+
    geom_text_repel(aes(x=TIM3, y=LAG3, label=patient), size=5)
  
  # BLCA
  g4_1 = ggplot(df.BLCA[df.BLCA$type == "CD8", ], aes(x=PDCD1, y=TIM3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("PD1 vs. TIM3")+
    xlab("PD1 expression")+ylab("TIM3 expression")+
    geom_text_repel(aes(x=PDCD1, y=TIM3, label=patient), size=5)
  
  g4_2 = ggplot(df.BLCA[df.BLCA$type == "CD8", ], aes(x=PDCD1, y=LAG3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("PD1 vs. LAG3")+
    xlab("PD1 expression")+ylab("LAG3 expression")+
    geom_text_repel(aes(x=PDCD1, y=LAG3, label=patient), size=5)
  
  g4_3 = ggplot(df.BLCA[df.BLCA$type == "CD8", ], aes(x=TIM3, y=LAG3, shape=tissue, col=status))+
    geom_point(size=5)+
    gtheme.GENERIC()+gt+gt2+
    ggtitle("TIM3 vs. LAG3")+
    xlab("TIM3 expression")+ylab("LAG3 expression")+
    geom_text_repel(aes(x=TIM3, y=LAG3, label=patient), size=5)
  
  # ==============================================================================
  
  png("../fig/scatter.png", width=1200, height=1600)
  tryCatch({
    
      grid.arrange(g1_1, g1_2, g1_3,
                   g2_1, g2_2, g2_3,
                   g3_1, g3_2, g3_3,
                   g4_1, g4_2, g4_3,
                   ncol=3)
      
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  
  # apply EM
  
  genes = c("PDCD1", "HAVCR2", "LAG3")
  gene_labels = c("PD1", "TIM3", "LAG3")
  
  df.list = lapply(mat.list, function(mat) data.frame(t(log_fn(mat[genes, ]))))
  
  df.list.cd8 = lapply(1:length(mat.list), function(j){
    mat = mat.list[[j]]
    temp_pheno = pheno.list[[j]]
    data.frame(t(log_fn(mat[genes, temp_pheno$type == "CD8"])))
  })
  
  em.list = lapply(df.list, function(df)
    lapply(df, Mclust, G=2:3))
  names(em.list) = names(mat.list)
  
  em.list$BLCA$PDCD1 = Mclust(df.list$BLCA$PDCD1, G=2)
  
  source("../code/plotEMdistr.R")
  
  plot.all = function(exprDat, emOut, pheno, ...){
    
    out = plotEMdistr(exprDat=exprDat, emOut=emOut, ...)
    
    split.point = out[1]
    
    barplot(exprDat[order(exprDat)]-split.point)
    
    out
  }
  
  # ==============================================================================
  
  png("../fig/density_EM.png", width=800, height=600)
  tryCatch({
    
    par(mfrow=c(4, 3))
    for(i in 1:4){
      for(j in 1:3){
        plotEMdistr(df.list[[i]][[j]], emOut=em.list[[i]][[j]], 
                    name=sprintf("%s, %s", names(df.list)[i], gene_labels[j]))
      }
    }
    
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  
  plotWaterfalls = function(exprDat, emOut, gene, pheno, j=NULL, title, cohort, ...){
    
    out = plotEMdistr(exprDat=exprDat, emOut=emOut, name=gene, plot=FALSE, ...)
    if(is.null(j))
      j = length(out)
    split.point = out[j]
    
    temp_df = data.frame(x=exprDat-split.point, pheno)
    temp_df = temp_df[which(temp_df$type == "CD8"), ]
    temp_df = temp_df[order(temp_df$x, decreasing=TRUE), ]
    temp_df$sample = factor(temp_df$sample, levels=as.character(temp_df$sample))
    temp_df$cohort = cohort
    temp_df$title = title
    temp_df$cell_type = paste(temp_df$tissue, temp_df$status)
    
    if("PBMC EXP" %in% temp_df$cell_type){
      temp_df$cell_type = factor(temp_df$cell_type,
                                 levels=c("PBMC ACTIVE", "PBMC NAIVE", "TUMOR EXP", "PBMC EXP"))
    }else{
      temp_df$cell_type = factor(temp_df$cell_type,
                                 levels=c("PBMC ACTIVE", "PBMC NAIVE", "TUMOR EXP"))
    }
    
    g = ggplot(temp_df, aes(y=x, x=sample, fill=cell_type))+
      geom_col(size=2)+
      scale_fill_manual(name="", values=pal.default)+
      ggtitle(cohort)+ylab(sprintf("expression - %.2f", split.point))+
      gtheme.GENERIC()+
      facet_wrap(~title, scale="free", ncol=1)
    
    temp_df$group = factor(ifelse(temp_df$x > 0, "HIGH", "LOW"), levels=c("LOW", "HIGH"))
    
    list(out=out, split.point=split.point, df=temp_df, g=g)
  }
  
  emResults = utils.lapply_i(df.list, function(x, i, y){
    utils.lapply_i(x, function(z, j, a){
      plotWaterfalls(z, em.list[[y]][[a]], gene=a, pheno=pheno.list[[y]], title=gene_labels[j], cohort=y)
    })
  })
  
  save(emResults, file="../obj/emResults.rda")
  
  gt4 = theme(axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"))
  gt5 = theme(axis.text.y = element_text(colour="grey20",size=10,angle=0,hjust=1,vjust=0,face="plain"))
  
  water.list.GBM = lapply(emResults$GBM, function(x) x$g+coord_flip()+gt2+gt4+gt0)
  water.list.PRAD = lapply(emResults$PRAD, function(x) x$g+coord_flip()+gt2+gt5+gt0)
  water.list.RCC = lapply(emResults$RCC, function(x) x$g+coord_flip()+gt2+gt4+gt0)
  water.list.BLCA = lapply(emResults$BLCA, function(x) x$g+coord_flip()+gt2+gt4+gt0)
  grobs.list = c(water.list.GBM, water.list.PRAD, water.list.RCC, water.list.BLCA)
  
  # ==============================================================================
  
  png("../fig/waterfall.png", width=1200, height=1800)
  tryCatch({
    
    grid.arrange(grobs=grobs.list, ncol=3)
    
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  
  # summary
  
  df.count = data.frame(Reduce(rbind, lapply(1:length(emResults), function(j){
    x = emResults[[j]]
    z = t(sapply(1:length(x), function(i){
      y = table(x[[i]]$df$group[which(x[[i]]$df$type == "CD8")])
      w = c(names(x)[i], y, 
            paste(round(y[1] * 100/sum(y), 1), "%", sep=""),
            paste(round(y[2] * 100/sum(y), 1), "%", sep=""))
    }))
    cbind(rep(names(emResults)[j], nrow(z)), z)
  })))
  colnames(df.count) = c("Tissue", "Gene", "Negative Samples", "Positive Samples", "Negative %", "Positive %")
  df.count$Gene = factor(df.count$Gene, levels=c("PDCD1", "HAVCR2", "LAG3"), labels=c("PD1", "TIM3", "LAG3"))
  
  u = lapply(emResults, function(z){ 
    out = lapply(z, function(x){ 
      out.2 = x$df$group[which(x$df$type == "CD8")] == "HIGH"
      names(out.2) = rownames(x$df)[which(x$df$type == "CD8")]
      out.2
    })
    out.3 = lapply(out, function(out.2){
      out.2[names(out[[1]])]
    })
    out.3
  })
  u.2 = lapply(u, function(out){
    out.2 = Reduce(cbind, out)
    colnames(out.2) = names(out)
    out.2
  })
  v = t(sapply(u.2, function(x) sapply(0:3, function(y) sum(rowSums(x) == y))))
  w = data.frame(v)
  colnames(w) = c("All Negative", "Single Positive", "Double Positive", "Triple Positive")
  w = data.frame(cohort=rownames(w), w)
  

  a = t(apply(v, 1, function(x) paste(round(x * 100 / sum(x), 2), "%", sep="") ))
  b = data.frame(a)
  colnames(b) = c("All Negative", "Single Positive", "Double Positive", "Triple Positive")
  b = data.frame(cohort=rownames(b), b)
  
  WriteXLS::WriteXLS(list(df.count, w, b), "../tables/table1.xls", SheetNames=c("A", "B", "C"))
  
  
}, error = function(e){
  print(e)
})

print(sessionInfo())

rm(list=ls())
gc()

sink()


