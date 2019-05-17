

sink("log/4.log.txt", split=TRUE)

tryCatch({

  set.seed(1)
  
  library(knitr)
  library(ggplot2)
  library(gridExtra)
  library(ggrepel)
  library(rutils)
  library(RColorBrewer)
  library(gplots)
  library(data.tree)
  library(limma)
  library(edgeR)
  library(WriteXLS)
  
  pal.default = c("#F8766D", "#619CFF", "#00BA38", "#C77CFF")
  log_fn = function(x) log2(x+1) 
  
  l1 = load("../data/data.rda")
  colnames(mat.list$PRAD) = pheno.list$PRAD[colnames(mat.list$PRAD), "sample"]
  
  l2 = load("../obj/emResults.rda")
  
  cat("Running differential expression analysis..\n")
  
  groups.list = lapply(emResults, function(x){
    
    df.list = lapply(x, function(y){
      z = y$df[, c("class", "sample", "group", "patient")]
      rownames(z) = z$sample
      z
    })
    samples = as.character(df.list$PDCD1$sample)
    
    df.all = sapply(df.list, function(y) y[samples, c("group")])
    df.all = data.frame(sample=samples, 
                        group=df.list[[1]][samples, "class"], 
                        patient=df.list[[1]][samples, "patient"], 
                        df.all)
    colnames(df.all)[4:5] = c("PD1", "TIM3")
    
    df.tree = df.all
    df.tree$PD1 = paste("PD1", ifelse(df.tree$PD1 == "HIGH", "+", "-"), sep="")
    df.tree$LAG3 = paste("LAG3", ifelse(df.tree$LAG3 == "HIGH", "+", "-"), sep="")
    df.tree$TIM3 = paste("TIM3", ifelse(df.tree$TIM3 == "HIGH", "+", "-"), sep="")
    
    df.tree.2 = df.tree
    
    df.tree$pathString = paste("*", df.tree$PD1, df.tree$LAG3, df.tree$TIM3, sep="/")
    df.node = as.Node(df.tree)
    
    df.tree.2$pathString = paste("*", df.tree$PD1, df.tree$LAG3, df.tree$TIM3, df.tree$group, sep="/")
    
    z = t(sapply(unique(df.tree.2$pathString), function(y){
      pd = df.tree.2$PD1[which(df.tree.2$pathString == y)[1]]
      pt = df.tree.2$patient[which(df.tree.2$pathString == y)]
      c(pd, paste(y, "[", paste(pt, collapse=","), "]"))
    }))
    
    w = t(sapply(unique(df.tree.2$pathString), function(y){
      pd = df.tree.2$PD1[which(df.tree.2$pathString == y)[1]]
      pt = df.tree.2$patient[which(df.tree.2$pathString == y)]
      c(pd, paste(y, paste(pt, collapse="\n"), sep="/"))
    }))
    
    w = data.frame(w)
    colnames(w) = c("PD1", "pathString")
    w$pathString = sapply(w$pathString, function(x) gsub("_", "\n", x))
    w.node.low = as.Node(w[which(w$PD1 == "PD1-"),])
    w.node.high = as.Node(w[which(w$PD1 == "PD1+"),])
    

    df.tree.3 = data.frame(z)
    colnames(df.tree.3) = c("PD1", "pathString")
    df.tree.3.list.low = ToListExplicit(as.Node(df.tree.3[which(df.tree.3$PD1 == "PD1-"), ]), 
                                        unname = TRUE)
    df.tree.3.list.high = ToListExplicit(as.Node(df.tree.3[which(df.tree.3$PD1 == "PD1+"), ]), 
                                         unname = TRUE)

    list(df.all=df.all,
         df.tree=df.tree,
         df.tree.2=df.tree.2,
         df.tree.3=df.tree.3,
         df.node=df.node,
         df.tree.3.list.low=df.tree.3.list.low,
         df.tree.3.list.high=df.tree.3.list.high,
         w.node.low=w.node.low,
         w.node.high=w.node.high)
  })
  
  
  
  
  pal = c("maroon", "steelblue", "cyan", "magenta")
  
  heatmat.list = lapply(groups.list, function(gl){
    
    df = gl$df.all
    
    u = grep("CD8", df$group)
    df = df[u, ]
    rownames(df) = df$sample
    
    df$sample = as.character(df$sample)
    
    ordered_groups = c("PBMC_CD8_ACTIVE", "PBMC_CD8_NAIVE", "TUMOR_CD8_EXP", "PBMC_CD8_EXP")
    ordered_groups = ordered_groups[ordered_groups %in% grep("CD8", levels(df$group), value=TRUE)]
    
    u = lapply(unique(df$patient), function(y){
      
      Reduce(rbind, lapply(ordered_groups, function(x){  
        
        z = as.vector(df[which(df$group == x & df$patient == y), ])
        sample = as.character(z["sample"])
        
        v = t(sapply(c("PD1", "LAG3", "TIM3"), 
                     function(w) c(sample, x, y, w, ifelse(z[w] == "HIGH", 1, 0)) ))
        rownames(v) = paste(v[, 2], v[, 4], sep="_")
        
        v
        
      }))
      
    })
    names(u) = unique(df$patient)
    
    v = sapply(u, function(x){
      
      y = as.numeric(x[, 5])
      names(y) = rownames(x)
      y
      
    })
    
    groups = u[[1]][, 2]
    genes = u[[1]][, 4]
    cols = as.character(factor(groups, levels=ordered_groups, 
                               labels=pal.default[1:length(ordered_groups)]))
    
    genes2 = genes

    
    list(u=u, v=v, groups=groups, genes=genes2, 
         cols=cols, ordered_groups=gsub("_CD8_", " ", ordered_groups))
    
  })
  
  
  # ==============================================================================
  
  png("../fig/heat_GBM.png", width=800, height=600)
  tryCatch({

    h = heatmat.list[[1]]
    
    heatmap.2(h$v, 
            xlab="patient",
            Rowv='none',
            trace='non',
            dendrogram='column',
            labRow=h$genes,
            col=c("gray20", "orange"),
            colsep=1:ncol(h$v),
            rowsep=1:nrow(h$v),
            sepcolor="black",
            sepwidth=c(0.01, 0.01),
            RowSideColors=h$cols,
            key=FALSE,
            margins=c(4, 5),
            main=names(heatmat.list)[1]
    )
    legend("bottomleft", legend=h$ordered_groups, fill=pal.default, cex=.7)
  
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  # ==============================================================================
  
  png("../fig/heat_PRAD.png", width=800, height=600)
  tryCatch({
    
    h = heatmat.list[[2]]
    
    heatmap.2(h$v, 
              xlab="patient",
              Rowv='none',
              trace='non',
              dendrogram='column',
              labRow=h$genes,
              col=c("gray20", "orange"),
              colsep=1:ncol(h$v),
              rowsep=1:nrow(h$v),
              sepcolor="black",
              sepwidth=c(0.01, 0.01),
              RowSideColors=h$cols,
              key=FALSE,
              margins=c(4, 5),
              main=names(heatmat.list)[2]
    )
    legend("bottomleft", legend=h$ordered_groups, fill=pal.default, cex=.7)
    
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  # ==============================================================================
  
  png("../fig/heat_RCC.png", width=800, height=600)
  tryCatch({
    
    h = heatmat.list[[3]]
    
    heatmap.2(h$v, 
              xlab="patient",
              Rowv='none',
              trace='non',
              dendrogram='column',
              labRow=h$genes,
              col=c("gray20", "orange"),
              colsep=1:ncol(h$v),
              rowsep=1:nrow(h$v),
              sepcolor="black",
              sepwidth=c(0.01, 0.01),
              RowSideColors=h$cols,
              key=FALSE,
              margins=c(4, 5),
              main=names(heatmat.list)[3]
    )
    legend("bottomleft", legend=h$ordered_groups, fill=pal.default, cex=.7)
    
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  # ==============================================================================
  
  png("../fig/heat_BLCA.png", width=800, height=600)
  tryCatch({
    
    h = heatmat.list[[4]]
    
    heatmap.2(h$v, 
              xlab="patient",
              Rowv='none',
              trace='non',
              dendrogram='column',
              labRow=h$genes,
              col=c("gray20", "orange"),
              colsep=1:ncol(h$v),
              rowsep=1:nrow(h$v),
              sepcolor="black",
              sepwidth=c(0.01, 0.01),
              RowSideColors=h$cols,
              key=FALSE,
              margins=c(4, 5),
              main=names(heatmat.list)[4]
    )
    legend("bottomleft", legend=h$ordered_groups, fill=pal.default, cex=.7)
    
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  

  
  # ==============================================================================
  x = groups.list$GBM
  
  png("../fig/tree_GBM_1.png", width=600, height=300)
  plot(x$df.node)
  dev.off()
  
  png("../fig/tree_GBM_2.png", width=600, height=300)
  plot(x$w.node.high)
  dev.off()
  
  png("../fig/tree_GBM_3.png", width=600, height=300)
  plot(x$w.node.low)
  dev.off()
  # ==============================================================================
  
  # differential expression analysis 1
  
  source("../code/util.R")
  
  table.list = lapply(names(groups.list), function(x){
    
    y = groups.list[[x]]$df.all
    y$subgroup = paste(
      y$group,
      "(",
      ifelse(y$PD1=="HIGH", "PD1+", "PD1-"), ",",
      ifelse(y$LAG3=="HIGH", "LAG3+", "LAG3-"), ",",
      ifelse(y$TIM3=="HIGH", "TIM3+", "TIM3-"),
      ")",
      sep=""
    )
    
    sel = grep("CD8", names(table(y$subgroup)), value=TRUE)
    
    y = y[which(y$subgroup %in% sel), ]
    
    sel.samples = as.character(y$sample)
    
    mat.exp = mat.list[[x]][, sel.samples]
    
    tables = lapply(sel, function(sel.group){
      
      z = y[, c("sample", "patient", "subgroup")]
      rownames(z) = as.character(z$sample)
      z$group = factor(z$subgroup == sel.group)
      
      mod = model.matrix(~ 0 + group, z)
      
      con = makeContrasts(group=groupTRUE-groupFALSE, levels=mod)
      
      dge = DGEList(counts=mat.exp, group=z$group, remove.zeros=TRUE, genes=rownames(mat.exp))
      keep = rowSums(cpm(dge) > 1) >= 1
      dge = dge[keep, ]
      dge$samples$lib.size = colSums(dge$counts)
      
      dge_n = calcNormFactors(dge)
      
      V = voom(dge_n, mod, plot=FALSE)
      lmfit_voom = lmFit(V, mod)
      lmfit_con = contrasts.fit(lmfit_voom, contrasts=con)
      
      top_table = topTable(eBayes(lmfit_con), coef=1, number=Inf, adjust.method="bon")
      
      top_table = top_table[order(top_table$adj.P.Val, decreasing=FALSE), ]
      
      top_table
    })
    
    names(tables) = sel
    
    tables
    
  })
  names(table.list) = names(groups.list)
  
  save(groups.list, table.list, file="../obj/marker.groups.rda")
  
  table.list.2 = lapply(table.list, function(x){
    y = x[["TUMOR_CD8_EXP(PD1+,LAG3+,TIM3+)"]]
    y
  })
  
  table.list.2 = global_p_adj(table.list.2, method="BH")
  table.list.2 = lapply(table.list.2, function(x) x[order(x$adj.P.Val), ])
  
  kable.list = utils.lapply_i(table.list.2, function(x, i, y){
    temp.table = x[1:20, ]
    rownames(temp.table) = 1:nrow(temp.table)
    list(table=temp.table, caption=sprintf("%s TUMOR_CD8_EXP(PD1+,LAG3+,TIM3+)", y))
  })
  
  threshold = 0.01
  
  tlist = lapply(table.list.2, function(x) x[which(x$adj.P.Val <= threshold), ])
  tlist.genes = lapply(tlist, function(x) rownames(x))
  
  
  # ==============================================================================
  
  png("../fig/venn_1.png", width=800, height=800)
  tryCatch({
    
    venn(tlist.genes)
    
  }, error = function(e){print(e)})
  dev.off()  
  
  # ==============================================================================
  
  utils.kable_vector(sort(Reduce(intersect, tlist.genes)), 10)
  
  WriteXLS(table.list.2, ExcelFileName="../tables/TriplePositive_CD8_TUMORvsAllOther.xls", 
           SheetNames=names(table.list.2))
  
  
  tb = t(sapply(table.list.2, function(x) c(rank=which(rownames(x) == "METRNL"), round(x["METRNL", -1], 6))))
  tb = data.frame(cohort=rownames(tb), tb)
  
  met.list = list()
  met.list$A = tb
  
  # differential expression analysis 2
  
  sel1 = "TUMOR_CD8_EXP"
  sel2 = "PBMC_CD8_ACTIVE"
  
  samples.common = lapply(groups.list, function(y){
    
    df = y$df.all
    sel.high = which(df$PD1 == "HIGH" & df$TIM3 == "HIGH" & df$LAG3 == "HIGH")
    
    z = lapply(c(sel1, sel2), 
               function(x) sort(as.numeric(as.character(df[intersect(sel.high, which(df$group == x)), "patient"]))))
    names(z) = c(sel1, sel2)
    z$common = intersect(z[[1]], z[[2]])
    
    z
    
  })
  

  table.list.common = lapply(1:length(groups.list), function(i){
    
    x = names(groups.list)[i]
    
    df = groups.list[[x]]$df.all
    rownames(df) = df$sample
    mat = mat.list[[x]]
    
    samples = samples.common[[x]]
    
    samples.ids = lapply(1:2, function(j) paste(samples$common, names(samples)[j], sep="_"))
    
    # sapply(samples.ids, function(y) all(y %in% colnames(mat)) )
    
    # sapply(samples.ids, function(y) all(y %in% df$sample) )
    
    samples.ids.v = unlist(samples.ids)
    
    mat = mat[, samples.ids.v]
    df = df[samples.ids.v, ]
    df$group = factor(df$group, levels=c(sel2, sel1))
    df$patient = factor(df$patient, levels=samples$common)
    
    mod = model.matrix(~ 0 + patient + group, df)
    
    colnames(mod) = gsub("group", "", colnames(mod))
    
    con = makeContrasts(TUMOR_CD8_EXP, levels=mod)
    
    dge = DGEList(counts=mat, group=df$group, remove.zeros=TRUE, genes=rownames(mat))
    keep = rowSums(cpm(dge) > 1) >= 1
    
    #print(table(keep))
    
    dge = dge[keep, ]
    dge$samples$lib.size = colSums(dge$counts)
    
    dge_n = calcNormFactors(dge)
    
    V = voom(dge_n, mod, plot=FALSE)
    lmfit_voom = lmFit(V, mod)
    lmfit_con = contrasts.fit(lmfit_voom, contrasts=con)
    
    top_table = topTable(eBayes(lmfit_con), coef=1, number=Inf, adjust.method="bon")
    
    top_table = top_table[order(top_table$adj.P.Val, decreasing=FALSE), ]
    
    top_table  
  })
  names(table.list.common) = names(groups.list)
  
  table.list.common = global_p_adj(table.list.common, method="BH")
  table.list.common = lapply(table.list.common, function(x) x[order(x$adj.P.Val), ])
  
  kable.list = utils.lapply_i(table.list.common, function(x, i, y){
    temp.table = x[1:20, ]
    rownames(temp.table) = 1:nrow(temp.table)
    list(table=temp.table, caption=sprintf("%s PD1+,LAG3+,TIM3+ TUMOR_CD8_EXP vs. PBMC_CD8_ACTIVE with patient", y))
  })
  
  tb = t(sapply(table.list.common, function(x) c(rank=which(rownames(x) == "METRNL"), round(x["METRNL", -1], 6))))
  tb = data.frame(cohort=rownames(tb), tb)
  
  kable(tb, caption="METRNL in PD1+,LAG3+,TIM3+ TUMOR_CD8_EXP vs. PBMC_CD8_ACTIVE with patient", row.names=FALSE)

  met.list$B = tb
  
  WriteXLS(met.list, ExcelFileName="../tables/METRNL_ranks.xls", SheetNames=names(met.list))
  
  WriteXLS(table.list.common, ExcelFileName="../tables/TriplePositive_CD8_TUMORvsPBMC.xls", SheetNames=names(table.list.common))
  
  
  threshold = 0.0001

  tlist = lapply(table.list.common, function(x) x[which(x$adj.P.Val <= threshold), ])
  tlist.genes = lapply(tlist, function(x) rownames(x))
  
  # ==============================================================================
  
  png("../fig/venn_2.png", width=800, height=800)
  tryCatch({
    
    venn(tlist.genes)
    
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


