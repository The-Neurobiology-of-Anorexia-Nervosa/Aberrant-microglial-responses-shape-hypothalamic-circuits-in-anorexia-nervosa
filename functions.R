## Fisher test

fisher_enrich <- function(gene_set, markers, background) {
  
  res <- lapply(names(markers), function(ct) {
    
    m <- markers[[ct]]
    
    gene_set <- intersect(gene_set, background)
    m <- intersect(m, background)
    
    overlap_genes <- intersect(gene_set, m)
    
    a <- length(overlap_genes)
    b <- length(m) - a
    c <- length(gene_set) - a
    d <- length(background) - (a + b + c)
    
    mat <- matrix(c(d, b, c, a), nrow = 2)
    
    ft <- fisher.test(mat, alternative = "greater")
    
    data.frame(
      celltype = ct,
      overlap = a,
      overlap_genes = paste(overlap_genes, collapse = ", "),
      unique_in_marker = b,
      unique_in_deg = c,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value
    )
  })
  
  res <- do.call(rbind, res)
  
  res$p_adj <- p.adjust(res$p_value, method = "BH")
  
  return(res)
}


## Visualize fisher tests

fisher_enrich_plot=function(df,subtitle){
    library(viridis)
    library(ggplot2)
    library(ggrepel)
  plot=
    ggplot(df,
           aes(x = -log10(p_adj),
               y = reorder(celltype, -log10(p_adj)),
               color = odds_ratio)) +
    scale_color_gradient(low = "skyblue", high = "red")+
    #scale_color_viridis_c(option = "plasma")+
    geom_point(size = 2) +
    geom_vline(xintercept = -log10(0.05),
               linetype = "dashed",
               linewidth = 0.3,
               color="grey10") +
    geom_text_repel(aes(label = overlap_genes),size = 2.5)+
    theme_bw() +
    facet_wrap(~source,ncol=2)+
    labs(
      title="Fisher enrichment test of DEGs with Microglia state markers",
      subtitle=paste(subtitle),
      x = "-log10(FDR)",
      y = "Cell Type",
      color = "Odds Ratio"
    )
  
  return(plot)
}

## fGSEA

run_fgsea <- function(deg_df, markers_list,
                      gene_col = "homo_gene",
                      stat_col = "logFC") {
  
  library(fgsea)

  # remove NA genes
  deg_df <- deg_df[!is.na(deg_df[[gene_col]]), ]
  
  # keep strongest signal per gene
  deg_df <- deg_df[order(abs(deg_df[[stat_col]]), decreasing = TRUE), ]
  deg_df <- deg_df[!duplicated(deg_df[[gene_col]]), ]
  
  # build ranked vector
  ranks <- deg_df[[stat_col]]
  names(ranks) <- deg_df[[gene_col]]
  ranks <- sort(ranks, decreasing = TRUE)
  
  # run fgsea
  fgsea::fgsea(
    pathways = markers_list,
    stats = ranks,
    nPermSimple = 10000
  )
}

## Visualize fGSEA

fgsea_enrich_plot=function(df,subtitle){
  
    library(viridis)
    library(ggplot2)
    library(ggrepel)
  plot=
    ggplot(df,
           aes(x = -log10(padj),
               y = reorder(pathway, -log10(padj)),
               color = NES)) +
    scale_color_gradient(low = "skyblue", high = "red")+
    #scale_color_viridis_c(option = "plasma")+
    geom_point(size = 2) +
    geom_vline(xintercept = -log10(0.05),
               linetype = "dashed",
               linewidth = 0.3,
               color="grey10") +
    geom_text_repel(
      aes(label = sapply(
        leadingEdge,
        function(x) paste(head(x, 10), collapse = ", ")
      )),
      size = 2
    )+
    theme_bw() +
    facet_wrap(~source,ncol=2)+
    labs(
      title="FGSEA between ranked DESeq fold changes and Microglia state markers",
      subtitle=paste(subtitle),
      x = "-log10(FDR)",
      y = "Microglia state",
      color = "NES"
    )
  
  return(plot)
}

## small function used to join multiple fisher / fgsea results in a single df in visualize scripts
AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}
