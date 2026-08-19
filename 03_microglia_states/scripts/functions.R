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
    nPermSimple = 1000000
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
    scale_color_gradient2(
      low = "blue",
      mid = "grey80",
      high = "red",
      midpoint = 0
    )+           
    #scale_color_gradient(low = "skyblue", high = "red")+
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
      size = 2,
      nudge_y = 0.15,
      point.padding = 0.3,
      box.padding = 0.4,
      min.segment.length = 0
    ) +
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
