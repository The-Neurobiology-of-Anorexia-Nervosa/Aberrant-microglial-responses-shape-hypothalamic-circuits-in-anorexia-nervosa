
# Arc
arc_fgsea_dolan=run_fgsea(deg_df = arc_homo,
                             markers_list = dolan_markers)
arc_fgsea_mccarrol=run_fgsea(deg_df = arc_homo,
                             markers_list = mus_markers)
arc_fgsea_marco=run_fgsea(deg_df = arc_homo,
                          markers_list = marco_markers)
arc_fgsea_chhatbar=run_fgsea(deg_df = arc_homo,
                          markers_list = chhatbar_markers)


# dmh
dmh_fgsea_dolan=run_fgsea(deg_df = dmh_homo,
                          markers_list = dolan_markers)
dmh_fgsea_mccarrol=run_fgsea(deg_df = dmh_homo,
                             markers_list = mus_markers)
dmh_fgsea_marco=run_fgsea(deg_df = dmh_homo,
                          markers_list = marco_markers)
dmh_fgsea_chhatbar=run_fgsea(deg_df = dmh_homo,
                             markers_list = chhatbar_markers)

#deg
deg_fgsea_dolan=run_fgsea(deg_df = deg_homo,
                          markers_list = dolan_markers,
                          gene_col = "external_gene_name",
                          stat_col = "log2FoldChange")

deg_fgsea_mccarrol=run_fgsea(deg_df = deg_homo,
                             markers_list = mus_markers,
                             gene_col = "external_gene_name",
                             stat_col = "log2FoldChange")

deg_fgsea_marco=run_fgsea(deg_df = deg_homo,
                          markers_list = marco_markers,
                          gene_col = "external_gene_name",
                          stat_col = "log2FoldChange")

deg_fgsea_chhatbar=run_fgsea(deg_df = deg_homo,
                          markers_list = chhatbar_markers,
                          gene_col = "external_gene_name",
                          stat_col = "log2FoldChange")


#htsyn
htsyn_fgsea_dolan=run_fgsea(deg_df = htsyn_homo,
                            markers_list = dolan_markers,
                            gene_col = "external_gene_name",
                            stat_col = "log2FoldChange")

htsyn_fgsea_mccarrol=run_fgsea(deg_df = htsyn_homo,
                               markers_list = mus_markers,
                               gene_col = "external_gene_name",
                               stat_col = "log2FoldChange")

htsyn_fgsea_marco=run_fgsea(deg_df = htsyn_homo,
                          markers_list = marco_markers,
                          gene_col = "external_gene_name",
                          stat_col = "log2FoldChange")

htsyn_fgsea_chhatbar=run_fgsea(deg_df = htsyn_homo,
                            markers_list = chhatbar_markers,
                            gene_col = "external_gene_name",
                            stat_col = "log2FoldChange")



#plot
library(ggplot2)



# FDR on X, color NES
ggplot(arc_fgsea_mccarrol,
       aes(x = -log10(padj),
           y = reorder(pathway, -log10(padj)),
           color = NES)) +
  geom_vline(xintercept = -log10(0.05),
             linetype = "dashed",
             linewidth = 0.5) +
  scale_color_gradient(low = "skyblue", high = "red")+
  geom_point(size = 2) +
  geom_text_repel(
    aes(label = sapply(
      leadingEdge,
      function(x) paste(head(x, 10), collapse = ", ")
    )),
    size = 2
  )+
  theme_bw() +
  labs(
    x = "-log10(FDR)",
    y = "Microglia State",
    color = "NES"
  )
