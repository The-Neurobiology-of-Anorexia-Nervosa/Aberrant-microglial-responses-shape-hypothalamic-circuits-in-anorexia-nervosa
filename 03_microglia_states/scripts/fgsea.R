
#deg
deg_fgsea_dolan=run_fgsea(deg_df = deg_homo,
                          markers_list = dolan_markers,
                          gene_col = "external_gene_name",
                          stat_col = "log2FoldChange")



#htsyn
htsyn_fgsea_dolan=run_fgsea(deg_df = htsyn_homo,
                            markers_list = dolan_markers,
                            gene_col = "external_gene_name",
                            stat_col = "log2FoldChange")


# Arc

arc_fgsea_chhatbar=run_fgsea(deg_df = arc_homo,
                          markers_list = chhatbar_markers)


# dmh
dmh_fgsea_chhatbar=run_fgsea(deg_df = dmh_homo,
                             markers_list = chhatbar_markers)





