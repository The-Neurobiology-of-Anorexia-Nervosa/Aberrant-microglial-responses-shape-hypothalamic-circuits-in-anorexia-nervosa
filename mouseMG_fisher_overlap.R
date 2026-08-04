##To check DEG from anx/anx mouse ARC/DMH MG vs ctrl with ipsc & mouse_McCarroll

arc=read.csv("shared_by_xingxing/DEGs_IBA1_ARC_GeoMX_20250617.csv")
arc_up=arc[arc$logFC > 0 & arc$adj.P.Val < 0.1,1]
arc_down=arc[arc$logFC < 0 & arc$adj.P.Val < 0.1,1]


res_arc <- gorth(
  query = unique(arc$X),
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)
head(res)
map_arc <- res_arc[, c("input", "ortholog_name")]
colnames(map_arc) <- c("mouse", "human")
map_arc <- na.omit(map_arc)


arc$homo_gene=map_arc$human[match(arc$X,map_arc$mouse)]

arc_homo=arc[!is.na(arc$homo_gene),]

arc_up=arc_homo[arc_homo$logFC > 0 & arc_homo$adj.P.Val < 0.1,"homo_gene"]
arc_down=arc_homo[arc_homo$logFC < 0 & arc_homo$adj.P.Val < 0.1,"homo_gene"]

##checking enrichmwnt

# with mouse MG from McCarroll

arc_fisher_up_mccarol <- fisher_enrich(arc_up, mus_markers, arc_homo$homo_gene)
arc_fisher_down_mccarol <- fisher_enrich(arc_down, mus_markers, arc_homo$homo_gene)


#dolan
arc_fisher_up_dolan <- fisher_enrich(arc_up, dolan_markers, arc_homo$homo_gene)
arc_fisher_down_dolan <- fisher_enrich(arc_down, dolan_markers, arc_homo$homo_gene)

#marco
arc_fisher_up_marco <- fisher_enrich(arc_up, marco_markers, arc_homo$homo_gene)
arc_fisher_down_marco <- fisher_enrich(arc_down, marco_markers, arc_homo$homo_gene)

#chhatbar
arc_fisher_up_chhatbar <- fisher_enrich(arc_up, chhatbar_markers, arc_homo$homo_gene)
arc_fisher_down_chhatbar <- fisher_enrich(arc_down, chhatbar_markers, arc_homo$homo_gene)


## same stuff for DMH
dmh=read.csv("shared_by_xingxing/DEGS_IBA1_DMH_GeoMX_20250711.csv")

dmh_up=dmh[dmh$logFC > 0 & dmh$adj.P.Val < 0.1,1]
dmh_down=dmh[dmh$logFC < 0 & dmh$adj.P.Val < 0.1,1]


res_dmh <- gorth(
  query = unique(dmh$X),
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)

map_dmh <- res_dmh[, c("input", "ortholog_name")]
colnames(map_dmh) <- c("mouse", "human")
map_dmh <- na.omit(map_dmh)


dmh$homo_gene=map_dmh$human[match(dmh$X,map_dmh$mouse)]

dmh_homo=dmh[!is.na(dmh$homo_gene),]

dmh_up=dmh_homo[dmh_homo$logFC > 0 & dmh_homo$adj.P.Val < 0.1,"homo_gene"]
dmh_down=dmh_homo[dmh_homo$logFC < 0 & dmh_homo$adj.P.Val < 0.2,"homo_gene"]

##checking enrichmwnt

# with mouse MG from McCarroll

dmh_fisher_up_mccarol <- fisher_enrich(dmh_up, mus_markers, dmh_homo$homo_gene)
dmh_fisher_down_mccarol <- fisher_enrich(dmh_down, mus_markers, dmh_homo$homo_gene)

#dolan
dmh_fisher_up_dolan <- fisher_enrich(dmh_up, dolan_markers, dmh_homo$homo_gene)
dmh_fisher_down_dolan <- fisher_enrich(dmh_down, dolan_markers, dmh_homo$homo_gene)

#marco
dmh_fisher_up_marco <- fisher_enrich(dmh_up, marco_markers, dmh_homo$homo_gene)
dmh_fisher_down_marco <- fisher_enrich(dmh_down, marco_markers, dmh_homo$homo_gene)

#chhatbar
dmh_fisher_up_chhatbar <- fisher_enrich(dmh_up, chhatbar_markers, dmh_homo$homo_gene)
dmh_fisher_down_chhatbar <- fisher_enrich(dmh_down, chhatbar_markers, dmh_homo$homo_gene)




## between human deg vs ARV/ DMH mouse deg from xing xing data

dmh_fisher_up_ipsc <- fisher_enrich(dmh_up, list("deg_up_AN"=deg_up), dmh_homo$homo_gene)
dmh_fisher_down_ipsc <- fisher_enrich(dmh_down, list("deg_down_AN"=deg_down), dmh_homo$homo_gene)

arc_fisher_up_ipsc <- fisher_enrich(arc_up, list("deg_up_AN"=deg_up), arc_homo$homo_gene)
arc_fisher_down_ipsc <- fisher_enrich(arc_down, list("deg_down_AN"=deg_down), arc_homo$homo_gene)

## if up ones are in down
dmh_fisher_up_ipsc_down <- fisher_enrich(dmh_up, list("deg_down_AN"=deg_down), dmh_homo$homo_gene)
dmh_fisher_down_ipsc_up <- fisher_enrich(dmh_down, list("deg_up_AN"=deg_up), dmh_homo$homo_gene)

arc_fisher_up_ipsc_down <- fisher_enrich(arc_up, list("deg_down_AN"=deg_down), arc_homo$homo_gene)
arc_fisher_down_ipsc_up <- fisher_enrich(arc_down, list("deg_up_AN"=deg_up), arc_homo$homo_gene)
