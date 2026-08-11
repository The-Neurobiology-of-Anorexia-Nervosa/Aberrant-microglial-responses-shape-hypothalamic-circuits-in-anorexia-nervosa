## to check overlap of ctrl vs AN (human ipsc MG) with HuMicAtlas


deg=read.csv('deg ipsc AN vs CTRL only deseq2.csv')

deg_homo=deg[deg$external_gene_name!="",]
deg_homo=deg_homo[!is.na(deg_homo$padj),]


## To be done in synapse fed iMG

htsyn=read.csv('deg ipsc AN htsyn vs base deseq2.csv')

htsyn_homo=htsyn[htsyn$external_gene_name!="",]
htsyn_homo=htsyn_homo[!is.na(htsyn_homo$padj),]

##To check DEG from anx/anx mouse ARC/DMH MG vs ctrl with ipsc & mouse_McCarroll

arc=read.csv("DEGs_IBA1_ARC_GeoMX_20250617.csv")
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


## same stuff for DMH
dmh=read.csv("DEGS_IBA1_DMH_GeoMX_20250711.csv")

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
