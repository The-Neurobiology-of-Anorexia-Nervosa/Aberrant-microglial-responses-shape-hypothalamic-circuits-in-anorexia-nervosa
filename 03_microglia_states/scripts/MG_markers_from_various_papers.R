##create marker lists for fisher
# finding the human symbols for mouse genes using gprofiler
library(gprofiler2)



## Dolan (iPSC Microglia in vitro)
## https://www.nature.com/articles/s41590-023-01558-2
## below file is supplementry table 2
## file link: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41590-023-01558-2/MediaObjects/41590_2023_1558_MOESM3_ESM.csv

dolan=read.csv("dolan_iMG_invitro_markers.csv")
dolan_up=dolan[dolan$avg_logFC > 0, ]
dolan_down=dolan[dolan$avg_logFC < 0, ]

dolan_markers=split(dolan_up$gene,dolan_up$cluster)
names(dolan_markers) <- paste0("MG", names(dolan_markers))
names(dolan_markers)=
  c("Transition Homeo to DAM",
    "Neurodegenerative DAM 1",
    "Antigen Presenting 1",
    "Antigen Presenting 2",
    "Homeostatic",
    "Proliferative 1",
    "Antigen Presenting 3",
    "Neurodegenerative DAM 2",
    "Proliferative 2",
    "Proliferative 3",
    "IFN responsive")


### New paper
### https://doi.org/10.1038/s41590-026-02472-z
### seurat markers from file: https://zenodo.org/records/16938034/files/markers_mouse_myeloid_cells_clusterwise_top20.csv?download=1

# mouse MG marker set from seurat table (fig 5)
##below mapping was created looking at the cluster plot in paper figure 5a

chhatbar_df=read.csv("prinz_mouse_top20.txt",sep="\t")

chhatbar_df[chhatbar_df$cluster %in% 1  | 
              chhatbar_df$cluster %in% 5:13 ,
           "supercluster"] = "Surveillance"

chhatbar_df[  chhatbar_df$cluster %in% 14:18  | 
                chhatbar_df$cluster %in% 55 |
                chhatbar_df$cluster %in% 61:69 ,
              "supercluster"] = "Neuroprotection"

chhatbar_df[  chhatbar_df$cluster %in% 2  | 
                chhatbar_df$cluster %in% 19:33 |
                chhatbar_df$cluster %in% 70 ,
              "supercluster"] = "Phagocytosis"
              
chhatbar_df[  chhatbar_df$cluster %in% 34:38  | 
                chhatbar_df$cluster %in% 56:59 |
                chhatbar_df$cluster %in% 71 ,
              "supercluster"] = "Inflammation"

chhatbar_df[  chhatbar_df$cluster %in% 39,
              "supercluster"] = "Cytokine_production"

chhatbar_df[  chhatbar_df$cluster %in% 40:41  | 
                chhatbar_df$cluster %in% 72:73 ,
              "supercluster"] = "Antigen_presentation"

chhatbar_df[  chhatbar_df$cluster %in% 42:49  | 
                chhatbar_df$cluster %in% 60 |
                chhatbar_df$cluster %in% 74:76 ,
              "supercluster"] = "Interferon_signature"

chhatbar_df[  chhatbar_df$cluster %in% 3:4  | 
                chhatbar_df$cluster %in% 50:54 |
                chhatbar_df$cluster %in% 77:78 ,
              "supercluster"] = "Proliferation"

chhatbar_df$phenotype=paste(chhatbar_df$celltype,chhatbar_df$condition,chhatbar_df$supercluster,sep = "_")
unique(chhatbar_df$phenotype)

res_chhatbar <- gorth(
  query = chhatbar_df$gene,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)

map_chhatbar <- res_chhatbar[, c("input", "ortholog_name")]
colnames(map_chhatbar) <- c("mouse", "human")
map_chhatbar <- na.omit(map_chhatbar)

#adding human symbol to mouse marker list
chhatbar_df$homo_gene=map_chhatbar$human[match(chhatbar_df$gene,map_chhatbar$mouse)]


#selecting only the ones that mapped to human symbol
#microglia only both homeo and non-homeo
chhatbar_mg=chhatbar_df[chhatbar_df$celltype == "microglia", ]
chhatbar_mg=chhatbar_mg[!is.na(chhatbar_mg$homo_gene),]

chhatbar_markers=split(chhatbar_mg$homo_gene,chhatbar_mg$supercluster)
chhatbar_markers <- lapply(chhatbar_markers, unique)


