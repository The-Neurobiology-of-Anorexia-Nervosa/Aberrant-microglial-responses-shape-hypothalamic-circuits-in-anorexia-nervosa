##create marker lists for fisher

library(readxl)



## mouse paper markers
## https://doi.org/10.1016/j.immuni.2018.11.004
mus=read.csv("mouse_MG_scrnaseq_lifespan_McCarrol/Cluster_markers_scRNAseq_MG_mouse_lifespan.txt",
             sep="\t")

mouse_genes <- unique(mus$Gene)


# finding the human symbols for mouse genes using gprofiler
library(gprofiler2)

res <- gorth(
  query = mouse_genes,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)

map <- res[, c("input", "ortholog_name")]
colnames(map) <- c("mouse", "human")
map <- na.omit(map)

#adding human symbol to mouse marker list
#selecting only the ones that mapped to human symbol
mus$homo_gene=map$human[match(mus$Gene,map$mouse)]
mus_homo=mus[!is.na(mus$homo_gene),]

#create mourse marker lists from mccarroll
mus_markers <- split(mus_homo$homo_gene, mus_homo$cluster)
names(mus_markers) <- paste0("MG", names(mus_markers))

names(mus_markers)=
  c("MG1_Arg1+",
    "Proliferative_2a",
    "Proliferative_2b",
    "Proliferative_2c",
    "Metabolically Active",
    "Axon Tract Associated",
    "MG5_Hmox1",
    "Brain Border Macrophage",
    "Homeostatic_7a",
    "Homeostatic_7b",
    "Homeostatic_7c",
    "Inflammatory_IFN_responsive",
    "Injury responsive")




## Dolan (iPSC Microglia in vitro)
## https://www.nature.com/articles/s41590-023-01558-2
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
### different marker sets mouse and human from same paper

# human MG marker set from seurat table 
marco_df=read.csv("prinz_human_top20.txt",sep="\t")

##below mapping was created looking at the cluster plot in paper figure 1b

marco_df[  marco_df$cluster %in% 1:11 | 
             marco_df$cluster %in% 48:63,
           "supercluster"] = "Surveillance"


marco_df[  marco_df$cluster %in% 64 ,
           "supercluster"] = "Neuroprotection"

marco_df[  marco_df$cluster %in% 12:24 | 
             marco_df$cluster %in% 65:94 |
             marco_df$cluster %in% 153 |
             marco_df$cluster %in% 160:163 |
             marco_df$cluster %in% 166:177,
           "supercluster"] = "Phagocytosis"

marco_df[  marco_df$cluster %in% 95:108 | 
             marco_df$cluster %in% 150:151 |
             marco_df$cluster %in% 154:157 |
             marco_df$cluster %in% 178:183 ,
           "supercluster"] = "Inflammation"

marco_df[  marco_df$cluster %in% 25:29 | 
             marco_df$cluster %in% 109:120 |
             marco_df$cluster %in% 158:159 |
             marco_df$cluster %in% 184:188 ,
           "supercluster"] = "Cytokine_production"


marco_df[  marco_df$cluster %in% 30:38 | 
             marco_df$cluster %in% 121:122 |
             marco_df$cluster %in% 152 |
             marco_df$cluster %in% 164:165 |
             marco_df$cluster %in% 189,
           "supercluster"] = "Antigen_presentation"


marco_df[  marco_df$cluster %in% 39:42 | 
             marco_df$cluster %in% 123:130 |
             marco_df$cluster %in% 190,
           "supercluster"] = "Interferon_signature"


marco_df[  marco_df$cluster %in% 43:47 | 
             marco_df$cluster %in% 131:149 |
             marco_df$cluster %in% 191:192,
           "supercluster"] = "Proliferation"


marco_df$phenotype=paste(marco_df$celltype,marco_df$condition,marco_df$supercluster,sep = "_")
unique(marco_df$phenotype)

#non homeostatis microglia
marco_mg_nonhomeo=marco_df[marco_df$celltype == "microglia" &
                             marco_df$condition =="non-homeostasis", ]

marco_markers=split(marco_mg_nonhomeo$gene,marco_mg_nonhomeo$supercluster)
marco_markers <- lapply(marco_markers, unique)

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


