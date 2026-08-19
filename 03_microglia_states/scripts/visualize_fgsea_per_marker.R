## final hd fgsea plots

library(stringr)
library(ggplot2)
library(ggrepel)


# mouse for marco prinz

mouse_final=AppendMe(c("arc_fgsea_chhatbar","dmh_fgsea_chhatbar"))
mouse_final$source=str_remove(string = mouse_final$source,
                                 pattern = "_fgsea_chhatbar")
mouse_final$source=str_replace(mouse_final$source,pattern = "arc","ARC")
mouse_final$source=str_replace(mouse_final$source,pattern = "dmh","DMH")
mouse_final$source=factor(mouse_final$source,levels = c("ARC","DMH"))

p=fgsea_enrich_plot(mouse_final,"A transcriptomic microglia taxonomy across mouse pathologies")

p

ggsave("mouse_hd2_modified_label_subtitle.tiff",
       plot = p,
       device = "tiff",
       width = 11,
       height = 7,
       units = "in",
       dpi = 600,
       compression = "lzw")


# iMG for dolan

iMG_final=AppendMe(c("deg_fgsea_dolan","htsyn_fgsea_dolan"))
iMG_final$source=str_remove(string = iMG_final$source,
                            pattern = "_fgsea_dolan")
table(iMG_final$source)
iMG_final[iMG_final$source=="deg","source"]="AN vs HC Microglia"
iMG_final[iMG_final$source=="htsyn","source"]="AN vs HC Microglia after phagocytosis of hypothalamic synaptosomes"

#iMG_final$source=factor(iMG_final$source,levels = c("deg","dmh"))
p= fgsea_enrich_plot(iMG_final,"iPSC derived human microglia exposed to brain substrates")

p

ggsave("iMG_hd2_modified_label.tiff",
       plot = p,
       device = "tiff",
       width = 11,
       height = 7,
       units = "in",
       dpi = 600,
       compression = "lzw")
