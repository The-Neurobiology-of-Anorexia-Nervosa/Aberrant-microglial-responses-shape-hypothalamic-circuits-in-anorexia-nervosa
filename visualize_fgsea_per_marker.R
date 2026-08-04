
## Plot fGSEA results for each of the DEGs (iMG and mouse)

# dolan ipsc
dolan_fgsea=AppendMe(c(ls(pattern = "_fgsea_dolan")))
dolan_fgsea$source=str_remove(string = dolan_fgsea$source,
                              pattern = "_fgsea_dolan")
dolan_fgsea$source=factor(dolan_fgsea$source,levels = c("deg","htsyn","arc","dmh"))

fgsea_enrich_plot(dolan_fgsea,"iPSC-derived human microglia states in vitro")

#mccarrol mouse
mccarrol_fgsea=AppendMe(c(ls(pattern = "_fgsea_mccarrol")))
mccarrol_fgsea$source=str_remove(string = mccarrol_fgsea$source,
                                 pattern = "_fgsea_mccarrol")
mccarrol_fgsea$source=factor(mccarrol_fgsea$source,levels = c("deg","htsyn","arc","dmh"))

fgsea_enrich_plot(mccarrol_fgsea,"Microglia from Mouse Lifespan and in the Injured Brain")

#marco
marco_fgsea=AppendMe(c(ls(pattern = "_fgsea_marco")))

marco_fgsea$source=str_remove(string = marco_fgsea$source,
                              pattern = "_fgsea_marco")
marco_fgsea$source=factor(marco_fgsea$source,levels = c("deg","htsyn","arc","dmh"))

fgsea_enrich_plot(marco_fgsea,"microglia taxonomy across mouse and human pathologies")

#chhatbar
chhatbar_fgsea=AppendMe(c(ls(pattern = "_fgsea_chhatbar")))

chhatbar_fgsea$source=str_remove(string = chhatbar_fgsea$source,
                                 pattern = "_fgsea_chhatbar")
chhatbar_fgsea$source=factor(chhatbar_fgsea$source,levels = c("deg","htsyn","arc","dmh"))

fgsea_enrich_plot(chhatbar_fgsea,"Marco Prinz mouse microglia seurat markers")
