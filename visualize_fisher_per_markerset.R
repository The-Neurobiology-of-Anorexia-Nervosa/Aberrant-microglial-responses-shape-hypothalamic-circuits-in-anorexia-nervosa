## Priority markers:
## McCarroll -> mouse microglia injury / lifespan
## Dolan -> iPSC Microglia
## Marco Prinz -> mouse and human both

library(stringr)

## For upregulated genes

## dolan ipsc
dolan_fis_up=AppendMe(c(ls(pattern = "_fisher_up_dolan")))
dolan_fis_up$source=str_remove(string = dolan_fis_up$source,
                               pattern = "_fisher_up_dolan")
dolan_fis_up$source=factor(dolan_fis_up$source,levels = c("deg","htsyn","arc","dmh"))

fisher_enrich_plot(dolan_fis_up,"iPSC-derived human microglia states in vitro")


## mccarol MOUSE
mccarol_fis_up=AppendMe(c(ls(pattern = "_fisher_up_mccarol")))
mccarol_fis_up$source=str_remove(string = mccarol_fis_up$source,
                                 pattern = "_fisher_up_mccarol")
mccarol_fis_up$source=factor(mccarol_fis_up$source,levels = c("deg","htsyn","arc","dmh"))

fisher_enrich_plot(mccarol_fis_up,"Microglia from Mouse Lifespan and in the Injured Brain")


## marco
marco_fis_up=AppendMe(c(ls(pattern = "_fisher_up_marco")))
marco_fis_up$source=str_remove(string = marco_fis_up$source,
                               pattern = "_fisher_up_marco")
marco_fis_up$source=factor(marco_fis_up$source,levels = c("deg","htsyn","arc","dmh"))

fisher_enrich_plot(marco_fis_up,"upregulated genes in human microglia markers from seurat results")


#chhatbar

chhatbar_fis_up=AppendMe(c(ls(pattern = "_fisher_up_chhatbar")))
chhatbar_fis_up$source=str_remove(string = chhatbar_fis_up$source,
                                  pattern = "_fisher_up_chhatbar")
chhatbar_fis_up$source=factor(chhatbar_fis_up$source,levels = c("deg","htsyn","arc","dmh"))

fisher_enrich_plot(chhatbar_fis_up,"upregulated genes vs Marco Prinz mouse microglia seurat markers")



### For downregulated DEGs


## dolan ipsc
dolan_fis_down=AppendMe(c(ls(pattern = "_fisher_down_dolan")))
dolan_fis_down$source=str_remove(string = dolan_fis_down$source,
                                 pattern = "_fisher_down_dolan")
dolan_fis_down$source=factor(dolan_fis_down$source,levels = c("deg","htsyn","arc","dmh"))

fisher_enrich_plot(dolan_fis_down,"iPSC-derived human microglia states in vitro")


## mccarol MOUSE
mccarol_fis_down=AppendMe(c(ls(pattern = "_fisher_down_mccarol")))
mccarol_fis_down$source=str_remove(string = mccarol_fis_down$source,
                                   pattern = "_fisher_down_mccarol")
mccarol_fis_down$source=factor(mccarol_fis_down$source,levels = c("deg","htsyn","arc","dmh"))

fisher_enrich_plot(mccarol_fis_down,"Microglia from Mouse Lifespan and in the Injured Brain")

## marco
marco_fis_down=AppendMe(c(ls(pattern = "_fisher_down_marco")))
marco_fis_down$source=str_remove(string = marco_fis_down$source,
                                 pattern = "_fisher_down_marco")
marco_fis_down$source=factor(marco_fis_down$source,levels = c("deg","htsyn","arc","dmh"))

fisher_enrich_plot(marco_fis_down,"downregulated genes in human microglia markers from seurat results")


#chhatbar

chhatbar_fis_down=AppendMe(c(ls(pattern = "_fisher_down_chhatbar")))
chhatbar_fis_down$source=str_remove(string = chhatbar_fis_down$source,
                                    pattern = "_fisher_down_chhatbar")
chhatbar_fis_down$source=factor(chhatbar_fis_down$source,levels = c("deg","htsyn","arc","dmh"))

fisher_enrich_plot(chhatbar_fis_down,"downregulated genes vs Marco Prinz mouse microglia seurat markers")