## to check overlap of ctrl vs AN (human ipsc MG) with HuMicAtlas

setwd("/Users/mridul.chaudhary/work/xingxing_microglia/")

deg=read.csv('/Users/mridul.chaudhary/work/xingxing_microglia/shared_by_xingxing/deg ipsc AN vs CTRL only deseq2.csv')

deg_homo=deg[deg$external_gene_name!="",]
deg_homo=deg_homo[!is.na(deg_homo$padj),]


lfc_cutoff=0.25
padj_cutoff=0.05

# inputs
deg_up <- deg_homo[deg_homo$log2FoldChange > 0 & deg_homo$padj < 0.05, "external_gene_name"]
deg_up <- na.omit(deg_up)
deg_up=as.character(deg_up)
deg_up=deg_up[deg_up != ""]

deg_down <- deg_homo[deg_homo$log2FoldChange < 0 & deg_homo$padj < 0.05, "external_gene_name"]
deg_down <- na.omit(deg_down)
deg_down=as.character(deg_down)
deg_down=deg_down[deg_down != ""]



## Checking deg from iMG in other markers

# with mouse MG from McCarroll

deg_fisher_up_mccarol <- fisher_enrich(deg_up, mus_markers, deg_homo$external_gene_name)
deg_fisher_down_mccarol <- fisher_enrich(deg_down, mus_markers, deg_homo$external_gene_name)


#dolan
deg_fisher_up_dolan <- fisher_enrich(deg_up, dolan_markers, deg_homo$external_gene_name)
deg_fisher_down_dolan <- fisher_enrich(deg_down, dolan_markers, deg_homo$external_gene_name)

#marco
deg_fisher_up_marco <- fisher_enrich(deg_up, marco_markers, deg_homo$external_gene_name)
deg_fisher_down_marco <- fisher_enrich(deg_down, marco_markers, deg_homo$external_gene_name)

#chhatbar
deg_fisher_up_chhatbar <- fisher_enrich(deg_up, chhatbar_markers, deg_homo$external_gene_name)
deg_fisher_down_chhatbar <- fisher_enrich(deg_down, chhatbar_markers, deg_homo$external_gene_name)


## To be done in synapse fed iMG

htsyn=read.csv('/Users/mridul.chaudhary/work/xingxing_microglia/shared_by_xingxing/deg ipsc AN htsyn vs base deseq2.csv')

htsyn_homo=htsyn[htsyn$external_gene_name!="",]
htsyn_homo=htsyn_homo[!is.na(htsyn_homo$padj),]

# inputs
htsyn_up <- htsyn_homo[htsyn_homo$log2FoldChange > 0 & htsyn_homo$padj < 0.05, "external_gene_name"]
htsyn_up <- na.omit(htsyn_up)
htsyn_up=as.character(htsyn_up)
htsyn_up=htsyn_up[htsyn_up != ""]

htsyn_down <- htsyn_homo[htsyn_homo$log2FoldChange < 0 & htsyn_homo$padj < 0.05, "external_gene_name"]
htsyn_down <- na.omit(htsyn_down)
htsyn_down=as.character(htsyn_down)
htsyn_down=htsyn_down[htsyn_down != ""]


# with mouse MG from McCarroll

htsyn_fisher_up_mccarol <- fisher_enrich(htsyn_up, mus_markers, htsyn_homo$external_gene_name)
htsyn_fisher_down_mccarol <- fisher_enrich(htsyn_down, mus_markers, htsyn_homo$external_gene_name)

#dolan
htsyn_fisher_up_dolan <- fisher_enrich(htsyn_up, dolan_markers, htsyn_homo$external_gene_name)
htsyn_fisher_down_dolan <- fisher_enrich(htsyn_down, dolan_markers, htsyn_homo$external_gene_name)

#marco
htsyn_fisher_up_marco <- fisher_enrich(htsyn_up, marco_markers, htsyn_homo$external_gene_name)
htsyn_fisher_down_marco <- fisher_enrich(htsyn_down, marco_markers, htsyn_homo$external_gene_name)

#chhatbar
htsyn_fisher_up_chhatbar <- fisher_enrich(htsyn_up, chhatbar_markers, htsyn_homo$external_gene_name)
htsyn_fisher_down_chhatbar <- fisher_enrich(htsyn_down, chhatbar_markers, htsyn_homo$external_gene_name)
