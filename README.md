# Aberrant-microglial-responses-shape-hypothalamic-circuits-in-anorexia-nervosa
This repository includes code for the 2027 paper "Aberrant microglial responses shape hypothalamic circuits in anorexia nervosa".


### R code files workflow description


#### Step 1: MG_markers_from_various_papers.R
Here, we create gene markers sets for different Microglia states obtained from major MG subtype related publications. These gene sets are then later used for fisher overlap tests and for fGSEA tests

#### Step 2: functions.R
Key functions used in rest of the steps. fisher test, fgsea analysis and plotting results from both

#### Step 3: human_iMG_fisher_overlap.R
Cleaning up the DEG result files from human induced microglia dataset and checking the overlap of microglia baseline genes (deg set) and synaptosomes fed microglia (htsyn set) with microglia state markers

#### Step 4: mouseMG_fisher_overlap.R
Cleaning up the DEG result files from mouse in-vivo microglia data from hypothalamus and checking the overlap of arc region DEGs and dmh region DEGs with microglia state markers

#### Step 5: visualize_fisher_per_markerset.R
Plotting fisher results

#### Step 6: fgsea.R
fast gene set enrichment analysis by looking at enrichment of MG state markers against ranked log2FC from the entire DEG data. This takes into account all markers as we check against entire DE table instead of taking overlap with significant DEGs only

#### Step 7: visualize_fgsea_per_marker.R
Plotting fGSEA results
