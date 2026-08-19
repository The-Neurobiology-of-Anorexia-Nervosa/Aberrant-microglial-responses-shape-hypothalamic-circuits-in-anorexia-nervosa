# Microglial-state enrichment analyses

This module contains FGSEA analyses comparing
study-derived results with published human and mouse microglial-state
markers.

### R code files workflow description


#### Step 1: MG_markers_from_various_papers.R
Here, we create gene markers sets for different Microglia states obtained from major MG subtype related publications. These gene sets are then later used for for fGSEA tests

#### Step 2: functions.R
Key functions used in rest of the steps. fgsea analysis and plotting results from both

#### Step 3: preparing_DEG_input_files.R
Cleaning up the DEG result files from human induced microglia dataset and checking the overlap of microglia baseline genes (deg set) and synaptosomes fed microglia (htsyn set) with microglia state markers


#### Step 4: fgsea.R
fast gene set enrichment analysis by looking at enrichment of MG state markers against ranked log2FC from the entire DEG data. This takes into account all markers as we check against entire DE table instead of taking overlap with significant DEGs only

#### Step 5: visualize_fgsea_per_marker.R
Plotting fGSEA results


## Directory contents

- `scripts/`: data preparation, Fisher-test, FGSEA and plotting scripts
- `data/external/`: public marker tables used in MG_markers_from_various_papers.R
- `data/processed/`: study-derived DEG tables used as inputs

