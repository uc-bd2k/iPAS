# iPAS R package

## Overview

iPAS (integrative Pathway Activity Signatures) is a method for detecting changes in 
signaling pathway activity via differential expression signatures. Our method uses 
the LINCS (Library of Integrative Network-Based Cellular Signatures) library of 
perturbation signatures to construct transcriptomic signatures that can detect 
alterations in activity of signaling pathways. We integrate transcriptomic response 
across a panel of 12 diverse cell lines using machine learning methods to create pathway signatures 
that are strong predictors of pathway activity in a variety of contexts.

## Installation

``` r
install.packages("remotes")
remotes::install_github("NicholasClark/iPAS")
```

## Usage

See the [vignette](https://nicholasclark.github.io/iPAS/index.html) for an example of usage.

The main function, *iPAS_enrich*, takes in a matrix of query signatures (column names should 
be condition/sample names, row names should be entrez gene ids) and outputs an object with 
z-scores and empirical p-values (computed by a permutation test) for each pathway.

There are functions to view the output as bar charts (*iPAS_bar*) of top pathways, heatmaps 
(*iPAS_heatmap*) of all pathway scores, or a density plot (*iPAS_density* and *iPAS_density_facet*) 
of the z-score for a particular pathway versus the null distribution of permutation scores for that pathway.

## Publication

Our manuscript is in preparation:

<b>Clark et al., <i>Integrative signatures of signaling pathway response 
increase robustness and accuracy of pathway predictions</i></b>

## Related work

This work builds on previous work published in Bioinformatics (2020):

<b>Ren et al. <i>Predicting mechanism of action of cellular perturbations with pathway activity signatures"</i></b>
https://doi.org/10.1093/bioinformatics/btaa590
