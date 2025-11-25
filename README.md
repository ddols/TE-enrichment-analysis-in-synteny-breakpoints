# Transposable Elements (TE) enrichment analysis in synteny breakpoints

# Summary 

As the name suggests, this repository will walk you through the steps towards performing a transposable elements (**TEs**) 
enrichment analysis in a macrosynteny context by reusing [`GENESPACE`](https://github.com/jtlovell/GENESPACE) and [`RepeatMasker`](https://www.repeatmasker.org/) outputs.
It is acknowledged that [TEs may play an important role](https://doi.org/10.1016/j.gene.2012.07.042) in the emergence of chromosomal rearrangements (**CRs**). 
The [outcomes of CRs](https://doi.org/10.1016/j.tree.2010.07.008) are diverse, and they may not always necessarily imply huge rearrangement events such as
fusion or fission events. 

The scripts shared in this repository are meant to perform an exploratory analysis to find candidate TE families that may be enriched or depleted across a distribution of **n*-bp flanking windows surrounding synteny breakpoints as identified by 
[GENESPACE] in a macrosynteny analysis. By comparing the observed number of TEs within these breakpoint regions to a null distribution generated via a permutation test, the analysis identifies TE families that are significantly associated with breaks in genomic conservation.

# Prerequisites

-  `R 3.5.0` or newer. The R script has been tested in `R 4.5.2`.
-  Access to the `RepeatMasker` outputs of the genomes we want to analyze. In particular, to `*.out.gff` and `*.out` files.
-  
