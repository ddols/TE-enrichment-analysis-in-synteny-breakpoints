# Transposable Elements (TE) enrichment analysis in synteny breakpoints

# Summary 

As the name suggests, this repository will walk you through the steps towards performing a transposable elements (**TEs**) 
enrichment analysis in a macrosynteny context by reusing [`GENESPACE`](https://github.com/jtlovell/GENESPACE) and [`RepeatMasker`](https://www.repeatmasker.org/) outputs.
It is acknowledged that [TEs may play an important role](https://doi.org/10.1016/j.gene.2012.07.042) in the emergence of chromosomal rearrangements (**CRs**). 
The [outcomes of CRs](https://doi.org/10.1016/j.tree.2010.07.008) are diverse, and they may not always necessarily imply huge rearrangement events such as
fusion or fission events. 

The scripts shared in this repository are meant to perform an exploratory analysis to find candidate TE families that may be enriched or depleted across a distribution of **n*-bp flanking windows surrounding synteny breakpoints as identified by 
`GENESPACE` in a macrosynteny analysis on chromosomes of interest. By comparing the observed number of TEs within these breakpoint regions to a null distribution generated via a permutation test, the analysis identifies TE families that are significantly associated with breaks in genomic conservation.

# Prerequisites

-  `R 3.5.0` or newer. The R script has been tested in `R 4.5.2`.
     -  Install [`BiocManager`](https://www.bioconductor.org/install/), [`GenomicRanges`](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), [`regioneR`](https://bioconductor.org/packages/release/bioc/html/regioneR.html), and `data.table` (`install.packages("data.table")`).
-  Access to the `RepeatMasker` (>= 4.1.2) outputs of the genomes we want to analyze. In particular, to `*.out.gff` and `*.out` files.
-  Access to the `syntenyBlock_coordinates.txt` output file from `GENESPACE`. 

# Step 1 - Preprocessing `RepeatMasker` outputs for the analysis

Contrary to what one could expect, some of `RepeatMasker` output files appear to be tab-delimited, but are not. This can be troublesome. You will find within the `/scripts` directory two scripts labeled by order of usage
to preprocess the output files of interest. I would recommend simplifying the names of the scaffolds to avoid problems with the `R` script that you will run to perform the enrichment analysis. If you do, **make sure** to do the same
for the `syntenicBlock_coordinates.txt` file. Otherwise, the `R` script will fail.

First, we will run `1_format_RM.py` on the `*.out` file of `RepeatMasker`. This script will yield two outputs:
-  A formatted version of the `*.out` file, which is tab-delimited.
-  A list of all motifs (families) in the input `*.out` file and their associated class as defined by `RepeatMasker`.

To run this script, type the following command. The first argument is the `*.out` file of `RepeatMasker`; the second argument, the name of the tab-delimited output; and the third argument, the name of the list file with the families and classes:

```
python 1_format_RM.py <input_RM.out> <output_RM.tab.out> <output_families.out>
```

Then, we will run `2_add_TEclass.py`, which will take the `*.out.gff` file from `RepeatMasker` and the `<output_families.out>` file as inputs, and output the final `gff` table required to run the analysis.
To run the script, type:

```
python 2_add_TEclass.py <input_RM.out.gff> <out_families.out> <final_out.gff>
```

# Step 2 - Run `breakpoints_analysis.R` 

The script that runs the analysis runs pretty automatically. First, we will have to place it in a folder alongside the required input files
that we will set as the working directory. Then, to run the script, we will type the following in the `R` console:

``` r
setwd("path/to/the/folder/containig_the_input_files")

source("breakpoints_analysis.R")
```
However, we will have to specify a series of parameters beforehand. Let's have a look at them:

``` r
# Step 1: Set Parameters
WINDOW_SIZE <- 20000        # Adjust length (in base pairs) of the flanking windows
N_ITERATIONS <- 1000        # Adjust number of permutations
AVG_COUNT_THRESHOLD <- 10   # Threshold to account only those families observed >= 10 times in a window for testing
P_VALUE_ADJ_METHOD <- "BH"  # Apply Benjamini-Hochberg correction for FDR
genespace_file <- "syntenicBlock_coordinates.txt" # Your GENESPACE file
te_gff_file <- "<final_out.gff>"      # Your GFF file after running 2_add_TEclass.py of the sp of interest
target_genome <- "your_Sp_ID"           # Indicate the name of the sp of interest as indicated in GENESPACE file column genome1
TARGET_CHROMOSOME <- "chr_ID"  # Indicate the name of the target chromosome you want to inspect as it appears in GENESPACE file column chr1
set.seed(42)                   # Set a seed to ensure reproducibility
```

As you can see, we can adjust different parameters to specify the length of the flanking regions around the detected synteny breakpoints in `syntenicBloc_coordinates.txt`, the number of permutations to establish our null distribution
of TEs in those regions, or set a threshold to avoid testing exceedingly rare elements. The most important thing is, however, to make sure that the `chr_ID` is the same 
in the `<final_out.gff>` and the `syntenicBlock_coordinates.txt` files. Otherwise, the script will not be able to parse the files and perform the analysis.

Once the analysis is finished, the script will output a `.csv` file with the results of the analysis labeled `breakpoint_enrichment_full_results_[chr_ID]_families.csv`
in the same directory were all the input files are.
