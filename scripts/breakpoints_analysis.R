# ---
# TE FAMILY/SUPERFAMILY Enrichment at Synteny Breakpoints Analysis
# FINAL VERSION: Reads a 10-column GFF, parses full motif names, and adds superfamily info to the output.
# ---

# Step 0: Load ALL Required Packages
cat("Loading required libraries...\n")
library(GenomicRanges) # For GRanges objects and core genomics functions
library(data.table)    # For the fast fread() function
library(regioneR)      # FOR THE PERMUTATION TEST (contains createRandomRegions)
cat("Libraries loaded successfully.\n\n")

# Step 1: Set Parameters
WINDOW_SIZE <- 20000        # Adjust length (in base pairs) of the flanking windows
N_ITERATIONS <- 1000        # Adjust number of permutations
AVG_COUNT_THRESHOLD <- 10   # Threshold to account only those families observed >= 10 times in a window for testing
P_VALUE_ADJ_METHOD <- "BH"  # Apply Benjamini-Hochberg correction for FDR
genespace_file <- "syntenicBlock_coordinates.txt" # Your GENESPACE file
te_gff_file <- "<final_out.gff>"      # Your GFF file after running 2_add_TEclass.py of the sp of interest
target_genome <- "your_sp_ID"   # Indicate the name of the sp of interest as indicated in GENESPACE file column genome1
TARGET_CHROMOSOME <- "chr_ID"  # Indicate the name of the target chromosome you want to inspect as it appears in GENESPACE file column chr1
set.seed(42)

# ---
# Step 2: Load and Prepare Input Data with New Parsing Logic
# ---
cat(sprintf("Loading and preparing data for chromosome %s...\n", TARGET_CHROMOSOME))
synteny_blocks_df <- fread(genespace_file)
synteny_blocks_df$genome1 <- trimws(synteny_blocks_df$genome1)
synteny_blocks_df$chr1 <- trimws(synteny_blocks_df$chr1)
dhepta_blocks_df <- na.omit(synteny_blocks_df[genome1 == target_genome & chr1 == TARGET_CHROMOSOME], cols = c("startBp1", "endBp1"))
cat(sprintf("Found %d synteny blocks for genome '%s' on chromosome '%s' after removing any with NA coordinates.\n", nrow(dhepta_blocks_df), target_genome, TARGET_CHROMOSOME))
if (nrow(dhepta_blocks_df) == 0) {
  stop("Analysis halted: No valid synteny blocks were found matching your criteria.")
}
synteny_gr <- GRanges(seqnames = dhepta_blocks_df$chr1, ranges = IRanges(start = dhepta_blocks_df$startBp1, end = dhepta_blocks_df$endBp1))

# --- Read the 10-column GFF file ---
te_df <- fread(te_gff_file, header = FALSE, sep = "\t", 
               col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "superfamily"))

te_df$seqid <- trimws(te_df$seqid)
te_df_filtered <- te_df[seqid == TARGET_CHROMOSOME]
if (nrow(te_df_filtered) == 0) {
  stop(sprintf("Analysis halted: No TE annotations found for chromosome '%s'.", TARGET_CHROMOSOME))
}
original_te_count <- nrow(te_df_filtered)
te_df_filtered <- na.omit(te_df_filtered, cols = c("start", "end"))
if(nrow(te_df_filtered) < original_te_count) {
  cat(sprintf("INFO: Removed %d TE annotations with NA coordinates.\n", original_te_count - nrow(te_df_filtered)))
}

# --- Capture the full motif name ---

cat("Parsing full TE family names from attributes...\n")
has_family_match <- grepl('Motif:[^"]+', te_df_filtered$attributes)
te_df_filtered$te_family <- "unknown_family"
family_names <- regmatches(
  te_df_filtered$attributes[has_family_match],
  regexpr('(?<=Motif:)[^"]+', te_df_filtered$attributes[has_family_match], perl = TRUE)
)
te_df_filtered$te_family[has_family_match] <- family_names

# --- Create a lookup table to map each family to its superfamily ---
cat("Creating a family-to-superfamily mapping...\n")
family_to_superfamily_map <- unique(te_df_filtered[, .(te_family, superfamily)])
# --- Handle cases where one family might be ambiguously assigned to multiple superfamilies ---
family_to_superfamily_map <- family_to_superfamily_map[, .(superfamily = paste(superfamily, collapse = ";")), by = te_family]

# --- The GRanges object only needs the te_family for the core analysis ---
te_annotations_gr <- GRanges(seqnames = te_df_filtered$seqid, ranges = IRanges(start = te_df_filtered$start, end = te_df_filtered$end), strand = te_df_filtered$strand, te_family = te_df_filtered$te_family)
cat("Successfully created TE GRanges object with parsed families.\n")

# ---
# Step 3: Define Windows and Genomic Universe
# ---
cat("Defining windows and genomic universe...\n")
breakpoint_windows_gr <- c(flank(synteny_gr, width = WINDOW_SIZE, start = TRUE), flank(synteny_gr, width = WINDOW_SIZE, start = FALSE))
chr_length <- max(te_df_filtered$end, na.rm = TRUE)
if (!is.finite(chr_length) || chr_length <= 0) {
  stop("Could not determine a valid chromosome length from the TE annotation file.")
}
genome_gr <- GRanges(seqnames = TARGET_CHROMOSOME, ranges = IRanges(start = 1, end = chr_length))
cat(sprintf("Created genome boundaries for permutation test. Chromosome '%s' length: %d\n", TARGET_CHROMOSOME, chr_length))

# ---
# Step 4 & 5: Permutation Test 
# ---
cat("Calculating observed TE counts at breakpoints...\n")
hits <- findOverlaps(breakpoint_windows_gr, te_annotations_gr)
observed_counts <- table(mcols(te_annotations_gr[subjectHits(hits)])$te_family)
cat("Starting permutation test...\n")
all_te_families <- unique(mcols(te_annotations_gr)$te_family)
random_counts_matrix <- matrix(0, nrow = length(all_te_families), ncol = N_ITERATIONS, dimnames = list(all_te_families, 1:N_ITERATIONS))
n_windows <- length(breakpoint_windows_gr)
cat(sprintf("Will generate %d random windows per iteration.\n", n_windows))
for (i in 1:N_ITERATIONS) {
  if (i %% 100 == 0) cat("  Iteration", i, "of", N_ITERATIONS, "...\n")
  random_windows_gr <- createRandomRegions(nregions = n_windows, length.mean = WINDOW_SIZE, genome = genome_gr)
  if (!is.null(random_windows_gr) && length(random_windows_gr) > 0) {
    random_hits <- findOverlaps(random_windows_gr, te_annotations_gr)
    random_counts_iter <- table(mcols(te_annotations_gr[subjectHits(random_hits)])$te_family)
    random_counts_matrix[names(random_counts_iter), i] <- random_counts_iter
  }
}
cat("Permutation test finished.\n")

# ---
# Steps 6 & 7: Analyze, Format, and Write Final Results
# ---
cat("Analyzing and formatting all results...\n")
results_df <- data.frame(te_family = rownames(random_counts_matrix), observed_count = as.numeric(observed_counts[rownames(random_counts_matrix)]), mean_random_count = rowMeans(random_counts_matrix), stringsAsFactors = FALSE)
results_df$observed_count[is.na(results_df$observed_count)] <- 0
results_to_test <- results_df[results_df$mean_random_count >= AVG_COUNT_THRESHOLD, ]
if (nrow(results_to_test) > 0) {
  cat("Calculating p-values for", nrow(results_to_test), "families that met the threshold.\n")
  results_to_test$p_value <- NA
  for (j in 1:nrow(results_to_test)) {
    family_name <- results_to_test$te_family[j]
    obs <- results_to_test$observed_count[j]
    random_dist <- random_counts_matrix[family_name, ]
    p_upper <- (sum(random_dist >= obs) + 1) / (N_ITERATIONS + 1)
    p_lower <- (sum(random_dist <= obs) + 1) / (N_ITERATIONS + 1)
    results_to_test$p_value[j] <- 2 * min(p_upper, p_lower)
  }
  results_to_test$fdr <- p.adjust(results_to_test$p_value, method = P_VALUE_ADJ_METHOD)
  stats_to_merge <- results_to_test[, c("te_family", "p_value", "fdr")]
  final_results_table <- merge(results_df, stats_to_merge, by = "te_family", all.x = TRUE)
} else {
  cat("No TE families met the average count threshold of", AVG_COUNT_THRESHOLD, "for p-value calculation.\n")
  final_results_table <- results_df
  final_results_table$p_value <- NA
  final_results_table$fdr <- NA
}

# --- Add the superfamily information to the final output table ---
cat("Merging superfamily information into the final results table...\n")
names(family_to_superfamily_map)[names(family_to_superfamily_map) == "te_family"] <- "family"
names(final_results_table)[names(final_results_table) == "te_family"] <- "family"
final_results_table <- merge(final_results_table, family_to_superfamily_map, by = "family", all.x = TRUE)

# --- Final formatting and column ordering ---
names(final_results_table)[names(final_results_table) == "observed_count"] <- "observed"
names(final_results_table)[names(final_results_table) == "mean_random_count"] <- "expected"
final_results_table$enrichment <- final_results_table$observed / (final_results_table$expected + 1e-9)

# --- Include the new 'superfamily' column in the final arrangement ---
final_results_table <- final_results_table[, c("family", "superfamily", "observed", "expected", "enrichment", "p_value", "fdr")]
final_results_table <- final_results_table[order(final_results_table$fdr, final_results_table$p_value, na.last = TRUE), ]

# --- Write the final output file ---
output_file <- paste0("breakpoint_enrichment_full_results_", TARGET_CHROMOSOME, "_families.csv")
write.csv(final_results_table, file = output_file, row.names = FALSE, na = "NA")

cat("\n--- Analysis Complete ---\n")
cat("A full table of results for all", nrow(final_results_table), "families has been saved to:\n", output_file, "\n\n")
cat("--- Top 20 families from the results file ---\n")
print(head(final_results_table, 20))