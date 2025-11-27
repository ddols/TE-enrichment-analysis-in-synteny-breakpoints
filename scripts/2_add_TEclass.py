import re
import argparse

def add_family_column(gff_input_file, family_list_file, output_file):
    """
    Reads a GFF, skips headers, and adds a final column with family information
    from a lookup file. Includes enhanced diagnostics for troubleshooting.

    Args:
        gff_input_file (str): Path to the input GFF file.
        family_list_file (str): Path to the 2-column, tab-delimited lookup file.
        output_file (str): Path for the new output file.
    """
    print("--- Starting Process ---")

    try:
        # --- Step 1: Create the lookup dictionary from the family list file ---
        print(f"Reading TE family map from '{family_list_file}'...")
        motif_to_family = {}
        with open(family_list_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    motif, family = line.strip().split('\t')
                    motif_to_family[motif] = family
                except ValueError:
                    print(f"Warning: Skipping malformed line in family list: '{line}'")
        
        if not motif_to_family:
            print("\nError: The family lookup dictionary is empty. Please check your family list file.")
            return

        print(f"Successfully mapped {len(motif_to_family)} motifs.")

        # --- Diagnostic variables ---
        found_count = 0
        unknown_count = 0
        unknown_motifs_found = set()
        motifs_in_gff = set()

        # --- Step 2: Process the GFF and write to the output file ---
        print(f"Processing '{gff_input_file}' and writing to '{output_file}'...")
        with open(gff_input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith("#"):
                    continue
                line_stripped = line.strip()
                if not line_stripped:
                    continue
                
                parts = line_stripped.split('\t')
                if len(parts) < 9:
                    continue
                
                attributes = parts[8]
                family_to_add = "UNKNOWN"

                # --- CORRECTED REGEX ---
                # This pattern finds 'Motif:' and then captures all characters
                # until it hits a quote or a space. It's much more robust.
                match = re.search(r'Motif:([^"\s]+)', attributes)
                
                if match:
                    motif = match.group(1)
                    motifs_in_gff.add(motif) # Add to our diagnostic set
                    
                    if motif in motif_to_family:
                        family_to_add = motif_to_family[motif]
                        found_count += 1
                    else:
                        unknown_count += 1
                        unknown_motifs_found.add(motif)
                else:
                    unknown_count += 1
                
                outfile.write(line_stripped + "\t" + family_to_add + "\n")

        # --- FINAL DIAGNOSTIC REPORT ---
        print("\n--- Process Finished ---")
        print(f"Output file saved to '{output_file}'")
        print("\n--- Diagnostic Summary ---")
        print(f"Successfully mapped motifs: {found_count}")
        print(f"Motifs not found (marked as UNKNOWN): {unknown_count}")
        
        if unknown_count > 0:
            print("\n--- Mismatch Analysis ---")
            print("To fix this, compare the names from the GFF with the names in your family list.")
            
            print("\nExamples of motifs found in your GFF file:")
            for i, motif in enumerate(list(motifs_in_gff)[:5]): # Show first 5 examples
                print(f"  - '{motif}'")
                
            print("\nExamples of motifs from your family list file (the lookup keys):")
            for i, (motif, family) in enumerate(motif_to_family.items()):
                if i >= 5: break # Show first 5 examples
                print(f"  - '{motif}'")
            
            print("\nCheck for subtle differences like parentheses, quotes, or extra characters.")

    except FileNotFoundError as e:
        print(f"Error: A required file was not found. Details: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# --- Argument parsing section ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Adds a final column to a GFF file based on a lookup table, skipping headers."
    )
    parser.add_argument("main_input_file", help="Path to the input GFF file.")
    parser.add_argument("te_family_input_list", help="Path to the 2-column, tab-delimited file mapping motifs to families.")
    parser.add_argument("output_file", help="Path for the new, annotated output file.")
    args = parser.parse_args()
    add_family_column(args.main_input_file, args.te_family_input_list, args.output_file)
