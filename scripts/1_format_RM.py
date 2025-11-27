import re
import argparse

def format_repeatmasker_output(input_file, main_output_file, extra_output_file):
    """
    Reads a RepeatMasker-style file, skips the header, converts all space
    separators to tabs, and generates two output files.

    Args:
        input_file (str): Path to the input RepeatMasker file.
        main_output_file (str): Path for the main tab-separated output file.
        extra_output_file (str): Path for the extra file with two specific columns.
    """
    try:
        # --- Step 1: Read all lines from the input file ---
        with open(input_file, 'r') as infile:
            all_lines = infile.readlines()

        # --- Step 2: Remove the first 3 lines ---
        data_lines = all_lines[3:]

        # --- Open output files for writing ---
        with open(main_output_file, 'w') as main_outfile, \
             open(extra_output_file, 'w') as extra_outfile:
            
            # Define and write the new header ONLY to the main output file
            header = [
                "SW_score", "perc_div", "perc_del", "perc_ins", "query_sequence",
                "begin", "end", "left", "strand", "matching_repeat", "repeat_class_family",
                "begin_in_repeat", "end_in_repeat", "left_in_repeat", "ID", "optional_star"
            ]
            main_outfile.write("\t".join(header) + "\n")

            # --- Step 3: Process each data line ---
            for line in data_lines:
                # Remove leading/trailing whitespace
                clean_line = line.strip()

                if not clean_line:
                    continue

                # Replace any sequence of one or more spaces with a single tab
                tab_delimited_line = re.sub(r'\s+', '\t', clean_line)
                
                # Write the fully converted line to the main output file
                main_outfile.write(tab_delimited_line + "\n")

                # Split the newly tabbed line to extract columns for the extra file
                columns = tab_delimited_line.split('\t')
                
                # Ensure the line has enough columns to avoid errors
                if len(columns) >= 11:
                    matching_repeat = columns[9]
                    repeat_class_family = columns[10]
                    extra_outfile.write(f"{matching_repeat}\t{repeat_class_family}\n")
                else:
                    # This helps debug if a line is malformed
                    print(f"Warning: Skipped a line with fewer than 11 columns: '{clean_line}'")

        print(f"Successfully processed '{input_file}'.")
        print(f"  - Main TSV saved to: '{main_output_file}'")
        print(f"  - Extra columns saved to: '{extra_output_file}'")

    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# --- Argument parsing section (no changes needed here) ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a space-aligned RepeatMasker-style file to tab-separated (TSV) format and extract specific columns."
    )

    parser.add_argument(
        "input_file",
        help="Path to the input text file that needs formatting."
    )
    
    parser.add_argument(
        "main_output_file",
        help="Path for the main tab-separated output file."
    )

    parser.add_argument(
        "extra_output_file",
        help="Path for the extra output file containing 'matching_repeat' and 'repeat_class_family'."
    )

    args = parser.parse_args()

    format_repeatmasker_output(args.input_file, args.main_output_file, args.extra_output_file)
