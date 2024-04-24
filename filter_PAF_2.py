import sys
import os

def filter_alignments(input_file, target_length, query_length, mapping_quality, alignment_length):
    # Read input file and process each line
    try:
        with open(input_file, 'r') as infile:
            for line in infile:
                line = line.strip().split('\t')
                if len(line) >= 3:
                    paf_file = line[0]
                    target_species = line[1]
                    query_species = line[2]
                    output_filename = f"{target_species}_{query_species}_filtered_alignment_q{mapping_quality}_l{alignment_length}.txt"
                    
                    # Process each PAF file specified in the input file
                    with open(paf_file, 'r') as paf, open(output_filename, 'w') as output:
                        # Write header with input thresholds
                        output.write(f"# Mapping Quality (MapQ) Threshold: {mapping_quality}\n")
                        output.write(f"# Alignment Length Threshold: {alignment_length}\n")
                        output.write("# Query sequence name\tQuery start position\tQuery end position\tTarget sequence name\tTarget start position\tTarget end position\tOrientation\tReference species ID\tTarget species ID\n")
                        
                        for paf_line in paf:
                            fields = paf_line.strip().split('\t')
                            if len(fields) >= 13:
                                query_name = fields[0]
                                query_len = int(fields[1])
                                target_name = fields[5]
                                target_len = int(fields[6])
                                quality = int(fields[11])
                                
                                if query_len < query_length and target_len < target_length and quality >= mapping_quality:
                                    orientation = "+" if fields[4] == "+" else "-"
                                    reference_id = 1  # Example reference species ID (you can modify this as needed)
                                    target_id = 2  # Example target species ID (you can modify this as needed)
                                    
                                    # Write filtered alignments to the output file
                                    output.write(f"{query_name}\t{fields[2]}\t{fields[3]}\t{target_name}\t{fields[7]}\t{fields[8]}\t{orientation}\t{reference_id}\t{target_id}\n")
                                    
                    print(f"Filtered alignments saved to {output_filename}")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 6 or sys.argv[1] == "-h":
        print("Usage: python filter_alignments.py <input_file.txt> <target_length> <query_length> <mapping_quality> <alignment_length>")
        sys.exit(1)

    input_file = sys.argv[1]
    target_length = int(sys.argv[2])
    query_length = int(sys.argv[3])
    mapping_quality = int(sys.argv[4])
    alignment_length = int(sys.argv[5])

    filter_alignments(input_file, target_length, query_length, mapping_quality, alignment_length)

