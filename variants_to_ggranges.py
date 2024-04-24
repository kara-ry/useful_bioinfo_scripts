from collections import defaultdict
import sys

# Function to parse VCF file and extract structural variant positions
def parse_vcf(vcf_file):
    sv_positions = defaultdict(list)
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            sv_type = fields[7].split(';')[0].split('=')[1]
            sv_length = int(fields[7].split(';')[1].split('=')[1])
            sv_positions[sv_type].append((chrom, pos, pos + abs(sv_length)))
    return sv_positions

# Function to write GRanges file
def write_granges_file(sv_positions, output_file):
    with open(output_file, 'w') as f:
        f.write("Chromosome\tStart\tEnd\tType\n")
        for sv_type, positions in sv_positions.items():
            for chrom, start, end in positions:
                f.write(f"{chrom}\t{start}\t{end}\t{sv_type}\n")

# Main function
def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_vcf> <output_granges>")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_file = sys.argv[2]
    
    sv_positions = parse_vcf(input_vcf)
    write_granges_file(sv_positions, output_file)

if __name__ == "__main__":
    main()

