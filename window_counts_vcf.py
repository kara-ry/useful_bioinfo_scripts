import argparse
import pysam

def parse_chromosome_lengths(lengths_file):
    chromosome_lengths = {}
    with open(lengths_file, 'r') as f:
        for line in f:
            chromosome, _, length = line.strip().split('\t')
            chromosome_lengths[chromosome] = int(length)
    return chromosome_lengths

def parse_bed_file(bed_file):
    bed_data = []
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 3:
                print(f"Error parsing line: {line.strip()}. Skipping...")
                continue
            chromosome, start, end = parts
            bed_data.append((chromosome, int(start), int(end)))
    return bed_data

def parse_vcf_file(vcf_file):
    variants = []
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            chromosome = record.chrom
            position = record.pos
            variants.append((chromosome, position))
    return variants

def count_variants_and_genes(chromosome_lengths, gene_locations, variant_locations, window_size):
    results = []
    for chromosome, length in chromosome_lengths.items():
        for window_start in range(1, length, window_size):
            window_end = min(window_start + window_size - 1, length)
            variant_count = sum(window_start <= variant[1] <= window_end for variant in variant_locations if variant[0] == chromosome)
            gene_count = sum(window_start <= gene[1] <= window_end for gene in gene_locations if gene[0] == chromosome)
            results.append((chromosome, window_start, window_end, variant_count, gene_count))
    return results

def write_output(results, output_file):
    with open(output_file, 'w') as f:
        f.write("Chromosome\tStart\tEnd\tVariant_Count\tGene_Count\n")
        for row in results:
            f.write('\t'.join(map(str, row)) + '\n')

def main(gene_bed_file, variant_vcf_file, lengths_file, output_file, window_size):
    chromosome_lengths = parse_chromosome_lengths(lengths_file)
    gene_locations = parse_bed_file(gene_bed_file)
    variant_locations = parse_vcf_file(variant_vcf_file)
    results = count_variants_and_genes(chromosome_lengths, gene_locations, variant_locations, window_size)
    write_output(results, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count variants and genes within windows on chromosomes")
    parser.add_argument("gene_bed_file", help="Path to the BED file containing gene locations")
    parser.add_argument("variant_vcf_file", help="Path to the VCF file containing variant locations")
    parser.add_argument("lengths_file", help="Path to the file containing chromosome lengths")
    parser.add_argument("output_file", help="Path to the output file")
    parser.add_argument("window_size", type=int, help="Size of the window")
    args = parser.parse_args()

    main(args.gene_bed_file, args.variant_vcf_file, args.lengths_file, args.output_file, args.window_size)

