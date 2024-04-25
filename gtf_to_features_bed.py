import sys

def parse_gtf(input_file):
    gene_regions = []
    exon_regions = []
    intron_regions = []

    gene_exon_counts = {}  # Track exon counts for each gene
    gene_intron_counts = {}  # Track intron counts for each gene

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            feature_type = fields[2]
            if feature_type == 'gene':
                gene_regions.append(fields)
            elif feature_type == 'exon':
                gene_id = fields[8].split('gene_id "')[1].split('"')[0]
                if gene_id not in gene_exon_counts:
                    gene_exon_counts[gene_id] = 1
                else:
                    gene_exon_counts[gene_id] += 1
                exon_number = gene_exon_counts[gene_id]
                fields[8] = f"{gene_id}.e{exon_number}"
                exon_regions.append(fields)
            elif feature_type == 'intron':
                gene_id = fields[8].split('gene_id "')[1].split('"')[0]
                if gene_id not in gene_intron_counts:
                    gene_intron_counts[gene_id] = 1
                else:
                    gene_intron_counts[gene_id] += 1
                intron_number = gene_intron_counts[gene_id]
                fields[8] = f"{gene_id}.i{intron_number}"
                intron_regions.append(fields)

    return gene_regions, exon_regions, intron_regions

def write_bed(output_file, regions):
    with open(output_file, 'w') as f:
        for region in regions:
            chrom = region[0]
            start = region[3]
            end = region[4]
            f.write(f'{chrom}\t{start}\t{end}\t{region[8]}\n')

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py input.gtf")
        sys.exit(1)

    input_file = sys.argv[1]
    gene_regions, exon_regions, intron_regions = parse_gtf(input_file)

    write_bed('gene_regions.bed', gene_regions)
    write_bed('exon_regions.bed', exon_regions)
    write_bed('intron_regions.bed', intron_regions)

    print("Bed files generated successfully!")

