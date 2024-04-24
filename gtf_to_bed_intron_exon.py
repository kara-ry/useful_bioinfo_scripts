import sys

def gtf_to_bed(gtf_file):
    bed_lines = []

    gene_info = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            feature_type = fields[2]
            if feature_type not in ('exon', 'intron', 'gene'):
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])
            gene_id = None
            transcript_id = None

            for info in fields[8].split(';'):
                if 'gene_id' in info:
                    gene_id = info.split('"')[1]
                elif 'transcript_id' in info:
                    transcript_id = info.split('"')[1]

            if gene_id is None or transcript_id is None:
                continue

            if gene_id not in gene_info:
                gene_info[gene_id] = {'chrom': chrom, 'start': start, 'end': end, 'exons': [], 'introns': []}
            else:
                if start < gene_info[gene_id]['start']:
                    gene_info[gene_id]['start'] = start
                if end > gene_info[gene_id]['end']:
                    gene_info[gene_id]['end'] = end

            if feature_type == 'exon':
                exon_id = f"{gene_id}.e{len(gene_info[gene_id]['exons']) + 1}"
                gene_info[gene_id]['exons'].append((start, end, exon_id))
            elif feature_type == 'intron':
                intron_id = f"{gene_id}.i{len(gene_info[gene_id]['introns']) + 1}"
                gene_info[gene_id]['introns'].append((start, end, intron_id))

    for gene_id, info in gene_info.items():
        bed_lines.append(f"{info['chrom']}\t{info['start']}\t{info['end']}\t{gene_id}\t0\t.\t{info['start']}\t{info['end']}\t0\t1\t{info['end'] - info['start']}\t0\n")
        for exon_start, exon_end, exon_id in info['exons']:
            bed_lines.append(f"{info['chrom']}\t{exon_start}\t{exon_end}\t{exon_id}\t0\t.\t{info['start']}\t{info['end']}\t0\t1\t{exon_end - exon_start}\t0\n")
        for intron_start, intron_end, intron_id in info['introns']:
            bed_lines.append(f"{info['chrom']}\t{intron_start}\t{intron_end}\t{intron_id}\t0\t.\t{info['start']}\t{info['end']}\t0\t1\t{intron_end - intron_start}\t0\n")

    return bed_lines

def write_bed_file(bed_lines, output_file):
    with open(output_file, 'w') as f:
        f.write("track name=features\n")
        for line in bed_lines:
            f.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python gtf_to_bed_intron_exon.py input.gtf output.bed")
        sys.exit(1)

    input_gtf = sys.argv[1]
    output_bed = sys.argv[2]

    bed_lines = gtf_to_bed(input_gtf)
    write_bed_file(bed_lines, output_bed)

    print(f"Conversion complete. BED file saved as {output_bed}")

