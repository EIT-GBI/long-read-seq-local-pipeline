#!/usr/bin/env python3
import argparse
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq

# Syn57 non-recoded codons (codons to exclude)
NON_RECODED = {"TCT", "TCC", "TCA", "TCG", "GCG", "GCA", "TAG"}


def get_recoded_codons_and_bed(genbank_file, bed_file=None):
    '''Extract all the recoded codon positions from genbank file and create a bed file.
    This assumes that the genbank file has misc_feature annotations for recoded codons.
    '''
    if bed_file is None:
        bed_file = genbank_file.replace('.gb', '.bed')

    # Read in the genbank file
    record = SeqIO.read(genbank_file, format='genbank')
    recoded_codons = []

    # Find all the CDS features
    cds_features = [feature for feature in record.features if feature.type == 'CDS']

    # Find recoded annotations (misc_feature with "to" in qualifiers)
    recoding_features = [feature for feature in record.features
                         if feature.type == 'misc_feature'
                         and 'to' in str(feature.qualifiers)]

    # iterate over CDS features and find recoded codons within them.
    # There must be a more efficient way to do this. Something like interval trees?
    for cds in cds_features:
        # iterate over recoding features
        for recoding_feature in recoding_features:
            # Check if the recoding feature is within the CDS
            # That is:
                # the start of the cds is less than or equal to the start of the recoding_feature
                # and
                # the end of the recoding_feature is less than or equal to the end of the cds
            if cds.location.start <= recoding_feature.location.start and recoding_feature.location.end <= cds.location.end:
                # Now determine the strand and sequence
                if cds.location.strand == 1:
                    strand = 'f'
                else:
                    strand = 'r'

                # add all in the codon info dictionary
                codon_info = {
                    'start': recoding_feature.location.start,
                    'end': recoding_feature.location.end,
                    'position': sorted(list(range(recoding_feature.location.start, recoding_feature.location.end))),
                    'strand': strand
                }
                recoded_codons.append(codon_info)

    recoded_codons = sorted(recoded_codons, key=lambda x: x['start'])
    print(f"Total recoded codons found: {len(recoded_codons)}")

    # Let's write all the recoded codons to a bed file
    print(f"Writing recoded codons to bed file: {bed_file}")
    with open(bed_file, 'w') as bed:
        for codon in recoded_codons:
            bed.write(f".\t{codon['start']}\t{codon['end']}\t{codon['strand']}\n")

    return recoded_codons


def parse_mpileup_bases(bases, ref_base):
    '''Let's parse the mpileup base string into a list of actual bases'''
    parsed = []
    i = 0
    while i < len(bases):
        base = bases[i]

        # if it's a reference match
        if base in '.,':
            parsed.append(ref_base)
            i += 1

        # if it's start of a read
        elif base == '^':
            i += 2  # skip the next character

        # if it's end of a read
        elif base == '$':
            i += 1  # just skip it

        # if it's an insertion
        elif base == '+':
            i += 1
            num_str = ''
            while i < len(bases) and bases[i].isdigit():
                num_str += bases[i]
                i += 1
            insert_length = int(num_str)
            i += insert_length  # skip the inserted bases

        # if it's a deletion
        elif base == '-':
            i += 1
            num_str = ''
            while i < len(bases) and bases[i].isdigit():
                num_str += bases[i]
                i += 1
            delete_length = int(num_str)
            i += delete_length  # skip the deleted bases

        # if it's a mismatch
        else:
            parsed.append(base.upper())
            i += 1
    return parsed


def analyze_all_codons_mpileup(bam_file, ref_fasta, bed_file, recoded_codons):
    '''Analyze all codons efficiently using samtools mpileup'''
    # Run mpileup with the bed file as input
    cmd = [
        'samtools', 'mpileup',
        '-f', ref_fasta,
        '-l', bed_file,
        '-Q', '0',  # no filtering for now
        '--output-QNAME',
        bam_file
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running samtools mpileup:", result.stderr)
        return []

    # Cheeky trick to map codon positions to their info
    # Let's create a lookup table -> position-to-codon start.
    # This will help us map each line in the mpileup output to the corresponding codon.
    pos_to_codon_start = {}
    for codon_info in recoded_codons:
        for pos in codon_info['position']:
            pos_to_codon_start[pos] = codon_info['start']

    # Let's create a mapping from codon start position to codon info
    codon_reads = {}
    for codon_info in recoded_codons:
        codon_key = codon_info['start']
        codon_reads[codon_key] = {
            'info': codon_info,
            'reads': {}  # read_name -> base
        }

    # Parse the mpileup output
    for line in result.stdout.strip().split('\n'):
        # ignore empty lines
        if not line:
            continue

        # split the line into parts
        parts = line.split('\t')

        # isolate relevant fields
        pos = int(parts[1]) - 1  # 0-based position
        ref_base = parts[2]
        bases = parts[4]
        qnames = parts[6].split(',')

        if pos in pos_to_codon_start:
            # extract the codon start position
            codon_start = pos_to_codon_start[pos]

            # parse the bases string into actual bases
            parsed_bases = parse_mpileup_bases(bases, ref_base)

            # map each read name to its base at this position
            for base, qname in zip(parsed_bases, qnames):
                # add the read_name to the codon entry if not already present
                if qname not in codon_reads[codon_start]['reads']:
                    codon_reads[codon_start]['reads'][qname] = {}
                codon_reads[codon_start]['reads'][qname][pos] = base

    # analyze each codon to see if it is recoded or not
    results = []
    for codon_key in sorted(codon_reads.keys()):
        codon_data = codon_reads[codon_key]
        codon_info = codon_data['info']
        positions = codon_info['position']
        strand = codon_info['strand']

        total_codons = 0
        non_recoded_codons = 0

        # check each read
        for read_bases in codon_data['reads'].values():
            # check if the read covers all positions of the codon
            if all(pos in read_bases for pos in positions):
                # valid read covering the codon, add to total
                total_codons += 1

                # extract the codon sequence from the read
                codon_seq = ''.join(read_bases[pos] for pos in positions)

                # check to see if it's negative strand
                if strand == 'r':
                    codon_seq = str(Seq(codon_seq).reverse_complement())

                # check if codon is non-recoded
                if codon_seq in NON_RECODED:
                    non_recoded_codons += 1

        frequency = round(non_recoded_codons / total_codons, 2) if total_codons > 0 else 0

        results.append({
            'position': codon_key + 1,
            'depth': total_codons,
            'non_recoded_codons': non_recoded_codons,
            'frequency': frequency
        })

    return results


def process_sample(bam_file, ref_fasta, bed_file, recoded_codons, output_dir, sample_name, genome_length):
    '''This processes a single bam file'''
    results = analyze_all_codons_mpileup(bam_file, ref_fasta, bed_file, recoded_codons)

    if not results:
        print(f"No results for sample {sample_name}")
        return None

    df = pd.DataFrame(results)

    if output_dir and sample_name:
        csv_dir = os.path.join(output_dir, 'csv')
        os.makedirs(csv_dir, exist_ok=True)
        csv_path = os.path.join(csv_dir, f"{sample_name}_recoding_analysis.csv")
        df.to_csv(csv_path, index=False)

    plot_sample(df, output_dir=output_dir, sample_name=sample_name, genome_length=genome_length)

    return df['frequency'].tolist()


def plot_sample(df, output_dir=None, sample_name=None, genome_length=None, window=1):
    """Plot the recoding landscape for a sample (rolling mean)"""
    df = df.sort_values("position").copy()

    plt.figure(figsize=(15, 5))
    sns.scatterplot(data=df, x='position', y='frequency', s=10, alpha=0.5)

    plt.title(f"Recoding Landscape for {sample_name}")
    plt.xlabel("Genome Position")
    plt.ylabel("Frequency of Non-Recoded Codons")
    plt.ylim(0, 1)

    if genome_length:
        plt.xlim(0, genome_length)

    if output_dir and sample_name:
        plots_dir = os.path.join(output_dir, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        output_file = os.path.join(
            plots_dir, f"{sample_name}_recoding_landscape.png"
        )
        plt.savefig(output_file, bbox_inches="tight")
    plt.close()


def plot_summary(positions, avg_frequency, genome_length, output_dir, window=200):
    '''Generate summary plot of average recoding frequency'''

    summary_path = os.path.join()

    plt.figure(figsize=(15, 5))
    df = pd.DataFrame({'position': positions, 'frequency': avg_frequency})
    #df['frequency_smooth'] = df['frequency'].rolling(window=window, center=True).mean()
    sns.scatterplot(x='position', y='frequency', data=df, s=10, alpha=0.5)
    plt.title("Average Recoding Landscape Across Samples")
    plt.xlabel("Genome Position")
    plt.ylabel("Average Frequency of Non-Recoded Codons")
    plt.ylim(0, 1)
    plt.xlim(0, genome_length)
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Generate recoding landscape plots from BAM files')
    parser.add_argument('--genbank', required=True, help='Path to GenBank file with recoded codon annotations')
    parser.add_argument('--bam-dir', help='Directory containing BAM files (default: OUTPUT_DIR/alignment/bam)')
    parser.add_argument('--output-dir', help='Output directory (overrides config)')
    args = parser.parse_args()

    bam_dir = args.bam_dir
    recoding_dir = args.output_dir 

    os.makedirs(recoding_dir, exist_ok=True)

    bed_file = args.genbank.replace('.gb', '.bed')
    print(f"Extracting recoded codons from {args.genbank}")
    recoded_codons = get_recoded_codons_and_bed(args.genbank, bed_file)

    # extract reference fasta from genbank file
    ref_fasta = args.genbank.replace('.gb', '.fasta')
    SeqIO.write(SeqIO.read(args.genbank, 'genbank'), ref_fasta, 'fasta')

    record = SeqIO.read(ref_fasta, 'fasta')
    genome_length = len(record.seq)

    bam_files = sorted([f for f in os.listdir(bam_dir) if f.endswith('.bam')])
    print(f"Processing {len(bam_files)} BAM files from {bam_dir}")

    all_frequencies = []
    positions = [c['start'] + 1 for c in recoded_codons]  # 1-based positions

    for bam_file in bam_files:
        sample_name = bam_file.replace('.sorted.bam', '').replace('.bam', '')
        print(f"  {sample_name}")
        freq = process_sample(
            os.path.join(bam_dir, bam_file), ref_fasta, bed_file, recoded_codons,
            recoding_dir, sample_name, genome_length
        )
        if freq:
            all_frequencies.append(freq)

    if all_frequencies:
        # Calculate average frequencies
        avg_frequencies = np.mean(all_frequencies, axis=0)

        summary_df = pd.DataFrame({
            'position': positions,
            'avg_frequency': avg_frequencies
        })
        summary_df.to_csv(os.path.join(recoding_dir, 'summary_recoding_analysis.csv'), index=False)

        summary_plot_file = os.path.join(recoding_dir, 'summary_recoding_landscape.png')
        plot_summary(positions, avg_frequencies, genome_length, output_dir==recoding_dir)
        print(f"Summary plot saved to {summary_plot_file}")


if __name__ == '__main__':
    main()
