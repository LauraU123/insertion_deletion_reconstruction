import argparse
from Bio import SeqIO
import pandas as pd

"""Keeping only insertions and deletions which have an effect on translation """

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Reconstructing all branch sequences, cutting out gene of interest and translating",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--alignedgene", required=True, help="aligned sequences from specific gene")
    parser.add_argument("--indels", required=True, help="tsv file with insertions and deletions")
    parser.add_argument("--output", help="tsv output file")
    args = parser.parse_args()

    tsv_input = pd.read_csv(args.indels, sep="\t")
    relevant_ids = []
    aligned_gene_sequences = SeqIO.parse(args.alignedgene, "fasta")
    for entry in aligned_gene_sequences:
        if 'X' in entry.seq:
            relevant_ids.append(entry.id)

    relevant_indels = tsv_input.loc[tsv_input['index'].isin(relevant_ids)]
    relevant_indels.to_csv(args.output, sep="\t")

