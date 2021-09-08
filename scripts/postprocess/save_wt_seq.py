"""
File: save_wt_seq.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Saves the WT sequence of fusions
"""
import pandas as pd
import argparse
from pyensembl import EnsemblRelease


def parse_arguments():
    info = 'Saves the WT sequence of fusions'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Merge afgusion results')
    parser.add_argument('-e', '--ensembl-release',
                        type=int, default=95,
                        help='Ensembl release version')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='File containing WT protein sequence')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # load ensembl db
    data = EnsemblRelease(opts['ensembl_release'])

    # read in fusion file
    df = pd.read_csv(opts['input'], sep='\t')

    output_list = []
    for ix, row in df.iterrows():
        # extract gene / tx
        gene5 = row["5'_gene"]
        gene3 = row["3'_gene"]
        tx_id5 = row["5'_transcript"]
        tx_id3 = row["3'_transcript"]

        # fetch prot sequence
        tx5 = data.transcript_by_id(tx_id5)
        tx3 = data.transcript_by_id(tx_id3)
        prot5 = tx5.protein_id
        prot3 = tx3.protein_id
        prot_seq5 = tx5.protein_sequence
        prot_seq3 = tx3.protein_sequence

        # append output
        output_list.append([gene5, tx_id5, prot5, prot_seq5])
        output_list.append([gene3, tx_id3, prot3, prot_seq3])

    # save output
    output_df = pd.DataFrame(output_list, columns=['gene', 'transcript_id', 'protein_id', 'protein_sequence'])
    output_df.drop_duplicates(subset=['gene', 'transcript_id', 'protein_id']).to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


