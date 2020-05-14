"""
File: make_snvbox_input.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Make SNVBox input file
"""
import pandas as pd
import numpy as np
import argparse

# global dictionary mapping codons to AA
codon_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
               'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
               'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
               'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
               'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
               'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
               'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
               'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
               'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
               'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
               'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*',
               'Splice_Site': 'Splice_Site'}

def parse_arguments():
    info = 'Make SNVBox input file'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Known degron file')
    parser.add_argument('-t', '--transcript',
                        type=str, required=True,
                        help='SNVBox transcript annotation')
    parser.add_argument('-u', '--uniprot',
                        type=str, required=True,
                        help='SNVBox uniprot annotation')
    parser.add_argument('-c', '--codon',
                        type=str, required=True,
                        help='SNVBox codon table annotation')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in data
    degron_df = pd.read_csv(opts['input'], sep='\t')
    degron_df = degron_df.dropna(subset=['standardized_uniprot_id'])
    tx_df = pd.read_csv(opts['transcript'], sep='\t')
    uniprot_df = pd.read_csv(opts['uniprot'], sep='\t')

    # filter to relevant data
    tx_df = tx_df.loc[tx_df['UniProtID'].isin(degron_df['standardized_uniprot_id'].unique()), :].copy()
    # drop extra transcripts
    tx_df = tx_df.sort_values('aaLen', ascending=False).dropna(subset=['RefseqP']).drop_duplicates('UniProtID')

    # filter uniprot pos mapping file
    is_uniprot_id = uniprot_df['Uniprot'].isin(degron_df['standardized_uniprot_id'].unique())
    is_correct_uid = uniprot_df['UID'].isin(tx_df['UID'].unique())
    uniprot_df = uniprot_df.loc[is_uniprot_id & is_correct_uid, :].dropna(subset=['UniprotPos']).copy()

    # create long-form degron data frame
    degron_list = []
    degron_df = degron_df[degron_df['standardized_uniprot_id'].isin(tx_df['UniProtID'].unique())]
    degron_df = pd.merge(degron_df, tx_df[['UniProtID', 'aaLen']],
                         left_on='standardized_uniprot_id', right_on='UniProtID', how='left')
    for ix, row in degron_df.iterrows():
        aa_len = row['aaLen']
        # create data for actually observed degrons
        start, end = row['start'], row['end']
        deg_id = row['gene'] + '_' + row['standardized_uniprot_id'] + '_{}-{}'.format(start, end)
        for pos in range(start, end+1):
            tmp = [deg_id, row['standardized_uniprot_id'], pos]
            degron_list.append(tmp)

    degron_long_df = pd.DataFrame(degron_list, columns=['DegronID', 'UniprotID', 'UniprotPos'])

    # unify col names
    tx_df = tx_df.rename(columns={'UniProtID' :'UniprotID'})
    uniprot_df = uniprot_df.rename(columns={'Uniprot' :'UniprotID'})

    # merge info
    merged = pd.merge(degron_long_df, uniprot_df, on=['UniprotID', 'UniprotPos'])
    merged = pd.merge(merged, tx_df[['UID', 'RefseqP', 'geneSymbol', 'aaLen']], on='UID', how='left')

    # read in codon table
    codon_df = pd.read_csv(opts['codon'], sep='\t')
    codon_df = codon_df.loc[codon_df['UID'].isin(merged['UID'].unique()), :]

    # merge codon information
    merged = pd.merge(merged, codon_df[['UID', 'Pos', 'bases']], on=['UID', 'Pos'], how='left')
    merged['aa'] = merged['bases'].apply(codon_table.get)

    # format mutation column
    merged['mutation'] = merged['aa'] + merged['Pos'].astype(str) + 'A'

    # save output
    merged.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

