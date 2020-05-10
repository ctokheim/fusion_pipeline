"""
File: convert_gene_to_transcript.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Convert gene symbols to canonical transcript
"""
import pandas as pd
import argparse
from pyensembl import ensembl_grch37, EnsemblRelease
import csv



def parse_arguments():
    info = 'Convert gene symbols to canonical transcript'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Agfusion input file with gene symbols')
    parser.add_argument('-t', '--transcript',
                        type=str, required=True,
                        help='Canonical transcript file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Agfusion file with transcript IDs')
    args = parser.parse_args()
    return vars(args)


def pick_transcript(data, gene, pos, mane_tx):
    """Use MANE canonical transcript if available, otherwise largest transcript by CDS."""
    # figure out all of the tx ids for a gene
    try:
        tx_ids = data.transcript_ids_of_gene_name(gene)
    except ValueError:
        return ''

    # iterate through each tx
    tx_dict = {}
    for t in tx_ids:
        # get transcript obj
        tx = data.transcript_by_id(t)

        # figure out whether the junction actually hits the tx
        offset_pos = None
        try:
            offset_pos = tx.spliced_offset(pos)
        except ValueError:
            pass
        if offset_pos is None: continue

        # add tx to the dict
        if tx.protein_sequence:
            tx_dict[t] = len(tx.protein_sequence)

    # return the appropriate tx
    if len(tx_dict) == 0:
        return ''
    elif mane_tx in tx_dict:
        return mane_tx
    else:
        sorted_txs = sorted(tx_dict.items(), key=lambda x: (x[1], x[0]), reverse=True)
        return sorted_txs[0][0]


def main(opts):
    # read in data
    df = pd.read_csv(opts['input'], sep='\t', header=None, names=['Gene1', 'Break1', 'Gene2', 'Break2'])
    tx = pd.read_csv(opts['transcript'], sep='\t')

    # merge in transcript
    rename_dict = {'symbol': 'Gene1', 'Ensembl_nuc': 'Transcript1'}
    df = pd.merge(df, tx.rename(columns=rename_dict)[['Gene1', 'Transcript1']],
                  on='Gene1', how='left')
    rename_dict = {'symbol': 'Gene2', 'Ensembl_nuc': 'Transcript2'}
    df = pd.merge(df, tx.rename(columns=rename_dict)[['Gene2', 'Transcript2']],
                  on='Gene2', how='left')

    # remove the transcript version number
    tmp1 = df['Transcript1'].str.split('.', expand=True)[0]
    tmp2 = df['Transcript2'].str.split('.', expand=True)[0]
    df['Transcript1'] = tmp1
    df['Transcript2'] = tmp2

    # Check transcripts
    data = EnsemblRelease(95)

    # figure out the appropriate canonical tx
    output_list = []
    for ix, row in df.iterrows():
        # get transcript
        mytx1 = pick_transcript(data, row['Gene1'], row['Break1'], row['Transcript1'])
        mytx2 = pick_transcript(data, row['Gene2'], row['Break2'], row['Transcript2'])

        # fill in whether to replace gene ID with transcript ID
        mygene1 = mytx1 if mytx1 else row['Gene1']
        mygene2 = mytx2 if mytx2 else row['Gene2']

        # append output
        output_list.append([mygene1, row['Break1'], mygene2, row['Break2']])

    # save output file
    with open(opts['output'], 'w') as whandle:
        mywriter = csv.writer(whandle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(output_list)

    # replace gene symbols
    #is_tx1_null = df['Transcript1'].isnull()
    #df.loc[~is_tx1_null, 'Gene1'] = df.loc[~is_tx1_null, 'Transcript1']
    #is_tx2_null = df['Transcript2'].isnull()
    #df.loc[~is_tx2_null, 'Gene2'] = df.loc[~is_tx2_null, 'Transcript2']
    # save file
    #df[['Gene1', 'Break1', 'Gene2', 'Break2']].to_csv(opts['output'], sep='\t', index=False, header=None)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


