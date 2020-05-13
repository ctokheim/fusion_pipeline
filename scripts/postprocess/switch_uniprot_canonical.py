"""
File: switch_uniprot_canonical.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Change isoform number if canonical
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Change isoform number if canonical'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Input')
    parser.add_argument('-c', '--canonical',
                        type=str, required=True,
                        help='Canonical uniprot protein')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in data
    df = pd.read_csv(opts['input'], sep='\t')
    with open(opts['canonical']) as handle:
        canonical = [l.strip() for l in handle]

    # drop isoform number for canonical
    df['standardized_uniprot_id'] = df['uniprot_id']
    is_canonical = df['uniprot_id'].isin(canonical)
    no_iso_num = df['uniprot_id'].str.split('-', expand=True)[0]
    df.loc[is_canonical, 'standardized_uniprot_id'] = no_iso_num[is_canonical]

    # save results
    df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
