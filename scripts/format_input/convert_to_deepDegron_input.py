"""
File: format_fusion_file.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Format the fusion calls into a format for deepDegron
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Format the fusion calls into a format for deepDegron'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='TCGA fusion calls')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Formatted varions for deepDegron')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    df = pd.read_csv(opts['input'], sep='\t')

    # Split the gene names
    df_split = df['Fusion'].str.split('--', expand=True)

    # rename columns
    df_split = df_split.rename(columns={0: 'first', 1: 'second'})

    # save results
    df_split.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


