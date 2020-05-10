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
                        help='Formatted varions for agfusion')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    df = pd.read_csv(opts['input'], sep='\t')

    # Split the gene names
    df_split = df['Fusion'].str.split('--', expand=True)

    # rename columns
    df_split = df_split.rename(columns={0: 'gene1', 1: 'gene2'})

    # add in break positions
    df_split['gene1 Junction'] = df['Breakpoint1'].str.extract(':([0-9]+):', expand=True)[0]
    df_split['gene2 Junction'] = df['Breakpoint2'].str.extract(':([0-9]+):', expand=True)[0]

    # remove lines with invalid gene symbols
    contains_slash = df_split['gene1'].str.contains('/') | df_split['gene2'].str.contains('/')
    df_split = df_split[~contains_slash]

    # save results
    out_cols = ['gene1', 'gene1 Junction', 'gene2', 'gene2 Junction']
    df_split[out_cols].to_csv(opts['output'], sep='\t', index=False, header=None)
    #import IPython ; IPython.embed() ; raise


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


