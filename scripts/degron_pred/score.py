"""
File: score.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Score motif hits for degron potential
"""
import pandas as pd
import argparse
import pickle


def parse_arguments():
    info = 'Score motif hits for degron potential'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Saved machine learning model in pickle format')
    parser.add_argument('-m', '--motif-info',
                        type=str, required=True,
                        help='Motif information that was used as input into snvbox')
    parser.add_argument('-s', '--snvbox-annot',
                        type=str, required=True,
                        help='Annotated features from SNVBox')
    parser.add_argument('-t', '--train',
                        type=str, required=True,
                        help='Training file for degron prediction')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='output')
    args = parser.parse_args()
    return vars(args)


def read_degron_info(mypath):
    """File containing info about protein length."""
    degron_info = pd.read_csv(mypath, sep='\t')
    degron_info['frac_protein_len'] = degron_info['Pos'] / degron_info['aaLen']
    degron_info['ID'] = degron_info['RefseqP'] + '_' + degron_info['mutation']
    return degron_info

def convert_id_to_zero_based(mylist):
    output_list = []
    for i in mylist:
        str_split = i.split('_')
        start, end = str_split[-1].split('-')
        new_str = '_'.join([str_split[0], str_split[-2]]) + '_' + '{}-{}'.format(int(start)-1, end)
        output_list.append(new_str)
    return output_list

def main(opts):
    # read in data
    df = pd.read_csv(opts['snvbox_annot'], sep='\t', na_values=['None'])
    degron_info = read_degron_info(opts['motif_info'])
    df = pd.merge(df, degron_info[['ID', 'frac_protein_len']], on='ID', how='left')

    # read in known degron instances from literature
    train = pd.read_csv(opts['train'], sep='\t')
    train_ids = train[train['class']=='real']['DegronID'].unique()
    train_ids = convert_id_to_zero_based(train_ids)

    # read in model
    with open(opts['input'], 'rb') as handle:
        model = pickle.load(handle)

    # prepare the feature dataframe and class labels
    feat_cols = [c for c in df.columns[2:]
                if 'DiffProb' not in c and 'ProbMut' not in c
                ]
    feat_df = df.groupby('UID')[feat_cols].mean()
    feat_df = feat_df.fillna(feat_df.mean())

    # score the motif examples for degron potential
    result = feat_df.copy()
    result['Random Forest score'] = model.predict_proba(feat_df)[:,1]
    # merge the score information back into the data frame
    result = pd.merge(result, degron_info[['DegronID', 'geneSymbol']].drop_duplicates('DegronID'),
                      left_index=True, right_on='DegronID', how='left')

    # degrons known in the literature
    result.loc[result['DegronID'].isin(train_ids), 'Random Forest score'] = 1

    # save the results
    result.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


