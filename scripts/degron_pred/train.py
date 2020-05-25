"""
File: train.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Train a model to predict degron potential
"""
import pandas as pd
import argparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, cross_val_predict, learning_curve
from sklearn.metrics import roc_auc_score, roc_curve
import pickle


def parse_arguments():
    info = 'Train a model to predict degron potential'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-m', '--motif-info',
                        type=str, required=True,
                        help='Motif information that was used as input into snvbox')
    parser.add_argument('-s', '--snvbox-annot',
                        type=str, required=True,
                        help='Annotated features from SNVBox')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Save the trained machine learning model as a pickle file')
    args = parser.parse_args()
    return vars(args)


def read_degron_info(mypath):
    """File containing info about protein length."""
    degron_info = pd.read_csv(mypath, sep='\t')
    degron_info['frac_protein_len'] = degron_info['Pos'] / degron_info['aaLen']
    degron_info['ID'] = degron_info['RefseqP'] + '_' + degron_info['mutation']
    return degron_info


def main(opts):
    # read in data
    df = pd.read_csv(opts['snvbox_annot'], sep='\t', na_values=['None'])
    degron_info = read_degron_info(opts['motif_info'])
    df = pd.merge(df, degron_info[['ID', 'frac_protein_len']], on='ID', how='left')

    # prepare the feature dataframe and class labels
    feat_cols = [c for c in df.columns[2:]
                if 'DiffProb' not in c and 'ProbMut' not in c
                ]
    feat_df = df.groupby('UID')[feat_cols].mean()
    feat_df = feat_df.fillna(feat_df.mean())
    y = [0 if ix.startswith('SIM_') else 1 for ix in feat_df.index.values]

    # train rf in cross-val to eval performance
    rf = RandomForestClassifier(n_estimators=1000, oob_score=True)
    myscores = cross_val_predict(rf, feat_df, y, cv=20, method='predict_proba')[:,1]
    y = pd.Series(y)
    print('auROC={:.2g}'.format(roc_auc_score(y, myscores)))

    # train on full data
    model = rf.fit(feat_df, y)

    # save model as a pickle file
    with open(opts['output'], 'wb') as output:  # Overwrites any existing file.
        pickle.dump(model, output, pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


