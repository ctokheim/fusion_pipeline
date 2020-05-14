"""
File: fusion_degron_impact.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Analyzes the impact of fusions on internal degrons
"""
import pandas as pd
import argparse
import csv


def parse_arguments():
    info = 'Analyzes the impact of fusions on internal degrons'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Fusion annotation file')
    parser.add_argument('-m', '--motif-info',
                        type=str, required=True,
                        help='Motif information from motif search')
    parser.add_argument('-d', '--degron-pred',
                        type=str, required=True,
                        help='Degron predictions')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def filter_fusions(fus_df):
    is_utr = fus_df['Fusion_effect'].str.contains('UTR')
    is_out_of_bounds = fus_df['Fusion']
    is_empty_seq = fusions['protein_sequence'].isnull()



def main(opts):
    # read data
    fusions = pd.read_csv(opts['input'], sep='\t')
    motif = pd.read_csv(opts['motif_info'], sep='\t')
    useful_cols = ['Random Forest score', 'geneSymbol', 'DegronID']
    degron = pd.read_csv(opts['degron_pred'], sep='\t', usecols=useful_cols)

    # merge degron info
    motif['DegronID'] = motif['gene']+'_'+motif['standardized_uniprot_id']+'_'+motif['start'].astype(str)+'-'+motif['end'].astype(str)
    degron = pd.merge(degron, motif[['DegronID', 'motif', 'start', 'end']], on='DegronID', how='left')
    # filter to only high scoring motifs
    degron = degron[degron['Random Forest score']>0.6]

    # filter to only fusions with prot sequence
    is_empty_seq = fusions['protein_sequence'].isnull()
    fusions = fusions[~is_empty_seq]

    # iterate over each fusion
    output_list = [['ID', 'type', 'gene', 'degron_id (loss)', 'motif (loss)',
                    'score (loss)', 'degron_id (all)', 'motif (all)', 'score (all)']]
    for ix, row in fusions.iterrows():
        g1, g2 = row["5'_gene"], row["3'_gene"]
        c1, c2 = row['CodonPos1'], row['CodonPos2']

        # get all degrons for those genes
        g1_degrons = degron.loc[degron['geneSymbol']==g1,:]
        g2_degrons = degron.loc[degron['geneSymbol']==g2,:]

        # figure out if fusion impacted degrons
        if len(g1_degrons):
            tmp_did_all = ','.join(g1_degrons['DegronID'].values)
            tmp_motif_all = ','.join(g1_degrons['motif'].values)
            tmp_scores_all = ','.join(g1_degrons['Random Forest score'].astype(str).values)

            # figure out if any degrons are lost
            g1_lost = g1_degrons.loc[g1_degrons['start']>=c1, :]
            if len(g1_lost):
                tmp_did = ','.join(g1_lost['DegronID'].values)
                tmp_motif = ','.join(g1_lost['motif'].values)
                tmp_scores = ','.join(g1_lost['Random Forest score'].astype(str).values)
                output_list.append([row['ID'], "5_prime", g1, tmp_did, tmp_motif, tmp_scores, tmp_did_all, tmp_motif_all, tmp_scores_all])
            else:
                output_list.append([row['ID'], "5_prime", g1, None, None, None, tmp_did_all, tmp_motif_all, tmp_scores_all])
        else:
            output_list.append([row['ID'], "5_prime", g1, None, None, None, None, None, None])

        if len(g2_degrons):
            tmp_did_all = ','.join(g2_degrons['DegronID'].values)
            tmp_motif_all = ','.join(g2_degrons['motif'].values)
            tmp_scores_all = ','.join(g2_degrons['Random Forest score'].astype(str).values)

            # figure out if any degrons are lost
            g2_lost = g2_degrons.loc[g2_degrons['start']>=c1, :]
            if len(g2_lost):
                tmp_did = ','.join(g2_lost['DegronID'].values)
                tmp_motif = ','.join(g2_lost['motif'].values)
                tmp_scores = ','.join(g2_lost['Random Forest score'].astype(str).values)
                output_list.append([row['ID'], "3_prime", g2, tmp_did, tmp_motif, tmp_scores, tmp_did_all, tmp_motif_all, tmp_scores_all])
            else:
                output_list.append([row['ID'], "3_prime", g2, None, None, None, tmp_did_all, tmp_motif_all, tmp_scores_all])
        else:
            output_list.append([row['ID'], "3_prime", g2, None, None, None, None, None, None])


        #else:
            #g1_lost = []
        #if len(g2_degrons):
            #g2_lost = g2_degrons.loc[g2_degrons['end']<=c2, :]
        #else:
            #g2_lost = []

        # append results
        #if len(g1_lost):
            #tmp_did = ','.join(g1_lost['DegronID'].values)
            #tmp_motif = ','.join(g1_lost['motif'].values)
            #tmp_scores = ','.join(g1_lost['Random Forest score'].astype(str).values)
            #output_list.append([row['ID'], "5_prime", g1, tmp_did, tmp_motif, tmp_scores])
        #else:
            #output_list.append([row['ID'], "5_prime", g1, None, None, None])
        #if len(g2_lost):
            #tmp_did = ','.join(g2_lost['DegronID'].values)
            #tmp_motif = ','.join(g2_lost['motif'].values)
            #tmp_scores = ','.join(g2_lost['Random Forest score'].astype(str).values)
            #output_list.append([row['ID'], "3_prime", g2, tmp_did, tmp_motif, tmp_scores])
        #else:
            #output_list.append([row['ID'], "3_prime", g2, None, None, None])

    # save output
    with open(opts['output'], 'w') as whandle:
        mywriter = csv.writer(whandle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(output_list)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


