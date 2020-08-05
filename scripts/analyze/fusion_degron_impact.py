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


def sum_rf_score(mydf):
    ivls = list(mydf[['start', 'end']].apply(lambda x: pd.Interval(x['start'], x['end'], closed='both'), axis=1).values)
    mysum = 0
    prev_intervals = []
    for ix, ivl in enumerate(ivls):
        has_overlap = any([ivl.overlaps(p) for p in prev_intervals])
        if not has_overlap:
            mysum += mydf['Random Forest score'].iloc[ix]
        prev_intervals.append(ivl)
    return mysum


def main(opts):
    # read data
    fusions = pd.read_csv(opts['input'], sep='\t')
    motif = pd.read_csv(opts['motif_info'], sep='\t')
    useful_cols = ['Random Forest score', 'geneSymbol', 'DegronID']
    degron = pd.read_csv(opts['degron_pred'], sep='\t', usecols=useful_cols)

    # process motif info to remove duplicates from similar motifs
    motif['DegronID'] = motif['gene']+'_'+motif['standardized_uniprot_id']+'_'+motif['start'].astype(str)+'-'+motif['end'].astype(str)
    motif = motif.drop_duplicates(subset=['DegronID'])

    # merge degron info
    degron = pd.merge(degron, motif[['DegronID', 'motif', 'start', 'end']], on='DegronID', how='left')
    # filter to only high scoring motifs
    degron = degron[degron['Random Forest score']>0.6]

    # filter to only fusions with prot sequence
    is_empty_seq = fusions['protein_sequence'].isnull()
    fusions = fusions[~is_empty_seq]

    # iterate over each fusion
    output_list = [['ID', 'type', 'gene', 'degron_id (loss)', 'motif (loss)',
                    'score (loss)', 'score sum (loss)', 'degron_id (retained)',
                    'motif (retained)', 'score (retained)', 'score sum (retained)',
                    'degron_id (all)', 'motif (all)', 'score (all)', 'score sum (all)']]
    for ix, row in fusions.iterrows():
        g1, g2 = row["5'_gene"], row["3'_gene"]
        c1, c2 = row['CodonPos1'], row['CodonPos2']

        # get all degrons for those genes
        g1_degrons = degron.loc[degron['geneSymbol']==g1,:].sort_values('Random Forest score', ascending=False)
        g2_degrons = degron.loc[degron['geneSymbol']==g2,:].sort_values('Random Forest score', ascending=False)

        # figure out if fusion impacted degrons
        if len(g1_degrons):
            tmp_did_all = ','.join(g1_degrons['DegronID'].values)
            tmp_motif_all = ','.join(g1_degrons['motif'].values)
            tmp_scores_all = ','.join(g1_degrons['Random Forest score'].astype(str).values)
            #tmp_scores_all_sum = g1_degrons['Random Forest score'].sum()
            tmp_scores_all_sum = sum_rf_score(g1_degrons)

            # figure out if any degrons are lost
            g1_lost = g1_degrons.loc[g1_degrons['end']>=c1, :]
            if len(g1_lost):
                tmp_did_lost = ','.join(g1_lost['DegronID'].values)
                tmp_motif_lost = ','.join(g1_lost['motif'].values)
                tmp_scores_lost = ','.join(g1_lost['Random Forest score'].astype(str).values)
                #tmp_scores_lost_sum = g1_lost['Random Forest score'].sum()
                tmp_scores_lost_sum = sum_rf_score(g1_lost)
            else:
                tmp_did_lost = None
                tmp_motif_lost = None
                tmp_scores_lost = None
                tmp_scores_lost_sum = None
                #output_list.append([row['ID'], "5_prime", g1, None, None, None, tmp_did_all, tmp_motif_all, tmp_scores_all])
            # figure out if any degrons are retained
            g1_retained = g1_degrons.loc[g1_degrons['start']<c1, :]
            if len(g1_retained):
                tmp_did_retained = ','.join(g1_retained['DegronID'].values)
                tmp_motif_retained = ','.join(g1_retained['motif'].values)
                tmp_scores_retained = ','.join(g1_retained['Random Forest score'].astype(str).values)
                #tmp_scores_retained_sum = g1_retained['Random Forest score'].sum()
                tmp_scores_retained_sum = sum_rf_score(g1_retained)
            else:
                tmp_did_retained = None
                tmp_motif_retained = None
                tmp_scores_retained = None
                tmp_scores_retained_sum = None
                #output_list.append([row['ID'], "5_prime", g1, None, None, None, tmp_did_all, tmp_motif_all, tmp_scores_all])
            output_list.append(
                [row['ID'], "5_prime", g1,
                 tmp_did_lost, tmp_motif_lost, tmp_scores_lost, tmp_scores_lost_sum,
                 tmp_did_retained, tmp_motif_retained, tmp_scores_retained, tmp_scores_retained_sum,
                 tmp_did_all, tmp_motif_all, tmp_scores_all, tmp_scores_all_sum,
            ])
        else:
            output_list.append([row['ID'], "5_prime", g1,
                                None, None, None, None,
                                None, None, None, None,
                                None, None, None, None])

        if len(g2_degrons):
            tmp_did_all = ','.join(g2_degrons['DegronID'].values)
            tmp_motif_all = ','.join(g2_degrons['motif'].values)
            tmp_scores_all = ','.join(g2_degrons['Random Forest score'].astype(str).values)
            #tmp_scores_all_sum = g2_degrons['Random Forest score'].sum()
            tmp_scores_all_sum = sum_rf_score(g2_degrons)

            # figure out if any degrons are lost
            g2_lost = g2_degrons.loc[g2_degrons['start']<=c2, :]
            if len(g2_lost):
                tmp_did_lost = ','.join(g2_lost['DegronID'].values)
                tmp_motif_lost = ','.join(g2_lost['motif'].values)
                tmp_scores_lost = ','.join(g2_lost['Random Forest score'].astype(str).values)
                #tmp_scores_lost_sum = g2_lost['Random Forest score'].sum()
                tmp_scores_lost_sum = sum_rf_score(g2_lost)
            else:
                tmp_did_lost = None
                tmp_motif_lost = None
                tmp_scores_lost = None
                tmp_scores_lost_sum = None
                #output_list.append([row['ID'], "3_prime", g2, None, None, None, tmp_did_all, tmp_motif_all, tmp_scores_all])
            # figure out if any degrons are retained
            g2_retained = g2_degrons.loc[g2_degrons['start']>c2, :]
            if len(g2_retained):
                tmp_did_retained = ','.join(g2_retained['DegronID'].values)
                tmp_motif_retained = ','.join(g2_retained['motif'].values)
                tmp_scores_retained = ','.join(g2_retained['Random Forest score'].astype(str).values)
                #tmp_scores_retained_sum = g2_retained['Random Forest score'].sum()
                tmp_scores_retained_sum = sum_rf_score(g2_retained)
            else:
                tmp_did_retained = None
                tmp_motif_retained = None
                tmp_scores_retained = None
                tmp_scores_retained_sum = None
                #output_list.append([row['ID'], "5_prime", g1, None, None, None, tmp_did_all, tmp_motif_all, tmp_scores_all])
            output_list.append(
                [row['ID'], "3_prime", g2,
                 tmp_did_lost, tmp_motif_lost, tmp_scores_lost, tmp_scores_lost_sum,
                 tmp_did_retained, tmp_motif_retained, tmp_scores_retained, tmp_scores_retained_sum,
                 tmp_did_all, tmp_motif_all, tmp_scores_all, tmp_scores_all_sum,
            ])
        else:
            output_list.append([row['ID'], "3_prime", g2,
                                None, None, None, None,
                                None, None, None, None,
                                None, None, None, None])

    # save output
    with open(opts['output'], 'w') as whandle:
        mywriter = csv.writer(whandle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(output_list)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


