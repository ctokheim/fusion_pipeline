"""
File: cterm_degron_loss_test.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Analyze whether internal degrons are preferentially lost in gene fusions
"""
import pandas as pd
import numpy as np
import utils
import argparse


def parse_arguments():
    info = 'Analyze whether degrons are preferentially retained in gene fusions'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Annotated fusion files')
    parser.add_argument('-d', '--domain',
                        action='store_true', default=False,
                        help='Flag on whether to analyze only fusions with an intact domain')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file prefix')

    args = parser.parse_args()
    return vars(args)


def main(opts):
    fusion_calls = pd.read_csv(opts['input'], sep='\t')
    fusion_calls['delta degron potential'] = -fusion_calls['delta degron potential']

    eoi = "5'"  # side of fusition that you're interested in
    other = "3'"

    # filter fusion events
    is_inframe = fusion_calls['Fusion_effect'].str.contains('in-frame').fillna(False)
    is_domain_retained = (~fusion_calls["{}_gene domain".format(eoi)].isnull())
    is_low_point_drivers = (fusion_calls['num_point_drivers'].fillna(0)<4)
    if opts["domain"]:
        myflag = is_inframe & is_low_point_drivers & is_domain_retained
    else:
        myflag = is_inframe & is_low_point_drivers
    tmp = fusion_calls[myflag].copy()
    gene_cts = tmp["{}_gene".format(eoi)].value_counts()
    minimal_genes = gene_cts[gene_cts>1].index.tolist()
    tmp = tmp[tmp["{}_gene".format(eoi)].isin(minimal_genes)].copy()

    prng = np.random.RandomState(101)
    random_result = []
    for i in range(10):
        rand_retention = tmp['delta degron potential'].fillna(0).sample(frac=1, replace=False, random_state=prng).values
        random_df = pd.DataFrame({'gene': tmp["{}_gene".format(eoi)], "score": rand_retention})
        tmp_random_result = random_df.groupby('gene')['score'].sum().tolist()
        random_result += tmp_random_result

    # null distribution
    null_pvals = utils.nullScore2pvalTable(random_result)
    real_scores = tmp.groupby("{}_gene".format(eoi))["delta degron potential"].sum().fillna(0)
    pvals = utils.compute_p_value(real_scores-0.0001, null_pvals.sort_index(ascending=False))
    myfdr = utils.bh_fdr(pvals)

    # figure out which genes are bonafide oncogenes
    #is_og = (fusion_calls["{}_Is Oncogene (OncoKB)".format(eoi)]=="Yes")
    #og_genes = fusion_calls[is_og]["{}_gene".format(eoi)].unique()
    is_og = (fusion_calls["{}_Is Oncogene (OncoKB)".format(eoi)]=='Yes') | (fusion_calls["{}_Is Oncogene (TCGA)".format(eoi)]=='Yes') | (fusion_calls["{}_Is Oncogene (CGC)".format(eoi)]=='Yes')
    is_not_tsg = (fusion_calls["{}_Is Tumor Suppressor Gene (OncoKB)".format(eoi)]!='Yes') & (fusion_calls["{}_Is Tumor Suppressor Gene (TCGA)".format(eoi)]!='Yes') & (fusion_calls["{}_Is Tumor Suppressor Gene (CGC)".format(eoi)]!='Yes')
    og_genes = fusion_calls[is_og & is_not_tsg]["{}_gene".format(eoi)].unique()
    # figure out which genes are actually bonafide TSGs
    is_tsg = (fusion_calls["{}_Is Tumor Suppressor Gene (OncoKB)".format(eoi)]=="Yes")
    tsg_genes = fusion_calls[is_tsg]["{}_gene".format(eoi)].unique()

    # count the in-frame proportion
    inframe_cts = fusion_calls[is_inframe.fillna(False)].groupby("{}_gene".format(eoi)).size()
    total_cts = fusion_calls[~fusion_calls['protein_sequence'].isnull()].groupby("{}_gene".format(eoi)).size()
    frame_df = pd.DataFrame({'in-frame': inframe_cts, 'all': total_cts})
    frame_df['in-frame'] = frame_df['in-frame'].fillna(0)
    frame_df['fraction in-frame'] = frame_df['in-frame'].div(frame_df['all'], fill_value=0)

    # save results
    result = pd.DataFrame({'p-value': pvals, 'q-value': myfdr, 'delta degron potential': real_scores})
    inframe_frac = frame_df.loc[result.index, 'fraction in-frame']
    result['inframe frac'] = inframe_frac
    result.to_csv(opts['output'], sep='\t')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
