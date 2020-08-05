"""
File: degron_retention_test.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Analyze whether degrons are preferentially retained in gene fusions
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
    parser.add_argument('-t', '--three-prime',
                        action='store_true', default=False,
                        help='Flag on whether to analyze three-pime gene')
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

    # decide whether to analyze the 5' gene or the 3' gene
    if opts['three_prime']:
        eoi = "3'"  # side of fusition that you're interested in
        other = "5'"
    else:
        eoi = "5'"
        other = "3'"

    # prep analysis
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
    tmp.loc[tmp['Fusion_effect'].str.contains('out').fillna(False), 'score sum (retained) {}'.format(other)] = 0

    # normalize by protein length
    """
    if eoi=="5''":
        len_mis = tmp['ProtLen2'] - tmp['CodonPos2']
        tmp['normalized score sum (retained) {}'.format(other)] = tmp['score sum (retained) {}'.format(other)].div(len_mis).fillna(0)
    else:
        len_mis = tmp['CodonPos1']
        tmp['normalized score sum (retained) {}'.format(other)] = tmp['score sum (retained) {}'.format(other)].div(len_mis).fillna(0)
    """

    # simulate random scores
    prng = np.random.RandomState(101)
    random_result = []
    for i in range(10):
        rand_retention = tmp['score sum (retained) {}'.format(other)].fillna(0).sample(frac=1, replace=False, random_state=prng).values
        #rand_retention = tmp['normalized score sum (retained) {}'.format(other)].fillna(0).sample(frac=1, replace=False, random_state=prng).values
        random_df = pd.DataFrame({'gene': tmp["{}_gene".format(eoi)], "score": rand_retention})
        tmp_random_result = random_df.groupby('gene')['score'].sum().tolist()
        random_result += tmp_random_result

    # calculate p-value
    null_pvals = utils.nullScore2pvalTable(random_result)
    real_scores = tmp.groupby("{}_gene".format(eoi))["score sum (retained) {}".format(other)].sum().fillna(0)
    #real_scores = tmp.groupby("{}_gene".format(eoi))["normalized score sum (retained) {}".format(other)].sum().fillna(0)
    pvals = utils.compute_p_value(real_scores-0.0001, null_pvals.sort_index(ascending=False))
    myfdr = utils.bh_fdr(pvals)

    # save results
    result = pd.DataFrame({'p-value': pvals, 'q-value': myfdr, 'score': real_scores})
    result.to_csv(opts['output'], sep='\t')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
