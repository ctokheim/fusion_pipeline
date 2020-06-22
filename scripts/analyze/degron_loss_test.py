"""
File: degron_loss_test.py
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
    #tmp.loc[tmp['Fusion_effect'].str.contains('out').fillna(False), 'score sum (loss) {}'.format(other)] = 0

    # simulate random scores
    prng = np.random.RandomState(101)
    random_result = []
    for i in range(10):
        rand_loss = tmp['score sum (loss) {}'.format(eoi)].fillna(0).sample(frac=1, replace=False, random_state=prng).values
        random_df = pd.DataFrame({'gene': tmp["{}_gene".format(eoi)], "score": rand_loss})
        tmp_random_result = random_df.groupby('gene')['score'].sum().tolist()
        random_result += tmp_random_result

    # calculate p-value
    null_pvals = utils.nullScore2pvalTable(random_result)
    real_scores = tmp.groupby("{}_gene".format(eoi))["score sum (loss) {}".format(eoi)].sum().fillna(0)
    pvals = utils.compute_p_value(real_scores-0.0001, null_pvals.sort_index(ascending=False))
    myfdr = utils.bh_fdr(pvals)

    # count the in-frame proportion
    inframe_cts = fusion_calls[is_inframe.fillna(False)].groupby("{}_gene".format(eoi)).size()
    total_cts = fusion_calls[~fusion_calls['protein_sequence'].isnull()].groupby("{}_gene".format(eoi)).size()
    frame_df = pd.DataFrame({'in-frame': inframe_cts, 'all': total_cts})
    frame_df['in-frame'] = frame_df['in-frame'].fillna(0)
    frame_df['fraction in-frame'] = frame_df['in-frame'].div(frame_df['all'], fill_value=0)

    # save results
    result_loss = pd.DataFrame({'p-value': pvals, 'q-value': myfdr,})
    inframe_frac = frame_df.loc[result_loss.index, 'fraction in-frame']
    result_loss['inframe frac'] = inframe_frac
    result_loss.to_csv(opts['output'], sep='\t')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
