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
    parser.add_argument('-n', '--num-simulations',
                        type=int, default=10000,
                        help='Number of simulations')
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


def calc_score_sum(mydf, gene_col, feat_col):
    myscores = mydf.groupby(gene_col)[feat_col].sum().fillna(0)
    return myscores


def statistical_test(tmp, eoi, num_sim=10000):
    feat_col_name = "score sum (loss) {}".format(eoi)
    real_scores = calc_score_sum(tmp, gene_col="{}_gene".format(eoi), feat_col=feat_col_name)

    # simulate
    prng = np.random.RandomState(101)
    perms = np.argsort(prng.rand(tmp[feat_col_name].fillna(0).values.shape[0], num_sim), axis=0)
    #permuted_ctype = np.hstack((ctype.values[:,np.newaxis], ctype.values[perms]))
    permuted_scores = tmp[feat_col_name].values[perms]

    # create a dataframe with random ctypes
    random_df = pd.DataFrame(permuted_scores, columns=map(str, range(num_sim)))
    random_df['gene'] = tmp["{}_gene".format(eoi)].values

    # calculate ctype counts
    rand_score_sums = random_df.groupby("gene").agg({x: 'sum' for x in map(str, range(num_sim))})

    # figure out the the p-value
    rand_score_sums = rand_score_sums.loc[real_scores.index]
    flags = rand_score_sums.apply(lambda x, y: x>=y, args=(real_scores,))
    null_cts = flags.sum(axis=1)
    # minimum of one count to indicate minimum p-value resolution
    null_cts.loc[null_cts==0] = 1
    # get p-values
    pvals = null_cts / num_sim

    # save results
    result_loss = pd.DataFrame({'p-value': pvals, 'score': real_scores})

    return result_loss


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
    #is_low_point_drivers = (fusion_calls['num_point_drivers'].fillna(0)<4)
    if opts["domain"]:
        myflag = is_domain_retained & is_inframe # & is_low_point_drivers
    else:
        myflag = is_inframe    # & is_low_point_drivers
    tmp = fusion_calls[myflag].copy()
    gene_cts = tmp[is_inframe]["{}_gene".format(eoi)].value_counts()
    minimal_genes = gene_cts[gene_cts>1].index.tolist()
    tmp = tmp[tmp["{}_gene".format(eoi)].isin(minimal_genes)].copy()

    # make sure non-inframe fusions don't contribute to scores
    tmp.loc[~is_inframe, 'score sum (loss) {}'.format(eoi)] = 0

    # figure out prot length col
    if eoi == "5'":
        prot_len_col = 'ProtLen1'
    else:
        prot_len_col = 'ProtLen2'

    # split small/large genes
    gcol = "{}_gene".format(eoi)
    is_large = (tmp[prot_len_col]>=1500)
    large_prot = tmp[is_large][gcol].unique()
    is_small = tmp[prot_len_col]<1500
    small_prot = tmp[is_small][gcol].unique()
    small_prot = set(small_prot) - set(large_prot)

    result1 = statistical_test(tmp[tmp[gcol].isin(small_prot)].copy(), eoi=eoi, num_sim=opts['num_simulations'])
    result2 = statistical_test(tmp[tmp[gcol].isin(large_prot)].copy(), eoi=eoi, num_sim=opts['num_simulations'])
    result_loss = pd.concat([result1, result2])
    result_loss['q-value'] = utils.bh_fdr(result_loss['p-value'])

    # count the in-frame proportion
    inframe_cts = fusion_calls[is_inframe.fillna(False)].groupby("{}_gene".format(eoi)).size()
    total_cts = fusion_calls[~fusion_calls['protein_sequence'].isnull()].groupby("{}_gene".format(eoi)).size()
    frame_df = pd.DataFrame({'in-frame': inframe_cts, 'all': total_cts})
    frame_df['in-frame'] = frame_df['in-frame'].fillna(0)
    frame_df['fraction in-frame'] = frame_df['in-frame'].div(frame_df['all'], fill_value=0)

    # save results
    inframe_frac = frame_df.loc[result_loss.index, 'fraction in-frame']
    result_loss['inframe frac'] = inframe_frac
    result_loss.to_csv(opts['output'], sep='\t')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
