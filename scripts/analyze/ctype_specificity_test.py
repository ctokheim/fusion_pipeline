"""
File: ctype_specificity_test.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Statistical test for cancer type specificity
"""
import numpy as np
import pandas as pd
import argparse


def parse_arguments():
    info = 'Statistical test for cancer type specificity'
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


def calc_entropy(mydf, gene_col, cancer_col="cancer_type"):
    ctype_cts = mydf.groupby(gene_col)[cancer_col].value_counts()
    gene_cts = ctype_cts.sum(level=0)
    ctype_frac = ctype_cts.div(gene_cts)
    entropy = ctype_frac.groupby(level=0).apply(lambda x: -np.sum(x*np.log2(x)))
    return entropy


def cummin(x):
    """A python implementation of the cummin function in R"""
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    """A python implementation of the Benjamani-Hochberg FDR method.
    This code should always give precisely the same answer as using
    p.adjust(pval, method="BH") in R.
    Parameters
    ----------
    pval : list or array
        list/array of p-values
    Returns
    -------
    pval_adj : np.array
        adjusted p-values according the benjamani-hochberg method
    """
    pval_array = np.array(pval)
    sorted_order = np.argsort(pval_array)
    original_order = np.argsort(sorted_order)
    pval_array = pval_array[sorted_order]

    # calculate the needed alpha
    n = float(len(pval))
    pval_adj = np.zeros(int(n))
    i = np.arange(1, int(n)+1, dtype=float)[::-1]  # largest to smallest
    pval_adj = np.minimum(1, cummin(n/i * pval_array[::-1]))[::-1]
    return pval_adj[original_order]


def main(opts):
    fusion_calls = pd.read_csv(opts['input'], sep='\t')

    # decide whether to analyze the 5' gene or the 3' gene
    if opts['three_prime']:
        eoi = "3'"  # side of fusition that you're interested in
    else:
        eoi = "5'"

    # filter fusion events
    is_inframe = fusion_calls['Fusion_effect'].str.contains('in-frame').fillna(False)
    is_domain_retained = (~fusion_calls["{}_gene domain".format(eoi)].isnull())
    is_low_point_drivers = (fusion_calls['num_point_drivers'].fillna(0)<4)
    if opts["domain"]:
        myflag = is_inframe & is_domain_retained # & is_low_point_drivers
    else:
        myflag = is_inframe # & is_low_point_drivers
    tmp = fusion_calls[myflag].copy()
    # only keep genes with at least 2 fusions
    gene_cts = tmp["{}_gene".format(eoi)].value_counts()
    minimal_genes = gene_cts[gene_cts>1].index.tolist()
    tmp = tmp[tmp["{}_gene".format(eoi)].isin(minimal_genes)].copy()

    # calculate entropy
    entropy = calc_entropy(tmp, gene_col="{}_gene".format(eoi))

    # simulate
    prng = np.random.RandomState(101)
    ctype = tmp['cancer_type']
    num_sim = opts['num_simulations']
    perms = np.argsort(prng.rand(ctype.values.shape[0], num_sim), axis=0)
    #permuted_ctype = np.hstack((ctype.values[:,np.newaxis], ctype.values[perms]))
    permuted_ctype = ctype.values[perms]

    # create a dataframe with random ctypes
    random_df = pd.DataFrame(permuted_ctype, columns=map(str, range(num_sim)))
    random_df['gene'] = tmp["{}_gene".format(eoi)].values

    # calculate ctype counts
    ctype_cts = random_df.groupby("gene").agg({x: 'value_counts' for x in map(str, range(num_sim))})
    gene_cts = ctype_cts["0"].sum(level=0)
    #ctype_frac = ctype_cts.div(gene_cts)
    ctype_frac = ctype_cts.apply(lambda x, y: x.div(y, level=0), args=(gene_cts,))
    rand_entropy = ctype_frac.groupby(level=0).apply(lambda x: -np.sum(x*np.log2(x)))

    # figure out the the p-value
    rand_entropy = rand_entropy.loc[entropy.index]
    flags = rand_entropy.apply(lambda x, y: x<=y, args=(entropy,))
    null_cts = flags.sum(axis=1)
    # minimum of one count to indicate minimum p-value resolution
    null_cts.loc[null_cts==0] = 1
    # get p-values
    null_pvals = null_cts / num_sim
    myfdr = bh_fdr(null_pvals)

    # figure out which genes are bonafide oncogenes
    #is_og = (fusion_calls["{}_Is Oncogene (OncoKB)".format(eoi)]=="Yes")
    is_og = (fusion_calls["{}_Is Oncogene (OncoKB)".format(eoi)]=='Yes') | (fusion_calls["{}_Is Oncogene (TCGA)".format(eoi)]=='Yes') | (fusion_calls["{}_Is Oncogene (CGC)".format(eoi)]=='Yes')
    is_not_tsg = (fusion_calls["{}_Is Tumor Suppressor Gene (OncoKB)".format(eoi)]!='Yes') & (fusion_calls["{}_Is Tumor Suppressor Gene (TCGA)".format(eoi)]!='Yes') & (fusion_calls["{}_Is Tumor Suppressor Gene (CGC)".format(eoi)]!='Yes')
    og_genes = fusion_calls[is_og & is_not_tsg]["{}_gene".format(eoi)].unique()

    # count the in-frame proportion
    inframe_cts = fusion_calls[is_inframe.fillna(False)].groupby("{}_gene".format(eoi)).size()
    total_cts = fusion_calls[~fusion_calls['protein_sequence'].isnull()].groupby("{}_gene".format(eoi)).size()
    frame_df = pd.DataFrame({'in-frame': inframe_cts, 'all': total_cts})
    frame_df['in-frame'] = frame_df['in-frame'].fillna(0)
    frame_df['fraction in-frame'] = frame_df['in-frame'].div(frame_df['all'], fill_value=0)

    # compile results
    result_loss = pd.DataFrame({'p-value': null_pvals, 'q-value': myfdr,})
    inframe_frac = frame_df.loc[result_loss.index, 'fraction in-frame']
    result_loss['inframe frac'] = inframe_frac
    result_loss_tmp = result_loss.loc[result_loss.index.intersection(og_genes)].dropna()
    result_loss_tmp['q-value2'] = bh_fdr(result_loss_tmp['p-value'])

    # save results
    mypath = '{}.txt'.format(opts['output'])
    result_loss.to_csv(mypath, sep='\t')
    mypath = '{}_og.txt'.format(opts['output'])
    result_loss_tmp.to_csv(mypath, sep='\t')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


