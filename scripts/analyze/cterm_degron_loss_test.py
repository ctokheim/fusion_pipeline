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
import re


def parse_arguments():
    info = 'Analyze whether degrons are preferentially retained in gene fusions'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Annotated fusion files')
    parser.add_argument('-d', '--domain',
                        action='store_true', default=False,
                        help='Flag on whether to analyze only fusions with an intact domain')
    parser.add_argument('-m', '--motif',
                        type=str, required=True,
                        help='C-terminal motifs')
    parser.add_argument('-w', '--wildtype-seq',
                        type=str, required=True,
                        help='Wildtype protein sequence')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file prefix')
    args = parser.parse_args()
    return vars(args)

def read_cterm_motifs(mypath):
    mydf = pd.read_csv(mypath, sep='\t')
    motif_series = mydf['motif'].str.replace('x', '.') + '$'
    return motif_series.values.tolist()

def regex_search(seq, motif_list):
    for m in motif_list:
        hit = re.search(m, str(seq))
        if hit:
            return 1
    return 0

def add_cterm_motif(fusion_calls, gps_motif_path, wt_prot_seq_path):
    """Add c-terminal degron annotation"""
    # read in c-terminal motifs and wt sequence
    cterm_motifs = read_cterm_motifs(gps_motif_path)
    wt_seq = pd.read_csv(wt_prot_seq_path, sep='\t').rename(columns={"transcript_id": "5'_transcript"})

    # search for c-terminal degron motifs
    fus_hits = fusion_calls['protein_sequence'].apply(regex_search, args=(cterm_motifs,))
    wt_hits = wt_seq['protein_sequence'].apply(regex_search, args=(cterm_motifs,))
    wt_seq['wt_cterm_deg'] = wt_hits

    # add results to fusion df
    fusion_calls = pd.merge(fusion_calls, wt_seq[["5'_transcript", 'wt_cterm_deg']], on="5'_transcript", how='left')
    fusion_calls['fusion_cterm_deg'] = fus_hits.values

    return fusion_calls

def calc_score_sum(mydf, gene_col, feat_col):
    myscores = mydf.groupby(gene_col)[feat_col].sum().fillna(0)
    return myscores

def statistical_test(tmp, num_sim=10000):
    feat_col_name = "delta degron potential"
    #feat_col_name = "second degron potential"
    real_scores = calc_score_sum(tmp, gene_col="5'_gene", feat_col=feat_col_name)
    wt_cterm_deg = tmp.groupby("5'_gene")['wt_cterm_deg'].max().fillna(0)
    real_scores2 = real_scores * 1 #* wt_cterm_deg

    # simulate
    prng = np.random.RandomState(101)
    second_deg_col = 'second degron potential'
    perms = np.argsort(prng.rand(tmp[feat_col_name].fillna(0).values.shape[0], num_sim), axis=0)
    #permuted_ctype = np.hstack((ctype.values[:,np.newaxis], ctype.values[perms]))
    permuted_scores = tmp[feat_col_name].values[perms]

    # create a dataframe with random ctypes
    random_df = pd.DataFrame(permuted_scores, columns=map(str, range(num_sim)))
    random_df = random_df#.multiply(tmp["wt_cterm_deg"].values, axis=0)
    random_df['gene'] = tmp["5'_gene"].values

    # calculate ctype counts
    rand_score_sums = random_df.groupby("gene").agg({x: 'sum' for x in map(str, range(num_sim))})

    # figure out the the p-value
    rand_score_sums = rand_score_sums.loc[real_scores.index]
    flags = rand_score_sums.apply(lambda x, y: x>= (y-0.001), args=(real_scores2,))
    null_cts = flags.sum(axis=1)
    # minimum of one count to indicate minimum p-value resolution
    null_cts.loc[null_cts==0] = 1
    # get p-values
    pvals = null_cts / num_sim

    # save results
    result_loss = pd.DataFrame({'p-value': pvals, 'score': real_scores2})
    result_loss['q-value'] = utils.bh_fdr(result_loss['p-value'])

    return result_loss


def main(opts):
    fusion_calls = pd.read_csv(opts['input'], sep='\t')
    fusion_calls['delta degron potential'] = -fusion_calls['delta degron potential']

    # test out trying to prevent negative degron evidence from making something
    # significant
    #fusion_calls['delta degron potential'] = (fusion_calls['first degron potential'] - fusion_calls['second degron potential'].apply(lambda x: max(0, x)))

    eoi = "5'"  # side of fusition that you're interested in
    other = "3'"

    # filter fusion events
    is_inframe = fusion_calls['Fusion_effect'].str.contains('in-frame').fillna(False)
    is_domain_retained = (~fusion_calls["{}_gene domain".format(eoi)].isnull())
    #is_low_point_drivers = (fusion_calls['num_point_drivers'].fillna(0)<4)
    if opts["domain"]:
        myflag = is_inframe & is_domain_retained # & is_low_point_drivers
    else:
        myflag = is_inframe # & is_low_point_drivers
    tmp = fusion_calls[myflag].copy()
    gene_cts = tmp["{}_gene".format(eoi)].value_counts()
    minimal_genes = gene_cts[gene_cts>1].index.tolist()
    tmp = tmp[tmp["{}_gene".format(eoi)].isin(minimal_genes)].copy()

    # add c-terminal degron motif info
    tmp = add_cterm_motif(tmp, opts['motif'], opts['wildtype_seq'])
    tmp = tmp[tmp["wt_cterm_deg"]==1].copy()
    result = statistical_test(tmp)
    #import IPython ; IPython.embed() ; raise

    """
    prng = np.random.RandomState(101)
    random_result = []
    for i in range(100):
        rand_loss = tmp['delta degron potential'].fillna(0).sample(frac=1, replace=False, random_state=prng).values
        #rand_second = tmp['second degron potential'].fillna(0).sample(frac=1, replace=False, random_state=prng).values
        random_df = pd.DataFrame({'gene': tmp["{}_gene".format(eoi)], "delta degron score": rand_loss, "motif": tmp['wt_cterm_deg']})
        #random_df = pd.DataFrame({'gene': tmp["{}_gene".format(eoi)], "delta degron score": -(tmp['first degron potential'].values - rand_second), "motif": tmp['wt_cterm_deg']})
        random_df["score"] = random_df["delta degron score"] * random_df['motif']
        tmp_random_result = random_df.groupby('gene')['score'].sum().tolist()
        random_result += tmp_random_result

    # null distribution
    null_pvals = utils.nullScore2pvalTable(random_result[:10000])
    real_scores = tmp.groupby("{}_gene".format(eoi))["delta degron potential"].sum().fillna(0)
    #tmp['score'] = tmp['delta degron potential'] * tmp['wt_cterm_deg']
    #real_scores = tmp.groupby("{}_gene".format(eoi))["score"].sum().fillna(0)
    pvals = utils.compute_p_value(real_scores-0.0001, null_pvals.sort_index(ascending=False))
    myfdr = utils.bh_fdr(pvals)
    """

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
    #result = pd.DataFrame({'p-value': pvals, 'q-value': myfdr, 'delta degron potential': real_scores})
    inframe_frac = frame_df.loc[result.index, 'fraction in-frame']
    result['inframe frac'] = inframe_frac
    result.to_csv(opts['output'], sep='\t')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
