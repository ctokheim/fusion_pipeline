import pandas as pd
import numpy as np
from collections import Counter
import bisect

def read_tcga_fusions(mypath):
    """Read reported TCGA fusions in a format more amenable for downstream analysis"""
    tmp = pd.read_csv(mypath, sep='\t')

    # split gene name out
    fusion_calls = tmp['Fusion'].str.split('--', expand=True)

    # rename columns
    fusion_calls = fusion_calls.rename(columns={0: "5'_gene", 1: "3'_gene"})

    # add in break positions
    fusion_calls['Break1'] = tmp['Breakpoint1'].str.extract(':([0-9]+):', expand=True)[0].astype(int)
    fusion_calls['Break2'] = tmp['Breakpoint2'].str.extract(':([0-9]+):', expand=True)[0].astype(int)

    # add aditional columns
    fusion_calls['cancer_type'] = tmp['Cancer']
    fusion_calls['sample_id'] = tmp['Sample']
    fusion_calls['patient_id'] = tmp['Sample'].str[:12]

    return fusion_calls

def read_degron_impact(mypath):
    tmp = pd.read_csv(mypath, sep='\t')

    # extract break positions
    breaks = tmp['ID'].str.split(':', expand=True)[1]
    breaks_split = breaks.str.split('-', expand=True).astype(int)

    # add break location to results
    tmp['break'] = -1
    is_5prime = tmp['type']=='5_prime'
    tmp.loc[is_5prime, 'break'] = breaks_split.loc[is_5prime, 0]
    is_3prime = tmp['type']=='3_prime'
    tmp.loc[is_3prime, 'break'] = breaks_split.loc[is_3prime, 1]

    return tmp

def read_drivers(mypath):
    tmp = pd.read_csv(mypath, sep='\t').rename(columns={'Hugo Symbol': 'gene'})
    return tmp

def myrename(mystr, suffix="5'"):
    if mystr=='gene':
        return "5'_gene" if suffix=="5'" else "3'_gene"
    elif mystr=='break':
        return 'Break1' if suffix=="5'" else "Break2"
    else:
        return mystr + ' ' + suffix

def read_fusion_domains(fus_path, wt_path, fusion_annot, min_domain_len=25, min_domain_frac=0.5):
    """Add domain annotation to fusion information"""
    # read in domains
    domains = pd.read_csv(fus_path, sep='\t')
    wt_domains = pd.read_csv(wt_path, sep='\t')

    # filter out transmembrane helix
    is_tm = domains["Domain_ID"].str.startswith('TMhelix')
    domains = domains[~is_tm].copy()

    # filter out short domains
    domain_len = domains['Protein_end'] - domains['Protein_start']
    domains = domains[domain_len>=min_domain_len].copy()

    # merge in info about where the fusion happened
    domains = pd.merge(domains, fusion_annot[['ID', 'CodonPos1', 'CodonPos2']].drop_duplicates('ID'), on='ID', how='left')
    domains['origin of domain'] = "3'_gene domain"
    domains.loc[domains['Protein_end']<=domains['CodonPos1'], 'origin of domain'] = "5'_gene domain"

    # separate out protein ID
    domains['ID_five'] = domains['ID'].str.split('-', expand=True)[0]
    domains['ID_three'] = domains['ID'].str.split('-', expand=True)[1].str.split(':', expand=True)[0]

    # check 5' domains
    domains_five = domains[domains["origin of domain"]=="5'_gene domain"].copy()
    is_good_domain = []
    for ix, row in domains_five.iterrows():
        domain_end = row['Protein_end']
        fus_codon_pos = row['CodonPos1']
        if fus_codon_pos==domain_end:
            pid = row['ID_five']
            is_pid = wt_domains['PROT_ID']==pid
            is_loc = (fus_codon_pos>=wt_domains['Protein_start']) & (fus_codon_pos<=wt_domains['Protein_end'])
            tmp_df = wt_domains.loc[is_pid & is_loc,:]
            if len(tmp_df)==0:
                is_good_domain.append(True)
                continue
            else:
                tmp = tmp_df.iloc[0]
            frac_domain = (fus_codon_pos - tmp['Protein_start']) / (tmp['Protein_end'] - tmp['Protein_start'])
            if frac_domain>=min_domain_frac:
                is_good_domain.append(True)
            else:
                is_good_domain.append(False)
        else:
            is_good_domain.append(True)
    domains_five = domains_five.loc[is_good_domain,:]
    # check 3' domains
    domains_three = domains[domains["origin of domain"]=="3'_gene domain"].copy()
    is_good_domain = []
    for ix, row in domains_three.iterrows():
        domain_start = row['Protein_start']
        fus_codon_pos = row['CodonPos1']
        if fus_codon_pos==domain_start:
            pid = row['ID_three']
            is_pid = wt_domains['PROT_ID']==pid
            is_loc = (fus_codon_pos>=wt_domains['Protein_start']) & (fus_codon_pos<=wt_domains['Protein_end'])
            tmp_df = wt_domains.loc[is_pid & is_loc,:]
            if len(tmp_df)==0:
                is_good_domain.append(True)
                continue
            else:
                tmp = tmp_df.iloc[0]
            frac_domain = (tmp['Protein_end']-fus_codon_pos) / (tmp['Protein_end'] - tmp['Protein_start'])
            if frac_domain>=min_domain_frac:
                is_good_domain.append(True)
            else:
                is_good_domain.append(False)
        else:
            is_good_domain.append(True)
    domains_three = domains_three.loc[is_good_domain,:]
    # append results
    domains = pd.concat([domains_five, domains_three])

    # aggregate all domains for one fusion into a single line
    domains_agg = domains.groupby(['ID', 'origin of domain'])['Domain_name'].agg(lambda x: ','.join(x))
    domains_agg = domains_agg.reset_index()
    domains_agg = domains_agg.pivot(index='ID', columns='origin of domain', values='Domain_name')
    domains_agg = domains_agg.reset_index()

    return domains_agg

def merge_lost_prot_domains(fusion_calls, wt_dom_path):
    # read in wildtype protein domains
    wt_domains = pd.read_csv(wt_dom_path, sep='\t')
    is_pfam = wt_domains['database']=='pfam'
    wt_domains = wt_domains[is_pfam]
    wt_domains_dict = wt_domains.groupby('PROT_ID')['Domain_description'].apply(list).to_dict()

    five_dom_lost, three_dom_lost = [], []
    for ix, row in fusion_calls.iterrows():
        # skip if no annotate prot sequence
        if 'ENSP' not in str(row['ID']):
            five_dom_lost.append(np.nan)
            three_dom_lost.append(np.nan)
            continue

        tmp_row = row.fillna('')

        # get the ensembl protein domains
        p1 = tmp_row['ID'].split('-')[0]
        p2 = tmp_row['ID'].split('-')[1].split(':')[0]

        # count the domains
        p1_fus_cts = Counter(tmp_row["5'_gene domain"].split(','))
        p1_wt_cts = Counter(wt_domains_dict.get(p1, []))
        p2_fus_cts = Counter(tmp_row["3'_gene domain"].split(','))
        p2_wt_cts = Counter(wt_domains_dict.get(p2, []))

        # figure out which domains are lost, if any
        p1_diff_fusion = p1_wt_cts - p1_fus_cts
        p1_domain_lost = ''.join([(x+',')*p1_diff_fusion[x] for x in p1_diff_fusion]).strip(',')
        p2_diff_fusion = p2_wt_cts - p2_fus_cts
        p2_domain_lost = ''.join([(x+',')*p2_diff_fusion[x] for x in p2_diff_fusion]).strip(',')

        # append result
        five_dom_lost.append(p1_domain_lost)
        three_dom_lost.append(p2_domain_lost)

    # add results to fusion annotation
    fusion_calls["5'_gene domain lost"] = five_dom_lost
    fusion_calls["3'_gene domain lost"] = three_dom_lost

    # remove any cases of empty strings
    fusion_calls.loc[fusion_calls["5'_gene domain lost"]=='', "5'_gene domain lost"] = np.nan
    fusion_calls.loc[fusion_calls["3'_gene domain lost"]=='', "3'_gene domain lost"] = np.nan

    return fusion_calls

def nullScore2pvalTable(scores):
    """Create a p-value table from null-distribution scores"""
    scores = pd.Series(scores)
    score_counts = scores.value_counts().sort_index(ascending=True)
    score_cum_counts = np.cumsum(score_counts.values[::-1])[::-1]
    count_sum = np.sum(score_counts)
    pval_array = score_cum_counts / count_sum
    pval_series = pd.Series(pval_array, index=score_counts.index)
    return pval_series

def compute_p_value(scores, null_p_values):
    """Get the p-value for each score by examining the list null distribution
    where scores are obtained by a certain probability.
    NOTE: uses score2pval function
    Parameters
    ----------
    scores : pd.Series
        series of observed scores
    null_p_values: pd.Series
        Empirical null distribution, index are scores and values are p values
    Returns
    -------
    pvals : pd.Series
        Series of p values for scores
    """
    num_scores = len(scores)
    pvals = pd.Series(np.zeros(num_scores))
    null_p_val_scores = list(reversed(null_p_values.index.tolist()))
    #null_p_values = null_p_values.ix[null_p_val_scores].copy()
    null_p_values.sort_values(inplace=True, ascending=False)
    pvals = scores.apply(lambda x: score2pval(x, null_p_val_scores, null_p_values))
    return pvals

def score2pval(score, null_scores, null_pvals):
    """Looks up the P value from the empirical null distribution based on the provided
    score.
    NOTE: null_scores and null_pvals should be sorted in ascending order.
    Parameters
    ----------
    score : float
        score to look up P value for
    null_scores : list
        list of scores that have a non-NA value
    null_pvals : pd.Series
        a series object with the P value for the scores found in null_scores
    Returns
    -------
    pval : float
        P value for requested score
    """
    # find position in simulated null distribution
    pos = bisect.bisect_right(null_scores, score)

    # if the score is beyond any simulated values, then report
    # a p-value of zero
    if pos == null_pvals.size and score > null_scores[-1]:
        return 0
    # condition needed to prevent an error
    # simply get last value, if it equals the last value
    elif pos == null_pvals.size:
        return null_pvals.iloc[pos-1]
    # normal case, just report the corresponding p-val from simulations
    else:
        return null_pvals.iloc[pos]

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

def process_cgc(path, return_dataframe=False, fusions=False):
    """Get the list of CGC genes with small somatic variants."""
    # read in data
    df = pd.read_table(path)

    # keep small somatic variants
    if not fusions:
        s = df['Mutation Types']
        is_small = s.str.contains('Mis|F|N|S').fillna(False)
        is_somatic = ~df['Tumour Types(Somatic)'].isnull()
        df = df[is_small & is_somatic].copy()

        # label oncogenes / TSG
        df['Is Oncogene (CGC)'] = 'No'
        df.loc[df['Role in Cancer'].fillna('').str.contains('oncogene'), 'Is Oncogene'] = 'Yes'
        df['Is Tumor Suppressor Gene (CGC)'] = 'No'
        df.loc[df['Role in Cancer'].fillna('').str.contains('TSG'), 'Is Tumor Suppressor Gene'] = 'Yes'
        df['Is Driver Gene (CGC)'] = 'Yes'

        # rename columns
        df = df.rename(columns={'Entrez GeneId': 'Entrez Gene ID', 'Gene Symbol': 'Hugo Symbol'})

        # get gene names
        if not return_dataframe:
            cgc_genes = df['Gene Symbol'].tolist()
        else:
            cgc_genes = df

        return cgc_genes
    else:
        # return fusion gene information
        has_fus_partner = ~df['Translocation Partner'].isnull()
        output_list = []
        for ix, row in df[has_fus_partner].iterrows():
            g1 = row["Gene Symbol"]
            for g2 in row['Translocation Partner'].split(', '):
                output_list.append([g1, g2])
        output_df = pd.DataFrame(output_list, columns=["Gene1", "Gene2"])
        output_df['GENE_ID'] = output_df['Gene1'] + '--' + output_df['Gene2']

        if not return_dataframe:
            cgc_genes = list(set(output_df["Gene1"].unique()) | set(output_df["Gene2"]))
        else:
            cgc_genes = output_df

        return cgc_genes

def read_oncokb_fusions(mypath):
    """Read in oncogenic fusions annotated from OncoKB."""
    tmp = pd.read_csv(mypath, sep='\t')
    is_fus = tmp['Alteration'].str.contains('Fusion')
    tmp = tmp[is_fus].copy()

    # fetch cases where both fusion partners specified
    oncokb_fusions_both = tmp["Alteration"].str.extract('(.+) Fusion', expand=False)
    oncokb_fusions_both = oncokb_fusions_both.str.replace('-', '--').dropna()

    # get cases where any fusion partner is allowed
    oncokb_fusions_one = tmp[tmp['Alteration']=='Fusions']['Gene'].unique()

    return oncokb_fusions_both, oncokb_fusions_one

def read_drugable_genes(mypath):
    """Read in druggable genes from DGIdb"""
    drug = pd.read_csv('data/drugability/interactions.tsv', sep='\t')
    drug_types = ['inhibitor', 'antagonist', 'blocker', 'channel blocker', 'gating inhibitor', 'negative modulator', 'channel blocker,gating inhibitor',
                  'allosteric modulator,antagonist', 'inhibitor allosteric modulator', 'suppressor', 'antagonist,inhibitor']
    drug = drug[drug.interaction_types.isin(drug_types)]
    cts = drug.groupby(['drug_name', 'interaction_claim_source'])['gene_name'].nunique().sort_values(ascending=True)
    single_drug_impact = cts[cts==1].copy().index.to_series().reset_index()['drug_name'].unique()
    drug = drug[drug['drug_name'].isin(single_drug_impact)].dropna(subset=['gene_name'])

    # create a table of gene 2 drug relationships
    gene2drug = drug.groupby('gene_name')['drug_name'].agg(lambda x: ','.join(list(set(x)))).reset_index()
    #gene_names = drug['gene_name'].unique()
    return gene2drug
