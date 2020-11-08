"""
File: add_fusion_annotations.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Adds annotation information to gene fusion data
"""
import pandas as pd
import utils
import argparse


def parse_arguments():
    info = 'Adds annotation information to gene fusion data'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Annotation for fusions')
    parser.add_argument('-d', '--degron-pred',
                        type=str, required=True,
                        help='Predicted impact on internal degron')
    parser.add_argument('-c', '--cterm-degron',
                        type=str, required=True,
                        help='Predicted impact on cterminal degrons')
    parser.add_argument('-fc', '--fusion-calls',
                        type=str, required=True,
                        help='original TCGA fusion calls')
    parser.add_argument('-wd', '--wildtype-domain',
                        type=str, required=True,
                        help='Wildtype protein domains')
    parser.add_argument('-fd', '--fusion-domain',
                        type=str, required=True,
                        help='Fusion protein domains')
    parser.add_argument('-cgc', '--cgc',
                        type=str, required=True,
                        help='Cancer gene census')
    parser.add_argument('-ot', '--og-tsg',
                        type=str, required=True,
                        help='Oncogene / TSG annotation file')
    parser.add_argument('-ok', '--oncokb',
                        type=str, required=True,
                        help='OncoKB annotations')
    parser.add_argument('-m', '--maf',
                        type=str, required=True,
                        help='MAF file of TCGA point mutations')
    parser.add_argument('-drug', '--druggability',
                        type=str, required=True,
                        help='Druggability annotation')
    parser.add_argument('-driver', '--driver',
                        type=str, required=True,
                        help='Driver annotation')
    parser.add_argument('-u', '--ubiquitin',
                        type=str, required=True,
                        help='Ubiquitin sites')
    parser.add_argument('-xref', '--cross-ref',
                        type=str, required=True,
                        help='Uniprot cross-ref file')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in data
    fusion_calls = utils.read_tcga_fusions(opts['fusion_calls'])
    fusion_annot = pd.read_csv(opts['input'], sep='\t')
    degron_impact = utils.read_degron_impact(opts['degron_pred'])
    mycols = degron_impact.columns[2:].values

    # add driver gene information
    drivers = utils.read_drivers(opts['og_tsg'])
    mycols_drivers = ['gene', 'Is Oncogene (OncoKB)', 'Is Tumor Suppressor Gene (OncoKB)',
                      'Is Oncogene (TCGA)', 'Is Tumor Suppressor Gene (TCGA)', 'Is Driver Gene (TCGA)',
                      'Is Oncogene (CGC)', 'Is Tumor Suppressor Gene (CGC)', 'Is Driver Gene (CGC)']
    rename_func = lambda x: "5'_" + x
    fusion_calls = pd.merge(fusion_calls, drivers[mycols_drivers].rename(columns=rename_func), on=["5'_gene"], how='left')
    rename_func = lambda x: "3'_" + x
    fusion_calls = pd.merge(fusion_calls, drivers[mycols_drivers].rename(columns=rename_func), on=["3'_gene"], how='left')

    # add agfusion annotations
    fusion_calls = pd.merge(fusion_calls, fusion_annot, on=["5'_gene", "3'_gene", "Break1", "Break2"], how='left')

    # add fusion-specific driver info from CGC
    cgc_fusions = utils.process_cgc(opts['cgc'], return_dataframe=True, fusions=True)
    fusion_calls['Is Fusion Driver (CGC)'] = 'No'
    tmp1 = fusion_calls["5'_gene"] + '--' + fusion_calls["3'_gene"]
    tmp2 = fusion_calls["3'_gene"] + '--' + fusion_calls["5'_gene"]
    is_tmp1 = tmp1.isin(cgc_fusions["GENE_ID"].fillna('').values)
    is_tmp2 = tmp2.isin(cgc_fusions["GENE_ID"].fillna('').values)
    fusion_calls.loc[is_tmp1 | is_tmp2, "Is Fusion Driver (CGC)"] = 'Yes'
    # add oncokb annotation
    oncokb_fusions, oncokb_fusions_one = utils.read_oncokb_fusions(opts['oncokb'])
    fusion_calls['Is Fusion Driver (OncoKB)'] = 'No'
    is_tmp1 = tmp1.isin(oncokb_fusions.fillna('').values)
    is_tmp2 = tmp2.isin(oncokb_fusions.fillna('').values)
    fusion_calls.loc[is_tmp1 | is_tmp2, "Is Fusion Driver (OncoKB)"] = 'Yes'
    is_tmp1 = fusion_calls["5'_gene"].isin(oncokb_fusions_one)
    is_tmp2 = fusion_calls["3'_gene"].isin(oncokb_fusions_one)
    is_inframe = fusion_calls["Fusion_effect"].fillna('').str.contains('in-frame')
    fusion_calls.loc[(is_tmp1 | is_tmp2) & is_inframe, "Is Fusion Driver (OncoKB)"] = 'Yes'

    # add drugability info
    dgenes = utils.read_drugable_genes(opts['druggability'])
    fusion_calls["Is 5' inhibitor?"] = 'No'
    fusion_calls.loc[fusion_calls["5'_gene"].isin(dgenes['gene_name']), "Is 5' inhibitor?"] = 'Yes'
    fusion_calls["Is 3' inhibitor?"] = 'No'
    fusion_calls.loc[fusion_calls["3'_gene"].isin(dgenes['gene_name']), "Is 3' inhibitor?"] = 'Yes'
    rename_dict = {'gene_name': "5'_gene", "drug_name": "5' inhibitor"}
    fusion_calls = pd.merge(fusion_calls, dgenes.rename(columns=rename_dict),
                            on="5'_gene", how='left')
    rename_dict = {'gene_name': "3'_gene", "drug_name": "3' inhibitor"}
    fusion_calls = pd.merge(fusion_calls, dgenes.rename(columns=rename_dict),
                            on="3'_gene", how='left')

    # add driver point mutation counts
    pancan_flags = pd.read_csv(opts['driver'], sep='\t', index_col=0)
    point_driver_cts = pancan_flags.sum(axis=1).reset_index(name='num_point_drivers').rename(columns={'ID': 'patient_id'})
    fusion_calls = pd.merge(fusion_calls, point_driver_cts, on='patient_id', how='left')

    # add tmb
    maf = pd.read_csv(opts['maf'], sep='\t')
    maf['patient_id'] = maf['Tumor_Sample_Barcode'].str[:12]
    tmb_cts = maf['patient_id'].value_counts().reset_index(name='TMB').rename(columns={'index': 'patient_id'})
    fusion_calls = pd.merge(fusion_calls, tmb_cts, on='patient_id', how='left')

    # Add degron impact info
    is_5prime = degron_impact['type']=='5_prime'
    fusion_calls = pd.merge(fusion_calls, degron_impact.loc[is_5prime, mycols].rename(columns=utils.myrename).drop_duplicates(["5'_gene", "Break1"]),
                            on=["5'_gene", 'Break1'], how='left')
    myrename2 = lambda x: utils.myrename(x, suffix="3'")
    is_3prime = degron_impact['type']=='3_prime'
    fusion_calls = pd.merge(fusion_calls, degron_impact.loc[is_3prime, mycols].rename(columns=myrename2).drop_duplicates(["3'_gene", "Break2"]),
                            on=["3'_gene", 'Break2'], how='left')

    # add protein domain info
    domains = utils.read_fusion_domains(opts['fusion_domain'], opts['wildtype_domain'], fusion_annot)
    fusion_calls = pd.merge(fusion_calls, domains, on='ID', how='left')

    # figure out whether domains were lost from the wt sequence
    fusion_calls = utils.merge_lost_prot_domains(fusion_calls, opts['wildtype_domain'])

    # add ubiquitin site data
    ubiq = pd.read_csv(opts['ubiquitin'], sep='\t')
    ubiq = ubiq[ubiq['ORGANISM']=='human']
    ubiq['POS'] = ubiq['MOD_RSD'].str.extract('([0-9]+)-ub', expand=True)[0].astype(int)
    # read in ensembl to uniprot cross-ref
    uniprot_xref = pd.read_csv(opts['cross_ref'], sep='\t',
                               header=None, names=['UniprotID', 'Type', 'XrefID'])
    uniprot_xref = uniprot_xref[uniprot_xref['Type']=='Ensembl_PRO']
    ensToUni = dict(uniprot_xref[['XrefID', 'UniprotID']].values)
    # get prot IDs from ID column
    prot_ids = fusion_calls['ID'].str.split(':', expand=True)[0]
    tmp_ids = prot_ids.str.split('-', expand=True)
    prot_id1, prot_id2 = tmp_ids[0], tmp_ids[1]
    num_ub1, num_ub2, total_ub1, total_ub2 = [], [], [], []
    # iterate through each fusion
    for ix, row in fusion_calls.iterrows():
        prot1, prot2 = prot_id1.iloc[ix], prot_id2.iloc[ix]
        if prot1=='' and prot2=='':
            num_ub1.append(0)
            num_ub2.append(0)
            total_ub1.append(0)
            total_ub2.append(0)
            continue
        uniprot1, uniprot2 = ensToUni.get(prot1, ''), ensToUni.get(prot2, '')
        if uniprot1:
            tmp_ub = ubiq[(ubiq['ACC_ID']==uniprot1)]
            total_ub1.append(len(tmp_ub))
            num_ub1.append(len(tmp_ub[(tmp_ub['POS']<=row['CodonPos1'])]))
        else:
            num_ub1.append(0)
            total_ub1.append(0)
        if uniprot2:
            tmp_ub = ubiq[(ubiq['ACC_ID']==uniprot2)]
            total_ub2.append(len(tmp_ub))
            num_ub2.append(len(tmp_ub[(tmp_ub['POS']>=row['CodonPos2'])]))
        else:
            num_ub2.append(0)
            total_ub2.append(0)
    fusion_calls['NumUbiqSites_Fusion1'] = num_ub1
    fusion_calls['NumUbiqSites_Fusion2'] = num_ub2
    fusion_calls['TotalUbiqSites_Prot1'] = total_ub1
    fusion_calls['TotalUbiqSites_Prot2'] = total_ub2

    # add translocation info
    fusion_calls['Is translocation?'] = ''
    is_null = fusion_calls['chrom1'].isnull() | fusion_calls['chrom2'].isnull()
    is_same_chrom = fusion_calls['chrom1'] == fusion_calls['chrom2']
    fusion_calls.loc[~is_null & is_same_chrom, 'Is translocation?'] = 'no'
    fusion_calls.loc[~is_null & ~is_same_chrom, 'Is translocation?'] = 'yes'

    # add cterm
    mycols = ['ID', 'fusion degron potential', 'first degron potential',
              'second degron potential', 'delta degron potential']
    cterm = pd.read_csv(opts['cterm_degron'], sep='\t')
    fusion_calls = pd.merge(fusion_calls, cterm[mycols], on='ID', how='left')

    # save file
    fusion_calls.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
