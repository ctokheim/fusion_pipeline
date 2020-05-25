"""
File: regex_degron_search.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Find COP1 degron motifs
"""
import argparse
import re
#from Bio import SeqIO
import csv
import pandas as pd

def parse_arguments():
    info = 'Find COP1 degron motifs'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Input file for analysis')
    parser.add_argument('-m', '--motif',
                        type=str, required=True,
                        help='File containing regex motifs for degrons')
    parser.add_argument('-u', '--uniprot-xref',
                        type=str, required=True,
                        help='Uniprot cross reference')
    parser.add_argument('-p', '--phospho',
                        type=str, required=True,
                        help='Phosphosites from phosphositeplusdb')
    parser.add_argument('-a', '--acetyl',
                        type=str, required=True,
                        help='Acetylation from phosphositeplusdb')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output results')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in sequence file
    df = pd.read_csv(opts['input'], sep='\t')
    df = df.dropna(subset=['protein_sequence'])
    df = df.rename(columns={'protein_id': 'ID'})
    # read uniprot to ensembl data
    uniprot_xref = pd.read_csv(opts['uniprot_xref'], sep='\t',
                               header=None, names=['UniprotID', 'Type', 'XrefID'])
    uniprot_xref = uniprot_xref[uniprot_xref['Type']=='Ensembl_PRO']
    ensToUni = dict(uniprot_xref[['XrefID', 'UniprotID']].values)
    # read motifs
    motif = pd.read_csv(opts['motif'], sep='\t')
    # read phosphorylation sites
    phospho = pd.read_csv(opts['phospho'], sep='\t')
    phospho = phospho[phospho['ORGANISM']=='human']
    phospho['POS'] = phospho['MOD_RSD'].str.extract('([0-9]+)-p', expand=True)[0].astype(int)
    # read acetylation sites
    acetyl = pd.read_csv(opts['acetyl'], sep='\t', skiprows=3)
    acetyl = acetyl[acetyl['ORGANISM']=='human']
    acetyl['POS'] = acetyl['MOD_RSD'].str.extract('([0-9]+)-ac', expand=True)[0].astype(int)

    # regex search
    output_list = [['motif', 'gene', 'transcript_id', 'ID', 'uniprot_id', 'start', 'end']]
    for seq_ix, seq_row in df.iterrows():
        gene, tx_id = seq_row['gene'], seq_row['transcript_id']
        prot_id = seq_row['ID']
        uniprot_id = ensToUni.get(prot_id, '')
        myseq = seq_row['protein_sequence']

        for ix, row in motif.iterrows():
            motif_id, regex = row['Name'], row['Degron']
            num_ptm = len(re.findall('\(', regex))
            for hit in re.finditer(regex, myseq):
                ptm_ct = 0
                # make sure phosphorylation matches, if necessary
                if num_ptm:
                    # extract appropriate ptm mark
                    if row['E3'] != 'CRBN':
                        tmp_ptm = phospho[phospho['ACC_ID']==uniprot_id]
                    else:
                        tmp_ptm = acetyl[acetyl['ACC_ID']==uniprot_id]

                    tmp_flag = False
                    for i in range(1, num_ptm+1):
                        tmp_pos = hit.span(i)[1]  # pick end so adjust for zero-to-one coords
                        if tmp_pos not in tmp_ptm['POS'].values: tmp_flag = True
                    if tmp_flag:
                        continue
                start, end = hit.start(), hit.end()
                output_list.append([motif_id, gene, tx_id, prot_id, uniprot_id, start, end, ])

    # save output
    with open(opts['output'], 'w') as handle:
        mywriter = csv.writer(handle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(output_list)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)



