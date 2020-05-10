"""
File: regex_degron_search.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Find COP1 degron motifs
"""
import argparse
import re
from Bio import SeqIO
import csv
import pandas as pd

def parse_arguments():
    info = 'Find COP1 degron motifs'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='FASTA file')
    parser.add_argument('-m', '--motif',
                        type=str, required=True,
                        help='File containing regex motifs for degrons')
    parser.add_argument('-p', '--phospho',
                        type=str, required=True,
                        help='Phosphosites from phosphositeplusdb')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output results')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read motifs
    motif = pd.read_csv(opts['motif'], sep='\t')
    motif = motif[motif['ELMIdentifier'].str.startswith('DEG')]
    #motif = motif[motif['ELMIdentifier']=='DEG_SCF_TRCP1_1']
    # read phosphorylation sites
    phospho = pd.read_csv(opts['phospho'], sep='\t')
    phospho = phospho[phospho['ORGANISM']=='human']
    phospho['POS'] = phospho['MOD_RSD'].str.extract('([0-9]+)-p', expand=True)[0].astype(int)

    # regex search
    output_list = [['MOTIF', 'NAME', 'ACC', 'START', 'END']]
    for record in SeqIO.parse(opts['input'], 'fasta'):
        uniprot_id = record.name.split('|')[1]
        name = record.name.split('|')[2]
        myseq = str(record.seq)
        for ix, row in motif.iterrows():
            motif_id, regex = row['ELMIdentifier'], row['Regex']
            num_phospho = len(re.findall('\(', regex))
            for hit in re.finditer(regex, myseq):
                phospho_ct = 0
                # make sure phosphorylation matches, if necessary
                if num_phospho:
                    tmp_phospho = phospho[phospho['ACC_ID']==uniprot_id]
                    tmp_flag = False
                    #if uniprot_id == 'P35222': import IPython ; IPython.embed() ; raise
                    for i in range(1, num_phospho+1):
                        tmp_pos = hit.span(i)[1]  # pick end so adjust for zero-to-one coords
                        if tmp_pos not in tmp_phospho['POS'].values: tmp_flag = True
                        #if tmp_pos in tmp_phospho['POS'].values: phospho_ct += 1
                    if tmp_flag:
                        continue
                start, end = hit.start(), hit.end()
                output_list.append([motif_id, name, uniprot_id, start, end, ])

    # save output
    with open(opts['output'], 'w') as handle:
        mywriter = csv.writer(handle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(output_list)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)



