"""
File: merge_agfusion_domains.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Merge protein domain info of fusions
"""
import csv
import argparse
import pandas as pd
import os
import glob


def parse_arguments():
    info = 'Merge protein domain info'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input-dir',
                        type=str, required=True,
                        help='Agfusion result directory')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Merged output results')
    args = parser.parse_args()
    return vars(args)


def read_fasta(path):
    """Read in fasta info"""
    seq = ''
    with open(path) as handle:
        # parse out meta info
        header = next(handle).strip()
        prot_ids = header[1:].split(' ')[0]
        tmp_split = header.split(',')[2]
        tx_ids = tmp_split.split(': ')[-1]

        # get full prot sequence
        for line in handle:
            seq += line.strip()

    return prot_ids, tx_ids, seq


def main(opts):
    output_list = []
    columns=["5'_gene","3'_gene", "5'_transcript", "3'_transcript",
             "5'_strand", "3'_strand", "Domain_ID", "Domain_name",
             "Domain_description", "Protein_start", "Protein_end", "ID"]
    column_list = None
    pattern = os.path.join(opts['input_dir'], '*')
    for d in glob.glob(pattern):
        # figure out break points
        mybase = os.path.basename(d)
        break1 = mybase.split('_')[0].split('-')[-1]
        break2 = mybase.split('_')[1].split('-')[-1]

        # read fasta
        prot_fa_paths = glob.glob(os.path.join(d, '*_protein.fa'))
        if prot_fa_paths:
            prot_fa_path = prot_fa_paths[0]
            p_id, tx_id, seq = read_fasta(prot_fa_path)
        else:
            p_id, tx_id, seq = '', '', ''

        # read in fusion info
        fus_info_path = glob.glob(os.path.join(d, '*protein_domains.csv'))[0]
        with open(fus_info_path) as handle:
            myreader = csv.reader(handle, delimiter=',')
            tmp_columns = next(myreader)
            output_list += [list(row)+[p_id+':'+break1+'-'+break2] for row in myreader]

    # fix terrible choice by agfusion to use csv
    output_list2 = []
    for l in output_list:
        mylen = len(l)
        if mylen==12:
            output_list2.append(l)
        elif mylen>12:
            gap = mylen - 12
            begin = l[:8]
            mid = ', '.join(l[8:8+gap+1])
            end = l[8+gap+1:]
            output_list2.append(begin+[mid]+end)

    output_df = pd.DataFrame(output_list2, columns=columns)

    # save results
    output_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

