"""
File: merge_agfusion_result.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Merge results from agfusion
"""
import pandas as pd
import argparse
import glob
import os
import csv
from pyensembl import ensembl_grch37, EnsemblRelease


def parse_arguments():
    info = 'Merge results from agfusion'
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


def get_relative_cds_pos(data, tx_id, pos):
    """Calculate the relative position along the CDS"""
    # get transcript obj
    tx = data.transcript_by_id(tx_id)

    # figure out whether the junction actually hits the tx
    offset_pos = None
    try:
        offset_pos = tx.spliced_offset(pos)
    except ValueError:
        pass
    if offset_pos is None:
        return None

    if tx.protein_sequence and tx.complete:
        # figure out the relative loc
        start_pos = tx.first_start_codon_spliced_offset
        last_pos = tx.last_stop_codon_spliced_offset
        mylen = last_pos - start_pos

        return (offset_pos - start_pos) / float(mylen)
    else:
        return None


def main(opts):
    # load ensembl db
    data = EnsemblRelease(95)

    output_list = []
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
        fus_info_path = glob.glob(os.path.join(d, '*fusion_transcripts.csv'))[0]
        with open(fus_info_path) as handle:
            myreader = csv.reader(handle, delimiter=',')
            tmp_columns = next(myreader)
            if column_list is None:
                column_list = tmp_columns
            tmp_list = next(myreader)
        #tmp_df = pd.read_csv(fus_info_path)
        #tmp_df['ID'] = p_id
        #tmp_df['TX_ID'] = tx_id
        #tmp_df['protein_sequence'] = seq

        # figure out the pos of the break
        relative_pos1 = get_relative_cds_pos(data, tmp_list[2], int(break1))
        relative_pos2 = get_relative_cds_pos(data, tmp_list[3], int(break2))

        # append results
        output_list.append(tmp_list + [p_id+':'+break1+"-"+break2, tx_id, seq, break1, break2, relative_pos1, relative_pos2])

    # merge results
    output_df = pd.DataFrame(output_list, columns=column_list+['ID', 'TX_ID', 'protein_sequence', 'Break1', 'Break2', 'RelativePos1', 'RelativePos2'])

    # add gene ID
    output_df['GENE_ID'] = output_df["5'_gene"] + '--' + output_df["3'_gene"]

    # save results
    output_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


