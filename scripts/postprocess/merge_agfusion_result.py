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
import math


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


def get_chrom(data, tx_id):
    """Fetch the chromosome for a transcript"""
    # get transcript obj
    tx = data.transcript_by_id(tx_id)
    #return chrom
    return tx.contig


def get_cds_pos(data, tx_id, pos):
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
        return None, None, None

    if tx.protein_sequence and tx.complete:
        # figure out pos of cds
        start_pos = tx.first_start_codon_spliced_offset
        last_pos = tx.last_stop_codon_spliced_offset
        mylen = last_pos - start_pos

        # figure out the pos relative to the CDS
        cds_offset = offset_pos - start_pos
        rel_pos = cds_offset / float(mylen)
        #if cds_offset<0 or cds_offset>mylen:
            #codon_pos = 0
        #else:
        codon_pos = math.ceil((cds_offset+1)/3)
        prot_len = math.ceil((mylen+1)/3)

        return codon_pos, rel_pos, prot_len
    else:
        return None, None, None


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

        # figure out the pos of the break
        codon_pos1, relative_pos1, prot_len1 = get_cds_pos(data, tmp_list[2], int(break1))
        codon_pos2, relative_pos2, prot_len2 = get_cds_pos(data, tmp_list[3], int(break2))
        # figure out the chromosome
        chrom1 = get_chrom(data, tmp_list[2])
        chrom2 = get_chrom(data, tmp_list[3])

        # append results
        output_list.append(tmp_list + [p_id+':'+break1+"-"+break2, tx_id, seq,
                                       chrom1, break1, chrom2, break2, codon_pos1, codon_pos2,
                                       relative_pos1, relative_pos2, prot_len1, prot_len2])

    # merge results
    mycols = column_list+['ID', 'TX_ID', 'protein_sequence', 'chrom1', 'Break1', 'chrom2', 'Break2', 'CodonPos1',
                          'CodonPos2', 'RelativePos1', 'RelativePos2', 'ProtLen1', 'ProtLen2']
    output_df = pd.DataFrame(output_list, columns=mycols)

    # add gene ID
    output_df['GENE_ID'] = output_df["5'_gene"] + '--' + output_df["3'_gene"]

    # save results
    output_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


