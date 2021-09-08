"""
File: save_wt_domains.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Save protien domain annotations to text file
"""
from agfusion.database import AGFusionDB
import csv
import argparse


def parse_arguments():
    info = 'Save protein domain annotations to text file'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-d', '--db',
                        type=str, required=True,
                        help='Agfusion database file')
    parser.add_argument('-e', '--ensembl-release',
                        type=int, default=95,
                        help='Ensembl release version')
    parser.add_argument('-p', '--protein-databases',
                        type=str, default='pfam,tmhmm',
                        help='Protein databases to use (separated by comma)')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in all prot domains
    mydb = AGFusionDB(opts['db'])
    prot_databases = opts['protein_databases'].split(',')
    output_list = [['translate_id', 'PROT_ID', 'Domain_ID', 'Protein_start', 'Protein_end', 'Domain_name', 'Domain_description', 'database']]
    for database in prot_databases:
        sqlite3_command = "SELECT * FROM homo_sapiens_{}_{}".format(opts['ensembl_release'], database)
        mydb.sqlite3_cursor.execute(sqlite3_command)
        output_list += [list(x) + [database] for x in mydb.sqlite3_cursor.fetchall()]

    # save output
    with open(opts['output'], 'w') as handle:
        mywriter = csv.writer(handle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(output_list)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

