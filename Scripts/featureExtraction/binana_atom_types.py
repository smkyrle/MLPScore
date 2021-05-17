'''
Script for obtaining all the possible binana features across multiple binana output text files
'''

import pandas as pd
import io
import os, sys
import concurrent.futures
from tqdm import tqdm

binana_data = pd.DataFrame(columns = ['Receptor-Ligand Complex', ])
binana_tally_types = dict([
    ('Atom-type pair counts within 2.5 angstroms:', list()),
    ('Atom-type pair counts within 4.0 angstroms:', list()),
    ('Ligand atom types:', list()),
    ('Summed electrostatic energy by atom-type pair, in J/mol:', list())
])

def parse_table(file):
    binana_file = open(file, 'r')
    binana_stats = binana_file.read().split('[- END INTRO -]')[1]

    headers = ['Atom-type pair counts within 2.5 angstroms:',
               'Atom-type pair counts within 4.0 angstroms:',
               'Ligand atom types:',
               'Summed electrostatic energy by atom-type pair, in J/mol:',
               'Number of rotatable bonds in the ligand:',
               'Active-site flexibility:',
               'Hydrogen bonds:',
               'Hydrophobic contacts (C-C):',
               'pi-pi stacking interactions:',
               'T-stacking (face-to-edge) interactions:',
               'Cation-pi interactions:',
               'Salt Bridges:']

    tables = {}
    for index, header in enumerate(headers):
        column_names = dict()
        if header == 'Number of rotatable bonds in the ligand:':
            table = binana_stats.split(header)[1].split(headers[index + 1])[0].strip().lstrip()

        elif index < (len(headers) - 1):
            table = binana_stats.split(header)[1].split(headers[index + 1])[0].strip().lstrip()
            table = table.split('Raw data')[0]
            table = io.StringIO(table)
            table = pd.read_table(table, sep="|")
            for column in list(table):
                column_names[column] = column.lstrip().rstrip()
                table = table.loc[~table[column].str.contains('----')]
            table.rename(columns = column_names, inplace = True)
        else:
            table = binana_stats.split(header)[1].strip().lstrip()
            table = table.split('Raw data')[0]
            table = io.StringIO(table)
            table = pd.read_table(table, sep="|")

            for column in list(table):
                column_names[column] = column.lstrip().rstrip()
                table = table.loc[~table[column].str.contains('----')]
            table.rename(columns = column_names, inplace = True)
        tables[header] = table
    return tables

def atom_tally_types(binana_file):

    try:
        tables = parse_table(binana_file)

        for key, value in binana_tally_types.items():
            if key == 'Ligand atom types:':
                for index in tables[key].index:
                    atom = tables[key]['Atom Type'][index].lstrip().rstrip()
                    atom = f'LA {atom}'
                    if atom in value:
                        pass
                    else:
                        value.append(atom)
                    #print(value)
            elif key == 'Atom-type pair counts within 2.5 angstroms:':

                for index in tables[key].index:
                    atom1 = tables[key].iloc[index-1,0].lstrip().rstrip()
                    atom2 = tables[key].iloc[index-1,1].lstrip().rstrip()
                    atom_pair = f'2.5 ({atom1}, {atom2})'
                    if atom_pair in value:
                        pass
                    else:
                        value.append(atom_pair)

            elif key == 'Atom-type pair counts within 4.0 angstroms:':

                for index in tables[key].index:
                    atom1 = tables[key].iloc[index-1,0].lstrip().rstrip()
                    atom2 = tables[key].iloc[index-1,1].lstrip().rstrip()
                    atom_pair = f'4.0 ({atom1}, {atom2})'
                    if atom_pair in value:
                        pass
                    else:
                        value.append(atom_pair)

            elif key == 'Summed electrostatic energy by atom-type pair, in J/mol:':

                for index in tables[key].index:
                    atom1 = tables[key].iloc[index-1,0].lstrip().rstrip()
                    atom2 = tables[key].iloc[index-1,1].lstrip().rstrip()
                    atom_pair = f'ElSum ({atom1}, {atom2})'
                    if atom_pair in value:
                        pass
                    else:
                        value.append(atom_pair)

    except:
        print(binana_file)

def parse_args(args):
    
    binana_path = args[args.index('-loc') + 1]

    return binana_path

def main():
    global binana_tally_types
    binana_path = parse_args(sys.argv)
    binana_files = [(binana_path + file) for file in os.listdir(binana_path)]
    for file in tqdm(binana_files):
        atom_tally_types(file)
    print(binana_tally_types)

if __name__ == '__main__':
    main()
