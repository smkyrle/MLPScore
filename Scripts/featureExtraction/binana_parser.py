'''
Parse binana tesxt files to a single csv file
'''

import pandas as pd
import os, io, sys
import concurrent.futures
from tqdm import tqdm

destination_path = None


def read_binana(binana_file):
    # read in binana output text file
    binana_file = open(binana_file, 'r')

    # cut out license and prose text
    binana_stats = binana_file.read().split('[- END INTRO -]')[1]

    # define headers for parsing
    file_headers = ['Atom-type pair counts within 2.5 angstroms:',
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

    # define dictionary for populating
    tables = dict()

    # extract markdown tables from text file and add to tables dict
    for index, header in enumerate(file_headers):
        column_names = dict()
        if header == 'Number of rotatable bonds in the ligand:':
            table = binana_stats.split(header)[1].split(file_headers[index + 1])[0].strip().lstrip()
        elif index < (len(file_headers) - 1):
            table = binana_stats.split(header)[1].split(file_headers[index + 1])[0].strip().lstrip()
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

def add_data(file):

    tables = read_binana(file)

    name = file.split('/')

    name = name[len(name)-1].replace('.txt','')

    df_headers = [
                'Name',
                '2.5 (HD, OA)', '2.5 (A, HD)', '2.5 (HD, NA)', '2.5 (C, HD)', '2.5 (OA, ZN)', '2.5 (HD, N)', '2.5 (HD, HD)', '2.5 (HD, ZN)', '2.5 (N, ZN)', '2.5 (OA, OA)', '2.5 (N, OA)', '2.5 (CO, HD)', '2.5 (NA, OA)', '2.5 (CL, HD)', '2.5 (MN, OA)', '2.5 (FE, OA)', '2.5 (HD, P)', '2.5 (F, HD)', '2.5 (HD, SA)', '2.5 (CA, HD)', '2.5 (CA, OA)', '2.5 (N, N)', '2.5 (C, OA)', '2.5 (C, NA)', '2.5 (C, C)', '2.5 (C, N)', '2.5 (A, C)', '2.5 (A, OA)', '2.5 (HD, MG)', '2.5 (MG, OA)', '2.5 (HD, NI)', '2.5 (C, ZN)', '2.5 (HD, MN)', '2.5 (CA, NA)', '2.5 (MN, N)', '2.5 (MG, N)', '2.5 (MG, NA)', '2.5 (A, N)', '2.5 (CL, N)', '2.5 (SA, SA)', '2.5 (SA, ZN)', '2.5 (NI, OA)', '2.5 (NA, ZN)', '2.5 (C, CA)', '2.5 (FE, HD)', '2.5 (A, CL)', '2.5 (A, SA)', '2.5 (C, P)', '2.5 (CL, ZN)', '2.5 (CO, OA)', '2.5 (CU, OA)', '2.5 (BR, HD)', '2.5 (F, OA)', '2.5 (CU, HD)', '2.5 (CO, N)', '2.5 (N, SA)', '2.5 (C, SA)', '2.5 (CL, OA)', '2.5 (C, MG)', '2.5 (C, CO)', '2.5 (A, ZN)', '2.5 (FE, NA)', '2.5 (A, A)', '2.5 (A, F)', '2.5 (OA, SA)', '2.5 (N, NA)', '2.5 (MN, NA)', '2.5 (CL, MG)', '2.5 (C, MN)', '2.5 (HD, S)', '2.5 (C, F)', '2.5 (C, CL)', '2.5 (CD, OA)', '2.5 (CU, N)', '2.5 (CU, NA)', '2.5 (F, N)', '2.5 (CU, F)', '2.5 (A, CU)', '2.5 (FE, N)', '2.5 (CD, NA)', '2.5 (C, CU)', '2.5 (A, MG)', '2.5 (C, NI)', '2.5 (C, CD)', '2.5 (N, NI)', '2.5 (K, OA)', '2.5 (C, K)', '2.5 (NA, NI)', '2.5 (CD, HD)', '2.5 (S, ZN)', '2.5 (F, ZN)', '2.5 (MG, P)', '2.5 (N, P)', '2.5 (F, SA)', '2.5 (OA, S)', '2.5 (HD, K)', '2.5 (OA, P)', '2.5 (C, FE)', '2.5 (CO, NA)', '2.5 (MG, SA)', '2.5 (HD, I)', '2.5 (CA, N)', '2.5 (F, MG)', '4.0 (C, C)', '4.0 (A, OA)', '4.0 (OA, SA)', '4.0 (HD, N)', '4.0 (C, OA)', '4.0 (HD, OA)', '4.0 (C, N)', '4.0 (N, N)', '4.0 (A, C)', '4.0 (N, OA)', '4.0 (OA, OA)', '4.0 (A, HD)', '4.0 (HD, HD)', '4.0 (C, HD)', '4.0 (C, F)', '4.0 (A, A)', '4.0 (A, N)', '4.0 (F, HD)', '4.0 (F, OA)', '4.0 (C, NA)', '4.0 (FE, NA)', '4.0 (N, NA)', '4.0 (FE, OA)', '4.0 (C, FE)', '4.0 (HD, NA)', '4.0 (C, P)', '4.0 (N, ZN)', '4.0 (F, N)', '4.0 (C, ZN)', '4.0 (HD, ZN)', '4.0 (OA, ZN)', '4.0 (HD, SA)', '4.0 (C, SA)', '4.0 (N, SA)', '4.0 (A, SA)', '4.0 (N, S)', '4.0 (A, ZN)', '4.0 (HD, S)', '4.0 (A, S)', '4.0 (OA, S)', '4.0 (S, ZN)', '4.0 (A, CL)', '4.0 (A, F)', '4.0 (A, NA)', '4.0 (NA, OA)', '4.0 (OA, P)', '4.0 (P, ZN)', '4.0 (HD, P)', '4.0 (C, S)', '4.0 (N, P)', '4.0 (F, SA)', '4.0 (F, ZN)', '4.0 (NA, ZN)', '4.0 (C, CO)', '4.0 (A, CO)', '4.0 (CO, HD)', '4.0 (CO, N)', '4.0 (CO, OA)', '4.0 (C, CL)', '4.0 (CL, HD)', '4.0 (CL, OA)', '4.0 (CL, N)', '4.0 (SA, SA)', '4.0 (C, MN)', '4.0 (MN, N)', '4.0 (MN, NA)', '4.0 (HD, MN)', '4.0 (A, FE)', '4.0 (BR, N)', '4.0 (BR, HD)', '4.0 (A, BR)', '4.0 (BR, OA)', '4.0 (MG, OA)', '4.0 (C, MG)', '4.0 (C, I)', '4.0 (I, SA)', '4.0 (A, I)', '4.0 (A, P)', '4.0 (MG, P)', '4.0 (MG, N)', '4.0 (HD, MG)', '4.0 (HD, NI)', '4.0 (NI, OA)', '4.0 (SA, ZN)', '4.0 (C, CA)', '4.0 (CA, HD)', '4.0 (CA, N)', '4.0 (BR, C)', '4.0 (S, SA)', '4.0 (CA, OA)', '4.0 (MN, OA)', '4.0 (CA, P)', '4.0 (NA, SA)', '4.0 (A, NI)', '4.0 (I, OA)', '4.0 (MG, SA)', '4.0 (A, MN)', '4.0 (CL, SA)', '4.0 (F, MG)', '4.0 (A, MG)', '4.0 (P, SA)', '4.0 (MN, P)', '4.0 (N, NI)', '4.0 (C, NI)', '4.0 (A, CA)', '4.0 (BR, SA)', '4.0 (FE, N)', '4.0 (FE, HD)', '4.0 (MG, NA)', '4.0 (CA, F)', '4.0 (FE, P)', '4.0 (CL, CU)', '4.0 (C, CU)', '4.0 (CO, P)', '4.0 (NA, NI)', '4.0 (CU, N)', '4.0 (CU, OA)', '4.0 (CU, HD)', '4.0 (CA, NA)', '4.0 (MG, S)', '4.0 (CL, ZN)', '4.0 (CU, NA)', '4.0 (HD, I)', '4.0 (I, N)', '4.0 (CU, P)', '4.0 (NA, P)', '4.0 (F, NA)', '4.0 (F, MN)', '4.0 (K, OA)', '4.0 (K, P)', '4.0 (C, K)', '4.0 (CD, HD)', '4.0 (CD, N)', '4.0 (C, CD)', '4.0 (CD, OA)', '4.0 (CA, SA)', '4.0 (HD, K)', '4.0 (K, N)', '4.0 (CL, MG)', '4.0 (CD, P)', '4.0 (CA, S)', '4.0 (A, CU)', '4.0 (NA, NA)', '4.0 (CU, F)', '4.0 (NI, SA)', '4.0 (CL, MN)', '4.0 (CL, NA)', '4.0 (CO, NA)', '4.0 (BR, NA)', '4.0 (A, CD)', '4.0 (NA, S)', '4.0 (CL, FE)', '4.0 (BR, FE)', '4.0 (K, NA)', '4.0 (BR, ZN)', '4.0 (CD, NA)', '4.0 (A, K)', '4.0 (CO, SA)', '4.0 (CA, CL)', '4.0 (FE, SA)', '4.0 (FE, S)', '4.0 (K, SA)', '4.0 (MN, SA)', '4.0 (MN, S)', '4.0 (BR, MN)', '4.0 (CD, S)', 'LA C', 'LA OA', 'LA HD', 'LA N', 'LA A', 'LA F', 'LA NA', 'LA P', 'LA SA', 'LA S', 'LA CL', 'LA BR', 'LA I', 'ElSum (C, C)', 'ElSum (A, OA)', 'ElSum (OA, SA)', 'ElSum (HD, N)', 'ElSum (C, OA)', 'ElSum (HD, OA)', 'ElSum (C, N)', 'ElSum (N, N)', 'ElSum (A, C)', 'ElSum (N, OA)', 'ElSum (OA, OA)', 'ElSum (A, HD)', 'ElSum (HD, HD)', 'ElSum (C, HD)', 'ElSum (C, F)', 'ElSum (A, A)', 'ElSum (A, N)', 'ElSum (F, HD)', 'ElSum (F, OA)', 'ElSum (C, NA)', 'ElSum (FE, NA)', 'ElSum (N, NA)', 'ElSum (FE, OA)', 'ElSum (C, FE)', 'ElSum (HD, NA)', 'ElSum (C, P)', 'ElSum (N, ZN)', 'ElSum (OA, ZN)', 'ElSum (F, N)', 'ElSum (C, ZN)', 'ElSum (HD, ZN)', 'ElSum (HD, SA)', 'ElSum (C, SA)', 'ElSum (N, SA)', 'ElSum (A, SA)', 'ElSum (N, S)', 'ElSum (A, ZN)', 'ElSum (HD, S)', 'ElSum (A, S)', 'ElSum (OA, S)', 'ElSum (S, ZN)', 'ElSum (A, CL)', 'ElSum (A, F)', 'ElSum (A, NA)', 'ElSum (NA, OA)', 'ElSum (OA, P)', 'ElSum (P, ZN)', 'ElSum (HD, P)', 'ElSum (C, S)', 'ElSum (N, P)', 'ElSum (F, SA)', 'ElSum (F, ZN)', 'ElSum (NA, ZN)', 'ElSum (C, CO)', 'ElSum (A, CO)', 'ElSum (CO, HD)', 'ElSum (CO, N)', 'ElSum (CO, OA)', 'ElSum (C, CL)', 'ElSum (CL, HD)', 'ElSum (CL, OA)', 'ElSum (CL, N)', 'ElSum (SA, SA)', 'ElSum (C, MN)', 'ElSum (MN, N)', 'ElSum (MN, OA)', 'ElSum (MN, NA)', 'ElSum (HD, MN)', 'ElSum (A, FE)', 'ElSum (BR, N)', 'ElSum (BR, HD)', 'ElSum (A, BR)', 'ElSum (BR, OA)', 'ElSum (MG, OA)', 'ElSum (C, MG)', 'ElSum (C, I)', 'ElSum (I, SA)', 'ElSum (A, I)', 'ElSum (A, P)', 'ElSum (MG, P)', 'ElSum (MG, N)', 'ElSum (HD, MG)', 'ElSum (HD, NI)', 'ElSum (NI, OA)', 'ElSum (SA, ZN)', 'ElSum (C, CA)', 'ElSum (CA, HD)', 'ElSum (CA, N)', 'ElSum (BR, C)', 'ElSum (S, SA)', 'ElSum (CA, OA)', 'ElSum (CA, P)', 'ElSum (NA, SA)', 'ElSum (A, NI)', 'ElSum (I, OA)', 'ElSum (MG, SA)', 'ElSum (A, MN)', 'ElSum (CA, NA)', 'ElSum (CL, SA)', 'ElSum (MG, NA)', 'ElSum (F, MG)', 'ElSum (A, MG)', 'ElSum (P, SA)', 'ElSum (MN, P)', 'ElSum (N, NI)', 'ElSum (C, NI)', 'ElSum (A, CA)', 'ElSum (BR, SA)', 'ElSum (FE, N)', 'ElSum (FE, HD)', 'ElSum (CA, F)', 'ElSum (FE, P)', 'ElSum (CL, CU)', 'ElSum (C, CU)', 'ElSum (CL, ZN)', 'ElSum (CO, P)', 'ElSum (NA, NI)', 'ElSum (CU, N)', 'ElSum (CU, OA)', 'ElSum (CU, HD)', 'ElSum (MG, S)', 'ElSum (CU, NA)', 'ElSum (HD, I)', 'ElSum (I, N)', 'ElSum (CU, P)', 'ElSum (NA, P)', 'ElSum (F, NA)', 'ElSum (F, MN)', 'ElSum (CL, MG)', 'ElSum (K, OA)', 'ElSum (K, P)', 'ElSum (C, K)', 'ElSum (CD, HD)', 'ElSum (CD, N)', 'ElSum (C, CD)', 'ElSum (CD, OA)', 'ElSum (CA, SA)', 'ElSum (HD, K)', 'ElSum (K, N)', 'ElSum (CD, P)', 'ElSum (CA, S)', 'ElSum (A, CU)', 'ElSum (NA, NA)', 'ElSum (CU, F)', 'ElSum (NI, SA)', 'ElSum (CD, NA)', 'ElSum (CL, MN)', 'ElSum (CL, NA)', 'ElSum (CO, NA)', 'ElSum (BR, NA)', 'ElSum (A, CD)', 'ElSum (NA, S)', 'ElSum (CL, FE)', 'ElSum (BR, FE)', 'ElSum (K, NA)', 'ElSum (BR, ZN)', 'ElSum (A, K)', 'ElSum (CO, SA)', 'ElSum (CA, CL)', 'ElSum (FE, SA)', 'ElSum (FE, S)', 'ElSum (K, SA)', 'ElSum (MN, SA)', 'ElSum (MN, S)', 'ElSum (BR, MN)', 'ElSum (CD, S)',
                'BPF ALPHA SIDECHAIN', 'BPF ALPHA BACKBONE', 'BPF BETA SIDECHAIN', 'BPF BETA BACKBONE', 'BPF OTHER SIDECHAIN', 'BPF OTHER BACKBONE',
                'HC ALPHA SIDECHAIN', 'HC ALPHA BACKBONE', 'HC BETA SIDECHAIN', 'HC BETA BACKBONE', 'HC OTHER SIDECHAIN', 'HC OTHER BACKBONE',
                'HB ALPHA SIDECHAIN LIGAND', 'HB ALPHA BACKBONE LIGAND', 'HB BETA SIDECHAIN LIGAND', 'HB BETA BACKBONE LIGAND', 'HB OTHER SIDECHAIN LIGAND', 'HB OTHER BACKBONE LIGAND',
                'HB ALPHA SIDECHAIN RECEPTOR', 'HB ALPHA BACKBONE RECEPTOR', 'HB BETA SIDECHAIN RECEPTOR', 'HB BETA BACKBONE RECEPTOR', 'HB OTHER SIDECHAIN RECEPTOR', 'HB OTHER BACKBONE RECEPTOR',
                'SB ALPHA', 'SB BETA', 'SB OTHER',
                'piStack ALPHA', 'piStack BETA', 'piStack OTHER',
                'tStack ALPHA', 'tStack BETA', 'tStack OTHER',
                'catPi ALPHA LIGAND', 'catPi BETA LIGAND', 'catPi OTHER LIGAND',
                'catPi ALPHA RECEPTOR', 'catPi BETA RECEPTOR', 'catPi OTHER RECEPTOR',
                'nRot'
            ]

    binana_data = pd.DataFrame(columns = df_headers)

    data = {'Name': str(name)}

    for key, value in tables.items():

        if key == 'Atom-type pair counts within 2.5 angstroms:':
            for index in value.index:
                atom1 = (value.iloc[index-1,0]).lstrip().rstrip()
                atom2 = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'2.5 ({atom1}, {atom2})'] = tally


        elif key == 'Atom-type pair counts within 4.0 angstroms:':
            for index in value.index:
                atom1 = (value.iloc[index-1,0]).lstrip().rstrip()
                atom2 = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'4.0 ({atom1}, {atom2})'] = tally

        elif key == 'Ligand atom types:':
            for index in value.index:
                atom = (value.iloc[index-1,0]).lstrip().rstrip()
                data[f'LA {atom}'] = 1

        elif key == 'Summed electrostatic energy by atom-type pair, in J/mol:':
            for index in value.index:
                atom1 = (value.iloc[index-1,0]).lstrip().rstrip()
                atom2 = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = float(value.iloc[index-1,2])
                data[f'ElSum ({atom1}, {atom2})'] = tally

        elif key == 'Active-site flexibility:':
            for index in value.index:
                chain = (value.iloc[index-1,0]).lstrip().rstrip()
                struc = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'BPF {struc} {chain}'] = tally

        elif key == 'Hydrogen bonds:':
            for index in value.index:
                donor = (value.iloc[index-1,0]).lstrip().rstrip()
                chain = (value.iloc[index-1,1]).lstrip().rstrip()
                struc = (value.iloc[index-1,2]).lstrip().rstrip()
                tally = int(value.iloc[index-1,3])
                data[f'HB {struc} {chain} {donor}'] = tally

        elif key == 'Hydrophobic contacts (C-C):':
            for index in value.index:
                chain = (value.iloc[index-1,0]).lstrip().rstrip()
                struc = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'HC {struc} {chain}'] = tally

        elif key == 'pi-pi stacking interactions:':
            for index in value.index:
                struc = (value.iloc[index-1,0]).lstrip().rstrip()
                tally = int(value.iloc[index-1,1])
                data[f'piStack {struc}'] = tally

        elif key == 'T-stacking (face-to-edge) interactions:':
            for index in value.index:
                struc = (value.iloc[index-1,0]).lstrip().rstrip()
                tally = int(value.iloc[index-1,1])
                data[f'tStack {struc}'] = tally

        elif key == 'Cation-pi interactions:':
            for index in value.index:
                char = (value.iloc[index-1,0]).lstrip().rstrip()
                struc = (value.iloc[index-1,1]).lstrip().rstrip()
                tally = int(value.iloc[index-1,2])
                data[f'catPi {struc} {char}'] = tally

        elif key == 'Salt Bridges:':
            for index in value.index:
                struc = (value.iloc[index-1,0]).lstrip().rstrip()
                tally = int(value.iloc[index-1,1])
                data[f'SB {struc}'] = tally

        elif key == 'Number of rotatable bonds in the ligand:':
            data['nRot'] = value

    binana_data = binana_data.append(data, ignore_index = True)

    binana_data = binana_data.fillna(0)

    #print(os.path.isfile(destination_path))
    if len(binana_data.columns) != 489:
        print('Col Num')
        print(len(binana_data.columns))
        for i in list(binana_data.columns):
            if i not in df_headers:
                print(i)

    if os.path.isfile(destination_path) == False:
        binana_data.to_csv(destination_path, index = False)

    else:
        binana_data.to_csv(destination_path, mode = 'a', index = False, header = False)

def parse_args(args): # parse CLI user inputs

    binana_files_path = args[args.index('-loc') + 1]

    destination_path = args[args.index('-csv') + 1]

    return binana_files_path, destination_path

def main():

    global destination_path

    binana_files_path, destination_path = parse_args(sys.argv)

    binana_files = [(binana_files_path + text_file) for text_file in os.listdir(binana_files_path)]

    for file in tqdm(binana_files):
        binana_data = add_data(file)


if __name__ == '__main__':
    main()
