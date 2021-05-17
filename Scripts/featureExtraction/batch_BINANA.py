'''
Script for running binana.py over multiple receptor-ligand pairs
'''

import os, shutil, subprocess, sys
from tqdm import tqdm
import concurrent.futures
from subprocess import PIPE

fatal_error_list = list()
errors = dict()
binana_path = None
destination_path = None
empty_dir = list()

def run_binana(file, ligand_path, receptor_path, ligand): # Running binana.py over a receptor-ligand pair

    try:
        output = subprocess.check_output(f'python {binana_path} -receptor {receptor_path} -ligand {ligand_path} > {destination_path}{ligand}.txt', shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        errors[ligand] = e.output
        fatal_error_list.append(ligand)
        if os.path.isfile('/home/sammk/Desktop/fatal_error_list.txt'):
            os.remove('/home/sammk/Desktop/fatal_error_list.txt')
        with open('/home/sammk/Desktop/fatal_error_list.txt','a+') as error_list:
            error_list.write(str(fatal_error_list) + '\n' + str(errors))
            error_list.close()
        print(f'FATAL PROBLEM WITH {file}: Added to list and skipping...')

def file_prep_and_run(file): # Check directory for receptor and ligand file
    x = True
    folder_name = file.split('/')
    folder_name = folder_name[len(folder_name) - 1]
    try:
        ligand_file = [filename for filename in os.listdir(file) if 'ligand.pdbqt' in filename][0]
        ligand_path = file + '/' + ligand_file
        ligand = ligand_file.replace('.pdbqt','')
    except:
        print('No ligand found for:')
        print(folder_name)
        empty_dir.append(folder_name)
        x = False

    try:
        receptor_file = [filename for filename in os.listdir(file) if 'receptor.pdbqt' in filename][0]
        receptor_path = file + '/' + receptor_file
    except:
        print('No receptor found for:')
        print(folder_name)
        empty_dir.append(folder_name)
        x = False

    if x == True:
        run_binana(folder_name, ligand_path, receptor_path, ligand)

def multithread_run(pdbqt_files,pdbqt_files_path ,threads): # Parallel processing for running multiple instances of binana at once
    threads = min(threads, len(pdbqt_files))
    print(f'Running BINANA on {pdbqt_files_path}...')
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        list(tqdm(executor.map(file_prep_and_run, pdbqt_files), total=len(pdbqt_files)))

def move_problematic_structures(fatal_error_list, pdbqt_files_path, problem_path): # Move erroneous structures to a new directory
    folders = os.listdir(pdbqt_files_path)
    fatal_error_list = list(set(list(fatal_error_list)))
    with tqdm(total=len(fatal_error_list)) as pbar:
        for file in fatal_error_list:
            folder = [folder for folder in folders if folder in file][0]
            shutil.copytree(f'{pdbqt_files_path}{folder}', f'{problem_path}{folder}')
            shutil.copy(f'{destination_path}{file}.txt',f'{problem_path}{folder}/{file}.txt')
            os.remove(f'{destination_path}{file}.txt')
            pbar.update(1)

def parse_args(args): # Parse user inputs

    binana_path = args[args.index('-binana') + 1]

    pdbqt_files_path = args[args.index('-loc') + 1]

    destination_path = args[args.index('-suc') + 1]

    problem_path = args[args.index('-prob') + 1]

    threads = int(args[args.index('-threads') + 1])

    return binana_path, pdbqt_files_path, destination_path, problem_path, threads

def main():

    global destination_path
    global binana_path

    binana_path, pdbqt_files_path, destination_path, problem_path, threads = parse_args(sys.argv)

    pdbqt_folders = [(pdbqt_files_path + folder) for folder in os.listdir(pdbqt_files_path)]

    multithread_run(pdbqt_folders, pdbqt_files_path, threads)

    move_problematic_structures(fatal_error_list, pdbqt_files_path, problem_path)


if __name__ == '__main__':
    main()
