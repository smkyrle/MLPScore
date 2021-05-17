'''
Batch parallel processing of train.py for faster model training
'''

from itertools import product
from tqdm import tqdm
import os
import subprocess
import sys
import concurrent.futures
from subprocess import PIPE

def run_train(layer):
    # Prepare layer for Unix command line
    num = layer[1] # Identifying number for model
    layer = str(layer[0])
    layer = layer.strip('[]')
    layer = layer.replace(' ','')

    try:
        # Run train.py
        output = subprocess.check_output(f'python {model_script} -h5 {h5_path} -layers {layer} -des {output_dir} -num {num}',shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(e.output)
        print(f'FATAL error with {layer}...')

def multithread(layers, num_threads): # Parallel processing for running multiple instances of train.py at once
    threads = min(int(num_threads), len(layers))

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        list(tqdm(executor.map(run_train, layers), total=len(layers)))

def parse_args(args): # Parse user inputs

    model_script = args[args.index('-train')+1]

    h5_path = args[args.index('-h5') + 1]

    output_dir = args[args.index('-des') + 1]

    layer = args[args.index('-layers') + 1]
    layer = list(map(int,layer))

    models = int(args[args.index('-models') + 1])

    num_threads = int(args[args.index('-threads') + 1])

    return model_script, h5_path, output_dir, layer, models, num_threads

def main():
    global output_dir
    global model_script
    global h5_path

    model_script, h5_path, output_dir, layer, models, num_threads = parse_args(sys.argv)

    # Create list of list of layers with a identifying numbers
    layers = [[layer, i] for i in range(models)]

    multithread(layers, num_threads)


if __name__ == '__main__':
    main()
