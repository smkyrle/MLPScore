'''
Script for training multiple architecture structures
'''
import os
# Configure GPU memory for tensorflow
os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
from tensorflow.compat.v1.keras.backend import set_session
config = tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True
sess = tf.compat.v1.Session(config=config)
set_session(sess)

from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.utils import class_weight
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.compat.v1.keras.backend import clear_session
from keras.activations import relu, sigmoid
from keras.regularizers import l2
import warnings
import random
import sys
from tensorflow.python.client import device_lib

sys.setrecursionlimit(10000)
pd.set_option('display.max_rows',None)
warnings.filterwarnings('ignore')
pd.set_option('display.max_columns', None)

def load_data_with_validation_sets(h5_path): # fetch x_train and y_train and validation sets from hdf5 file

    # load dataframe
    data = pd.read_hdf(h5_path, key='df', mode='r')

    # make label for weak crystal binders to mitigate any decoy bias
    data['WeakCrystal'] = np.where(((~data['PDBCode'].str.contains('decoy')) & (data['Label'] == 0)), 1, 0)

    # set aside a small validation set stratified by class label and weak crystal binders
    splitter = StratifiedShuffleSplit(n_splits = 1, test_size=0.1, random_state=12)
    for train_index, val_index in splitter.split(data, data[['Label','WeakCrystal']]):
        training_set = data.iloc[train_index]
        val_set = data.iloc[val_index]

    # split features and identifier columns
    x_train = training_set.drop(['PDBCode', 'BindingDataType', 'BindingValue', 'BindingUnits', 'Database', 'Binding_Data_in_uM', 'Label','WeakCrystal'], axis=1)
    x_val = val_set.drop(['PDBCode', 'BindingDataType', 'BindingValue', 'BindingUnits', 'Database', 'Binding_Data_in_uM', 'Label','WeakCrystal'], axis=1)

    # split labels and ensure integer type
    y_train = training_set['Label'].fillna(0).astype(int)
    y_val = val_set['Label'].fillna(0).astype(int)

    print(f'Training on {len(x_train)} instances')
    print(f'Validating on {len(x_val)} instances')

    x_train.to_numpy
    y_train.to_numpy
    x_val.to_numpy
    y_val.to_numpy

    return x_train, y_train, x_val, y_val

def build_model(nodes, num_inputs): # input list of hidden layers and return a compiled model

    # Create model
    model = Sequential()
    # Add input layer dropout
    model.add(Dropout(0.4, input_shape=(num_inputs,)))
    model.add(Dense(nodes[0], activation='relu'))

    if len(nodes) > 1:
        # Sequentially add hidden layers with dropout
        for i in nodes[1:]:
            model.add(Dropout(0.2))
            model.add(Dense(i, activation='relu'))

    model.add(Dropout(0.2))
    # Output layer
    model.add(Dense(1, activation = 'sigmoid'))

    METRICS = [tf.keras.metrics.AUC(),tf.keras.metrics.AUC(curve='PR',name='PR'),'accuracy', tf.keras.metrics.Precision(), tf.keras.metrics.Recall()]

    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=METRICS)

    return model


def make_nodes(layers, node): # Generate list of nodes to parse to build model

    nodes = list()
    for layer in range(layers):
        nodes.append(node)

    return nodes


def test_arch(layers, node, path, des): # Iteratively train different architectures
    # Load data
    x_train, y_train, x_val, y_val = load_data_with_validation_sets(path)

    for l in range(layers): # Loop through layer range

        if l > 0 :
            l=l+1 # Minimum number of layers is 2
            for n in node: # Loop though user input node options

                print(f'Testing precision recall of {l} layer(s) with {n} nodes...')
                print(f'With early stopping patience at 30...')

                # Create model
                es = EarlyStopping(monitor = 'val_PR', mode='max', verbose=1, patience=30)
                nodes = make_nodes(l, n)
                model = build_model(nodes, x_train.shape[1])
                history = model.fit(x_train, y_train, epochs=1000, batch_size=16, validation_data=(x_val, y_val), callbacks = [es], verbose = 1)

                # Save training metrics to a csv
                pr = max(history.history['val_PR'])
                print(f'Max precision recall for {l} layer(s) with {n} nodes: {pr}')
                h = pd.DataFrame(history.history)
                h.to_csv(f'{des}{pr}_max_{l}_layers_{n}_nodes_results.csv')
                clear_session()


def parse_args(args): # Parse user inputs

    h5_path = args[args.index('-h5') + 1]

    layers = int(args[args.index('-layers') + 1])

    nodes = args[args.index('-nodes') + 1].split(',')

    des = args[args.index('-des') + 1]

    return h5_path, layers, nodes, des


def main():

    h5_path, layers, nodes, des = parse_args(sys.argv)

    nodes = list(map(int,nodes))

    test_arch(layers, nodes, h5_path, des)

if __name__ == '__main__':

    main()
