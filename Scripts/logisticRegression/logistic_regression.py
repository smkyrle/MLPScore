'''
Script for fitting a regression model to validation set predictions and testing on test data
Regression model and results are saved to an output directory
'''
import os
# Configure GPU memory for Tensorflow
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
import keras
from keras.models import load_model
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, f1_score, auc, precision_score, recall_score
from sklearn.linear_model import LogisticRegression
import warnings
import random
import sys
from pprint import pprint
from tqdm import tqdm
import pickle
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

    return x_train, y_train, x_val, y_val

def load_data(h5_path): # fetch x_test and y_test sets from hdf5 file

    # load dataframe
    data = pd.read_hdf(h5_path, key='df', mode='r')

    # split label, features and details columns
    y_test = data[['Label']].fillna(0).astype(int).copy()
    x_test = data.drop(['Label','PDBCode','Binding_Data'], axis=1).copy()

    binding_data = data[['Binding_Data','PDBCode']].copy()
    codes = data['PDBCode'].to_list()

    return x_test, y_test, binding_data, codes

def list_models(csv_path, model_path):

    csvs = os.listdir(csv_path)
    csvs = sorted(csvs, reverse = True)
    csvs = csvs[1:]

    models = list()
    model_nums = list()

    for csv in csvs:
        csv = csv.split('_')
        models.append(f'{model_path}model_{csv[-2]}.h5')
        model_nums.append(f'model_{csv[-2]}.h5')

    return models, model_nums

def model_predictions(models, model_nums, X):

    predictions = pd.DataFrame()
    for model, num in zip(models,model_nums):
        model = load_model(model)
        X_predictions = model.predict(X)
        predictions[num] = list(X_predictions.flat)

    return predictions

def fit_regression(X, y):
    reg = LogisticRegression().fit(X, y)
    return reg

def regression_metrics(reg, X, y):
    probs = reg.predict_proba(X)[:,1]
    yhat = reg.predict(X)
    p_score = precision_score(y,probs.round())
    r_score = recall_score(y,probs.round())
    precision, recall, thresholds = precision_recall_curve(y, probs)
    f1, auc_pr = f1_score(y, yhat), auc(recall, precision)
    auc_roc = roc_auc_score(y, probs)
    fpr, tpr, _ = roc_curve(y, probs)

    return precision, recall, f1, auc_pr, fpr, tpr, auc_roc, probs, p_score, r_score

def parse_args(args): # Parse user inputs

    des = args[args.index('-des') + 1] # Path to save test predictions and regression model

    csv_path = args[args.index('-csv') + 1] # Path to training csv files

    model_path = args[args.index('-models') + 1]) # Path to saved models

    train = args[args.index('-train') + 1] # Path to training set

    test = args[args.index('-test') + 1] # Path to testing set

    num_models = args[args.index('-num_models') + 1] # Number of models to use for logistic regression

    return des, csv_path, model_path, train, test, num_models

if __name__ == '__main__':

    des, csv_path, model_path, train, test, num_models = parse_args(sys.argv)

    x_train, y_train, x_val, y_val = load_data_with_validation_sets(train_path)
    models, model_nums = list_models(csv_path, model_path)

    models = models[:num_models]
    model_nums = model_nums[:num_models]
    predictions = model_predictions(models, model_nums, x_val)
    regression = fit_regression(predictions, y_val)

    X, y, bd, codes = load_data(path)

    X = model_predictions(models, model_nums, X)
    results = regression_metrics(regression, X, y)

    pickle.dump(results, open(f'{des}results.pkl'))
    pickle.dump(regression, open(f'{des}regression_model.pkl'))
