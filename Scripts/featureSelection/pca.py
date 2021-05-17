'''
Script for performing PCA on the training set
'''

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import joblib
import sys
pd.set_option('display.max_rows',None)

def pca(data, output_dir, n):

    print('Fitting data...')
    pca = PCA(n_components = n).fit(data)
    print('Saving PCA parameters...')
    joblib.dump(pca, f'{output_dir}pca_params.save')
    variance = pca.explained_variance_
    components = pca.components_
    print('Transforming data...')
    transformed = pca.transform(data)
    components_df = pd.DataFrame(pca.components_,columns=data.columns)

    headers = [f'PC-{y+1}' for y in range(transformed.shape[1])]

    transformed_df = pd.DataFrame(transformed,columns=headers)

    return transformed_df, components_df, var_perc, cumulative_var_perc, headers

def parse_args(args):

    training_data = args[args.index('-loc') + 1]

    output_dir = args[args.index('-des') + 1]

    rfecv_csv = args[args.index('-csv') + 1]

    n = int(args[args.index('-n') + 1])

    n = n/100

    return training_data, output_dir, rfecv_csv, n # Parse User arguments

def main():

    training_data, output_dir, rfecv_csv, n = parse_args(sys.argv)
    labels = ['PDBCode','BindingDataType', 'BindingValue','BindingUnits','Label', 'Database', 'Binding_Data_in_uM']
    rfecv_df = pd.read_csv(rfecv_csv)
    rfecv_df = rfecv_df["Feature"].loc[rfecv_df["Choices"] == True]
    rfecv_list = list(rfecv_df)


    data = pd.read_hdf(training_data, key="df", mode='r')
    k = data[labels].copy()
    keys = data[labels].copy().reset_index()
    df = data[rfecv_list].copy()

    transformed_df, feature_contribution_df, var_perc, cumulative_var_perc, headers = pca(df, output_dir, n)

    selected_df[labels] = keys[labels]
    transformed_df[labels] = keys[labels]

    selected_df.to_hdf(f'{output_dir}selected_components.h5', key="df", mode="w")
    transformed_df.to_hdf(f'{output_dir}all_components.h5', key="df", mode="w")
    feature_contribution_df.to_csv(f'{output_dir}feature_contribution.csv')



if __name__ == '__main__':
    main()
