'''
This script strips features from the dataset that are below
a user defined variance and above a user defined pairwise
pearson correlation, then scales the data using the sklearn
MaxAbsScaler
'''

import pandas as pd
import numpy as np
from scipy import stats
import sys
from sklearn.preprocessing import MaxAbsScaler
import joblib

pd.set_option('display.max_rows',None)

def find_low_variance_features(feature_df, var_thresh): # drop columns from dataframe with variance equal to or below threshold

    # calculate variance of features columnwise
    var_df = feature_df.var()

    # make column from column name index
    var_df = var_df.reset_index()

    # get list of columns with variance equal to or below threshold
    no_variance_df = var_df.loc[var_df[0] <= var_thresh]
    to_drop = list(no_variance_df['index'])

    return to_drop

def remove_co_correlated_features(varied_features_df, corr_thresh, y): # drop co-correlated features based on pairwise pearson and tiebreak with pearson to target variable

    # calculate pearson correlations for all features to label
    print('Populating target pearson dictionary...')
    lc_df = varied_features_df.copy()
    lc_df['Label'] = y
    lc_corr = lc_df.corr(method='pearson').abs()
    lc_corr_df = lc_corr.reset_index()

    # store pearson correlations in dictionary
    corr_dict = dict(zip(list(lc_corr_df['index']), list(lc_corr_df['Label'])))

    # calculate pairwise correlations for all features in dataset
    print('Calculating pairwise pearson correlations...')
    corr_mtx = varied_features_df.corr(method="pearson").abs()

    # remove duplicate correlations by taking upper triangle of correlation matrix
    upper = corr_mtx.where(np.triu(np.ones(corr_mtx.shape), k=1).astype(np.bool))

    # find pairs of co-correlated features
    print('Identifying co-correlated features over threshold...')
    pairwise_dict = dict()
    for col in upper.columns:
        sub_df = upper.loc[upper[col]>corr_thresh]
        pairs = list(sub_df.index)
        pairwise_dict[col] = pairs

    keep = list()
    drop = list()
    for key, value in pairwise_dict.items():
        if len(value) > 1:
            print(f'Key = {key}, Value = {value}')
        if len(value) == 1:
            k_corr = corr_dict.get(key)
            v_corr = corr_dict.get(value[0])
            if k_corr > v_corr:
                keep.append(key)
                drop.append(value[0])
            else:
                keep.append(value[0])
                drop.append(key)

    final_drops = [f for f in drop if f not in keep]

    return final_drops

def max_abs_scale_data(clean_df, output_dir): # scale the data with MaxAbsScaler maintaining sparsity

    # fit scaler to data
    scaler = MaxAbsScaler().fit(clean_df)

    # save fitted scaler for use on test data
    joblib.dump(scaler,f'{output_dir}maxabs_scaler_params.save')

    # scale the dataframe
    clean_df[clean_df.columns] = scaler.transform(clean_df)

    return clean_df

def parse_args(args): # parse CLI user inputs

    h5_path = args[args.index('-h5') + 1]

    var_thresh = float(args[args.index('-var_thresh') + 1])

    corr_thresh = float(args[args.index('-corr_thresh') + 1])

    output_dir = args[args.index('-out') + 1]

    scale = [True if '-scale' in args else False][0]

    return h5_path, var_thresh, corr_thresh, output_dir, scale

def main(): # run script using CLI

    h5_path, var_thresh, corr_thresh, output_dir, scale = parse_args(sys.argv)

    # load dataframe and drop label and non-feature columns
    df = pd.read_hdf(h5_path)
    print(f'Full Set: {len(list(df))}')

    feature_df = df.drop(['PDBCode','BindingDataType', 'BindingValue','BindingUnits','Label', 'Database', 'Binding_Data_in_uM'], axis=1)
    print(f'Feature Set: {len(list(feature_df))}')

    # define label column as separate variable
    y = df['Label'].copy()

    # find features with variance equal to or below threshold
    low_variance_to_drop = find_low_variance_features(feature_df, var_thresh)

    # drop low variance columns and return new dataframe
    varied_features_df = df.drop(low_variance_to_drop, axis=1)
    print(f'High Variance Feature Set: {len(list(varied_features_df))}')

    # re-remove label and non-feature columns
    varied_features_df = varied_features_df.drop(['PDBCode','BindingDataType', 'BindingValue','BindingUnits','Label', 'Database', 'Binding_Data_in_uM'], axis=1)

    # find co-correlated features to drop
    co_correlated_to_drop = remove_co_correlated_features(varied_features_df, corr_thresh, y)

    # drop low variance columns from original dataframe
    clean_df = df.drop(low_variance_to_drop, axis=1)
    print(f'Clean High Variance Feature Set: {len(list(clean_df))}')

    # drop co-correlated columns from original dataframe
    clean_df = clean_df.drop(co_correlated_to_drop, axis=1)
    print(f'Clean High Variance Low Correlation Feature Set: {len(list(clean_df))}')

    if scale: # scale the data then save

        # fetch features only again to pass to scaler
        clean_df_for_scaler = clean_df.copy()
        clean_df_for_scaler = clean_df_for_scaler.drop(['PDBCode','BindingDataType', 'BindingValue','BindingUnits','Label', 'Database', 'Binding_Data_in_uM'], axis=1)
        print(f'Data For Scaling Function: {len(list(clean_df_for_scaler))}')

        # scale the data
        print('Scaling data...')
        scaled_df = max_abs_scale_data(clean_df_for_scaler, output_dir)
        print(f'Scaling Function Output Data: {len(list(scaled_df))}')

        # return info columns to dataframe
        scaled_df[['PDBCode','BindingDataType', 'BindingValue','BindingUnits','Label', 'Database', 'Binding_Data_in_uM']] = clean_df[['PDBCode','BindingDataType', 'BindingValue','BindingUnits','Label', 'Database', 'Binding_Data_in_uM']]
        print(f'Scaling Function Output Data With Info Columns: {len(list(scaled_df))}')

        # save new dataframe of remaining features
        scaled_df.to_hdf(f'{output_dir}scaled_clean_df.h5',key="df",mode="w")

    else: # save the data

        # save new dataframe of remaining features
        clean_df.to_hdf(f'{output_dir}clean_df.h5',key="df",mode="w")


if __name__ == '__main__':
    main()
