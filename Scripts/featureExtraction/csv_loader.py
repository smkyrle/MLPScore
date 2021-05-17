'''
Script from lab book for combining dataframes into one hdf file
'''


import pandas as pd
pd.set_option('display.max_rows', None)

binana = pd.read_csv("/home/sammk/Research_Project/Data/Final/all_binana_data.csv")

ecif = pd.read_csv("/home/sammk/Research_Project/Data/Final/ecif_data.csv")

binding_data = pd.read_csv("/home/sammk/Research_Project/Data/Final/convertedBindingData.csv")

binana['PDBCode'] = binana['Name'].apply(lambda x: x.replace('_ligand',''))


final_df = ecif.merge(binana,how="left",on=["PDBCode"])
final_df = final_df.merge(binding_data,how="left",on=["PDBCode"])
print(len(list(final_df)))
final_df = final_df.drop(["Name"],axis=1)
final_df["Label"].fillna(0,inplace=True)
final_df["BindingValue"].fillna(1e+6,inplace=True)
final_df["Binding_Data_in_uM"].fillna(1e+6,inplace=True)
final_df["BindingDataType"].fillna("Decoy",inplace=True)
final_df["BindingUnits"].fillna("Decoy",inplace=True)
final_df.sort_values(by="PDBCode",inplace=True)
final_df["Database"].ffill(inplace=True)
print(len(list(final_df)))
print(len(final_df))

#final_df.to_hdf("/home/sammk/Research_Project/Data/Final/final_df.h5",key="df",mode="w")
