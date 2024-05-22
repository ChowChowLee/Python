import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict

file_location = "C:\\Users\\leetseng\\TWtest"
# input_file_path = file_location+'\\input\\sludgeDatasetMergeCalculated.tsv'
input_file_path_1 = file_location+'\\input\\sludge_Original_raw.tsv'    #'input/sludgeWithSmiles.tsv'    "Your concatenate the different iterations, so your index is not continuous!!!!"
input_file_path_2 = file_location+'\\input\\sludge_Leo_raw.tsv'
input_file_path_3 = file_location+'\\input\\sludge_Rich_raw.tsv'
output_file_path_full = file_location+'\\output\\sludge_raw_bay3_check.tsv'  #sludgeDatasetMergeCalculated.tsv
# output_file_path_full = file_location+'\\output\\sludgeDatasetMergeCalculated.tsv'

data1 = pd.read_csv(input_file_path_1, sep='\t')
data2 = pd.read_csv(input_file_path_2, sep='\t')
data3 = pd.read_csv(input_file_path_3, sep='\t')
# data1 = data1.reset_index(drop=True)
# data2 = data2.reset_index(drop=True)
# data3 = data3.reset_index(drop=True)
data_merge = pd.concat([data1, data2, data3])
data_merge.index = np.arange(1, len(data_merge) + 1)
data_merge = data_merge.drop('index', axis=1)  #old column
data_merge.index.name = 'index'
# data_merge = data_merge.set_index('index')
# print(data_merge.head(2))
# print(data_merge.columns.tolist())


# plt.figure(figsize=(35, 5))
# sns.barplot(data=data0, x='compound_name', y=cpd, palette="Greens_d")
# plt.xticks(fontsize=8, rotation=70)
# plt.savefig(file_location+'\\output\\figures\\std_sludge_Merge_test.pdf')
# plt.close()


def init_k_data(data_merge: pd.DataFrame) -> list:
    """
    Collect all the rate constant(k) and biomass corrected rate constant(k_biomass) first. Then the user can use
    the collected k and k_biomass for the subsequent transformation between rate constant and half-live.
    :param data:
    :return:
    """
    k_list = []
    k_biomass_list = []
    for index, row in data_merge.iterrows():
        print(row['scenario_id'] + '\n')
        k = get_k(row)
        k_biomass = get_k_biomass(row, k)
        k_list.append(k)
        k_biomass_list.append(k_biomass)
    raise NotImplementedError

def init_k_list(df_merge: pd.DataFrame) -> list:
    """
    Collect all the rate constant(k) first. Then the user can use
    the collected k for the subsequent transformation between k and half-live.
    :param data:
    :return:
    """
    k_list = []
    for index, row in df_merge.iterrows():
        k = get_k(row)
        k_list.append(k)
    return k_list


def init_k_biomass_list(df_merge: pd.DataFrame) -> list:
    """
    Collect all the biomass corrected rate constant(k_biomass) first. Then the user can use
    the k_biomass for the subsequent transformation between k_biomass and biomass corrected half-live.
    :param data:
    :return:
    """
    k_biomass_list = []
    for index, row in df_merge.iterrows():
        k_biomass = get_k_biomass(row, k)
        k_biomass_list.append(k_biomass)
    return k_biomass_list




def init_hl_data(data: pd.DataFrame) -> None:
    hl_list = []
    hl_biomass_list = []
    hl_biomass_list_2 = []
    for index, row in data_merge.iterrows():
        print(row['scenario_id'] + '\n')
        DT50 = get_DT50(row, k)  # k_combined does not yet exist in table
        DT50_biomass = get_DT50_biomass(row, DT50, k_biomass)
        hl_list.append(DT50)
        hl_biomass_list.append(DT50_biomass)
        DT50_biomass_2 = get_DT50_biomass_double_check(row, DT50, k_biomass)
        hl_biomass_list_2.append(DT50_biomass_2)
        print(DT50, DT50_biomass, DT50_biomass_2)

def main():

    data_merge['k_combined'] = init_k_list(data_merge)  # k_combined = k given + k calculated from halflife
    data_merge['k_biomass_corrected'] = init_k_biomass_list(data_merge)

    data_merge['halflife'] = hl_list
    data_merge['hl_biomass_corrected'] = hl_biomass_list
    data_merge['hl_biomass_corrected_2'] = hl_biomass_list_2

    data_merge['log_k_combined'] = np.log10(data_merge['k_combined'])
    data_merge['log_k_biomass_corrected'] = np.log10(data_merge['k_biomass_corrected'])
    data_merge['halflife_log'] = np.log10(data_merge['halflife'])
    data_merge['log_hl_biomass_corrected'] = np.log10(data_merge['hl_biomass_corrected'])
    data_merge.to_csv(output_file_path_full, mode='w', sep="\t")



#before the calculation of hl, check if you have the rateconstant

def legal_value(value: float) -> bool:
    raise NotImplementedError



def get_k(row):
    k_given = row['rateconstant']
    k_unit = row['rateconstant_unit']
    k_true = np.NaN
    TSS = row['total_suspended_solids_concentration_start']
    hl = row['halflife_raw']
    order = row['halflife_model']
    if not np.isnan(k_given):
        if k_given != 0 and k_unit == '1 / day' and not np.isnan(TSS):
            k_true = k_given
        elif k_given != 0 and k_unit == 'L / (g TSS * day)' and not np.isnan(TSS):
            k_true = k_given * TSS
        elif k_given != 0 and k_unit == '㎍ / (g TSS * day)' and not np.isnan(TSS):   ################
            pass
        else:
            if np.isnan(TSS):
                print('Problem: no TSS')
            elif k_given == 0:
                print('Problem: given rate constant is 0')
    elif not np.isnan(hl):    #elif not np.isnan(hl):
        if order == 'Zero order':
            k_true = TSS/(2 * hl)
        elif order == 'First order':
            k_true = np.log(2)/hl
        elif order == 'Pseudo first order': # it's a biomass corrected hl
            real_hl = hl / TSS
            k_true = np.log(2)/real_hl
        else:   #By default, using the 1st order reaction formula
            k_true = np.log(2)/hl
    return k_true

def get_k_biomass(row, k):
    k_given = row['rateconstant']
    k_unit = row['rateconstant_unit']
    TSS = row['total_suspended_solids_concentration_start']
    hl = row['halflife_raw']
    k_biomass = np.NaN
    if not np.isnan(k_given) and k_given != 0:
        if k_unit == '1 / day':
            k_biomass = k / TSS           ################# should be k_given / TSS
        elif k_unit == 'L / (g TSS * day)':
            k_biomass = k_given
        else:
            if k_given == 0:
                print('Error: rate constant is 0')
    elif np.isnan(k_given) and not np.isnan(hl):     #add this conditional expressions for the Rich's dataset
        k_biomass = k / TSS
    return k_biomass

def get_DT50(row, k):
    hl_given = row['halflife_raw']
    hl = np.NaN
    if np.isnan(hl_given):
        if not np.isnan(k): # removed k_combined, does not yet exist at this point
            hl = np.log(2)/k
        else:
            hl = np.NaN
    elif hl_given != 0:   #check yourself
        hl = hl_given
    else:
        print('Error: half-life == 0')
    return hl

def get_DT50_biomass(row, hl, k_biomass):   #can generate in two ways   1. take hl list/TSS  2. ln2 / k_biomass need to be in consistent. just safety check.
    TSS = row['total_suspended_solids_concentration_start']
    hl_biomass = np.NaN
    if not np.isnan(hl):
        hl_biomass = hl/TSS
    elif not np.isnan(k_biomass):
        hl_biomass = np.log(2)/k_biomass
    return hl_biomass

# We use an alternative way to see if the outcome of DT50 biomass is consistent in different ways.
def get_DT50_biomass_double_check(row, hl_list, k_biomass):
    TSS = row['total_suspended_solids_concentration_start'] = row['halflife']
    hl_biomass_2 = np.NaN
    if not np.isnan(hl):
        hl_biomass = hl_list / TSS
    elif not np.isnan(k_biomass):
        hl_biomass = np.log(2) / k_biomass
    return hl_biomass_2




if __name__ == '__main__':
    main()

#create the set of SMILES
# list_of_canonicalize_smiles = data_merge['canonicalize_smiles'].values.tolist()
# set_of_canonicalize_smiles = set(list_of_canonicalize_smiles)
# print(set_of_canonicalize_smiles)

