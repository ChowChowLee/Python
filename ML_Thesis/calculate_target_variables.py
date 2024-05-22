import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gmean
# from Bayesian import *

output_file_path = "C:\\Users\\Lee\\PycharmProjects\\pythonProjectAlpacaBaby\\Python\\ML_Thesis\\output_data\\"
input_file_path = "C:\\Users\\Lee\\PycharmProjects\\pythonProjectAlpacaBaby\\Python\\ML_Thesis\\input_data\\"
df_raw = pd.read_csv(input_file_path + "sludgeDatasetMerge.tsv", sep='\t')  # sludgeDatasetMergeCalculated.tsv  #sludgeDatasetMergeCalculatedAvglogHL.tsv ### sludge_raw_use_this_for_baymean_test.tsv
CURATE_DATA_POINTS = False

def main() -> None:
    print(df_raw.describe())



#
# def main():
#     df = calculate_target_variable(df_raw)
#     df1 = describe_dropna_halflife(df)
#     # df1.to_csv(output_file_path+"sludge_biomass_bay3.tsv", sep='\t', index=False)
#     df_ = calculate_bay_mean_std(df1)
#     df_.to_csv(output_file_path+'sludge_bay_PriorMuStd_bay3_biomass_BI.tsv', sep='\t')

    # plot_distribution(output_path=output_file_path + 'Distribution_comparison_PriorMuStd_2.pdf')

def calculate_target_variable(df):
    df['hl_log_gmean'] = np.log10(get_geometric_mean(df, 'halflife'))
    df['hl_log_median'] = np.log10(get_median(df, 'halflife'))
    df['hl_log_std'] = np.log10(get_std(df, 'halflife'))    # you might need to change the column name!!
    df['hl_log_spread'] = get_hl_spread(df)
    df['biomass_hl_log_gmean'] = np.log10(get_geometric_mean(df, 'hl_biomass_corrected'))
    df['biomass_hl_log_median'] = np.log10(get_median(df, 'hl_biomass_corrected'))
    df['biomass_hl_log_std'] = np.log10(get_std(df, 'hl_biomass_corrected'))   #
    df['biomass_hl_log_spread'] = get_biomass_hl_spread(df)
    df['acidity_std'] = get_std(df, 'acidity')
    df['temperature_std'] = get_std(df, 'temperature')
    df['biomass_log_std'] = get_std(df, 'total_suspended_solids_concentration_start')
    return df

def describe_dropna_halflife(df):
    df.dropna(subset=['halflife', 'halflife_log'], inplace=True)  # Remove the NaN in halflife column, otherwise you get ValueError in bmean, bstd calculation.
    # df.to_csv(output_file_path + 'sludge_calculated_test_for_baycalculation.tsv', sep='\t')   #'sludge_calculated.tsv'
    df_ = df.copy()
    description = df_.describe()
    print("Summary of loaded data:\n------------------\n", description)
    description.to_csv(output_file_path + 'sludge_calculated_test_for_baycalculation_describe.tsv', sep='\t')
    return df

def calculate_bay_mean_std(df):
    bmean, bstd = get_bayesian_stats(df)
    df['biomass_hl_log_bayesian_mean'] = bmean
    df['biomass_hl_log_bayesian_std'] = bstd
    # df.to_csv(output_file_path+'sludge_calculated_test_for_baycalculation.tsv', sep='\t')
    return df


def get_std(df, column):
    new = []
    for index, row in df.iterrows():
        this = df.loc[df['reduced_smiles'] == row['reduced_smiles']]
        std = np.nanstd(this[column])
        if std == 0:
            new.append(np.NaN)
        else:
            new.append(std)
    return new

def get_mean(df, column):
    new = []
    for index, row in df.iterrows():
        this = df.loc[df['reduced_smiles'] == row['reduced_smiles']]
        std = np.nanmean(this[column])
        new.append(std)
    return new

def get_median(df, column):
    new = []
    for index, row in df.iterrows():
        this = df.loc[df['reduced_smiles'] == row['reduced_smiles']]
        std = np.nanmedian(this[column])
        new.append(std)
    return new

def get_hl_spread(df):
    new = []
    for index, row in df.iterrows():
        this = df.loc[df['reduced_smiles'] == row['reduced_smiles']]
        spread = max(this['halflife_log']) - min(this['halflife_log'])
        new.append(spread)
    return new

def get_biomass_hl_spread(df):
    new = []
    for index, row in df.iterrows():
        this = df.loc[df['reduced_smiles'] == row['reduced_smiles']]
        spread = max(np.log10(this['hl_biomass_corrected'])) - min(np.log10(this['hl_biomass_corrected']))
        new.append(spread)
    return new


def legal(value, name):
    if np.isnan(value):
        print(f"Problem: no {name}")
        return False
    elif value == 0:
        print(f"Problem: {name} is 0")
        return False
    return True

def g_mean(x):
    a = np.log(x)
    return np.exp(a.mean())

def get_geometric_mean(df, column):
    new = []
    for index, row in df.iterrows():
        this = df.loc[df['reduced_smiles'] == row['reduced_smiles']]
        gmean = g_mean(this[column])
        new.append(gmean)
    return new

# def get_gmean(df, column):
#     new = []
#     for index, row in df.iterrows():
#         this = df.loc[df['reduced_smiles'] == row['reduced_smiles']
#         v = this[column]
#         if not legal(v, f"{column}"):
#             return np.NaN
#
#         if legal(v, f"{column}"):
#             std = gmean(v)
#             new.append(std)
#     return new

def process_comment_list(comment_list):
    new_list = []
    for comment in comment_list:
        if type(comment) == float:
            new_list.append('')
        elif '<' in comment:
            new_list.append('<')
        elif '>' in comment:
            new_list.append('>')
        else:
            new_list.append('')
    return new_list


# If you want to calculate the hl_log_bayesian_mean with biomass corrected values,
# you need to modify the line 156

def get_bayesian_stats(df):
    mean_list = []
    std_list = []
    results = {} # {'index': (mean, std)}
    for index, row in df.iterrows():
        if row['reduced_smiles'] in results.keys():
            mean, std = results[row['reduced_smiles']]
        else:
            this = df.loc[df['reduced_smiles'] == row['reduced_smiles']]
            comment_list_raw = process_comment_list(this["halflife_comment"])
            # y_raw = np.array(this['halflife_log'])   # this is for calculation of hl_log_bayesian_mean
            y_raw = np.array(this['log_hl_biomass_corrected'])     # # this is for calculation of hl_log_bayesian_mean biomass corrected
            if CURATE_DATA_POINTS == True:
                pass
            else:
                y = y_raw
                comment_list = comment_list_raw
            print("\nCOMPOUND reduced_smiles {}".format(row['reduced_smiles']))
            print("Compute bayes for {} with comments {}".format(y, comment_list))
            bayesian = Bayesian(y=y, comment_list=comment_list)
            bayesian.set_prior_mu(mean=-0.08, std=2)     #(Original: mean=1.5, std=2) Set prior_mu_std as 2
            bayesian.set_prior_sigma(mean=0.23, std=1.0) #(Original: mean=0.2, std=0.5)
            mean, std, posteriorMuStd = bayesian.get_posterior_distribution()   # I add the third variable "posteriorMuStd since the unpacked error expected 2"
            results[row['reduced_smiles']] = (mean, std)
            print('mean: {}, std: {}'.format(mean, std))
            bayesian.set_path_to_output_folder(output_file_path+'figures\\')
            # bayesian.plot_distribution(output_path=output_file_path + 'figures\\Distribution_comparison_PriorMuStd_2_{}.pdf'.format(row['compound_name']))  ## I add this line to plot the data distribution.
            bayesian.plot_distribution('Distribution_comparison_PriorMuStd_bay3_biomass_{}.pdf'.format(row['compound_name']))
        mean_list.append(round(mean, 2))
        std_list.append(round(std, 2))

    return mean_list, std_list

def print_y_stats(y):
    print('Original mean:', np.mean(y))
    print('Original median:', np.median(y))
    print('Original std:', np.std(y))

def run_bayes_simple(y, comment=[]):
    print_y_stats(y)
    bayesian = Bayesian(y=y, comment_list=comment)
    bayesian.set_save_backend(True)
    bayesian.set_iterations(2000)
    bayesian.set_prior_mu(mean=1.5, std=2)
    bayesian.set_prior_sigma(mean=0.2, std=0.5)
    bayesian.plot_emcee_chain()
    bayesian.plot_distribution()
    bayesian.log_prob_distribution_plot()
    bayesian.corner_plot()

    new_mean, new_std = bayesian.get_posterior_distribution()
    print('\nPrior mean: {}'.format(bayesian.get_prior_mu()))
    print('Prior std: {}'.format(bayesian.get_prior_sigma()))
    print('\nBayesian inference mean:', new_mean)
    print('Bayesian inference std:', new_std)

    # column_groupby_smiles_sludge(df_raw, 'log_hl_biomass_corrected')
# def rename_column(df):
#     df_rename = df.rename(columns={
#                        f'mean_log_hl_combined': 'hl_log_mean',
#                        f'median_log_hl_combined': 'hl_log_median',
#                        f'geomean_log_hl_combined': 'hl_log_gmean',
#                        f'std_log_hl_combined': 'hl_log_std',
#                        f'mean_log_hl_biomass_corrected': 'biomass_hl_log_mean',
#                        f'median_log_hl_biomass_corrected': 'biomass_hl_log_median',
#                        f'geomean_log_hl_biomass_corrected': 'biomass_hl_log_gmean',
#                        f'std_log_hl_biomass_corrected': 'biomass_hl_log_std'
#                        })
#     return df_rename

def column_groupby_smiles_sludge(df, column):
    mean_column_list = []
    median_column_list = []
    geomean_column_list = []
    std_column_list = []
    for index, row in df.iterrows():
    #     if not legal(row[column], f"row[{column}]"):
    #         return np.NaN
    #
    #     elif legal(row[column], f"row[{column}]"):

        x = df.loc[df['reduced_smiles'] == row['reduced_smiles']]
        mean_ = np.nanmean(x[column])
        median_ = np.nanmedian(x[column])
        g_ = gmean(x['hl_biomass_corrected'])
        geomean_ = np.log10(g_)
        std_ = np.nanstd(x[column])
        mean_column_list.append(mean_)
        median_column_list.append(median_)
        geomean_column_list.append(geomean_)
        std_column_list.append(std_)
    df[f'mean_{column}'] = mean_column_list
    df[f'median_{column}'] = median_column_list
    df[f'geomean_{column}'] = geomean_column_list
    df[f'std_{column}'] = std_column_list
    df.rename(columns={
                        'mean_log_hl_combined': 'hl_log_mean',
                        'median_log_hl_combined': 'hl_log_median',
                        'geomean_log_hl_combined': 'hl_log_gmean',
                        'std_log_hl_combined': 'hl_log_std',
                        'mean_log_hl_biomass_corrected': 'biomass_hl_log_mean',
                        'median_log_hl_biomass_corrected': 'biomass_hl_log_median',
                        'geomean_log_hl_biomass_corrected': 'biomass_hl_log_gmean',
                        'std_log_hl_biomass_corrected': 'biomass_hl_log_std'
                        }, inplace=True)
    df_new = df.to_csv(output_file_path + "sludgeDatasetMergeCalculatedAvglogHLBiomass_Munich_09032023.tsv", sep='\t')

    return df_new

if __name__ == '__main__':
    main()

"""    
The following code is my initial script for fetching the avg values of logDT50 and logDT50' 

def column_groupby_smiles_sludge(df, column):
    mean_column_list = []
    median_column_list = []
    geomean_column_list = []
    std_column_list = []
    for index, row in df.iterrows():
    #     if not legal(row[column], f"row[{column}]"):
    #         return np.NaN
    #
    #     elif legal(row[column], f"row[{column}]"):

        x = df.loc[df['smiles'] == row['smiles']]
        mean_ = np.nanmean(x[column])
        median_ = np.nanmedian(x[column])
        g_ = gmean(x['hl_biomass_corrected'])
        geomean_ = np.log10(g_)
        std_ = np.nanstd(x[column])
        mean_column_list.append(mean_)
        median_column_list.append(median_)
        geomean_column_list.append(geomean_)
        std_column_list.append(std_)
    df[f'mean_{column}'] = mean_column_list
    df[f'median_{column}'] = median_column_list
    df[f'geomean_{column}'] = geomean_column_list
    df[f'std_{column}'] = std_column_list
    df.rename(columns={
                        'mean_log_hl_combined': 'hl_log_mean',
                        'median_log_hl_combined': 'hl_log_median',
                        'geomean_log_hl_combined': 'hl_log_gmean',
                        'std_log_hl_combined': 'hl_log_std',
                        'mean_log_hl_biomass_corrected': 'biomass_hl_log_mean',
                        'median_log_hl_biomass_corrected': 'biomass_hl_log_median',
                        'geomean_log_hl_biomass_corrected': 'biomass_hl_log_gmean',
                        'std_log_hl_biomass_corrected': 'biomass_hl_log_std'
                        }, inplace=True)
    df_new = df.to_csv(output_file_path + "sludgeDatasetMergeCalculatedAvglogHLBiomass.tsv", sep='\t')

    return df_new
"""