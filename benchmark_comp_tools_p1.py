# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 08:36:22 2024

@author: ProtoQSAR,S.L.
"""

############################### CLEAN EVERYTHING #############################
for i in list(globals().keys()):
    if not i.startswith('_'):
        exec('del ' + i)
##############################################################################
#################################### IMPORTS #################################
import os
from pathlib import Path
import pandas as pd
from scipy.stats import zscore
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import itertools
import re
##############################################################################
################################ INITIAL VARIABLES ###########################

parent = Path(__file__).resolve().parent
os.chdir(parent)
print(f'Working on: \n {os.getcwd()}')
data_folder =  '.' + os.path.sep + 'input'

results_folder = '.' + os.path.sep + 'results'

if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    print("Output directory created: ../results/")
else:
    print("Output directory already exists: ../results/")

##############################################################################


###############################  FUNCTIONS ###################################
def order_input_df(df,name, desired_units):

    standard_columns = ['ID', 'orig_ID', 'SMILES', 'y', 'CAS', 'NAME', 'UNITS', 'smiles_source']
    additional_columns = set(list(df.columns)).difference(set(standard_columns))


    df.rename(columns = {'ID':f'ID_{name}', 'orig_ID':f'orig_ID_{name}', 'y':f'yorig_{name}', 'CAS': f'CAS_{name}', 'NAME' : f'NAME_{name}',
                         'UNITS': f'units_{name}', 'smiles_source': f'smilessource_{name}'}, inplace = True)

    if len(additional_columns) !=0:

        for add_col in additional_columns:
            df.rename(columns = {add_col: f'{add_col}_{name}'}, inplace = True)

    df.set_index('SMILES', inplace = True)
    df[f'SMILES_{name}'] = df.index

    units = list(set(df[f'units_{name}']))[0]

    df[f'desired_units_{name}'] = desired_units
    df[f'y_{name}'] = df.apply(lambda row: conversion(row[f'yorig_{name}'], units, desired_units), axis=1)


    return df, additional_columns



def conversion(row, units, desired_units):



    if units == 'K' and desired_units == 'C':
        return row - 273.15

    elif units == '%PBR' and desired_units == 'ratioFU' :
        return 1- row/100

    elif units == 'LogkPa' and desired_units == 'LogmmHg':

        return np.log10((10**row)*0.00750082)

    elif units == '%' and desired_units == 'classHIA':

        if row >= 30: return 1
        else: return 0

    elif units == '%' and desired_units == 'classF30':

        if row >= 30: return 1
        else: return 0

    elif units == '[10-6]cm/s' and desired_units == 'logcm/s':

        return np.log10(row*10**(-6))


    else:
        return row


def check_duplicated_molecules(df,name):

    unique_smiles = []
    nonunique_smiles = []
    for smi in df_ordered[f'SMILES_{name}']:
        if smi not in unique_smiles:
            unique_smiles.append(smi)
        else:
            nonunique_smiles.append(smi)
    if len(nonunique_smiles) != 0:
        print(f'\t\tWARNING!!!!!!!!! you have duplicated molecules in {name}')



def representation_zscore(info, threshold_zscore):
    number_datasets = len(info.keys())

    number_of_fig_rows = 4
    fig, axs = plt.subplots(nrows = number_of_fig_rows, ncols=number_datasets,figsize=(number_datasets*8,number_of_fig_rows*8))

    if number_datasets !=1:

        for idx, element in enumerate(info.keys()):

            # representation of zscore with all molecules
            sns.scatterplot(x="ID_all", y=f"zscore_{element}", data=all_merged,ax=axs[0,idx]).set_title(f"zscore_{element}" )


            # representation of zscore threshold analysis with all molecules
            for plot in range (len(info.keys())):
                ax=axs[0,plot]
                ax.axhline(y=threshold_zscore, linewidth=2, color='red', ls=':')
                ax.axhline(y=-threshold_zscore, linewidth=2, color='red', ls=':')
                ax.set_ylim(-5,5)


            sns.scatterplot(x="ID_all", y=f"zscore_{element}_thres", data=all_merged,ax=axs[1,idx]).set_title(f"zscore_{element}_thres")
            for plot in range (len(info.keys())):
                ax=axs[1,plot]
                ax.set_ylim(-1.2,1.2)

            # representation of zscore with filtered molecules
            sns.scatterplot(x="ID_all", y=f'zscore_{element}_nointraoutliers', data=all_merged,ax=axs[2,idx]).set_title(f'zscore_{element}_nointraoutliers' )
            for plot in range (len(info.keys())):
                ax=axs[2,plot]
                ax.axhline(y=threshold_zscore, linewidth=2, color='red', ls=':')
                ax.axhline(y=-threshold_zscore, linewidth=2, color='red', ls=':')
                ax.set_ylim(-5,5)

            # representation of zscore threshold analysis with filtered molecules
            sns.scatterplot(x="ID_all", y=f'zscore_{element}_nointraoutliers_thres', data=all_merged,ax=axs[3,idx]).set_title(f'zscore_{element}_nointraoutliers_thres')
            for plot in range (len(info.keys())):
                ax=axs[3,plot]
                ax.set_ylim(-1.2,1.2)

    else:


        for idx, element in enumerate(info.keys()):

            # representation of zscore with all molecules
            sns.scatterplot(x="ID_all", y=f"zscore_{element}", data=all_merged,ax=axs[0]).set_title(f"zscore_{element}" )


            # representation of zscore threshold analysis with all molecules

            ax=axs[0]
            ax.axhline(y=threshold_zscore, linewidth=2, color='red', ls=':')
            ax.axhline(y=-threshold_zscore, linewidth=2, color='red', ls=':')
            ax.set_ylim(-5,5)


            sns.scatterplot(x="ID_all", y=f"zscore_{element}_thres", data=all_merged,ax=axs[1]).set_title(f"zscore_{element}_thres")
            ax=axs[1]
            ax.set_ylim(-1.2,1.2)

            # representation of zscore with filtered molecules
            sns.scatterplot(x="ID_all", y=f'zscore_{element}_nointraoutliers', data=all_merged,ax=axs[2]).set_title(f'zscore_{element}_nointraoutliers' )

            ax=axs[2]
            ax.axhline(y=threshold_zscore, linewidth=2, color='red', ls=':')
            ax.axhline(y=-threshold_zscore, linewidth=2, color='red', ls=':')
            ax.set_ylim(-5,5)

            # representation of zscore threshold analysis with filtered molecules
            sns.scatterplot(x="ID_all", y=f'zscore_{element}_nointraoutliers_thres', data=all_merged,ax=axs[3]).set_title(f'zscore_{element}_nointraoutliers_thres')

            ax=axs[3]
            ax.set_ylim(-1.2,1.2)



    fig.suptitle(f'{propertie} intra analysis', fontsize=20)
    output_file = f'{propertie}_zscore_analysis.png'
    fig.savefig(results_folder + os.path.sep + output_file)
    print(f'\t\t[+++] Output file created: {output_file}')



def filter_row(ref_col, row):

    if (row[cols_stdv] > 0.2).any():
        return np.nan
    else:
        return ref_col





def compare_stds(comparisons, all_merged, tag):

    for pair_to_compare in comparisons:

        dataset1 = pair_to_compare[0]
        col_df1 = f'y_{dataset1}_{tag}'
        dataset2 = pair_to_compare[1]
        col_df2 = f'y_{dataset2}_{tag}'
        output_col_std = f'stand_desv_{dataset1}_{dataset2}_{tag}'
        cols_stdv.append(output_col_std)

        comparison_label =  dataset1 + '_' + dataset2

        all_merged.insert(all_merged.shape[1], output_col_std, list(abs(all_merged[[col_df1,col_df2]].std(axis=1, numeric_only=True)/all_merged[[col_df1,col_df2]].mean(axis=1, numeric_only=True))))

        mask = ~np.isnan(all_merged[col_df1]) & ~np.isnan(all_merged[col_df2])


        perc_mols = round( mask[mask==True].shape[0]*100/ all_merged.shape[0],2)

        if mask[mask==True].shape[0] != 0   :
            slope, intercept, r_value, p_value, std_err = stats.linregress(all_merged[col_df1][mask],all_merged[col_df2][mask])
            label = f"y={slope:.3f}x+{intercept:.3f}. R^2: {r_value:.3f}, {perc_mols} % of mols common"
        else:
            slope = intercept = r_value = p_value = std_err = np.nan
            label = "No common molecules"
        if comparison_label not in compare_dict.keys():
            compare_dict[comparison_label] = {}
        compare_dict[comparison_label][f'label_{tag}'] = label
        compare_dict[comparison_label][f'r_value_{tag}'] = r_value




def representation_absdesvest(compare_dict, all_merged, threshold_stand_dev):

    number_comparisons = len(compare_dict.keys())

    number_of_fig_rows = 4
    fig, axs = plt.subplots(nrows = number_of_fig_rows, ncols=number_comparisons,figsize=(number_comparisons*8,number_of_fig_rows*8))

    for idx, (key, element) in enumerate(compare_dict.items()):

        comparison1_element = key.split('_')[0]
        comparison1_column = f'y_{comparison1_element}_nointraoutliers'
        comparison1_column_filtered = f'y_{comparison1_element}_nointraoutliers_nointeroutliers'

        comparison2_element = key.split('_')[1]
        comparison2_column = f'y_{comparison2_element}_nointraoutliers'
        comparison2_column_filtered = f'y_{comparison2_element}_nointraoutliers_nointeroutliers'

        STD_column = f'stand_desv_{comparison1_element}_{comparison2_element}_nointraoutliers'
        STD_column_filtered = f'stand_desv_{comparison1_element}_{comparison2_element}_nointraoutliers_nointeroutliers'


        if number_comparisons !=1:


            list_no_nans = [x for x in list(all_merged[comparison1_column]) if pd.notnull(x) ]
            xmin_lim = min(list_no_nans) - np.mean(list_no_nans)/10
            xmax_lim = max(list_no_nans) + np.mean(list_no_nans)/10

           # representation of stdved with all molecules
            sns.regplot(x=comparison1_column, y=comparison2_column, data=all_merged[[comparison1_column,comparison2_column]],ax=axs[0,idx], line_kws=dict(color="r"),
                        label=compare_dict[key]['label_nointraoutliers']).legend(loc="best")

            for plot,name in enumerate(compare_dict.keys()):
                ax=axs[0,plot]
                ax.set_title(f"std_{name}" )


            sns.scatterplot(x=comparison1_column, y=STD_column, data=all_merged[[comparison1_column,STD_column]], ax=axs[1,idx]).set_title(f"std_{key}_thres")

            # representation of stdved threshold analysis with all molecules
            for plot in range (len(compare_dict.keys())):
                ax=axs[1,plot]
                ax.axhline(y=threshold_stand_dev, linewidth=2, color='red', ls=':')
                ax.set_ylim(-0.1,1)


            # representation of stdved with molecules over threshold

            sns.regplot(x=comparison1_column_filtered, y=comparison2_column_filtered, data=all_merged[[comparison1_column_filtered,comparison2_column_filtered]],ax=axs[2,idx], line_kws=dict(color="r"),
                        label=compare_dict[key]['label_nointraoutliers_nointeroutliers']).legend(loc="best")

            for plot,name in enumerate(compare_dict.keys()):
                ax=axs[2,plot]
                ax.set_title(f"std_{name}" )


            sns.scatterplot(x=comparison1_column_filtered, y=STD_column_filtered, data=all_merged[[comparison1_column_filtered, STD_column_filtered]], ax=axs[3,idx]).set_title(f"std_{key}_thres_filtered")

            # representation of stdved threshold analysis with molecules over threshold
            for plot in range (len(compare_dict.keys())):
                ax=axs[3,plot]
                ax.axhline(y=threshold_stand_dev, linewidth=2, color='red', ls=':')
                ax.set_ylim(-0.1,1)


        else:

            list_no_nans2 = [x for x in list(all_merged[comparison1_column]) if pd.notnull(x) ]
            xmin_lim2 = min(list_no_nans2) - np.mean(list_no_nans2)/10
            xmax_lim2 = max(list_no_nans2) + np.mean(list_no_nans2)/10

           # representation of stdved with all molecules
            sns.regplot(x=comparison1_column, y=comparison2_column, data=all_merged[[comparison1_column,comparison2_column]],ax=axs[0], line_kws=dict(color="r"),
                        label=compare_dict[key]['label_nointraoutliers']).legend(loc="best")



            ax=axs[0]
            ax.set_title(f"std_{key}" )
            ax.set_xlim(xmin_lim2,xmax_lim2)

            sns.scatterplot(x=comparison1_column, y=STD_column, data=all_merged[[comparison1_column,STD_column]], ax=axs[1]).set_title(f"std_{key}_thres")

            # representation of stdved threshold analysis with all molecules
            # for plot in range (len(compare_dict.keys())):
            ax=axs[1]
            ax.axhline(y=threshold_stand_dev, linewidth=2, color='red', ls=':')
            ax.set_ylim(-0.1,1)
            ax.set_xlim(xmin_lim2,xmax_lim2)

            # representation of stdved with molecules over threshold

            sns.regplot(x=comparison1_column_filtered, y=comparison2_column_filtered, data=all_merged[[comparison1_column_filtered,comparison2_column_filtered]],ax=axs[2], line_kws=dict(color="r"),
                        label=compare_dict[key]['label_nointraoutliers_nointeroutliers']).legend(loc="best")


            ax=axs[2]
            ax.set_title(f"std_{key}" )
            ax.set_xlim(xmin_lim2,xmax_lim2)


            sns.scatterplot(x=comparison1_column_filtered, y=STD_column_filtered, data=all_merged[[comparison1_column_filtered, STD_column_filtered]], ax=axs[3]).set_title(f"std_{key}_thres_filtered")

            # representation of stdved threshold analysis with all molecules

            ax=axs[3]
            ax.axhline(y=threshold_stand_dev, linewidth=2, color='red', ls=':')
            ax.set_ylim(-0.1,1)
            ax.set_xlim(xmin_lim2,xmax_lim2)

    fig.suptitle(f'{propertie} inter analysis', fontsize=20)
    output_file = f'{propertie}_AbsSTD_analysis.png'
    fig.savefig(results_folder + os.path.sep + output_file)
    print(f'\t\t[+++] Output file created: {output_file}')

def characterize_binary(val1,val2):


    if pd.isna(val1) and pd.isna(val2):
        return np.nan
    if val1 == val2 == 1:
        return 1
    elif val1 == val2 == 0:
        return 0
    elif val1 == 1 and pd.isna(val2) or val2 == 1 and pd.isna(val1):
        return 1
    elif val1 == 0 and pd.isna(val2) or val2 == 0 and pd.isna(val1):
        return 0

    else:
        return -1


def filter_binary(val1,val2):


    if pd.isna(val1) and pd.isna(val2):
        return np.nan
    if val1 == val2 == 1:
        return 1
    elif val1 == val2 == 0:
        return 0
    elif val1 == 1 and pd.isna(val2) or val2 == 1 and pd.isna(val1):
        return 1
    elif val1 == 0 and pd.isna(val2) or val2 == 0 and pd.isna(val1):
        return 0

    else:
        return np.nan

def representation_quantitative(compare_dict, all_merged):

    number_comparisons = len(compare_dict.keys())
    number_of_fig_rows = 2

    fig, axs = plt.subplots(nrows = number_of_fig_rows, ncols=number_comparisons,figsize=(number_comparisons*8,number_of_fig_rows*8))

    for idx, (key, element) in enumerate(compare_dict.items()):
        comparison1_element = key.split('_')[0]
        comparison2_element = key.split('_')[1]

        output_col = f'comparison_{comparison1_element}_{comparison2_element}'
        output_col_filtered = f'onlyequal_{comparison1_element}_{comparison2_element}'
        label = element["label"]
        cols_comparisons.append(output_col)



        if number_comparisons ==1:

            sns.scatterplot(x="ID_all", y=output_col, data=all_merged,ax=axs[0], label= label).legend(loc="best")

            ax=axs[0]
            ax.set_ylim(-1.5,1.5)
            sns.scatterplot(x="ID_all", y=output_col_filtered, data=all_merged,ax=axs[1]).set_title(f"{key}_filtered" )
            ax=axs[1]
            ax.set_ylim(-1.5,1.5)

        else:


            sns.scatterplot(x="ID_all", y=output_col, data=all_merged,ax=axs[0, idx], label= label).legend(loc="best")
            for plot in range (len(compare_dict.keys())):
                ax=axs[0,plot]
                ax.set_ylim(-1.5,1.5)

            sns.scatterplot(x="ID_all", y=output_col_filtered, data=all_merged,ax=axs[1, idx]).set_title(f"{key}_filtered" )
            for plot in range (len(compare_dict.keys())):
                ax=axs[1,plot]
                ax.set_ylim(-1.5,1.5)

    fig.suptitle(f'{propertie} inter analysis \n[-1 value means ambiguous results]', fontsize=20)
    output_file = f'{propertie}_comparison_analysis.png'
    fig.savefig(results_folder + os.path.sep + output_file)
    print(f'\t\t[+++] Output file created: {output_file}')


def filter_row_quanty(ref_col, row, cols_comparisons):

    if (row[cols_comparisons] == (-1)).any():
        return np.nan
    else:
        return ref_col

def filter_row_quanty_several(row, cols_comparisons):

    if (row[cols_comparisons] == (-1)).any():
        return np.nan
    else:
        values = [x for x in list(row[cols_comparisons])if pd.notnull(x)]

        return values[0]



##############################################################################


#############################  INPUT DATASETS ################################
'''
legend for dictios_dfs

dictios_dfs= {'property' : {'name_dataset': {'file' : file name wo extension}

input datasets should be previously curated and without duplicated molecules
'''

#dictios_for_WP3_ONTOX
dictio_properties = {
    'PhCh_BP':      {'desired_units' : 'C',
                     'classreg' : 'reg',},
    'PhCh_logD':    {'desired_units' : 'logD',
                     'classreg' : 'reg'},
    'PhCh_logH':    {'desired_units' : 'logH',
                     'classreg' : 'reg'},
    'PhCh_logP':    {'desired_units' : 'logKow',
                     'classreg' : 'reg'},
    'PhCh_WS':      {'desired_units' : 'LogMolar',
                     'classreg' : 'reg'},
    'PhCh_VP':      {'desired_units' : 'LogmmHg',
                     'classreg' : 'reg'},
    'PhCh_MP':      {'desired_units' : 'C',
                     'classreg' : 'reg'},
    'PhCh_pKa':     {'desired_units' : 'pKa',
                     'classreg' : 'reg'},
    'PhCh_pKaA':     {'desired_units' : 'pKa',
                     'classreg' : 'reg'},
    'PhCh_pKaB':     {'desired_units' : 'pKa',
                     'classreg' : 'reg'},
    'TK_Caco2':     {'desired_units' : 'logcm/s',
                     'classreg' : 'reg'},
    'TK_FU':        {'desired_units' : 'ratioFU',
                     'classreg' : 'reg'},
    'TK_logKp':     {'desired_units' : 'logcm/h',
                     'classreg' : 'reg'},
    'TK_BBB':       {'desired_units' : 'classBBB',
                     'classreg' : 'class'},
    'TK_F30':       {'desired_units' : 'classF30',
                     'classreg' : 'class'},
    'TK_HIA':       {'desired_units' : 'classHIA',
                     'classreg' : 'class'},
    'TK_Psub' :     {'desired_units' : 'classPsub',
                     'classreg' : 'class'},
    'TK_Pinh':      {'desired_units' : 'classPinh',
                     'classreg' : 'class'},

    }


dictios_dfs = {
    'PhCh_BP' :
                {'Hall':        {'file' : 'PhCh_BP_Hall'},
                 'Katritzky':   {'file' : 'PhCh_BP_Katritzky'},
                 'PhyspropBP' : {'file' : 'PhCh_BP_PHYSPROP'},
                },
    'PhCh_logD' :
                {'Liu':         {'file' : 'PhCh_logD_Liu'},
                 'Wu':          {'file' : 'PhCh_logD_Wu'},
                },
    'PhCh_logH' :
                {'Moldarresi':  {'file' : 'PhCh_logH_Moldarresi'},
                 'Yao2002':         {'file' : 'PhCh_logH_Yao2002'},
                },
    'PhCh_logP' :
                {'PhyspropLogP':{'file' : 'PhCh_logP_PHYSPROP'},
                 'Hughes':      {'file' : 'PhCh_logP_Hughes'},

                 'Martelv2':      {'file' : 'PhCh_logP_Martelv2'},
                },
    'PhCh_WS' :
                {'PhyspropWS':  {'file' : 'PhCh_logS_PHYSPROP'},
                 'Hughes':      {'file' : 'PhCh_logS_Hughes'},
                },
    'PhCh_VP' :
                {'PhyspropVP':     {'file' : 'PhCh_VP_PHYSPROP'},
                 'Katritzky':     {'file' : 'PhCh_logVP_Katritzky'},
                },

    'PhCh_MP' :
                {'PhyspropMP':  {'file' : 'PhCh_MP_PHYSPROP'},
                 'Katrizy':     {'file' : 'PhCh_MP_Katrizy'},
                 'Habibi':      {'file' : 'PhCh_MP_Habibi'},
                 'Hughes':      {'file' : 'PhCh_MP_Hughes'},
                },
    'PhCh_pKaA' :
                {'avdeefAcid':      {'file' : 'PhCh_pKaA_avdeef'},
                 'sanderAcid':      {'file' : 'PhCh_pKaA_sander'},
                 'liaoAcid':        {'file' : 'PhCh_pKaA_liao'},
                 'settimoAcid':     {'file' : 'PhCh_pKaA_settimo'},
                },
    'PhCh_pKaB' :
                {'avdeefBasic':      {'file' : 'PhCh_pKaB_avdeef'},
                 'sanderBasic':      {'file' : 'PhCh_pKaB_sander'},
                 'liaoBasic':        {'file' : 'PhCh_pKaB_liao'},
                 'settimoBasic':     {'file' : 'PhCh_pKaB_settimo'},
                },
    'TK_Caco2' :
                {'PhamThe':     {'file' : 'TK_Caco2_PhamThe'},
                 'Wang':        {'file' : 'TK_Caco2_Wang'},
                },
    'TK_FU' :
                {'Tonnelier':   {'file' : 'TK_FUB_Tonnelier'},
                 'Yamazaki':    {'file' : 'TK_FUB_Yamazaki'},
                 'Lombardo':    {'file' : 'TK_FUB_Lombardo'},
                 'Riley':       {'file' : 'TK_FUB_Riley'},
                 'Votano':      {'file' : 'TK_FUB_Votano'},
                 'cran':        {'file' : 'TK_FUB_cran'},
                 'Zhu':         {'file' : 'TK_FUB_Zhu'},
                },
    'TK_logKp' :
                {'Potts':       {'file' : 'TK_logKp_Potts'},
                 'Khajeh':      {'file' : 'TK_logKp_Khajeh'},
                },
    'TK_BBB' :
                {'Wang':        {'file' : 'TK_BBB_Wang'},
                },
    'TK_F30' :
                {'Fagerholm':   {'file' : 'TK_F30_Fagerholm'},
                 'Kim':        {'file' : 'TK_F30_Kim'},
                },
    'TK_HIA' :
                {'PhamThe':     {'file' : 'TK_HIA_PhamThe'},
                 'Wang':        {'file' : 'TK_HIA_Wang'},
                  },
    'TK_Psub' :
                {'Wang' :       {'file' : 'TK_Psub_Wang'},
                  'Li' :         {'file' : 'TK_Psub_Li'},
                },
    'TK_Pinh' :
                {'Broccatelli': {'file' : 'TK_Pinh_Broccatelli'},
                },

    }





##############################################################################

##########  PROCESS INITIAL DATASETS AND MERGE BY PROPERTY ###################



for propertie, info in dictios_dfs.items():

    print(f'\n[+] Analysing {len(info.keys())} datasets for {propertie}. {dictio_properties[propertie]["classreg"]}')

    propertie_type = dictio_properties[propertie]['classreg']

    all_dfs = []

    names = []
    add_cols = []

    for name, details in info.items():

        df = pd.read_csv(data_folder + os.path.sep + details['file'] + '_final.csv', sep = ';')

        tag = f"df_{propertie}_{name}"
        details['tag'] = tag

        desired_units = dictio_properties[propertie]['desired_units']

        df_ordered, additional_columns = order_input_df(df,name,desired_units )


        add_cols.append(additional_columns)

        details['df_clean'] = df_ordered

        all_dfs.append(df_ordered)
        duplic_mols_check = check_duplicated_molecules(df,name)

        names.append(name)


    all_merged = pd.concat(all_dfs, axis = 1, join="outer")
    all_merged.reset_index(inplace = True)
    all_merged.index.name = 'ID_all'
    all_merged.reset_index(inplace = True)




    ##################### internal outliers analysis ##########################
    if propertie_type == 'reg':
        print('\t[++] Continuos data: performing analysis of internal outliers')
        threshold_zscore = 3

        cols_prop = []

        for element in info.keys():

            output_col_ynooutliers = f'y_{element}_nointraoutliers'

            cols_prop.append(output_col_ynooutliers)

            all_merged[f'zscore_{element}'] = all_merged[[f'y_{element}']].apply(zscore, nan_policy='omit')
            all_merged[f'zscore_{element}_thres']  = all_merged.apply(lambda row: -1 if pd.isna(row[f'zscore_{element}']) else (1 if abs(row[f'zscore_{element}'])>threshold_zscore else 0), axis=1)

            all_merged[output_col_ynooutliers]  = all_merged.apply(lambda row: row[f'y_{element}'] if row[f'zscore_{element}_thres'] == 0 else np.nan, axis=1)
            all_merged[f'zscore_{element}_nointraoutliers'] = all_merged[[f'y_{element}_nointraoutliers']].apply(zscore, nan_policy='omit')
            all_merged[f'zscore_{element}_nointraoutliers_thres']  = all_merged.apply(lambda row: -1 if pd.isna(row[f'zscore_{element}_nointraoutliers']) else (1 if abs(row[f'zscore_{element}_nointraoutliers'])>threshold_zscore else 0), axis=1)

        representation_zscore(info, threshold_zscore)





    else:
        print('\t[++] Discrete data: performing analysis of internal outliers cannot be performed')


    ##################### inter-databases outliers analysis ###################

    comparisons = list(itertools.combinations(info.keys(), 2))

    if len(comparisons) == 0:
        print('\t[++] Only one dataset: performing analysis of interdataset outliers cannot be performed')

        if propertie_type == 'reg':


            all_merged['y_consensus'] = all_merged[f'y_{name}_nointraoutliers']

            all_merged['valid/outlier']  = all_merged.apply(lambda row: 'outlier' if pd.isna(row['y_consensus']) else 'valid', axis=1)

        else:


            all_merged['y_consensus'] = all_merged[f'y_{name}']

            all_merged['valid/outlier']  = all_merged.apply(lambda row: 'outlier' if pd.isna(row['y_consensus']) else 'valid', axis=1)



    else:
        print('\t[++] Performing analysis of interdataset outliers')

        compare_dict = {}


        if propertie_type == 'reg':

            cols_stdv = []

            compare_stds(comparisons, all_merged, tag = "nointraoutliers")

            for ref_col in cols_prop:
                all_merged[f'{ref_col}_nointeroutliers'] = all_merged.apply(lambda row: filter_row(row[ref_col], row), axis=1)

            compare_stds(comparisons, all_merged, tag = "nointraoutliers_nointeroutliers")


            for keycomparison, infocomparison in compare_dict.items():

                print('\t\t[+++] Correlations:')
                print(f'\t\t {keycomparison} :')
                print(f'\t\t\t Before filtering outliers: {infocomparison["label_nointraoutliers"]}' )

                if infocomparison["r_value_nointraoutliers"] < 0.8 : print('\t\t\t\tWARNING!!!!!') ### sttablish threshold!



                print(f'\t\t\t After filtering outliers : {infocomparison["label_nointraoutliers_nointeroutliers"]}' )
            threshold_stand_dev = 0.2
            representation_absdesvest(compare_dict, all_merged, threshold_stand_dev)

            y_columns = [f'y_{x}_nointraoutliers_nointeroutliers' for x in names if '*' not in name]


            all_merged['y_consensus'] = all_merged[y_columns].mean(axis=1)


            all_merged['y_consensus_sd'] = all_merged[y_columns].std(axis=1)

            all_merged['valid/outlier']  = all_merged.apply(lambda row: 'outlier' if pd.isna(row['y_consensus']) else 'valid', axis=1)

            cols_final_y = []

            for name, addit_col in zip(names, add_cols):
                cols_final_y.append(f'ID_{name}')
                cols_final_y.append(f'NAME_{name}')
                cols_final_y.append(f'CAS_{name}')
                cols_final_y.append(f'y_{name}_nointraoutliers_nointeroutliers')
                if len(addit_col) !=0:
                    cols_final_y.append(f'{list(additional_columns)[0]}_{name}')


            final_cols = ['ID_all', 'SMILES', 'y_consensus', 'y_consensus_sd', 'valid/outlier'] + cols_final_y



        else:

            cols_prop = set()

            cols_comparisons = []

            for pair_to_compare in comparisons:

                dataset1 = pair_to_compare[0]
                col_df1 = f'y_{dataset1}'
                dataset2 = pair_to_compare[1]
                col_df2 = f'y_{dataset2}'
                output_col = f'comparison_{dataset1}_{dataset2}'
                output_col_filtered = f'onlyequal_{dataset1}_{dataset2}'
                cols_prop.add(col_df1)
                cols_prop.add(col_df2)



                comparison_label =  dataset1 + '_' + dataset2


                all_merged[output_col] = all_merged.apply(lambda row: characterize_binary(row[col_df1],row[col_df2]), axis=1)


                all_merged[output_col_filtered] = all_merged.apply(lambda row: filter_binary(row[col_df1],row[col_df2]), axis=1)


                output_col_nonans = [x for x in list(all_merged[output_col]) if pd.notna(x)]


                nonequal_number = output_col_nonans.count(-1)
                nonequal = nonequal_number/len(output_col_nonans)
                percent = nonequal*100



                label = f"{percent:.3f} % ({nonequal_number}/{len(output_col_nonans)})"

                if comparison_label not in compare_dict.keys():
                    compare_dict[comparison_label] = {}
                compare_dict[comparison_label]['label'] = label
                compare_dict[comparison_label]['percent'] = percent

            for keycomparison, infocomparison in compare_dict.items():

                label = infocomparison['label']
                print('\t\t[+++] Correlations:')
                print(f'\t\t {keycomparison}:')
                print(f'\t\t\t {label}' )




            representation_quantitative(compare_dict, all_merged)

            for ref_col in list(cols_prop):
                all_merged[f'{ref_col}_nonambiguous'] = all_merged.apply(lambda row: filter_row_quanty(row[ref_col], row, cols_comparisons), axis=1)


            if len(comparisons) == 2:

                all_merged['y_consensus'] = all_merged[f'onlyequal_{names[0]}_{names[1]}']

                all_merged['valid/outlier']  = all_merged.apply(lambda row: 'outlier' if pd.isna(row['y_consensus']) else 'valid', axis=1)


            else:

                all_merged['y_consensus'] = all_merged.apply(lambda row: filter_row_quanty_several(row, cols_comparisons), axis=1)
                all_merged['valid/outlier']  = all_merged.apply(lambda row: 'outlier' if pd.isna(row['y_consensus']) else 'valid', axis=1)



            cols_final_y = []

            for name, addit_col in zip(names, add_cols):
                cols_final_y.append(f'ID_{name}')
                cols_final_y.append(f'NAME_{name}')
                cols_final_y.append(f'CAS_{name}')
                cols_final_y.append(f'y_{name}_nonambiguous')
                if len(addit_col) !=0:
                    cols_final_y.append(f'{list(additional_columns)[0]}_{name}')




    output_file_df = f'all_merged_{propertie}.csv'
    print(f'\t[++] All merged dataframe file created: {output_file_df}')
    all_merged.to_csv(results_folder + os.path.sep + output_file_df, sep = ';', index = False)
