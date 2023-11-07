# -*- coding: utf-8 -*-

"""
Created on Tue Aug 30 10:59:09 2022

This module takes the output of MCSM simulations for characterizing all 
metabolic interactions taking place within the microbial community 
(iteration 5 at Ginatt et. al.) to construct a trophic network linking the 
growing species. 

Secretions and uptake data are gathered from the MCSM and directly from the 
models, based on their exchanges, respectively. Compounds are those present 
in the intial rhizosphere environemnt along with secreted compounds along 
iterations. 

Secretion exchanges are filtered to have only organic compounds (excluding
sink compounds and inorganic compounds) from MCSM output. 

The output of this module is an array of tuples representing the 
potential directed metabolic interactions between community members in terms 
of uptake and secretion of metabolites.


@author: ginatta
"""

# =============================================================================
# librarys
# =============================================================================

import pandas as pd
import cobra 
import os
from collections import defaultdict
from itertools import chain
from operator import methodcaller


# =============================================================================
# paths and directories
# =============================================================================

base_dir = ''
models_dir = base_dir + 'models/' # directory of all GSMMs xml files
secretions_dir = base_dir + 'target/secretions/' # all bacterial secretions from MCSM are stored within the target directory (prior stage)
medium_dir = base_dir + 'target/media/' # all media derived from MCSM are stored within the target directory (prior stage)
final_medium_file = 'medium_5.csv' # the final medium composed of compounds secreted along iterations MCSM + the intial medium medium
target_path = base_dir + 'target/network/' # path for saving data


# =============================================================================
# functions
# =============================================================================


def read_all_secretions(secretions_dir):
    df_list = []
    for sec in os.listdir(secretions_dir):
        if sec.startswith("."): continue
        df = pd.read_csv(secretions_dir + sec)
        df = df.set_index('Unnamed: 0', drop=True)
        df_list.append(df)
    return df_list


def turn_secretion_df_to_dict(secretions_dir):
    dict_list = []
    secretions_df_list = read_all_secretions(secretions_dir)
    for df in secretions_df_list:
        secretions_dict = {}
        for col in df.columns:
            s = df[col]
            s = s[s.values < 0]
            secretions_dict[col] = s.index.tolist()
        dict_list.append(secretions_dict)
    return dict_list


def create_default_dict(secretions_dir):
    dict_list = turn_secretion_df_to_dict(secretions_dir)
    dd = defaultdict(list)
    dict_items = map(methodcaller('items'), [(d) for d in dict_list])
    for k, v in chain.from_iterable(dict_items):
        dd[k].extend(v)
    return dd


def generate_final_secretions(secretions_dir):
    dd = create_default_dict(secretions_dir)
    res = {}
    for k, v in dd.items():
        v_mod = list(set(v))
        res[k] = v_mod
    return res


def get_filtered_secretions(secretions_dir):
    organic_compounds_df = pd.read_csv(base_dir + 'organic_metabolites_formulas.csv', index_col=0)
    final_secretions = generate_final_secretions(secretions_dir)
    organic_dict = dict(organic_compounds_df.values)
    filtered_secretions= {}
    for k, v in final_secretions.items():
        l = []
        for i in v:
            if i in organic_dict.keys():
                l.append(i)
                filtered_secretions[k] = l
    return filtered_secretions
                

def create_models(models_dir: str, sceretions_dir):
    secretions = generate_final_secretions(secretions_dir)
    secreting_gsmms = sorted([i for i in list(secretions.keys())])
    models = []
    for model_name in secreting_gsmms:
        print(model_name)
        try:
            model = cobra.io.read_sbml_model(models_dir + model_name + '.xml')
        except:
            continue
        models.append(model)
    return models


def get_exchanges(model):
    sep = ' '
    model_data = []
    data = model.exchanges
    for datum in data:
        datum = str(datum)
        datum = datum.split(sep,1)[0]
        datum = datum.replace(':', '')
        model_data.append(datum)
    return model_data


def get_model_specific_medium_mod(model, medium_dict):
    exchanges = get_exchanges(model)
    temp_li = []
    for k, v in medium_dict.items():
        if k in exchanges:
            temp_li.append(k)
    return temp_li


def get_final_medium(medium_dir, final_medium_file):
    medium = pd.read_csv(medium_dir + final_medium_file, index_col=0)
    medium = dict(medium.values)
    return medium


def get_models_uptakes(secretions_dir, models_dir, medium_dir, final_medium_file):
    models = create_models(models_dir, secretions_dir)
    uptake_dict = {}
    medium = get_final_medium(medium_dir, final_medium_file)
    for model in models:
        medium_specific = get_model_specific_medium_mod(model, medium)
        uptake_dict[model.id] = medium_specific
    return uptake_dict


def set_uptake_tuples(secretions_dir, models_dir, medium_dir, final_medium_file):
    uptake_dict = get_models_uptakes(secretions_dir, models_dir, medium_dir, final_medium_file)
    models_uptakes = []
    for k, v in uptake_dict.items():
        for exc in v:
            tup = (exc, k)
            models_uptakes.append(tup)
    return models_uptakes
        

def set_secretion_tuples(secretions_dir):
    final_secretions = get_filtered_secretions(secretions_dir)
    models_secretion = []
    for k, v in final_secretions.items():
        for exc in v:
            tup = (k, exc)
            models_secretion.append(tup)
    return models_secretion


def main(secretions_dir, models_dir, medium_dir, final_medium_file, target_path=None):
    models_uptakes = set_uptake_tuples(secretions_dir, models_dir, medium_dir, final_medium_file)
    models_secretions = set_secretion_tuples(secretions_dir)
    data = models_uptakes + models_secretions
    return data


if __name__ == '__main__':
    tup_li = main(secretions_dir, models_dir, medium_dir, final_medium_file, target_path=None)
    df = pd.DataFrame(tup_li, columns=['from', 'to'])
    df.to_csv(target_path + 'network_tuples_df.csv')

