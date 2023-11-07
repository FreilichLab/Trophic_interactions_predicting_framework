# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 11:12:54 2022

This module inspects the growth and secretion profiles of Genome-Scale Metabolic 
Models (GSMMs) derived from a specific microbial community. The algorithm 
implemented here gathers the Flux varibility analysis (FVA) solutions of all 
models growth simulations in the rhizosphere environment, extracts their 
secretion profiles, and iteratively adds the secreted compounds to the 
environment (medium) for further growth simulations.

Community growth simulations begin with medium 0 - the intial rhizosphere 
enviornment medium.

Outputs of this module are csv files generated after each iteration, 
indicating the models growing, their secretion profiles and the media generated.

@author: ginatta
"""

# =============================================================================
# libraries, paths and constants
# =============================================================================

import cobra
import pandas as pd
import os

base_dir = ''
media_dir = base_dir + 'media/' # directory of all media
models_dir = base_dir + 'models/' # directory of all models
target_dir = base_dir + 'target/' # target dir in which output files will be saved

# =============================================================================
# functions
# =============================================================================

def read_models(path):
    '''
    A function for creating model objects using COBRA
    which represent different microbial species

    Parameters
    ----------
    path : str
        the path of all models

    Returns
    -------
    models : list
        array of COBRA-model objects.
    '''
    models = []
    for model_name in os.listdir(path):
        model = cobra.io.read_sbml_model(path + '/' + model_name)
        print(model.id)
        models.append(model)
    return models


def get_exchanges(model):
    '''
    A function to gather all exchanges a model holds

    Parameters
    ----------
    model : COBRA model object

    Returns
    -------
    model_data : list
        array of exchanges.

    '''
    sep = ' '
    model_exchanges = []
    data = model.exchanges
    for datum in data:
        datum = str(datum)
        datum = datum.split(sep,1)[0]
        datum = datum.replace(':', '')
        model_exchanges.append(datum)
    return model_exchanges


def get_model_specific_medium(model, medium_dict):
    '''
    A function for specifically modifying the medium dict, to suit a specific model.
    This function makes sure that the environment medium used by a model, would be
    composed of exchanges found within the model.

    Parameters
    ----------
    model : COBRA-model object

    medium_dict : dict
        A dictionary representing the medium in which the model is intended to grow

    Returns
    -------
    specific_medium : dict
        A suited specific medium to be used by the model

    '''
    exchanges = get_exchanges(model)
    temp_li = []
    for k, v in medium_dict.items():
        if k in exchanges:
            temp_li.append(k)
    specific_medium = {i: 1000 for i in temp_li}
    return specific_medium


def get_secretion_profile(model):
    '''
    A function for gathering the FVA secretion profile of a model grown in a
    specific medium.Importantly, this function assumes that the desired medium
    which the model should be testedin, has already been assigned to the model
    via the model.medium attribute.

    Parameters
    ----------
    model : COBRA-model object

    Returns
    -------
    solution : pandas DataFrame
        A dataframe representing the secretion profile of a model

    '''
    print(model)
    # assuming desired medium has already been assigned to model as medium
    solution = model.summary(fva=0.9).secretion_flux # this is a df
    solution.drop(['reaction', 'metabolite', 'flux', 'maximum'], axis=1, inplace=True)
    solution = solution.loc[~(solution==0).all(axis=1)]
    return solution


def create_secretions_df(secretion_df_list ,relevant_models):
    '''
    A function for gathering the secretion profiles of all models
    who grew in a specific growth iteration

    Parameters
    ----------
    secretion_df_list : list
        An array of secretion profiles dataframes, of models who grew in that
        specific iteration
    relevant_models : list
        an array of models who grew in the environment in that specific
        iteration

    Returns
    -------
    secretion_df : pandas DataFrame
        a matrix represnting the secretion profiles of all models
        who grew in that specific iteration.

    '''
    df_dict = {}
    for model, df_ in zip(relevant_models, secretion_df_list):
        df_dict[model.id] = df_
    mod_df_dict = {k: v for k, v in df_dict.items() if v.empty != True}
    dfs = list(mod_df_dict.values())
    names = [name for name in list(mod_df_dict.keys())]
    secretion_df = pd.concat(dfs, axis=1)
    secretion_df.columns = names
    return secretion_df


def iterate_growth_and_secretion(models, medium):
    '''
    An operative function using all the above functions to simulate the growth
    of community GSMMs in a specific enviornment, as well as obtaining both
    the secretion profile of the community, and consequently the medium for
    the next iteration to be executed.

    Parameters
    ----------
    models : list
        array of all GSMM model objects.

    medium : dict
        A dict representing the medium used in that specific iteration.

    Returns
    -------
    A tuple holding:

        1. next_round_medium
        2. growths_dict
        3. secretions_combined_df:

        next_round_medium : dict
            a dictionary holding all compounds secreted by the bacteria grew in
            the last round, added to the compounds already present from the last round

        growths_dict : dict
            dictionary holding models as keys and their growth as value for each round

        secretions_combined_df : df
            .df holding all secretion values for each itereations, added to the one before.

    '''
    growths_dict = {}
    secretion_dfs = []
    relevant_models = []
    for model in models:
        medium_ = get_model_specific_medium(model, medium) # assign the correct medium as input
        model.medium = medium_
        growth = model.slim_optimize()
        if growth > 0:
            relevant_models.append(model)
            growths_dict[model.id] = growth
            secretion_profile = get_secretion_profile(model)
            secretion_dfs.append(secretion_profile)

    secretions_combined_df = create_secretions_df(secretion_dfs, relevant_models)
    secretions_combined = secretions_combined_df.index.tolist()
    secretions_combined = {i : 1000 for i in secretions_combined}
    next_round_medium = {**medium, **secretions_combined}
    return next_round_medium, growths_dict, secretions_combined_df



def medium_dict_to_df(medium, target_path, medium_name: str):
    '''
    A function to save the medium dict after each iteration as a dataframe.

    Parameters
    ----------
    medium : dict
        a dictionary of the new medium composed of secretions at current
        iteration and the former iteration's medium.

    target_path : str
        A path for the directory in which the medium is saved

    medium_name : str
        The name of the given medium

    Returns
    -------
    None.

    '''
    df_ = pd.DataFrame(medium.items(), columns=['exchange', 'flux'])
    df_.to_csv(target_path + 'media/' + medium_name + '.csv')

def save_secretion_and_growths_dfs(growths, secretions_df, number: str):
    '''
    A function to save the secretions and growths dataframes.

    Parameters
    ----------
    growths : dict
        A dictionary holding the growth values of all GSMMs who grew in
        the current iteration.
    secretions_df : pandas DataFrame
        A dataframe representing the secretions of commnity members at
        that specific iteration.
    number : str
        a string of a number representing the iteration in which the growth
        and secretion data were gathered

    Returns
    -------
    None.

    '''
    growth_df = pd.DataFrame(growths.items(), columns=['GSMM', 'growth'])
    growth_df.to_csv(target_dir + 'growths/' + 'growths_' + number + '.csv')
    secretions_df.to_csv(target_dir + 'secretions/' + 'secretion_' + number + '.csv')



# =============================================================================
# =============================================================================
# ITERATION RUNS ARE SEPARATED SINCE RUNTIME IS LONG AND DATA ACQUISITION IS 
# NECESSARY 
# =============================================================================
# =============================================================================

# =============================================================================
# preliminary runs of the iterations.
# =============================================================================

models = read_models(models_dir)

medium_0 = pd.read_csv(media_dir + 'initial_root_environment_medium_for_iterations.csv', index_col=0)
medium_0 = dict(medium_0.values)

# =============================================================================
# iteration 1
# =============================================================================

medium_1, growths_1, secretions_1 = iterate_growth_and_secretion(models, medium_0)

medium_dict_to_df(medium_1, target_dir, 'medium_1')
save_secretion_and_growths_dfs(growths_1, secretions_1, '1')

# =============================================================================
# iteration 2
# =============================================================================


medium_2, growths_2, secretions_2 = iterate_growth_and_secretion(models, medium_1)

medium_dict_to_df(medium_2, target_dir, 'medium_2')
save_secretion_and_growths_dfs(growths_2, secretions_2, '2')

# =============================================================================
# iteration 3
# =============================================================================

medium_3, growths_3, secretions_3 = iterate_growth_and_secretion(models, medium_2)

medium_dict_to_df(medium_3, target_dir, 'medium_3')
save_secretion_and_growths_dfs(growths_3, secretions_3, '3')

# =============================================================================
# iteration 4
# =============================================================================

P_compounds = pd.read_csv(media_dir + 'P_compounds_from_minimal_media.csv', index_col=0)
P_compounds = dict(P_compounds.values) # an addition of more compounds to medium

comb_medium = {**P_compounds, **medium_3}

medium_4, growths_4, secretions_4 = iterate_growth_and_secretion(models, comb_medium)

medium_dict_to_df(medium_4, target_dir, 'medium_4')
save_secretion_and_growths_dfs(growths_4, secretions_4, '4')

# =============================================================================
# iteration 5
# =============================================================================

medium_5, growths_5, secretions_5 = iterate_growth_and_secretion(models, medium_4)

medium_dict_to_df(medium_5, target_dir, 'medium_5')
save_secretion_and_growths_dfs(growths_5, secretions_5, '5')
