# -*- coding: utf-8 -*-

"""

This module deconstructs the communal trophic network to sub-networks
sustaining the rule of stemming from an exudate, via GSMMs as proxies, 
to other metabolites. 
This module yields sub networks of various lengths (3-6), but only 
sub-networks in lengths of 3 (exudate-GSMM-final_metabolite), and 5 
(exudate-GSMM1-intermediate_metabolite-GSMM2-final_metabolite) were kept 
and used in order to maintain biological relevance. 

Sub-networks were functionally classified based on the differential abundance 
of bacteria participating in them, allowing their characterization in terms of 
microbial metabolic activity effect over soil health.

Written by Einam Castel, 18.4.23

"""

import networkx as nx
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import os

# =============================================================================
# parameters for the functions: to be filled by user.
# =============================================================================

base_dir = ''
network_dir = base_dir + 'target/network/'  # path to the directory holding network data
net_filename = 'network_tuples_df.csv'  # name of network csv file - a pandas dataframe of tuples (two columns, each holding one tuple)
exudates_dir = base_dir + 'media/'  # path to the directory holding exudates data
ex_filename = 'exudates.csv'  # name of the exudates csv file - a pandas dataframe holding exudate identifiers names.
edge_course_list = ['BjSA', 'NTC', 'NA']  # a list holding names (strings) of differential abundance classifications.
da_file_name = 'DA_classification.csv' # a string of the name of the dataframe holding classification of GSMMs to different treatments
target_dir = base_dir + 'target/paths/' # a directory to hold all files generated throughout the analyses


# =============================================================================
# utility functions - for creating variables
# =============================================================================

def get_DA(base_dir, da_file_name):
    '''
  returns dataframe with GSMM and their differential abundance (DA) classification.

    '''
    DA = pd.read_csv(base_dir + da_file_name)
    return DA

def get_exudates(exudates_dir, ex_filename):
    """
    Returns exudates list from exudates csv file.
    """
    exudates = pd.read_csv(exudates_dir + ex_filename)['metabolite'].to_list()
    return exudates

def build_G(net_filename, network_dir):
    """
   builds the network object, using networkx library

  Parameters
  ----------
  PATH : str
      A path leading for the network data (list of tuples representing edges)

  Returns
  -------
  G : networkx directed grpah object

  """
    G = nx.DiGraph()
    df = pd.read_csv(network_dir + net_filename, index_col=0)
    l = list(df.itertuples(index=False, name=None))
    G.add_edges_from(l)
    return G


def create_parquet_dir(base_dir):
    '''
    creates a directory to hold all parquest tables 

    Parameters
    ----------
    base_dir : str
        path to base dir.

    Returns
    -------
    None.

    '''
    directory = 'parquet'
    path = os.path.join(base_dir, directory)
    os.mkdir(path)
    return path + '/'
    

def make_parquet_tables(G, EXUDATES, praquet_dir_path):
    """
      Generates a table of exchange paths for each exudate, and saves it in a parquet table in the designated directory.

      Parameters
      ----------
      net_filename : str
          A directory of the network data

      Returns
      -------
      None.

  """

    min_len = 3
    nodes = list(G.nodes())
    for start in EXUDATES:
        dic = {start: []}
        for end in nodes:
            try:
                for path in nx.all_shortest_paths(G, start, end, weight=None, method='dijkstra'):
                    if len(path) >= min_len:
                        dic[start].append(path)
                        if len(path) < min_len:
                            print()
            except:
                continue
        table = pa.table(dic)
        pq.write_table(table, praquet_dir_path + '/' + start + '.parquet', use_dictionary=True)
    

# =============================================================================
# global variables creation
# =============================================================================

parquet_dir_path = create_parquet_dir(base_dir)
DA = get_DA(base_dir, da_file_name)
EXUDATES = get_exudates(exudates_dir, ex_filename)
G = build_G(net_filename, network_dir)

# =============================================================================
# Data generation functions
# =============================================================================

def get_classification(edge_course_list):
    '''
  Given a single path, returns a tuple of GSMM types (in length of 2).

  Parameters
  ----------
  edge_course_list : list
      List of items in a single edge_course_list.

  Returns
  -------
  TYPE
      Tuple containing all DA classifications in a path (either length 1 or 2).

  '''

    GSMM_list = []
    for item in edge_course_list:
        if item[0] == "G":
            GSMM_list.append(DA['DA_final_score'][DA.index[DA['GSMM'] == item]].values[0])
            res = ['NA' if pd.isnull(i) else i for i in GSMM_list]
    return tuple(res)


def build_PMM_edge_course_matrix(parquet_dir, target_dir, length=5):
    '''

  This function builds a CSV file for edge_courses, with emphasis on the
  attrbiute in each position. length=5 for PMMs.

  Parameters
  ----------
  praquet_dir : str
      directoriy location of all paquets tables.
  length : int, optional
      Length of edge course path. Default is 5.

  Returns
  -------
  pairs_dict: dict
      a dictionary holds all PMM course-edge path with respect to position and attributes.

  '''
    pairs_dict = {'exudate': [],
                  'GSMM1': [],
                  'metabolite1': [],
                  'GSMM2': [],
                  'metabolite2': [],
                  'classification': []}

    for exudate in EXUDATES:
        table = pq.read_table(parquet_dir + exudate + '.parquet').to_pandas()
        for edge in table[exudate]:
            if len(edge) == length and len(edge) > 0:                
                pairs_dict['exudate'].append(edge[0])
                pairs_dict['GSMM1'].append(edge[1])
                pairs_dict['metabolite1'].append(edge[2])
                pairs_dict['GSMM2'].append(edge[3])
                pairs_dict['metabolite2'].append(edge[4])
                pairs_dict['classification'].append(get_classification(edge))

    df = pd.DataFrame.from_dict(pairs_dict)
    # df = df.fillna('NA')
    df.to_csv(target_dir + 'PMM_edge_courses_df.csv')
    return df


def build_PM_edge_course_matrix(praquet_dir, target_dir, length=3):
    """

  This function builds a CSV file for edge_courses, with emphasis on the
  attrbiute in each position. length=3 for PMs.


  Parameters
  ----------
  parquet_dir : str
      directory location.
  length : int, optional
      DESCRIPTION. The default is 3.

  Returns
  -------
  pairs_dict: dict
      a dictionary holds all PM course-edge path with respect to position and attributes.

  """
    pairs_dict = {'exudate': [],
                  'GSMM1': [],
                  'metabolite1': [],
                  'classification': []}

    for exudate in EXUDATES:

        table = pq.read_table(praquet_dir + exudate + '.parquet').to_pandas()

        for edge in table[exudate]:
            if len(edge) == length and len(edge) > 0:
                pairs_dict['exudate'].append(edge[0])
                pairs_dict['GSMM1'].append(edge[1])
                pairs_dict['metabolite1'].append(edge[2])
                pairs_dict['classification'].append(get_classification(edge))

    df = pd.DataFrame.from_dict(pairs_dict)
    df['classification'] = [x[0] for x in df['classification']]
    df = df.fillna('NA')
    df.to_csv(target_dir + 'PM_edge_courses_df.csv')
    return df

if __name__ == "__main__":
    
    make_parquet_tables(G, EXUDATES, parquet_dir_path)
    df_PMM = build_PMM_edge_course_matrix(parquet_dir_path, target_dir, length=5)
    df_PM = build_PM_edge_course_matrix(parquet_dir_path, target_dir, length=3)
