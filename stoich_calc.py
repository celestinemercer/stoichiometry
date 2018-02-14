# stoich_calc.py

import numpy as np
import pandas as pd
from mineral_system import MineralSystem
import util

# anhydrous_silicates_stoich function.
def anhydrous_silicates_stoich(dataset,sram_lib,min_systems):
  # Prompt user to select mineral system; get active columns.
  asms = prompt_min_sys(min_systems,False,True)
  active_cols = find_active_columns(dataset,asms)
  # Perform stoichiometry calculations.
  formula_cations = calc_cations_per_formula_unit(dataset,active_cols,asms,sram_lib)
  # Report success.
  print('Anhydrous stoichiometric calculations complete.\n')
  return formula_cations

# hydrous_silicates_stoich function.
def hydrous_silicates_stoich(dataset,sram_lib,min_systems):
  # Prompt user to select mineral system; get active columns.
  hsms = prompt_min_sys(min_systems,True,True)
  active_cols = find_active_columns(dataset,hsms)
  # Normalize anhydrous oxides.
  normalized_dataset = normalize_anhydrous_oxides(dataset,active_cols,hsms)
  # Perform stoichiometry calculations.
  formula_cations = calc_cations_per_formula_unit(normalized_dataset,active_cols,hsms,sram_lib)
  # Report success.
  print('Hydrous stoichiometric calculations complete.\n')
  return formula_cations

# --------------------------------------------------------------------------------
# Common functions.
# --------------------------------------------------------------------------------

# prompt_min_sys function.
def prompt_min_sys(min_systems,hydrous,silicate):
  '''
  Convenience function to prompt the user to select a mineral system. This function
  has the same arguments as the util.filter_mineral_systems function.
  '''
  filtered_systems = util.filter_mineral_systems(min_systems,hydrous,silicate)
  syslist = list(filtered_systems.keys())
  idx = util.prompt_options('Please select a mineral system',syslist)
  sms = filtered_systems[syslist[idx]]
  print('You chose the {:} mineral system.'.format(sms.getSystemName()))
  print('The number of anions per formula unit is: {:}'.format(sms.getAnionsPerFormulaUnit()))
  return sms

# find_active_columns function.
def find_active_columns(dataset,mineral_sys):
  '''
  Returns a map of the column names (keys) in the specified dataset that contain
  the oxides (values) in the specified MineralSystem.
  '''
  # Gather list of oxides that exist in the dataset, and in the selected mineral system.
  cols = list(dataset.columns)
  active_cols = {}
  # Outer for loop; traverse list of oxides.
  for ox in mineral_sys.getOxideList():
    # Inner for loop; traverse remaining columns in dataset.
    for cname in cols:
      if ox in cname:
        active_cols[cname] = ox
        cols.remove(cname)
        break # Break out of the inner loop.
  return active_cols

# new_results_dataframe function.
def new_results_dataframe(dataset,active_cols,min_sys,intermediate=True):
  '''
  Creates a results dataframe for the specified dataset, active columns, and
  MineralSystem. If 'intermediate' is set to 'True', the column names will be
  those from the dataset; otherwise, they will be renamed and sorted according
  to the MineralSystem.
  '''
  # Get row and column names from dataset; create empty dataset.
  rnames = dataset.index
  cnames = list(active_cols.keys())
  data = np.empty((len(rnames),len(cnames)))
  data.fill(np.nan)
  # Modify column names if not an intermediate structure.
  if not intermediate:
    for cni, colname in enumerate(cnames):
      ecounts = util.parse_compound(active_cols[colname])
      elems = list(ecounts.keys())
      cnames[cni] = '{:}/{:} {:}'.format(elems[0],min_sys.getAnionsPerFormulaUnit(),elems[1])
    # Add 'Total' column.
    cnames.append('Total')
    data = np.hstack([data,np.zeros((len(rnames),1))])
  # Construct and return new dataframe.
  return pd.DataFrame(data,rnames,cnames)

# normalize_anhydrous_oxides function.
def normalize_anhydrous_oxides(dataset,active_cols,min_sys):
  '''
  Normalizes the anhydrous oxides in the specified dataset. Returns a copy of the
  specified dataset with the anhydrous oxides normalized to 100.
  '''
  normalized_dataset = dataset.copy()
  # Collect column indices.
  col_indices = []
  for colname in active_cols.keys():
    col_indices.append(list(normalized_dataset.columns).index(colname))
  # Normalize one analysis at a time.
  for r in range(len(normalized_dataset)):
    temp_sum = 0.0
    for cidx in col_indices:
      if normalized_dataset.iloc[r,cidx] > 0.0: # This makes sure the avoid empty values.
        temp_sum += normalized_dataset.iloc[r,cidx]
    if temp_sum > 0.0:
      norm_factor = 100.0/temp_sum
    else: 
      print('WARNING: No data for analysis: {:}'.format(normalized_dataset.index[r]))
    for cidx in col_indices:
      normalized_dataset.iloc[r,cidx] *= norm_factor
  return normalized_dataset

# calc_cations_per_formula_unit function.
def calc_cations_per_formula_unit(dataset,active_cols,min_sys,sram_lib):
  '''
  Computes the number of cations per formula unit for the specified dataset
  and MineralSystem. Returns a dataframe of results.
  '''
  # Prepare data structures.
  molecular_props = new_results_dataframe(dataset,active_cols,min_sys)
  anion_props = new_results_dataframe(dataset,active_cols,min_sys)
  cation_props = new_results_dataframe(dataset,active_cols,min_sys)
  formula_cations = new_results_dataframe(dataset,active_cols,min_sys,False)
  # Compute cation and anion proportions.
  for colname in active_cols.keys():
    # Calculate molecular weight and parse compound.
    weight = util.molecular_weight(active_cols[colname],sram_lib)
    element_counts = util.parse_compound(active_cols[colname])
    elements = list(element_counts.keys())
    # Calculate molecular proportions.
    molecular_props[colname] = dataset[colname]/weight
    # Calculate atomic proportion of oxygen from each molecule. (Anion proportion)
    anion_props[colname] = molecular_props[colname]*element_counts[elements[1]]
    # Calculate atomic proportion of cation from each molecule. (Cation proportion)
    cation_props[colname] = molecular_props[colname]*element_counts[elements[0]]
  # Collect column indices.
  col_indices = []
  for colname in active_cols.keys():
    col_indices.append(list(anion_props.columns).index(colname))
  # Determine conversion factors from anions.
  factors = []
  for r in range(len(anion_props)):
    temp_sum = 0.0
    for cidx in col_indices:
      if anion_props.iloc[r,cidx] > 0.0: # This makes sure to avoid empty values.
        temp_sum += anion_props.iloc[r,cidx]
    factors.append(min_sys.getAnionsPerFormulaUnit()/temp_sum)
  # Determine number of cations (cation proportions * factors).
  for r in range(len(cation_props)):
    temp_sum = 0.0
    for cidx in col_indices:
      element_counts = util.parse_compound(active_cols[colname])
      elements = list(element_counts.keys())
      if cation_props.iloc[r,cidx] > 0: # This makes sure to avoid empty values.
        formula_cations.iloc[r,cidx] = cation_props.iloc[r,cidx]*factors[r]
        temp_sum += formula_cations.iloc[r,cidx]
    formula_cations.iloc[r,list(formula_cations.columns).index('Total')] = temp_sum
  return formula_cations
