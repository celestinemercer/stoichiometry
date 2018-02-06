from mineral_system import MineralSystem
import util

# anhydrous_silicates_stoich function.
def anhydrous_silicates_stoich(dataset,sram_lib,min_systems):
  # Prompt user to select mineral system.
  anhydrous_silicate_systems = util.filter_mineral_systems(min_systems,False,True)
  syslist = list(anhydrous_silicate_systems.keys())
  idx = util.prompt_options('Please select a mineral system',syslist)
  asms = anhydrous_silicate_systems[syslist[idx]]
  print('You chose the {:} mineral system.'.format(asms.getSystemName()))
  print('The number of anions per formula unit is: {:}\n'.format(asms.getAnionsPerFormulaUnit()))
  # Gather list of oxides that exist in the dataset, and in the selected mineral system.
  cols = list(dataset.columns)
  active_cols = {}
  # outer for loop
  for ox in asms.getOxideList():
    # inner for loop
    for cname in cols:
      if ox in cname:
        active_cols[cname] = ox
        cols.remove(cname)
        break # break out of the inner loop
  # Perform stoichiometry calculations.
  molecular_props = dataset.copy()
  anion_props = dataset.copy()
  cation_props = dataset.copy()
  formula_cations = dataset.copy()
  for colname in active_cols.keys():
    # Calculate molecular weight and parse compound.
    weight = util.molecular_weight(active_cols[colname],sram_lib)
    element_counts = util.parse_compound(active_cols[colname])
    elements = list(element_counts.keys())
    # Calculate molecular proportions.
    molecular_props[colname] = molecular_props[colname]/weight
    # Calculate atomic proportion of oxygen from each molecule. (Anion proportion)
    anion_props[colname] = molecular_props[colname]*element_counts[elements[1]]
    # Calculate atomic proportion of cation from each molecule. (Cation proportion)
    cation_props[colname] = molecular_props[colname]*element_counts[elements[0]]
  # Determine conversion factors from anions.
  factors = []
  for r in range(len(anion_props)):
    temp_sum = 0.0
    for colname in active_cols.keys():
      cidx = list(anion_props.columns).index(colname)
      if anion_props.iloc[r,cidx] > 0: # This makes sure to avoid empty values.
        temp_sum += anion_props.iloc[r,cidx]
    factors.append(asms.getAnionsPerFormulaUnit()/temp_sum)
  # Determine number of cations (cation proportions * factors).
  new_colnames = list(formula_cations.columns)
  for r in range(len(cation_props)):
    temp_sum = 0.0
    for colname in active_cols.keys():
      cidx = list(cation_props.columns).index(colname)
      element_counts = util.parse_compound(active_cols[colname])
      elements = list(element_counts.keys())
      if cation_props.iloc[r,cidx] > 0: # This makes sure to avoid empty values.
        formula_cations.iloc[r,cidx] = cation_props.iloc[r,cidx]*factors[r]
        temp_sum += formula_cations.iloc[r,cidx]
        new_colnames[cidx] = '{:}/{:} Anions'.format(elements[0],asms.getAnionsPerFormulaUnit())
    formula_cations.iloc[r,list(formula_cations.columns).index('Total')] = temp_sum
  formula_cations.columns = new_colnames
  # Report success.
  print('Stoichiometric calculations are complete.\n')
  return formula_cations

# hydrous_silicates_stoick function.

# non_silicates_stoich function.
