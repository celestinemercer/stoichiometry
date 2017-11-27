# Stoichometry
# util.py
#
# This file contains utility funcitons for the stoichometry module.
#
# Author(s):
# - Cameron M. Mercer
#
# Copyright (C), 2017
#

sram_types = {'unknown':'standard relative atomic mass not known',
              'interval':'standard relative atomic mass varies widely in natural materials, so an interval is published',
              'quantity':'standard relative atomic mass given as a value with decisional uncertainties',
              'most_stable':'standard relative atomic mass is that of the most stable isotope of the element of interest'}

def load_atomic_weights(filepath):
  '''
  Loads atomic weights from the specified NIST linearized ASCII file,
  and returns a dictionary of elements with their atomic weights.

  Arguments:
  - filepath - the path to the NIST datafile.
  '''
  # Create dict container for elements data.
  elems = dict()
  # Read file and parse the contents.
  with open(filepath,'r') as f:
    lines = f.readlines()
  # Define variables.
  cz, lz = 0, 0     # Current and last Z values (proton number).
  ce = ''           # Current element symbol.
  mse, msa = [], [] # Most stable nuclide search lists (element, mass number A).
  for line in lines:
    if line == "":
      # Skip blank lines
      continue
    else:
      # Parse relevant information from each line.
      items = line.split(' = ')
      for i in range(len(items)):
        items[i] = items[i].strip()
      if items[0] == 'Atomic Number' and cz == 0:
        cz = int(items[1])
      elif items[0] == 'Atomic Symbol' and cz > 0 and cz != lz:
        ce = items[1]
        elems[ce] = {'Z': cz}
      elif items[0] == 'Mass Number' and cz > 0 and cz != lz:
        ca = int(items[1])
        elems[ce]['A'] = ca
        elems[ce]['N'] = elems[ce]['A'] - elems[ce]['Z']
      elif items[0] in ['Relative Atomic Mass','Isotopic Composition']:
        # Skip these lines of information.
        continue
      elif items[0] == 'Standard Atomic Weight' and cz > 0 and cz != lz:
        # Parse the Standard Relative Atomic Mass, a.k.a. Standard Atomic Weight, if possible.
        if not items[1]:
          elems[ce]['SRAM'] = {'Type':'unknown'}
        else:
          sram = items[1]
          if '[' in sram and ']' in sram:
            if ',' in sram:
              # The sram is given as an interval.
              bounds = sram.replace('[','').replace(']','').split(',')
              elems[ce]['SRAM'] = {'Type':'interval','lower':float(bounds[0]),'upper':float(bounds[1])}
            else:
              # The sram indicates the mass number of the most stable isotope.
              elems[ce]['A'] = int(sram.replace('[','').replace(']',''))
              elems[ce]['N'] = elems[ce]['A'] - elems[ce]['Z'] # Recalculate N.
              mse.append(ce)
              msa.append(elems[ce]['A'])
          else:
            # The sram is given as a value with an uncertainty.
            nums = sram.replace(')','').split('(')
            ndec = len(nums[0]) - nums[0].index('.') - 1
            uncert = '0.' + '0'*(ndec - 1) + nums[1]
            elems[ce]['SRAM'] = {'Type':'quantity','value':float(nums[0]),'uncertainty':float(uncert)}
      elif items[0] == 'Notes' and cz > 0:
        if cz != lz:
          elems[ce]['Notes'] = items[1]
        # Set cz to 0 to indicate we are done parsing data for this element.
        lz = cz
        cz = 0
  # Traverse the file again looking for RAM values of most stable nuclides.
  ca = 0
  for line in lines:
    if line == "":
      # Skip blank lines
      continue
    else:
      # Parse relevant information from each line.
      items = line.strip().split(' = ')
      if items[0] == 'Atomic Symbol':
        ce = items[1]
      elif items[0] == 'Mass Number':
        ca = int(items[1])
      elif items[0] == 'Relative Atomic Mass':
        if ce in mse and ca == msa[mse.index(ce)]:
          nums = items[1].replace(')','').split('(')
          ndec = len(nums[0]) - nums[0].index('.') - 1
          uncert = '0.' + '0'*(ndec - len(nums[1])) + nums[1]
          elems[ce]['SRAM'] = {'Type':'most_stable','value':float(nums[0]),'uncertainty':float(uncert)}
  # Return the element dictionary.
  return elems
  
