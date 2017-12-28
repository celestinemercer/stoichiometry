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

import re # Import module for regular expressions.

# Variables and methods for loading atomic weight information.
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

def update_sram_lib(subject,patch):
  '''
  Replaces element information stored in the subject sram_lib with that of the
  patch sram_lib, only for those elements that are present in the patch lib.
  The updated subject lib is returned.
  '''
  for (k,v) in patch.items():
    subject[k] = v
  return subject

# Regular expressions and methods for parsing chemical formulas.
re_elem = re.compile(r'([A-Z][a-z]?)([0-9]+)?') # elem, with repeats
re_enor = re.compile(r'([A-Z][a-z]?)')          # elem, no repeats
re_lpar = re.compile(r'\(')                     # opening left paren
re_rpar = re.compile(r'\)([0-9]+)?')            # closing right paren, with repeats

def parse_compound(comp):
  '''
  Parses the specified chemical compound string and returns a dict of the
  constituent elements, with the numbers of each element that are present
  in the compound formula. Calls the rpc helper method.
  '''
  return rpc(comp)[0]

def rpc(comp,stack=[],nolp=0):
  '''
  Recursive helper function for parse_compound() method.
  
  Arguments:
  - comp  - The chemical compound to parse.
  - stack - List of dicts containing elements and how many of them there are.
  - nolp  - The number of opening left parentheses that have not been closed.
  '''
  # Check for bad input.
  if not comp:
    print('Error - cannot parse empty chemical formula.')
    return
  # Declare variables for formula fragments. 
  tail, head, num = [], [], 1
  # Initialize the stack (if stack == []), or use existing one (e.g., if stack == [{}]).
  stack = stack or [{}] 
  # Test chemical formula for RE matches.
  em = re_elem.match(comp)
  lm = re_lpar.match(comp)
  rm = re_rpar.match(comp)
  if em:
    head = em.group()
    tail = comp[len(head):]
    elem = re_enor.match(head).group()
    num = 1
    if len(head) > len(elem):
      num = int(head[len(elem):])
    if elem not in stack[-1]:
      # First time adding this element to the dict.
      stack[-1][elem] = num
    else:
      # Increment number for this element (head).
      stack[-1][elem] += num
  elif lm:
    tail = comp[1:]
    nolp += 1
    stack.append({})
  elif rm:
    head = rm.group()
    tail = comp[len(head):]
    num = 1
    if len(head) > 1:
      num = int(head[1:])
    nolp -= 1
    if nolp < 0:
      print('Error - unmatched right parenthesis in formula: {:}'.format(comp))
    for (k,v) in stack.pop().items():
      if k not in stack[-1]:
        # First time adding this element to dict.
        stack[-1][k] = num*v
      else:
        # Increment number for this element.
        stack[-1][k] += num*v
  else:
    print('Error - {:} is an invalid chemical formula.')
    return
  # Terminate or continue recursion.
  if len(tail) > 0:
    stack = rpc(tail,stack,nolp)
    return stack
  else:
    if nolp > 0:
      print('Error - unmatched left parenthesis in formula: {:}'.format(comp))
      return
    return stack

# Method for getting molecular weight.
def molecular_weight(comp,sram_lib):
  '''
  Computes the molecular weight for the specified chemical compound, using
  the specified library of standard relative atomic masses. The sram_lib
  should be imported with the load_atomic_weights method, and modified as
  needed with the update_sram_lib method.
  '''
  # Parse the compound.
  atoms = parse_compound(comp)
  # Loop over constituent atoms, summing up the molecular weight.
  mw = 0.0
  for (k,v) in atoms.items():
    if k in sram_lib:
      if sram_lib[k]['SRAM']['Type'] in ['quantity','most_stable']:
        mw += sram_lib[k]['SRAM']['value']*v
      else:
        print('Standard relative atomic mass for {:} is unknown or published as an interval.'.format(k))
    else:
      print('Standard relative atomic mass library missing element: {:}.'.format(k))
      return
  return mw
