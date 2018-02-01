# mineral_system.py
# Class definition file for dealing with mineral systems in the stoichiometry program.

import util

# Class definition.
class MineralSystem(object):
  '''
  This class contains information specific to a particular mineral system for
  use in stoichiometry calculations.
  '''
  # Instance variables (fields) to add in the future: 
  # - list of validity check parameters (possibly for each oxide; enhancement for later?)
  # - possible enhancement: store lists of charges for elements in each oxide... e.g., Fe+2 and Fe+3, O-2, etc.
  
  # Constructor method.
  def __init__(self,filepath):
    # Set instance variables to None, for safety checks later.
    self.system_name = None
    self.hydrous = None
    self.silicate = None
    self.oxides = None
    self.anions = None
    self.endmembers = None
    # Read definition file, parse the contents, and overwrite instance variables.
    with open(filepath,'r') as f:
      lines = f.readlines()
    for line in lines:
      # Split line and trim away whitespace.
      items = line.split('=')
      for i in range(len(items)):
        items[i] = items[i].strip()
      # Store information to instance variables.
      if items[0] == 'System Name':
        self.system_name = items[1]
      elif items[0] == 'Hydrous':
        self.hydrous = items[1].lower() == 'true'
      elif items[0] == 'Silicate':
        self.silicate = items[1].lower() == 'true'
      elif items[0] == 'Oxides':
        temp_oxides = items[1].split(',')
        for i in range(len(temp_oxides)):
          temp_oxides[i] = temp_oxides[i].strip()
        self.oxides = temp_oxides
      elif items[0] == 'Anions':
        temp_anions = items[1].split(',')
        for i in range(len(temp_anions)):
          temp_anions[i] = temp_anions[i].strip()
        self.anions = temp_anions
      elif items[0] == 'Endmembers':
        temp_endmembers = items[1].split(',')
        for i in range(len(temp_endmembers)):
          temp_endmembers[i] = temp_endmembers[i].strip()
        self.endmembers = temp_endmembers
    # Check that all instance variables have been initialized.
    if None in [self.system_name,self.hydrous,self.silicate,self.oxides,self.anions,self.endmembers]:
      raise RuntimeException('Failed to load mineral system from: {:}'.format(filepath))
  
  # Instance methods.
  
  # getSystemName method.
  def getSystemName(self):
    return self.system_name
  
  # isHydrous method.
  def isHydrous(self):
    return self.hydrous
  
  # isSilicate method.
  def isSilicate(self):
    return self.silicate
  
  # getOxideList method.
  def getOxideList(self):
    return self.oxides
  
  # getAnions method.
  def getAnionList(self):
    return self.anions
  
  # getAnionBasis method.
  def getAnionBasis(self):
    # Count anions present in first endmember. Assume that all endmembers have the same basis for now.
    elems = util.parse_compound(self.endmembers[0])
    anion_count = 0
    for ai in self.anions:
      if ai in elems.keys():
        anion_count += elems[ai]
    return anion_count
  
  # getEndmembers method.
  def getEndmemberList(self):
    return self.endmembers
  
  
  