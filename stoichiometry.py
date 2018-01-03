# Stoichometry                                                                                                                                           
# stoichiometry.py                                                                                                                                                
#                                                                                                                                                        
# This file contains funcitons to enable interactive stoichometry calculations.                                                                                      
#                                                                                                                                                        
# Author(s):
# - Celeste N. Mercer
# - Cameron M. Mercer                                                                                                                               
#                                                                                                                                                        
# Copyright (C), 2018                                                                                                                                   
#

import util
import pandas as pd

datasets = {}

# start function.
def start():
  '''
  The start of the main program loop for performing stoichiometry.
  '''
  # Define menu options, prepare to loop.
  opts = ['import dataset', 'list datasets', 'exit']
  keep_going = True
  while keep_going == True:
    # Present user with menu options.
    choice = util.prompt_options('choose your own adventure', opts)
    # Perform requested action.
    if choice == 0:
      # Load dataset.
      filename, dataset = util.import_dataset()
      datasets[filename] = dataset
    elif choice == 1:
      print(datasets.keys())
    elif choice == 2:
      # Exit program.
      keep_going = False
  # Do any cleanup tasks here, after while loop terminates.
  
  

  
  
  

