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
  # Print welcome message.
  print('-'*80)
  print('Welcome to Stoichiometry Calculator 5000!')
  print('')
  # Define menu options, prepare to loop.
  opts = ['Import dataset', 'Remove dataset', 'Manual', 'Auto', 'Preferences', 'Exit']
  keep_going = True
  while keep_going == True:
    # Present user with main menu options.
    choice = util.prompt_options('Main Menu options', opts)
    # Perform requested action.
    if choice == 0:
      # Load dataset.
      filename, dataset = util.import_dataset()
      datasets[filename] = dataset
    elif choice == 1:
      # Remove dataset.
      print('Removal functionality coming soon... Here are the datasets you have loaded:')
    elif choice == 2:
      # Manual calculation options.
      manual_calc()
    elif choice == 3:
      # Auto calculation options.
      print('Auto calculation options coming soon...')
    elif choice == 4:
      # Edit preferences.
      print('Editing preferences will come later...')
    elif choice == 5:
      # Exit program.
      keep_going = False
  # Do any cleanup tasks here, after while loop terminates.
  print('Thank you, come again soon!')
  print('-'*80)
  
# manual_calc function. 
def manual_calc():
  # Select active dataset.
  activeKey = None
  if len(datasets) == 0:
    print('There are no datasets loaded, please import one.\n')
    return
  elif len(datasets) == 1:
    # Use the dataset that is loaded.
    activeKey = list(datasets.keys())[0]
  else:
    # Prompt user to select dataset.
    idx = util.prompt_options('Please select a dataset to process',list(datasets.keys()))
    activeKey = list(datasets.keys())[idx]
  active = datasets[activeKey]
  # Prepare to show manual calculation options.
  manual_opts = ['Anhydrous silicates', 'Hydrous silicates','Non-silicates', 'Return to Main Menu']
  keep_going = True
  while keep_going == True:
    # Show manual calculation options menu.
    choice = util.prompt_options('Manual Calculation options for: {:}'.format(activeKey), manual_opts)
    # Perform requested action.
    if choice == 0:
      # Anhydrous silicates.
      print('Anhydrous silicates options coming soon...\n')
    elif choice == 1:
      # Hydrous silicates.
      print('Hydrous silicates options coming soon...\n')
    elif choice == 2:
      # Non-silicates.
      print('Non-silicates options coming soon...\n')
    elif choice == 3:
      # Return to Main Menu.
      keep_going = False
  print('Rock on! Returning to main menu.\n')