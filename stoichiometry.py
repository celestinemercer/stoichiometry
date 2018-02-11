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
import os
import stoich_calc

# Define global variables.
def_prefs_path = 'resources/stoichiometry.prefs'
prefs = {}
def_sram_nist_path = 'resources/AtomicWeights_IsotopicCompositions_NIST_4.1.txt'
def_sram_patch_path = 'resources/SelectedGeologicAtomicWeights.txt'
def_minsys_dir = 'resources/mineral_systems'
sram_lib = {}
datasets = {}
results = {}
min_systems = {}

# start function.
def start(prefs_path=def_prefs_path,sram_nist=def_sram_nist_path,
          sram_patch=def_sram_patch_path,minsys_dir=def_minsys_dir):
  '''
  The start of the main program loop for performing stoichiometry calculations.
  '''
  # Print welcome message.
  print('-'*80); print('/\\'*40); print('\\/'*40); print('-'*80)
  print('\nWelcome to Stoichiometry Calculator 5000!\n')
  # Perform startup tasks; get access to global variables.
  global prefs, datasets, results
  run_startup_tasks(prefs_path,sram_nist,sram_patch,minsys_dir)
  # Define menu options, prepare to loop.
  opts = ['Import dataset', 'Remove dataset', 'Export results', 'Manual', 'Auto', 'Edit preferences', 'Exit']
  keep_going = True
  while keep_going == True:
    # Present user with main menu options.
    choice = util.prompt_options('Main Menu options', opts)
    # Perform requested action.
    if choice == 0:
      # Load dataset. Oxides must be in mixed case format for now.
      filename, dataset = util.import_dataset(prefs['wdir'],prefs['delimiter'])
      if filename is not None:
        print('Dataset imported: {:}\n'.format(filename))
        datasets[filename] = dataset
      else:
        print('Import dataset aborted.\n')
    elif choice == 1:
      # Remove dataset.
      print('Removal functionality coming soon...\n')
    elif choice == 2:
      # Export results.
      results_list = list(results.keys())
      selection = util.prompt_options('Specify results to export', results_list)
      res_name = results_list[selection]
      target = os.path.join(prefs['wdir'], res_name + '.stoichres')
      ok2save = False
      if os.path.isfile(target):
        print('Warning: there is already a file named {:} in the working directory.'.format(res_name + '.stoichres'))
        yesno = input('Would you like to overwrite the file? (y/n) >> ')
        if yesno == 'y':
          ok2save = True
      else:
        ok2save = True
      if ok2save:
        results[res_name].to_csv(target,prefs['delimiter'])
    elif choice == 3:
      # Manual calculation options.
      print('Manual stoichiometry calculations:')
      manual_calc()
    elif choice == 4:
      # Auto calculation options.
      print('Auto calculation options coming soon...\n')
    elif choice == 5:
      # Edit preferences.
      edit_prefs_main(prefs_path)
    elif choice == 6:
      # Exit program.
      keep_going = False
  # Do any cleanup tasks here, after while loop terminates.
  run_shutdown_tasks(prefs_path)
  print('Thank you, come again soon!')
  print('-'*80); print('/\\'*40); print('\\/'*40); print('-'*80)
  
# manual_calc function. 
def manual_calc():
  # Get access to global variables.
  global prefs, datasets, results, min_systems
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
  manual_opts = ['Anhydrous silicates', 'Hydrous silicates','Non-silicates', 'Return to: Main Menu']
  keep_going = True
  while keep_going == True:
    # Show manual calculation options menu.
    choice = util.prompt_options('Manual Calculation options for: {:}'.format(activeKey), manual_opts)
    # Perform requested action.
    if choice == 0:
      # Anhydrous silicates.
      results[activeKey] = stoich_calc.anhydrous_silicates_stoich(active, sram_lib, min_systems)
    elif choice == 1:
      # Hydrous silicates.
      print('Hydrous silicates options coming very soon...\n')
#      results[activeKey] = stoich_calc.hydrous_silicates_stoich(active, sram_lib, min_systems)
    elif choice == 2:
      # Non-silicates.
      print('Non-silicates options coming soon...\n')
    elif choice == 3:
      # Return to Main Menu.
      keep_going = False

# auto_calc function.



# run_startup_tasks function.
def run_startup_tasks(prefs_path,nist,patch,minsys_dir):
  # Get access to global variables; clear them on startup.
  global prefs, sram_lib, datasets, min_systems
  prefs, sram_lib, datasets, min_systems = {}, {}, {}, {}
  # Initialize/load preferences.
  prefs = init_prefs(prefs_path)
  print('Preferences loaded.')
  # Load and patch library of standard relative atomic weights (SRAMs).
  sram_lib = util.load_atomic_weights(nist)
  patch = util.load_atomic_weights(patch)
  sram_lib = util.update_sram_lib(sram_lib,patch)
  print('Standard relative atomic mass library loaded.')
  # Load mineral systems.
  min_systems = util.load_mineral_systems(minsys_dir)
  if min_systems is not None and len(min_systems) > 0:
    print('Mineral systems loaded (N = {:}).'.format(len(min_systems)))
  else:
    print('Warning: no mineral systems were loaded.')
  print('')

# run_shutdown_tasks function.
def run_shutdown_tasks(prefs_path):
  # Get access to global variables.
  global prefs
  # Autosave preferences if possible.
  if prefs['autosave_prefs']:
    util.write_dict(prefs,os.path.abspath(prefs_path))
    print('Preferences autosaved.')
  print('')

# init_prefs function.
def init_prefs(prefs_path):
  '''
  Initializes the default preferences for the stoichiometry program, and
  loads any preferences that have been saved to the specified file.

  Returns a dictionary of preferences.
  '''
  # Define default preferences.
  dprefs = {'wdir':os.path.abspath('.'),
            'autosave_prefs':True,
            'delimiter':'\t'}
  # Check for preferences file; load it if it exists.
  lprefs = {}
  if os.path.isfile(prefs_path):
    lprefs = util.load_dict(prefs_path)
  # Supply default preferences that are missing.
  for key in dprefs.keys():
    if key not in lprefs:
      lprefs[key] = dprefs[key]
  # Return preferences.
  return lprefs

# edit_prefs_main function.
def edit_prefs_main(prefs_path):
  '''
  Displays the main preference editing menu to the user.
  '''
  # Prepare to show main preference editing menu.
  main_edit_opts = ['General preferences','Manual preferences','Auto preferences','Save preferences',
                    'Return to: Main Menu']
  keep_going = True
  while keep_going == True:
    # Show main preference editing menu.
    choice = util.prompt_options('Edit Preferences - Main Menu',main_edit_opts)
    # Handle user selection.
    if choice == 0:
      # General prefs.
      edit_prefs_general()
    elif choice == 1:
      # Manual stoichiometry prefs.
      print('Manual stoichiometry preferences coming later...')
    elif choice == 2:
      # Auto stoichiometry prefs.
      print('Auto stoichiometry preferences coming later...')
    elif choice == 3:
      # Save preferences manually.
      util.write_dict(prefs,os.path.abspath(prefs_path))
      print('Preferences saved.\n')
    elif choice == 4:
      # Return to Main Menu.
      keep_going = False
  print('Rock on! Returning to main menu.\n')

# edit_prefs_general function.
def edit_prefs_general():
  '''
  Allows the user to edit general preferences for the stoichiometry program.
  '''
  # Get access to global variables.
  global prefs
  # Prepare to show menu.
  dels = {'\t':'tab',',':'comma'}
  keep_going = True
  while keep_going == True:
    # Show main preference editing menu.
    general_edit_opts = ['Autosave preferences = {:} (toggle value)'.format(prefs['autosave_prefs']),
                         'Working directory = {:} (edit)'.format(prefs['wdir']),
                         'File delimiter = {:} (choose)'.format(dels[prefs['delimiter']]),
                         'Return to: Edit Preferences - Main Menu']
    choice = util.prompt_options('Edit General Preferences',general_edit_opts)
    # Handle user selection.
    if choice == 0:
      # Invert boolean value.
      prefs['autosave_prefs'] = not prefs['autosave_prefs']
    elif choice == 1:
      # Track line count.
      lc = 0
      # Prompt for new wdir path.
      new_wdir = None
      invalid = True
      while invalid:
        print('Current working directory: {:}'.format(prefs['wdir']))
        new_wdir = input('\nEnter new working directory: ')
        lc += 3
        if os.path.isdir(new_wdir):
          invalid = False
          util.clear_stdout_lines(lc)
        else:
          print('Invalid directory path; please try again.\n')
          lc += 2
      # Save new wdir.
      prefs['wdir'] = new_wdir
    elif choice == 2:
      # Choose file delimiter.
      idx = util.prompt_options('Select file delimiter',list(dels.values()))
      prefs['delimiter'] = list(dels.keys())[idx]
    elif choice == 3:
      # Return to Main Menu.
      keep_going = False
  print('Rock on! Returning to main preferences editor menu.\n')
