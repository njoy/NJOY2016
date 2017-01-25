#!/usr/bin/env python

import difflib
import filecmp
import fileinput
import glob
import os
import re
import subprocess

retained_tapes = set( glob.glob('tape*') )
reference_tapes = glob.glob('referenceTape*')

with open('input', 'r') as i, \
     open('output', 'w') as o, \
     open('error', 'w') as e :
    njoy = os.path.join('..', '..', 'njoy')
    child = subprocess.Popen(njoy, stdin=i, stdout=o, stderr=e)
    child.communicate()
    if ( child.poll() ):
        print("Error enountered while running NJOY!")
        exit( child.poll() )

    for reference_tape in reference_tapes:
        trial_tape = 'tape' + reference_tape[-2:]
        if not filecmp.cmp(reference_tape, trial_tape):
            with open(reference_tape, 'r') as reference_file, \
                 open(trial_tape, 'r') as trial_file, \
                 open(trial_tape + '_diff', 'w') as diff_file :
                should_exit = False
                reference_lines = reference_file.readlines()
                trial_lines = trial_file.readlines()
                reference_lines = [re.sub(r'\d{2}/\d{2}/\d{2}',
                                          r'XX/XX/XX', line)
                                   for line in reference_lines]
                trial_lines = [re.sub(r'\d{2}/\d{2}/\d{2}',
                                      r'XX/XX/XX', line)
                               for line in trial_lines]

                if reference_lines != trial_lines:
                    should_exit = True
                    for line in difflib.context_diff(reference_lines,
                                                     trial_lines,
                                                     fromfile=reference_tape,
                                                     tofile=trial_tape, n=0):
                        diff_file.write(line)
                      
            if should_exit: exit(99)

    removed_tapes = list( set( glob.glob('tape*') ) - retained_tapes )
    for tape in removed_tapes:
        os.remove( tape )
        
    for diff in glob.glob('*_diff'):
        os.remove( diff )
        
    os.remove('output')
    os.remove('error')
