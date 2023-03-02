#!/usr/bin/env python

import filecmp
import glob
import os
import re
import subprocess
from math import isclose


floatPattern = re.compile("""(
                                (
                                    [-+]?				# sign (optional)
                                    \d+					# digit(s)
                                )
                                (
                                    \.              # period
                                    \d+  		    # digit(s)
                                )?                  # period+digits  (optional)
                                (
                                    ([eE])?			# designation (optional)
                                    [-+]			# sign (required)
                                    \d+				# digit(s) (required)
                                )? 					# exponent (optional)
                            )
                          """, re.VERBOSE)
datePattern = re.compile(r'\d{2}/\d{2}/\d{2}')


def lineEquivalence(ref, trial, n, relativeError=1E-9, absoluteError=1E-10):
    """
    Look to see if two lines (ref and trial) are identical. Returns a boolean.
    """
    equivalent = True
    if ref != trial:

        # Look for numbers
        refFloats = floatPattern.findall(ref)
        if refFloats:
            floatEquivalence = True
            trialFloats = floatPattern.findall(trial)

            if len(refFloats) != len(trialFloats):
                print("Found wrong number of floats on line: {}".format(n))
                floatEquivalence = False
            else:
                # Look at all the floats on the lines
                for rF, tF in zip(refFloats, trialFloats):
                    equal = isclose( makeFloat(rF), makeFloat(tF),
                                     rel_tol = relativeError,
                                     abs_tol = absoluteError )

                    if not equal:
                        print("{} and {} are not equal".format(
                            rF[0], tF[0]))
                        floatEquivalence = False

            if not floatEquivalence:
                equivalent = False
            else:
                refNoNumbers = floatPattern.sub("", ref)
                trialNoNumbers = floatPattern.sub("", trial)
                if refNoNumbers != trialNoNumbers:
                    equivalent = False

        # Lines contain only text (no numbers)
        else:
            equivalent = False

    return equivalent


def identicalLines(refLines, trialLines, diffFile,
                   relativeError=1E-9, absoluteError=1E-10):
    """
    Look at two lists of strings (presumably lines from two files) and compare
    them; that is, determine if they are the same within the given tolerance.
    Returns True if everything is the same, returns False if anythin differs.
    """
    if len(refLines) != len(trialLines):
        print("Reference file and trial file have different number of lines!")
        return False

    fullEquivalence = True
    for n, (ref, trial) in enumerate(zip(refLines, trialLines)):
        equivalent = lineEquivalence(ref, trial, n,
                                     relativeError, absoluteError)
        if not equivalent:
            fullEquivalence = False
            writeDiff(diffFile, ref, trial, n)

    return fullEquivalence


def makeFloat(D):
    """
    """
    try:
        return float(D[0])
    except (TypeError, ValueError):
        if not D[-1]:
            d = "{}{}E{}".format(D[1], D[2], D[3])
            f = float(d)
            return f


def writeDiff(diff_file, refLine, trialLine, lineNumber):
    """
    Write the differences to the diff file.
    """
    diff_file.write("***************\n")
    diff_file.write("*** {} ***\n".format(lineNumber+1))
    diff_file.write("!{}".format(refLine))
    diff_file.write("--- {} ---\n".format(lineNumber+1))
    diff_file.write("!{}".format(trialLine))

# Diff checking for original format files
#retained_tapes = set(glob.glob('SNL-NJOY-2016_Test_Tape*'))
reference_tapes = glob.glob('SNL-NJOY-2016_referenceTape*')

for reference_tape in reference_tapes:
  trial_tape = 'SNL-NJOY-2016_Test_Tape' + reference_tape[-2:]
  if not filecmp.cmp(reference_tape, trial_tape):
    with open(reference_tape, 'r') as reference_file, \
         open(trial_tape, 'r') as trial_file, \
         open(trial_tape + '_diff', 'w') as diff_file:
      should_exit = False
      reference_lines = reference_file.readlines()
      trial_lines = trial_file.readlines()
      reference_lines = [datePattern.sub(r'XX/XX/XX', line)
                         for line in reference_lines]
      trial_lines = [datePattern.sub(r'XX/XX/XX', line)
                     for line in trial_lines]

      diff_file.write("*** {} ***\n".format(reference_tape))
      diff_file.write("--- {} ---\n".format(trial_tape))

      identical = identicalLines(reference_lines, trial_lines, diff_file, 1E-5, 1E-5)

    if not identical:
      print("Diff found in NJOY test, exiting")
      exit(99)
      

# Diff checking for SNL-MANIPULATE series
#retained_tapes = set(glob.glob('SNL-MANIP_Test_Tape*'))
reference_tapes = glob.glob('SNL-MANIP_referenceTape*')

for reference_tape in reference_tapes:
  trial_tape = 'SNL-MANIP_Test_Tape' + reference_tape[-2:]
  if not filecmp.cmp(reference_tape, trial_tape):
    with open(reference_tape, 'r') as reference_file, \
         open(trial_tape, 'r') as trial_file, \
         open(trial_tape + '_diff', 'w') as diff_file:
      should_exit = False
      reference_lines = reference_file.readlines()
      trial_lines = trial_file.readlines()
      reference_lines = [datePattern.sub(r'XX/XX/XX', line)
                         for line in reference_lines]
      trial_lines = [datePattern.sub(r'XX/XX/XX', line)
                     for line in trial_lines]

      diff_file.write("*** {} ***\n".format(reference_tape))
      diff_file.write("--- {} ---\n".format(trial_tape))

      identical = identicalLines(reference_lines, trial_lines, diff_file, 1E-5, 1E-5)


    if not identical:
      print("Diff found in MANIPULATE-NJOY tests, exiting")
      exit(99)

print("No diff found")

#removed_tapes = list(set(glob.glob('tape*')) - retained_tapes)
#for tape in removed_tapes:
#  os.remove(tape)

for diff in glob.glob('*_diff'):
  os.remove(diff)
