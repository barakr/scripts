#!/usr/bin/python

import prepare
import Bio.PDB
import Bio.PDB.Polypeptide
import numpy
import itertools
from optparse import OptionParser
import re
import sys
import math
import os
import subprocess
import errno
from random import choice
sys.path.append(os.path.dirname(__file__)) # If next imports won't work
import get_interface
import transform_peptide

################# MAIN ############
'''
For each subfolder in current directory, prepare a run folder:
1) look for native.pdb, and for numerical subfolders with prepacked start structures,
2) put peptide in start structures away from the receptor interface
3) add flag files, fragment files, etc. - all that is needed for a run
'''
#    options, args = get_cmdline_options()
for dirname, pdb_dirnames, filenames in os.walk('.'):
    break
for pdb in pdb_dirnames:
    print pdb
    prepare.process_one_pdb( pdb )
