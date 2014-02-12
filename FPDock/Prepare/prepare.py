#!/usr/bin/python

import Bio.PDB
import Bio.PDB.Polypeptide
from Bio import SeqIO
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

ROSETTA=os.environ['ROSETTA']

flags_basic= '''
# io flags:
-s Run/start.pdb
-native Run/native.pdb
-out:file:silent_struct_type binary

# optimization flags:
-pep_refine
-ex1
-ex2aro
#-use_input_sc
#-unboundrot Run/start.pdb
-flexPepDocking:flexpep_score_only

#mute logging:
-mute protocols.moves.RigidBodyMover
-mute core.chemical
-mute core.scoring.etable
-mute protocols.evalution
-mute core.pack.rotamer_trials
-mute protocols.abinitio.FragmentMover
-mute core.fragment
-mute protocols.jd2.PDBJobInputter
'''

# flags_refine = flags_basic + '''
# -nstruct 1000
# '''

# flags_refine_prelow = flags_refine + '''
# -lowres_preoptimize
# '''

site_constraints='''

#constraint:
-constraints:cst_file Run/site.cst
-constraints.cst_fa_file Run/site.cst
-constraints:cst_fa_weight 1.0
-constraints:cst_weight 1.0
'''

flags_abinitio35= flags_basic + site_constraints + '''
-nstruct 50000

# optimization
-flexPepDocking:lowres_abinitio

#fragment picker flags:
-frag3 Run/Frags/frags.3mers.offset
-flexPepDocking:frag5 Run/Frags/frags.5mers.offset
-flexPepDocking:frag5_weight 0.25
'''

flags_abinitio359 = flags_abinitio35 + '''
-frag9 Run/Frags/frags.9mers.offset
-flexPepDocking:frag9_weight 0.1
'''

def get_pep_length(pdbfile):
    model = Bio.PDB.PDBParser().get_structure(pdbfile, pdbfile)[0]
    pep = model["B"]
    return sum(1 for r in pep.get_residues())

def switch_chain(pdbfile,chain_id="A"):
    structure = Bio.PDB.PDBParser().get_structure(pdbfile, pdbfile)
    chain = structure.get_chains().next()
    chain.id=chain_id
    io=Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdbfile)

def print_string_to_file(filename, string):
    FILE = open(filename,'w')
    print >>FILE, string
    FILE.close()

def soft_link_if_needed(filename, linkname):
    print filename
    print linkname
    if(not os.path.lexists(linkname)):
        os.symlink(filename, linkname)
    if(not os.path.exists(filename)):
        print "WARING: broken soft link %s to non-existing file %s" \
            % (filename, linkname)

def force_symlink(filename, linkname):
    try:
        os.symlink(filename, linkname)
    except OSError, e:
        if e.errno == errno.EEXIST:
            os.remove(linkname)
            os.symlink(filename, linkname)

def is_non_zero_file(fpath):
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 \
        else False

def choose_intable(L):
    '''
    choose a random entry from all
    entries in the list L that convert to int
    '''
    n_trials = 20
    for i in range(n_trials):
        c=choice(L)
        try:
            int(c)
            return c
        except ValueError:
            pass
    raise Exception('No valid integet found in list %s' % L)

def get_nres(pdb_fname,chain_id=0):
    ''' get nres of chain #chain_id (default, first chain)
    '''
    struct = Bio.PDB.PDBParser().get_structure("complex", pdb_fname)
    chain=struct.get_chains().next()
    return len(list(chain.get_residues()))

def process_one_pdb(pdb):
    '''
    Enter folder <pdb> that is assumed to have a native.pdb
    file, and anchor folders with <pdb>.ex.pdb start files.
    Create fragment files and start files for the complex
    for a random anchor, in folders Frags and Run resp.
    '''
    cwd=os.getcwd() # save to restore in the end
    try:
        os.chdir(pdb)
        if(not os.path.exists("native.pdb")):
            raise IOError("native.pdb not found in folder %s" % pdb)
        if not os.path.exists("Frags"):
            os.mkdir("Frags")
        if not is_non_zero_file("Frags/frags.3mers.offset"):
            os.chdir("Frags")
            subprocess.call(["../../Scripts/Prepare/make_frags.sh"])
            os.chdir("..")
            assert is_non_zero_file("Frags/frags.3mers.offset")
        # save  path to all subdirectories in dirnames
        # and choose random anchor dir
        for dirname, dirnames, filenames in os.walk('.'):
            break
        anchor_dir = choose_intable(dirnames)
        # prepare run flags files
        if(not os.path.exists("Run")):
            os.mkdir("Run")
        os.chdir("Run")
        force_symlink("../%s/%s.ex.pdb" % (anchor_dir, pdb), "start_raw.pdb")
        transform_peptide.transform_peptide("start_raw.pdb", "start.pdb") # move peptide away from interface
        soft_link_if_needed("../Frags", "Frags")
        soft_link_if_needed("../native.pdb", "native.pdb")
        soft_link_if_needed("../../Scripts/", "Scripts")
        if(get_pep_length("native.pdb") < 9):
            print_string_to_file("flags_abinitio", flags_abinitio35)
        else:
            print_string_to_file("flags_abinitio", flags_abinitio359)
        # Extract interface constraints from native
        # REPLACE WITH USER INTERFACE SITES CONSTRAINTS IN REAL SETTINGS
        with open("site.cst", "w") as sites_file_handle:
            get_interface.print_rosetta_sites("native.pdb", sites_file_handle)
        print "Finished with %s" % pdb
    except (Exception,ValueError,NameError,OSError) as e:
        print "Skipping %s" % pdb,
        print "Exception: ", e
    os.chdir(cwd) # restore initial folder


def create_fasta_from_coords(pdb, fasta_file, chain='B'):
    ''' create fasta file fasta_file from specified chain in pdb'''
    command = "%s/rosetta_tools/perl_tools/getFastaFromCoords.pl " + \
        " -pdbfile %s -chain %c > %s" % (ROSETTA, pdb, chain, fasta_file)
    subprocess.call(command, shell=True)

def create_fasta_from_seq(seq, fasta_file, header="pep"):
    FASTA=open(fasta_file,"w")
    assert( not "\n" in header )
    print >>FASTA, ">%s" % header
    print >>FASTA, seq

def make_frags(pep_fasta, receptor_pdb):
    ''' make frags for peptide sequence in pep_fasta with receptor pdb receptor_pdb '''
    if not os.path.exists("Frags"):
        os.mkdir("Frags")
    if is_non_zero_file("Frags/frags.3mers.offset"):
        print "NOTE: Frags files already prepared (at least 3mers - delete if want to redo)"
    os.chdir("Frags")
    if not os.path.isabs(receptor_pdb):
        receptor_pdb = "../" + receptor_pdb
    if not os.path.isabs(pep_fasta):
        pep_fasta = "../" + pep_fasta
    print "Receptor", receptor_pdb
    receptor_nres = get_nres(receptor_pdb, 0)
    pep_record = SeqIO.read(open(pep_fasta), "fasta")
    pep_nres = len(pep_record.seq)
    print "RECEPTOR LENGTH:", receptor_nres
    print "PEPTIDE LENGTH:", pep_nres
    print "ROSETTA", ROSETTA
    if pep_nres >= 9:
        frag_sizes="3,5,9"
    else:
        frag_sizes="3,5"
    command= ROSETTA + "/rosetta_tools/fragment_tools/make_fragments.pl" + \
        " -frag_sizes %s -n_frags 500 -id fragB %s" % (frag_sizes, pep_fasta)
    subprocess.call(command, shell=True)
    # shift frags by receptor length
    for fs in frag_sizes.split(","):
        command = "awk '{if ( substr ( $0,1,3 ) == \"pos\" ) " + \
            " {print substr ( $0,0,18 ) sprintf (\"%4d\",substr ( $0,19,4 ) + " + \
            str(receptor_nres) + " ) substr ( $0,23,1000 ) ; } else {print ; }}' " + \
            " fragB.500." + fs + "mers > frags." + fs + "mers.offset"
        subprocess.call(command, shell=True)
    os.chdir("..")


def prepare_from_receptor_pdb_and_pep_seq(receptor_pdb,
                                          pep_seq,
                                          native_pdb=None):
    '''
    prepare a run folder with frag lib, etc. from a receptor pdb file
    and a fasta file with the peptide sequence (receptor_pdb and pep_fasta,
    resp.)

    native_pdb - if specified, pdb with native complex
    '''
    cwd=os.getcwd() # save to restore in the end
    if native_pdb is None:
        native_pdb="start.pdb"
    else:
       if not os.path.exists(native_pdb):
           raise IOError("%s not found" % native_pdb)
    switch_chain(receptor_pdb, "A") # TODO: work on copy!!
    pep_fasta = "pep.fasta"
    pep_pdb = "pep.pdb"
    create_fasta_from_seq(pep_seq, pep_fasta)
    # prepare run flags files
    if(not os.path.exists("Input")):
        os.mkdir("Input")
    os.chdir("Input")
    command = ROSETTA + "/rosetta_source/bin/BuildPeptide.default.linuxgccrelease -database " + \
        ROSETTA + "/rosetta_database/ -in:file:fasta ../" + pep_fasta + \
        " -out:file:o " + pep_pdb
    subprocess.call(command, shell=True)
    switch_chain(pep_pdb, "B")
    subprocess.call("cat ../%s %s > start_raw.pdb" % (receptor_pdb, pep_pdb),
                    shell=True)
    transform_peptide.transform_peptide("start_raw.pdb", "start.pdb") # move peptide away from interface
    soft_link_if_needed("../Frags", "Frags")
    soft_link_if_needed(native_pdb, "native.pdb")
    soft_link_if_needed("../../Scripts/", "Scripts")
    if(get_pep_length("start.pdb") < 9):
        print_string_to_file("flags_abinitio", flags_abinitio35)
    else:
        print_string_to_file("flags_abinitio", flags_abinitio359)
    os.chdir(cwd)
    # Prepare frags
    make_frags(pep_fasta, receptor_pdb)
    assert is_non_zero_file("Frags/frags.3mers.offset")
    # Finito
    os.chdir(cwd) # restore initial folder
    print "Finished with %s" % receptor_pdb


prepare_from_receptor_pdb_and_pep_seq(sys.argv[1], sys.argv[2])
