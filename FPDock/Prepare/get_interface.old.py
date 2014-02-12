# Copyright 2007 Peter Cock, all rights reserved.
# Licenced under the GPL v2 (or later at your choice)
#
# Please see this website for details:
# http://www.warwick.ac.uk/go/peter_cock/python/protein_superposition/
import Bio.PDB
import numpy
import sys

def print_string_to_file(filename, string):
    FILE = open(filename,'w')
    print >>FILE, string
    FILE.close()

def get_atoms(chain):
    ret=[]
    for res in chain:
        for atom in res:
            ret.append(atom)
    return ret

def exit_usage():
   print "Usage: %s <complex pdb>" % sys.argv[0]
   print "Example: %s native.pdb" % sys.argv[0]
   exit(-1)

def print_rosetta_sites(pdb_fname, out_filehandle, receptor_chain_id="A", peptide_chain_id="B", interface_cutoff=5):
    '''
    print rosetta site constraints to file handle out_filehandle based
    on interface of pdb_fname complex between specified receptor and peptide
    chains, with specified interface cutoff
    '''
    struct = Bio.PDB.PDBParser().get_structure("complex", pdb_fname)
    receptor_chain=struct[0][receptor_chain_id]
    peptide_chain=struct[0][peptide_chain_id]
    receptor_atoms = get_atoms(receptor_chain)
    peptide_atoms = get_atoms(peptide_chain)
    ns = Bio.PDB.NeighborSearch(receptor_atoms)
    interface_set = set()
    for atom in peptide_atoms:
        neighbors = ns.search(atom.get_coord(),radius=interface_cutoff,level="R")
        interface_set.update(neighbors)
    for res in sorted(interface_set):
        res_id=res.id[1]
        ch_id=res.parent.id
        # penalty beyond 0-6 range, with std-dev 2.0
        atom_name="CB"
        if(res.resname=="GLY"):
            atom_name="CA"
        print >>out_filehandle, \
            "SiteConstraint %s %d%s %s BOUNDED 0 6 2.0 0.5 intrf" \
            % (atom_name, res_id, ch_id, peptide_chain_id)



######## MAIN #########
if (__name__ == "__main__"):
    if(len(sys.argv) <> 2):
        exit_usage()
    pdb_filename = sys.argv[1]
    receptor_chain_id="A"
    peptide_chain_id="B"
    interface_cutoff = 5 # in Angstroms
    print_rosetta_sites(pdb_filename, sys.stdout, receptor_chain_id, peptide_chain_id, interface_cutoff)
