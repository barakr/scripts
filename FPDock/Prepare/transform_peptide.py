import Bio.PDB
from Bio.PDB import Vector
import numpy
import sys
import math

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

def get_center_of_mass(chain):
    ''' get chain center of mass as x,y,z tuple '''
    coords=[]
    for res in chain:
        try:
            coords.append( res['CA'].get_coord() )
        except:
            pass
        try:
            coords.append( res['CB'].get_coord() )
        except:
            pass
    sum_x = sum(c[0] for c in coords)
    sum_y = sum(c[1] for c in coords)
    sum_z = sum(c[2] for c in coords)
    s = [sum_x,sum_y,sum_z]
    center_of_mass = [a/len(coords) for a in s]
    return center_of_mass

def exit_usage():
   print "Usage: %s <complex pdb>" % sys.argv[0]
   print "Example: %s native.pdb" % sys.argv[0]
   exit(-1)

def transform_peptide(pdb_fname, outpdb_fname,
                      receptor_chain_id="A", peptide_chain_id="B"):
    '''
    transform peptide chain in pdb_fname away from the receptor center of mass,
    and add a 90 degrees rotation. Save to outpdb_fname.
    '''
    struct = Bio.PDB.PDBParser().get_structure("complex", pdb_fname)
    receptor_chain=struct[0][receptor_chain_id]
    peptide_chain=struct[0][peptide_chain_id]
    receptor_center_of_mass = get_center_of_mass(receptor_chain)
    peptide_center_of_mass = get_center_of_mass(peptide_chain)
    rot_id=Bio.PDB.rotmat(Vector([1,0,0]),Vector([1,0,0]))
    trans_to_zero=numpy.array([-x for x in peptide_center_of_mass],'f')
    peptide_chain.transform(rot_id,trans_to_zero)
    trans_vec=[2*a-b for a,b in zip (peptide_center_of_mass, receptor_center_of_mass)]
    rotation=Bio.PDB.rotaxis(math.pi, Vector(trans_vec))
    translation = numpy.array(trans_vec,'f')
    peptide_chain.transform(rotation,translation)
    io=Bio.PDB.PDBIO()
    io.set_structure(struct)
    io.save(outpdb_fname)


###########################
#           MAIN          #
###########################
if (__name__ == "__main__"):
    ###### params #######
    receptor_chain_id="A"
    peptide_chain_id="B"
    #####################
    if(len(sys.argv) <> 2):
        exit_usage()
    pdb_fname = sys.argv[1]
    transform_peptide(pdb_fname, "start_rand.pdb",
                      receptor_chain_id, peptide_chain_id)
