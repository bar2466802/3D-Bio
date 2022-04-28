#################################################################
# FILE : my_align.py
# WRITER : Bar Melinarskiy & Rachel Ben Hamozeg
# EXERCISE : 3D-Bio ex1 2021
# DESCRIPTION: A program that aligns the given proteins structures.

# BONUS 9.1 IMPLEMENTED”.
#################################################################

# Load libraries
import sys

import numpy as numpy
from Bio.PDB import Superimposer, MMCIFIO
import os
from pathlib import Path
from Bio.PDB import Select
from Bio.PDB import PDBList, PDBParser
from typing import Text

# Constants
PDB = "pdb"
EXPECTED_ARGS_COUNT = 5
BONUS_EXPECTED_ARGS_COUNT = 6
PDB_FILE_FOLDER = "PdbFiles"
foldr = Path(PDB_FILE_FOLDER)

# Messages
NON_VALID_ARGS_ERROR = "ERROR: Non valid arguments were given"


# functions

#   Create a PDB object.
#     Arguments:
#      - name - string - name of the PDB
#      - chain - string - id of the desired chain inside the PDB
#     Return:
#       an Object with the PDB relevant info
def get_structure(pdb_name: Text):
    pdb_list = PDBList()
    filename = pdb_list.retrieve_pdb_file(pdb_name, pdir=PDB_FILE_FOLDER, file_format=PDB)
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    return parser.get_structure(pdb_name, filename)


# get C-alpha atoms from PDB file in the given chain
#     Arguments:
#       - pdb_object - Object - an Object with the PDB relevant info
def get_calpha_atoms(pdb_object, chain_id):
    channel = tuple(filter(lambda ch: ch.id == chain_id, pdb_object.get_chains()))
    if not channel:
        raise ValueError("Invalid channel id")
    # Iterate of all residues in each model in order to find proper atoms
    return [atom for atom in channel[0].get_atoms() if atom.name == "CA"]


class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered()


# Export the given PDB info after the alignment
#     Arguments:
#       - pdb_object - Object - an Object with the PDB relevant info
def create_output_file(pdb_object):
    io = MMCIFIO()
    io.set_structure(pdb_object)
    io.save(pdb_object.id + ".cif", select=NotDisordered())


# Calc the Bottleneck between the two given atoms lists -
# as taught in class Bottleneck =  max ||aki – bti||
#     Arguments:
#       - first_atoms_list - list - the atoms of the first structure
#       - second_atoms_list - list - the atoms of the second structure
def calc_bottleneck(first_atoms_list, second_atoms_list):
    max_dist = 0
    for atom1, atom2 in zip(first_atoms_list, second_atoms_list):
        vector1 = atom1.get_vector()
        vector2 = atom2.get_vector()
        dist = numpy.linalg.norm(vector1 - vector2)
        if max_dist < dist:
            max_dist = dist

    print("The achieved Bottleneck is: " + str(max_dist))


argumentsCounter = len(sys.argv)
print("Total arguments passed:", argumentsCounter)
if argumentsCounter not in [EXPECTED_ARGS_COUNT, BONUS_EXPECTED_ARGS_COUNT]:
    sys.exit(NON_VALID_ARGS_ERROR)

if not os.path.exists(PDB_FILE_FOLDER):
    os.makedirs(PDB_FILE_FOLDER)
pdbFoldr = Path(PDB_FILE_FOLDER)

firstPDB, secondPDB = get_structure(sys.argv[1]), get_structure(sys.argv[3])

firstPdbChain = get_calpha_atoms(firstPDB, sys.argv[2])
secondPdbChain = get_calpha_atoms(secondPDB, sys.argv[4])

use_bottleneck = argumentsCounter == BONUS_EXPECTED_ARGS_COUNT and sys.argv[5].lower() == 'true'

if len(firstPdbChain) == len(secondPdbChain):
    superImposer = Superimposer()
    superImposer.set_atoms(firstPdbChain, secondPdbChain)
    superImposer.apply(secondPDB.get_atoms())
    if use_bottleneck:  # Add an option to use your alternative measure for RMSD
        calc_bottleneck(firstPDB.get_atoms(), secondPDB.get_atoms())
    else:
        print("The achieved RMSD is: " + str(superImposer.rms))

else:
    raise ValueError("Atom's len don't match.")
create_output_file(firstPDB)
create_output_file(secondPDB)
