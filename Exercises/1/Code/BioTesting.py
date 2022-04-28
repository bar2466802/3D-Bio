# Load the Pkg
import Bio
import sys
from Bio.PDB import *
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO

# print(Bio.__version__)


# Check the Attributes, print with new-line between each attribute
# print("\n".join(dir(Seq)))

# Create a General DNA sequence
# my_seq = Seq("AGTACACTGGT")
# print("Original Seq: " + my_seq)
# print("Complement seq: " + my_seq.complement())
# print("Reverse complement seq: " + my_seq.reverse_complement())

# for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))


# for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))

# pdbl.retrieve_pdb_file("1FAT", file_format=PDB)

# imported libs
import sys
import os
import math
from Bio.PDB import *
from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import parse_pdb_header
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Constants
PDB = "pdb"
FASTA = "fasta"
EXPECTED_ARGS_COUNT = 5
PDB_FILE_FOLDER = "PdbFiles"
NON_VALID_ARGS_ERROR = "ERROR: Non valid arguments were given"


# functions
def createPdbObject(name, chain):
    pdbFileInfo = dict.𝑓𝑟𝑜𝑚𝑘𝑒𝑦𝑠(('name', 'filePath', 'fileName', 'chain', 'residue'))
    pdbFileInfo["name"] = name
    pdbFileInfo["chainId"] = chain.upper()
    fileName = "pdb" + name + ".ent"
    foldr = Path(PDB_FILE_FOLDER)
    filePath = foldr / fileName
    pdbFileInfo["fileName"] = fileName
    pdbFileInfo["filePath"] = filePath
    pdbFilesList.append(pdbFileInfo)

def getChinById(pdb):
    # Parsing structure from the PDB file
    # The PERMISSIVE flag indicates that a number of common problems associated with
    # PDB files will be ignored. (but note that some atoms and/or residues will be missing)
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(pdb['name'], pdb['filePath'])
    chains = structure.get_chains()
    for chain in chains:
        if chain.get_id() == pdb['chainId']:
            pdb['chain'] = chain
            return


# Main
# Check we got the 4 arguments as needed
argumentsCounter = len(sys.argv)
print("Total arguments passed:", argumentsCounter)
if argumentsCounter != EXPECTED_ARGS_COUNT:
    sys.exit(NON_VALID_ARGS_ERROR)

pdbFilesList = []
createPdbObject(sys.argv[1], sys.argv[2])
createPdbObject(sys.argv[3], sys.argv[4])

# Create folder for PDB files
if not os.path.exists(PDB_FILE_FOLDER):
    os.makedirs(PDB_FILE_FOLDER)
pdbFoldr = Path(PDB_FILE_FOLDER)
# Fetch the PDB files from RCSB if necessary
# try and except block to prevent re downloading the same file on multiple executions of the code
for pdb in pdbFilesList:
    try:
        if os.path.isfile(pdb['filePath']):
            print("pdb file: " + pdb['name'] + " was already fetched")
        else:
            # will download the file only if it doen not already exists
            pdbl = PDBList()
            filename = pdbl.retrieve_pdb_file(pdb['name'], pdir=pdbFoldr, file_format=PDB)
    except:
        print("an error occurred while trying to fetch the pdb file: " + pdb['name'])

    getChinById(pdb)
    # fasta_record_list = []
    # # Now that is clarified, let’s return to parsing the PDB header.
    # # The structure object has an attribute called header which is a Python dictionary
    # # that maps header records to their values.
    # resolution = structure.header["resolution"]
    # # The available keys are name, head, deposition_date, release_date, structure_method,
    # # resolution, structure_reference (which maps to a list of references),
    # # journal_reference, author, compound (which maps to a dictionary with various information
    # # about the crystallized compound), has_missing_residues, missing_residues, and astral
    # # (which maps to dictionary with additional information about the domain if present).
    # keywords = structure.header["keywords"]
    #
    # # The dictionary can also be created without creating a Structure object, ie.
    # # directly from the PDB file:
    #
    # with open(filePath, "r") as handle:
    #     header_dict = parse_pdb_header(handle)
    #
    # # Parsing primary structure from the PDB file
    # # record_list = SeqIO.parse(filePath, "pdb-atom")
    # # fasta_record_list = []
    # # for record in record_list:
    # #     header = ">" + record.id
    # #     sequence = record.seq
    # #     fasta_record_list.append(SeqRecord(id=header, seq=sequence))
    # #
    # # # writing the primary structures to fasta file
    # # fastaFileName = pdb['name'] + ".fasta"
    # # fastaPath = pdbFoldr/fastaFileName  # destination path for out fasta file
    # # SeqIO.write(fasta_record_list, fastaPath, format=FASTA)
    #
    # # The id of the Model object is an integer, which is derived from the position of the model
    # # in the parsed file (they are automatically numbered starting from 0).
    # # Crystal structures generally have only one model (with id 0),
    # # while NMR files usually have several models. Whereas many PDB parsers assume that there
    # # is only one model, the Structure class in Bio.
    # # PDB is designed such that it can easily handle PDB files with more than one model.
    # model = structure[0]  # The Model object stores a list of Chain children.
    #
    # # The id of a Chain object is derived from the chain identifier in the PDB/mmCIF file,
    # # and is a single character (typically a letter).
    # # Each Chain in a Model object has a unique id. As an example,
    # # to get the Chain object with given identifier from a Model object, use:
    # chain = model[pdb['chain']]  # The Chain object stores a list of Residue children.
    #
    # # a Residue object stores a set of Atom children. It also contains a string that specifies the
    # # residue name (e.g. “ASN”) and the segment identifier of the residue
    # hasResidue = chain.has_id(6888)
    # if hasResidue:
    #     # Full id
    #     residue = chain[(" ", 6888, " ")]
    #     # Shortcut id
    #     residueS = chain[6888]
    #
    #     name = residue.get_resname()  # returns the residue name, e.g. "ASN"
    #     is_disordered = residue.is_disordered()  # returns 1 if the residue has disordered atoms
    #     segid = residue.get_segid()  # returns the SEGID, e.g. "CHN1"
    #     hasAtoms = residue.has_id("N") & residue.has_id("CA")  # test if a residue has a certain atom
    #     # You can use is_aa(residue) to test if a Residue object is an amino acid.
    #     isAminoAcid = is_aa(residue)
    #
    #     # The Atom object stores the data associated with an atom, and has no children.
    #     # The id of an atom is its atom name (e.g. “OG” for the side chain oxygen of a Ser residue).
    #     # An Atom id needs to be unique in a Residue. Again, an exception is made for disordered atoms
    #     # In a PDB file, an atom name consists of 4 chars, typically with leading and trailing spaces.
    #     # Often these spaces can be removed for ease of use
    #     # (e.g. an amino acid C α atom is labeled “.CA.” in a PDB file, where the dots represent spaces)
    #     #
    #     # The atomic data stored includes the atom name, the atomic coordinates
    #     # (including standard deviation if present), the B factor
    #     # (including anisotropic B factors and standard deviation if present),
    #     # the altloc specifier and the full atom name including spaces.
    #     # Less used items like the atom element number or the atomic charge sometimes specified
    #     # in a PDB file are not stored.
    #     # To manipulate the atomic coordinates, use the transform method of the Atom object.
    #     # Use the set_coord method to specify the atomic coordinates directly.
    #     #
    #     # An Atom object has the following additional methods:
    #     if hasAtoms:
    #         atom = residue["N"]
    #         atomName = atom.get_name()  # atom name (spaces stripped, e.g. "CA")
    #         atomId = atom.get_id()  # id (equals atom name)
    #         atomCoordinates = atom.get_coord()  # atomic coordinates
    #         atomVector = atom.get_vector()  # atomic coordinates as Vector object
    #         atomBfactor = atom.get_bfactor()  # isotropic B factor
    #         atomOccupancy = atom.get_occupancy()  # occupancy
    #         altloc = atom.get_altloc()  # alternative location specifier
    #         sigatm = atom.get_sigatm()  # standard deviation of atomic parameters
    #         siguij = atom.get_siguij()  # standard deviation of anisotropic B factor
    #         anisou = atom.get_anisou()  # anisotropic B factor
    #         fullname = atom.get_fullname()  # atom name (with spaces, e.g. ".CA.")
    #
    #         # get atom coordinates as vectors
    #         n = residue["N"].get_vector()
    #         c = residue["C"].get_vector()
    #         ca = residue["CA"].get_vector()
    #         # center at origin
    #         n = n - ca
    #         c = c - ca
    #         # find rotation matrix that rotates n
    #         # -120 degrees along the ca-c vector
    #         rot = rotaxis(-math.pi * 120.0 / 180.0, c)
    #         # apply rotation to ca-n vector
    #         cb_at_origin = n.left_multiply(rot)
    #         # put on top of ca atom
    #         cb = cb_at_origin + ca
    # # Iterating through all atoms of a structure
    # for model in structure:
    #     for chain in model:
    #         for residue in chain:
    #             for atom in residue:
    #                 print(atom)
    #
    # # There is a shortcut if you want to iterate over all atoms in a structure:
    # atoms = structure.get_atoms()
    # for atom in atoms:
    #     print(atom)
    #
    # # Similarly, to iterate over all atoms in a chain, use
    # atoms = chain.get_atoms()
    # for atom in atoms:
    #     print(atom)
    #
    # # Iterating over all residues of a model
    # residues = model.get_residues()
    # for residue in residues:
    #     print(residue)
    #
    # # You can also use the Selection.unfold_entities function to get all residues from a structure:
    # res_list = Selection.unfold_entities(structure, "R")
    # # or to get all atoms from a chain:
    # atom_list = Selection.unfold_entities(chain, "A")
    # # Obviously, A=atom, R=residue, C=chain, M=model, S=structure. You can use this
    # # to go up in the hierarchy, e.g. to get a list of (unique) Residue or Chain parents
    # # from a list of Atoms:
    # residue_list = Selection.unfold_entities(atom_list, "R")
    # chain_list = Selection.unfold_entities(atom_list, "C")
    #
    # # Print out the coordinates of all CA atoms in a structure with B factor greater than 50
    # for model in structure.get_list():
    #     for chain in model.get_list():
    #         for residue in chain.get_list():
    #             if residue.has_id("CA"):
    #                 ca = residue["CA"]
    #                 if ca.get_bfactor() > 50.0:
    #                     print(ca.get_coord())
    # # Extracting polypeptides from a Structure object
    # ppb = PPBuilder()
    # model_nr = 1
    # polypeptide_list = ppb.build_peptides(structure, model_nr)
    # for polypeptide in polypeptide_list:
    #     print(polypeptide)
    #
    # # Note that in the this case only model 0 of the structure is considered by PolypeptideBuilder.
    # # However, it is possible to use PolypeptideBuilder to build Polypeptide objects from Model
    # # and Chain objects as well.
    #
    # # Using C-N
    # for pp in ppb.build_peptides(structure):
    #     print(pp.get_sequence())
    #
    # # Using CA-CA
    # ppb = CaPPBuilder()
    # for pp in ppb.build_peptides(structure):
    #     print(pp.get_sequence())

    # Writing mmCIF files, The MMCIFIO class can be used to write structures to the mmCIF file format
    # io = MMCIFIO()
    # io.set_structure(structure)
    # mmCIFFileName = pdb['name'] + ".cif"
    # mmCIFPath = pdbFoldr / mmCIFFileName  # destination path for out MMCIFIO file
    # io.save(mmCIFPath)

# Superimposing two structures
# Use a Superimposer object to superimpose two coordinate sets.
# This object calculates the rotation and translation matrix that rotates two lists of atoms
# on top of each other in such a way that their RMSD is minimized.
# Of course, the two lists need to contain the same number of atoms.
# The Superimposer object can also apply the rotation/translation to a list of atoms.
# The rotation and translation are stored as a tuple in the rotran attribute of
# the Superimposer object (note that the rotation is right multiplying!).
# The RMSD is stored in the rmsd attribute.
# The algorithm used by Superimposer comes from [Golub & Van Loan]
# and makes use of singular value decomposition
# (this is implemented in the general Bio.SVDSuperimposer module).
sup = Superimposer()
# Specify the atom lists
# 'fixed' and 'moving' are lists of Atom objects
# The moving atoms will be put on the fixed atoms
fixed = pdbFilesList[0]['chain'].get_atoms()
moving = pdbFilesList[1]['chain'].get_atoms()
sup.set_atoms(fixed, moving)
# Print rotation/translation/rmsd
print(sup.rotran)
print(sup.rms)
# Apply rotation/translation to the moving atoms
sup.apply(moving)


