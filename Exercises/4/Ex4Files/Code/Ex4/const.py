from tensorflow.keras.models import load_model
import argparse
import numpy as np
import utils
from Bio.PDB import *


# change these parameters

DIST_STD = 0
OMEGA_STD = 0
THETA_STD = 0
PHI_STD = 0


def write_const_dist(const_file, constraints, cdr3_s, seq):
    """
    writes the distance constraints into const_file in Rosetta format, uses harmonic function
    :param const_file: file to write the constraints into
    :param constraints: constraints matrix with no padding
    :param cdr3_s: start index of cdr3
    :param seq: nanobody sequence
    :return: None
    """
    length = len(constraints)
    for i in range(length):
        for j in range(i+1, length):  # symmetry
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s

            atom_i_type = 'CA' if seq[atom_i] == 'G' else 'CB'  # GLY
            atom_j_type = 'CA' if seq[atom_j] == 'G' else 'CB'  # GLY

            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1

            const_file.write("AtomPair {} {} {} {} HARMONIC {:.5f} {}\n".format(atom_i_type, atom_i, atom_j_type, atom_j, constraints[i,j], DIST_STD))


def write_const_omega(const_file, constraints, cdr3_s, seq):
    """
   writes the omega constraints into const_file in Rosetta format, uses circular harmonic function
   :param const_file: file to write the constraints into
   :param constraints: constraints matrix with no padding
   :param cdr3_s: start index of cdr3
   :param seq: nanobody sequence
   :return: None
   """
    length = len(constraints)
    for i in range(length):
        for j in range(i + 1, length):  # symmetry
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s
            if seq[atom_i] == 'G' or seq[atom_j] == 'G':  # GLY
                continue
            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1

            const_file.write("Dihedral CA {} CB {} CB {} CA {} CIRCULARHARMONIC {:.5f} {}\n".format(atom_i, atom_i, atom_j, atom_j, constraints[i, j], OMEGA_STD))


def write_const_theta(const_file, constraints, cdr3_s, seq):
    """
   writes the theta constraints into const_file in Rosetta format, uses circular harmonic function
   :param const_file: file to write the constraints into
   :param constraints: constraints matrix with no padding
   :param cdr3_s: start index of cdr3
   :param seq: nanobody sequence
   :return: None
   """
    length = len(constraints)
    for i in range(length):
        for j in range(length):
            if i == j:  # same atom...
                continue
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s
            if seq[atom_i] == 'G' or seq[atom_j] == 'G':  # GLY
                continue
            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1
            const_file.write("Dihedral N {} CA {} CB {} CB {} CIRCULARHARMONIC {:.5f} {}\n".format(atom_i, atom_i, atom_i, atom_j, constraints[i, j], THETA_STD))


def write_const_phi(const_file, constraints, cdr3_s, seq):
    """
   writes the phi constraints into const_file in Rosetta format, uses circular harmonic function
   :param const_file: file to write the constraints into
   :param constraints: constraints matrix with no padding
   :param cdr3_s: start index of cdr3
   :param seq: nanobody sequence
   :return: None
   """
    length = len(constraints)
    for i in range(length):
        for j in range(length):
            if i == j:  # same atom...
                continue
            atom_i = i + cdr3_s
            atom_j = j + cdr3_s
            if seq[atom_i] == 'G' or seq[atom_j] == 'G':  # GLY
                continue

            atom_i += 1  # pdb numbering starts from 1
            atom_j += 1  # pdb numbering starts from 1

            const_file.write("Angle CA {} CB {} CB {} CIRCULARHARMONIC {:.5f} {}\n".format(atom_i, atom_i, atom_j, constraints[i, j], PHI_STD))


def write_const_file(pdb_sequence, constraints_list, output_file):
    """
    writes the constraints file in Rosetta format.
    :param pdb_sequence: nanobody sequence (String)
    :param constraints_list: list of constraints obtained by using predict
    :param output_file: output file path
    :return: None
    """
    cdr_start, cdr_end = find_cdr3(pdb_sequence)
    cdr_len = (cdr_end + 1) - cdr_start

    # remove padding
    distance_const = remove_padding(constraints_list[0][0], cdr_len)
    omega_const = remove_padding(constraints_list[1][0], cdr_len)
    theta_const = remove_padding(constraints_list[2][0], cdr_len)
    phi_const = remove_padding(constraints_list[3][0], cdr_len)

    # retrieve the real distances and angles
    distance_const = distance_const[:, :, 0] * 10  # we divided by factor 10 in the processing of the data
    omega_const = np.arctan2(omega_const[:, :, 0], omega_const[:, :, 1])  # angle = arctan(sin, cos)
    theta_const = np.arctan2(theta_const[:, :, 0], theta_const[:, :, 1])  # angle = arctan(sin, cos)
    phi_const = np.arctan2(phi_const[:, :, 0], phi_const[:, :, 1])  # angle = arctan(sin, cos)

    with open(output_file, "w") as const_file:
        write_const_dist(const_file, distance_const, cdr_start, pdb_sequence)
        write_const_omega(const_file, omega_const, cdr_start, pdb_sequence)
        write_const_theta(const_file, theta_const, cdr_start, pdb_sequence)
        write_const_phi(const_file, phi_const, cdr_start, pdb_sequence)


if __name__ == '__main__':
    """
    this program recieves a path to a nanobody file and a trained network and creates a constraints file named constraints_file 
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file", help="nanobody pdb file")
    parser.add_argument("network", help="trained model")
    parser.add_argument("output_file", help="path of output file")
    args = parser.parse_args()

    # load the trained model
    network_model = load_model(args.network)

    # get the sequence of the nanobody (so we can skip GLY when needed)
    pdb_chain = PDBParser().get_structure(args.pdb_file, args.pdb_file)[0]["H"]
    sequence, aa_chains = get_seq_aa(pdb_chain)

    # generate the constraints using the network you built
    pdb_constraints = network_model.predict(np.array([generate_input(args.pdb_file)]))

    # write the constraints into a file, replace the empty list with your list.
    write_const_file(sequence, pdb_constraints, args.output_file)


