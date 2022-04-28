from HackatonUtils import *
import argparse
import os
import numpy as np
import tensorflow as tf

HEADER = "HEADER    IMMUNE SYSTEM - NANOBODY                           \nTITLE     COMPUTATIONAL MODELING     \nREMARK 777 MODEL GENERATED BY NANONET \n"
ATOM_LINE = "ATOM{}{}  CA  {} H{}{}{}{:.3f}{}{:.3f}{}{:.3f}  1.00 0.00           C\n"
END_LINE = "END"


def matrix_to_pdb(pdb_file, seq, coord_matrix):
    """
    translates a matrix of Ca x,y,z coordinates to PDB format
    :param pdb_file: pdb file to write to
    :param seq: Nb sequence
    :param coord_matrix: NanoNet output
    :return: None
    """
    pdb_file.write(HEADER)
    seq = pad_seq(seq)
    i = 1
    for aa in range(len(seq)):
        if seq[aa] != "-":
            first_space = (7 - len(str(i))) * " "
            second_space = (4 - len(str(i))) * " "
            third_space = (12 - len("{:.3f}".format(coord_matrix[aa][0]))) * " "
            forth_space = (8 - len("{:.3f}".format(coord_matrix[aa][1]))) * " "
            fifth_space = (8 - len("{:.3f}".format(coord_matrix[aa][2]))) * " "
            pdb_file.write(ATOM_LINE.format(first_space, i, Polypeptide.one_to_three(seq[aa]), second_space, i, third_space, coord_matrix[aa][0],forth_space, coord_matrix[aa][1],fifth_space, coord_matrix[aa][2]))
            i += 1
    pdb_file.write(END_LINE)


if __name__ == '__main__':
    """
    receives path to a Nb fasta file and a path to a trained neural network and creates a pdb file (Ca only) according to
    the network prediction. the output file name is: "<fasta file name>_nanonet_ca.pdb"
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Nb fasta file")
    parser.add_argument("network", help="nanonet trained model")

    args = parser.parse_args()
    nanonet = tf.keras.models.load_model(args.network)

    # nanobody sequence
    sequence = get_sequence(args.fasta)

    # ca coordinates
    ca_coords = nanonet.predict(np.array([generate_input(sequence, fasta=False)]))[0]

    # create ca pdb file
    ca_file_name = "{}_nanonet_ca.pdb".format(args.fasta.split(".")[0])

    with open(ca_file_name, "w") as ca_file:
        matrix_to_pdb(ca_file, sequence, ca_coords)
