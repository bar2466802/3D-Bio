import subprocess


def fast_rmsd(model1, model2):
    """
    calculates the RMSD between two PDB files
    :param model1: path for pdb file 1
    :param model2: path for pdb file 2
    :return: RMSD (float)
    """
    cmd = "rmsd -t {} {} | tail -n1 ".format(model1, model2)
    result = subprocess.run(cmd, shell=True, capture_output=True).stdout.strip()
    return float(result)
