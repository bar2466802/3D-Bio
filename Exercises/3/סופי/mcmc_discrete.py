import argparse
import numpy as np
import pprint


def get_cmdline_parser():
    parser = argparse.ArgumentParser(
        description='Run MCMC on a discrete n x n configuration space.')
    parser.add_argument('n', type=int, default=5, nargs='?',
                        help='number of rows/columns of grid')
    parser.add_argument('m', type=int, default=1000, nargs='?',
                        help='number of iterations of MCMC optimization')
    parser.add_argument('kT', type=float, default=1.0, nargs='?',
                        help='kT - the denominator for the metropolis criterion'
                             ' (Boltzmann constant times temperature)')
    return parser


def is_valid(c, n):
    ''' Return True if c is a valid 2-D coordinate on an n x n grid
    with 0-based indices '''
    if len(c) != 2:
        return False
    return (c[0] >= 0 and c[1] >= 0 and c[0] < n and c[1] < n)


def get_p_accept_metropolis(dE, kT, p_forward, p_backward):
    '''
    return the probability to accept the metropolis criteria
    for the specified conditions
    dE - change in energy from current to proposed configuration
    kT - the factor of Boltzmann constant (kB) and the temperature (T)
    p_forward - probability to propose a move from current to proposed configuration
    p_backward - probability to propose a move from proposed to current configuration
    '''
    p = np.exp(-dE / kT) * p_backward / p_forward
    return min(p, 1.0)


def E(c):
    assert (len(c) == 2)
    return 1.0 * c[0] + 0.5 * c[1]


def get_neighbours(c, n):
    ''' get up/down/left/right neighbours on an n x n grid with 0-based indices'''
    assert (is_valid(c, n))
    ret_value = []
    if c[0] > 0:
        ret_value.append((c[0] - 1, c[1]))
    if c[0] < n - 1:
        ret_value.append((c[0] + 1, c[1]))
    if c[1] > 0:
        ret_value.append((c[0], c[1] - 1))
    if c[1] < n - 1:
        ret_value.append((c[0], c[1] + 1))
    return ret_value


def init_conf(n):
    return np.zeros((n, n))


def get_random_coord(n):
    return np.random.randint(n)


def prob_to_go_to_dest(org, dest, n):
    neighbours = get_neighbours(org, n)
    assert dest in neighbours
    return 1 / len(neighbours)


def prob_accept_new_proposal(current_state, new_prop, kT, n):
    return get_p_accept_metropolis(
        E(new_prop) - E(current_state),
        kT,
        prob_to_go_to_dest(current_state, new_prop, n),
        prob_to_go_to_dest(new_prop, current_state, n)
    )
import seaborn
import matplotlib.pyplot as plt

def q2_print_heat_map(configuration):
    plt.figure()
    m = np.sum(configuration)
    seaborn.heatmap(configuration/m, annot=True, fmt=".2f", linewidths=.5)
    plt.show()


def q4_print_heat_map(configuration):
    plt.figure()
    seaborn.heatmap(np.max(np.log1p(configuration)) - np.log1p(configuration), annot=True, fmt=".2f", linewidths=.5)
    plt.show()


def q8_print_heat_map(configuration, kT):
    plt.figure()
    seaborn.heatmap(kT*(np.max(np.log1p(configuration)) - np.log1p(configuration)), annot=True, fmt=".2f", linewidths=.5)
    plt.show()


def print_config(config):
    for row in config:
        print(",".join(map(str, row)))
        print("\n")

def run_mcmc(n, m, kT):
    conf = init_conf(n)
    coord = get_random_coord(n), get_random_coord(n)
    for x in range(m):
        neighbours = get_neighbours(coord, n)
        suggest_coord = neighbours[np.random.randint(len(neighbours))]
        coord = suggest_coord if np.random.random() <= prob_accept_new_proposal(coord, suggest_coord, kT, n) else coord
        conf[coord[0], coord[1]] += 1
    assert np.sum(conf) == m
    print_config(conf)
    return conf

def main():
    params = get_cmdline_parser().parse_args()
    n, m, kT = params.n, params.m, params.kT
    conf = run_mcmc(n, m, kT)
    q2_print_heat_map(conf)
    q4_print_heat_map(conf)
    q8_print_heat_map(conf, kT)


if __name__ == "__main__":
    main()
