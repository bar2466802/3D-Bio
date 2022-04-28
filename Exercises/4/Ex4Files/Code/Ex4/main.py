
if __name__ == '__main__':
    #  you can make all the data for the network in this section. use picke dump to save all the 5 matrices.
    # this way you won't have to generate them each time you train a newtork.
    # you can save the matrices to your drive and load them in your google colab file later.

    input_matrix = []
    dist_matrix = []
    omega_matrix = []
    theta_matrix = []
    phi_matrix = []

    data_path = "/Ex4_data"  # TODO: change path if needed

    for pdb in tqdm(os.listdir(data_path)):
        dist, omega, theta, phi = generate_label(os.path.join(data_path, pdb))

        input_matrix.append(generate_input(os.path.join(data_path, pdb)))
        dist_matrix.append(dist)
        omega_matrix.append(omega)
        theta_matrix.append(theta)
        phi_matrix.append(phi)

    save_path = data_path + "/Ex4Out"  # TODO: change path if needed

    pickle.dump(np.array(dist_matrix), open(os.path.join(save_path, "train_dist.pkl"), "wb"))
    pickle.dump(np.array(omega_matrix), open(os.path.join(save_path, "train_omega.pkl"), "wb"))
    pickle.dump(np.array(theta_matrix), open(os.path.join(save_path, "train_theta.pkl"), "wb"))
    pickle.dump(np.array(phi_matrix), open(os.path.join(save_path, "train_phi.pkl"), "wb"))

    pickle.dump(np.array(input_matrix), open(os.path.join(save_path, "train_input.pkl"), "wb"))

    print("Number of samples: {}".format(len(input_matrix)))
