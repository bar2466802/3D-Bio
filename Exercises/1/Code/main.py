import gudhi.wasserstein
import numpy as np

dgm1 = np.array([[2.7, 3.7],[9.6, 14.],[34.2, 34.974]])
dgm2 = np.array([[2.8, 4.45],[9.5, 14.1]])
message = "Wasserstein distance value = " + '%.2f' % gudhi.wasserstein.wasserstein_distance(dgm1, dgm2, order=1., internal_p=2.)
print(message)

