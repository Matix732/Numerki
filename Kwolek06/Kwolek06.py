import numpy as np
np.set_printoptions(suppress=True, linewidth=1000)

A = np.array([[2, -1, 0, 0, 1],
            [-1, 2, 1, 0, 0],
            [0, 1, 1, 1, 0],
            [0, 0, 1, 2, -1],
            [1, 0, 0, -1, 2]])

EIGEN_VALUE = 0.38197


# Algorytm iteracyjny do znajdowania wektora własnego
def find_eigen_vector(A, eigen_value, iterations=1000, tol=1e-10):
    size = len(A)

    vec = np.random.rand(size)
    vec /= np.linalg.norm(vec)
    for j in range(iterations):
        vec_z = np.linalg.solve(A - eigen_value * np.eye(size), vec)
        vec_new = vec_z / np.linalg.norm(vec_z)

        if np.linalg.norm(vec_new - vec) < tol:
            break

        vec = vec_new

    return vec_new

eigen_vector = find_eigen_vector(A, EIGEN_VALUE)

# Wypisanie wyniku
print("Eigen vector for eigen value of:", EIGEN_VALUE, "is:\n", eigen_vector)