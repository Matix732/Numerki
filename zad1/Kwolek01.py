import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

A = csc_matrix([
    [4, 1, 0, 0, 0, 0, 1],
    [1, 4, 1, 0, 0, 0, 0],
    [0, 1, 4, 1, 0, 0, 0],
    [0, 0, 1, 4, 1, 0, 0],
    [0, 0, 0, 1, 4, 1, 0],
    [0, 0, 0, 0, 1, 4, 1],
    [1, 0, 0, 0, 0, 1, 4]
])

b = np.array([1, 2, 3, 4, 5, 6, 7])

# Do rozwiązania układu z zadaną macierzą rzadką wykorzystałem
# metodę eliminacji Gaussa-Crouta.

# Opcja #1 - bogactwo bibliotek Pythona
def solve_matrix_scipy(A,b):
    # Skorzystanie z biblioteki SciPy,
    # oraz funkcji spsolve przeznaczonej do rozwiązywania układów z macierzami rzadkimi.
    x = spsolve(A, b)
    return x

# Opcja #2 - własna implementacja algorytmu
def gauss_crout(A, b):
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.zeros((n, n))

    # Rozkład LU metodą Gaussa-Crouta
    for i in range(n):
        for j in range(i, n):
            L[j, i] = A[j, i] - sum(L[j, k] * U[k, i] for k in range(i))
        for j in range(i, n):
            if i == j:
                U[i, j] = 1
            else:
                U[i, j] = (A[i, j] - sum(L[i, k] * U[k, j] for k in range(i))) / L[i, i]

    # Rozwiązanie układu Ly = b
    y = np.zeros(n)
    for i in range(n):
        y[i] = (b[i] - sum(L[i, j] * y[j] for j in range(i))) / L[i, i]

    # Rozwiązanie układu Ux = y
    x = np.zeros(n)
    for i in range(n-1, -1, -1):
        x[i] = y[i] - sum(U[i, j] * x[j] for j in range(i+1, n))

    return x

# Wywołanie funkcji
C = gauss_crout(A, b)
C1 = solve_matrix_scipy(A, b)
print("Rozwiązanie układu równań:")
print("1) Rozwiązanie układu równań funkcją z biblioteki SciPy:\n", C1)
print("2) Rozkład LU metodą Gaussa-Crouta:\n", C)
