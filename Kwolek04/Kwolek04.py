import numpy as np
import copy
np.set_printoptions(suppress=True, linewidth=1000)

A = np.array([[19, 13, 10, 10, 13, -17],
              [13, 13, 10, 10, -11, 13],
              [10, 10, 10, -2, 10, 10],
              [10, 10, -2, 10, 10, 10],
              [13, -11, 10, 10, 13, 13],
              [-17, 13, 10, 10, 13, 19]]) / 12

# Trójdiagonalizacja macierzy metodą Householdera
def householder(A):
    size = len(A)
    R = copy.deepcopy(A)
    for i in range(size - 1):
        alpha = -np.sign(R[i + 1][i]) * np.linalg.norm(R[i + 1:, i])
        r = np.sqrt((alpha * alpha - R[i + 1, i] * alpha) / 2)
        v = np.zeros_like(R[i:, i])
        v[1] = (R[i + 1, i] - alpha) / (2 * r)
        v[2:] = (R[i + 2:, i] / (2 * r))
        P = np.eye(len(v)) - 2 * np.outer(v, v)
        R[i:, i:] = np.dot(np.dot(P, R[i:, i:]), P)
    return R

# Poszukiwanie wartości własnych macierzy trójdiagonalnej metodą rozkładu QR
def qr_algorithm(T, num_iter=1000, tol=1e-9):
    n = T.shape[0]
    A_k = np.copy(T)
    
    for _ in range(num_iter):
        Q, R = np.linalg.qr(A_k)
        A_k = R @ Q
        
        if np.allclose(A_k - np.diag(np.diag(A_k)), 0, atol=tol):
            break
    
    eigenvalues = np.diag(A_k)
    return eigenvalues

# Sprowadzenie macierzy A do postaci trójdiagonalnej
T = householder(A)

# Znalezienie wszystkich wartości własnych macierzy trójdiagonalnej
eigenvalues = qr_algorithm(T)

# Wyświetlenie wyników
print("Matrix A:\n", A)
print()
print("Tridiagonal matrix:\n", T)
print()
print("Eigenvalues:\n", eigenvalues)