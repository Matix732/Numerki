import numpy as np

A = np.array([
    [19/12, 13/12, 5/6, 5/6, 5/6, -17/12],
    [13/12, 13/12, 5/6, -1/6, 5/6, -11/12],
    [5/6, 5/6, 5/6, -1/6, 5/6, 5/6],
    [5/6, -1/6, -1/6, -1/6, 5/6, 5/6],
    [13/12, -11/12, 5/6, 5/6, 13/12, 13/12],
    [-17/12, 13/12, 5/6, 5/6, 13/12, 19/12]
])

# Trójdiagonalizacja mecierzy metodą Householdera
def householder_algorithm(A):
    n = A.shape[0]
    T = np.copy(A)
    Q = np.eye(n)
    
    for k in range(n-2):
        x = T[k+1:, k]
        e = np.zeros_like(x)
        e[0] = np.linalg.norm(x)
        u = x - e
        v = u / np.linalg.norm(u)
        
        Q_k = np.eye(n)
        Q_k[k+1:, k+1:] -= 2.0 * np.outer(v, v)
        
        T = Q_k @ T @ Q_k.T
        Q = Q @ Q_k
    
    return T, Q

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
T, Q = householder_algorithm(A)

# Znalezienie wszystkich wartości własnych macierzy trójdiagonalnej
eigenvalues = qr_algorithm(T)

# Wyświetlenie wyników
print("Matrix A:\n", A)
print("Tridiagonal matrix:\n", T)
print("Eigenvalues:\n", eigenvalues)
