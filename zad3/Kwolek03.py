import numpy as np

A = [
    [19/12, 13/12, 5/6, 5/6, 5/6, -17/12],
    [13/12, 13/12, 5/6, -1/6, 5/6, -11/12],
    [5/6, 5/6, 5/6, -1/6, 5/6, 5/6],
    [5/6, -1/6, -1/6, -1/6, 5/6, 5/6],
    [13/12, -11/12, 5/6, 5/6, 13/12, 13/12],
    [-17/12, 13/12, 5/6, 5/6, 13/12, 19/12]
]

# Metoda potęgowa (znalezienie jednej największej wartości)
def power_method(A, num_iter=1000, tol=1e-9):
    n = len(A)
    b_k = np.random.rand(n)
    
    for _ in range(num_iter):
        b_k1 = np.dot(A, b_k)
        b_k1_norm = np.linalg.norm(b_k1)
        b_k = b_k1 / b_k1_norm
        
        if np.linalg.norm(np.dot(A, b_k) - b_k1_norm * b_k) < tol:
            break
    
    eigenvalue = b_k1_norm
    eigenvector = b_k
    return eigenvalue, eigenvector

# Deflacja macierzy A względem wektora własnego i wartości własnej (znalezienie pozostałej wartości własnej)
def deflate(A, eigenvector, eigenvalue):
    eigenvector = eigenvector.reshape(-1, 1)
    A = A - eigenvalue * np.dot(eigenvector, eigenvector.T)
    return A

# Znajdywanie największej wartości własnej i odpowiadającego jej wektora własnego
eigenvalue1, eigenvector1 = power_method(A)

# Deflacja macierzy A
A_deflated = deflate(np.array(A), eigenvector1, eigenvalue1)

# Znajdywanie DRUGIEJ największej wartości własnej i odpowiadającego jej wektora własnego
eigenvalue2, eigenvector2 = power_method(A_deflated)

# Wyświetlenie wyników
print("Matrix A:\n", A)
print("Largest eigenvalue:", eigenvalue1)
print("Corresponding eigenvector:", eigenvector1)
print()
print("Second largest eigenvalue:", eigenvalue2)
print("Corresponding eigenvector:", eigenvector2)