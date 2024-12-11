import numpy as np

def gauss_seidel(A, b, tolerance=1e-10, max_iterations=1000):
    n = len(A)
    x = np.zeros_like(b, dtype=np.float64)

    for iteration in range(max_iterations):
        x_new = np.copy(x)
        for i in range(n):
            s1 = sum(A[i, j] * x_new[j] for j in range(i))
            s2 = sum(A[i, j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - s1 - s2) / A[i, i]

        if np.linalg.norm(x_new - x, ord=np.inf) < tolerance:
            return x_new

        x = x_new

    raise Exception("Gauss-Seidel method did not converge")

def conjugate_gradient(A, b, tolerance=1e-10, max_iterations=1000):
    x = np.zeros_like(b, dtype=np.float64)
    r = b - A @ x
    p = r.copy()
    rs_old = np.dot(r, r)

    for iteration in range(max_iterations):
        Ap = A @ p
        alpha = rs_old / np.dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rs_new = np.dot(r, r)

        if np.sqrt(rs_new) < tolerance:
            return x

        p = r + (rs_new / rs_old) * p
        rs_old = rs_new

    raise Exception("Conjugate Gradient method did not converge")

def solve_system_methods():
    n = 128
    A = np.zeros((n, n))

    # Tworzenie macierzy A zgodnie ze strukturą z zadania
    for i in range(n):
        A[i, i] = 4
        if i > 0:
            A[i, i - 1] = 1
        if i < n - 1:
            A[i, i + 1] = 1
        if i == 0 or i == n - 1:
            A[i, (i + n - 1) % n] = 1

    b = np.ones(n)  # Wektor e (wszystkie elementy równe 1)

    # Rozwiązanie metodą Gaussa-Seidela
    x_gauss_seidel = gauss_seidel(A, b)

    # Rozwiązanie metodą gradientów sprzężonych
    x_conjugate_gradient = conjugate_gradient(A, b)

    return x_gauss_seidel, x_conjugate_gradient

# Wywołanie funkcji i wypisanie rozwiązań
solution_gauss_seidel, solution_conjugate_gradient = solve_system_methods()
print("Rozwiązanie metodą Gaussa-Seidela:", solution_gauss_seidel)
print("Rozwiązanie metodą gradientów sprzężonych:", solution_conjugate_gradient)
