import numpy as np
import matplotlib.pyplot as plt

def gauss_seidel(A, b, tolerance=1e-10, max_iterations=1000):
    n = len(A)
    x = np.zeros_like(b, dtype=np.float64)
    norms = []

    for iteration in range(max_iterations):
        x_new = np.copy(x)
        for i in range(n):
            s1 = sum(A[i, j] * x_new[j] for j in range(i))
            s2 = sum(A[i, j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - s1 - s2) / A[i, i]

        norm = np.linalg.norm(x_new - x)
        norms.append(norm)

        if norm < tolerance:
            return x_new, norms

        x = x_new

    raise Exception("Gauss-Seidel method did not converge")

def conjugate_gradient(A, b, tolerance=1e-10, max_iterations=1000):
    x = np.zeros_like(b, dtype=np.float64)
    r = b - A @ x
    p = r.copy()
    rs_old = np.dot(r, r)
    norms = []

    for iteration in range(max_iterations):
        Ap = A @ p
        alpha = rs_old / np.dot(p, Ap)
        x += alpha * p
        r -= alpha * Ap
        rs_new = np.dot(r, r)

        norm = np.sqrt(rs_new)
        norms.append(norm)

        if norm < tolerance:
            return x, norms

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
    x_gauss_seidel, norms_gauss_seidel = gauss_seidel(A, b)

    # Rozwiązanie metodą gradientów sprzężonych
    x_conjugate_gradient, norms_conjugate_gradient = conjugate_gradient(A, b)

    return norms_gauss_seidel, norms_conjugate_gradient

def plot_convergence(norms_gauss_seidel, norms_conjugate_gradient):
    plt.figure(figsize=(10, 6))
    plt.semilogy(range(1, len(norms_gauss_seidel) + 1), norms_gauss_seidel, label="Gauss-Seidel", marker='o')
    plt.semilogy(range(1, len(norms_conjugate_gradient) + 1), norms_conjugate_gradient, label="Conjugate Gradient", marker='x')
    plt.xlabel("Iteration")
    plt.ylabel("Norm of difference (log scale)")
    plt.title("Convergence comparison")
    plt.legend()
    plt.grid()
    plt.show()

def compare_complexity():
    n = 128
    cholesky_complexity = (1 / 3) * n**3
    gauss_seidel_complexity = n**2  # Approximation per iteration
    conjugate_gradient_complexity = n**2  # Approximation per iteration

    print("Cholesky complexity (one-time):", cholesky_complexity)
    print("Gauss-Seidel complexity (per iteration):", gauss_seidel_complexity)
    print("Conjugate Gradient complexity (per iteration):", conjugate_gradient_complexity)

# Wywołanie funkcji i porównanie
norms_gauss_seidel, norms_conjugate_gradient = solve_system_methods()
plot_convergence(norms_gauss_seidel, norms_conjugate_gradient)
compare_complexity()
