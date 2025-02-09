import numpy as np

# Formowanie funkcji
def integrand(x):
    return np.sin(np.pi * (1 + np.sqrt(x)) / (1 + x**2)) * np.exp(-x)

# Metoda trapezów
def trapezoidal_rule(f, a, b, n):
    h = (b - a) / n
    integral = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        integral += f(a + i * h)
    integral *= h
    return integral

# Implementacja metody Romberga
def romberg_integration(f, a, b, tol=1e-10):
    R = [[trapezoidal_rule(f, a, b, 1)]]
    n = 1
    while True:
        h = (b - a) / 2**n
        sum_f = sum(f(a + (2*k - 1) * h) for k in range(1, 2**(n-1) + 1))
        R.append([0.5 * R[n-1][0] + h * sum_f])
        for m in range(1, n + 1):
            R[n].append(R[n][m-1] + (R[n][m-1] - R[n-1][m-1]) / (4**m - 1))
        if abs(R[n][n] - R[n-1][n-1]) < tol:
            break
        n += 1
    return R

# Parametry całkowania
a = 0
b = 10 
tol = 1e-7

# Obliczanie całki na przedziale [0, 10]
result_tab = romberg_integration(integrand, a, b, tol)

# Wypisanie tabelki wyników
for i, row in enumerate(result_tab):
    print(f"R[{i}]: {[float(x) for x in row]}")

# Wynik końcowy
result = float(result_tab[-1][-1])
print(f"\nCałka wynosi: {result}\n")