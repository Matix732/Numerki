import cmath

# Przykład użycia:
coeffs_a = [243, -486, 783, -990, 558, -28, -72, 16]
coeffs_b = [1, 1, 3, 2, -1, -3, -11, -8, -12, -4, -4]
coeffs_c = [1, 1j, -1, -1j, 1]

# obliczanie wartości wielomianu w punkcie z (z - przybliżenie pierwiastka)
def evaluate_polynomial(coeffs, z):
    return sum(coeff * (z ** (len(coeffs) - i - 1)) for i, coeff in enumerate(coeffs))

# oblicz wartość pierwszej pochodnej wielomianu w punkcie z
def evaluate_derivative(coeffs, z):
    return sum((len(coeffs) - i - 1) * coeff * (z ** (len(coeffs) - i - 2)) for i, coeff in enumerate(coeffs[:-1]))

# oblicz wartość drugiej pochodnej wielomianu w punkcie z
def evaluate_second_derivative(coeffs, z):
    return sum((len(coeffs) - i - 1) * (len(coeffs) - i - 2) * coeff * (z ** (len(coeffs) - i - 3)) for i, coeff in enumerate(coeffs[:-2]))

# metoda Laguerre'a do znajdowania pierwiastków wielomianu
def laguerre(coeffs, z, tol=1e-10, max_iter=100):
    n = len(coeffs) - 1
    for _ in range(max_iter):
        Pz = evaluate_polynomial(coeffs, z)
        if abs(Pz) < tol:
            return z
        G = evaluate_derivative(coeffs, z) / Pz
        H = G**2 - evaluate_second_derivative(coeffs, z) / Pz
        denom1 = G + cmath.sqrt((n - 1) * (n * H - G**2))
        denom2 = G - cmath.sqrt((n - 1) * (n * H - G**2))
        if abs(denom1) > abs(denom2):
            a = n / denom1
        else:
            a = n / denom2
        z = z - a
        if abs(a) < tol:
            return z
    return z

# deflacja wielomianu przez usunięcie znalezionego pierwiastka
def deflate_polynomial(coeffs, root):
    new_coeffs = [coeffs[0]]
    for i in range(1, len(coeffs)):
        new_coeffs.append(new_coeffs[-1] * root + coeffs[i])
    return new_coeffs[:-1]

# znajdowanie wszystkich pierwiastków wielomianu
def find_all_roots(coeffs, tol=1e-10, max_iter=100):
    roots = []
    current_coeffs = coeffs.copy()
    while len(current_coeffs) > 1:
        root = laguerre(current_coeffs, 0, tol, max_iter)
        roots.append(root)
        current_coeffs = deflate_polynomial(current_coeffs, root)
    return roots

# obliczanie pierwiastków równań
roots_a = find_all_roots(coeffs_a)
roots_b = find_all_roots(coeffs_b)
roots_c = find_all_roots(coeffs_c)

# wyświetlenie wyników
print("Pierwiastki równania 12a:\n", end=" ")
for root in roots_a:
    print(root, end=" ")
    print()
print("Pierwiastki równania 12b:\n", end=" ")
for root in roots_b:
    print(root, end=" ")
    print()
print("Pierwiastki równania 12c:\n", end=" ")
for root in roots_c:
    print(root, end=" ")
    print()