import cmath

# Współczynniki zadanych wielomianów
coeffs_12a = [243, -486, 783, -990, 558, -28, -72, 16]
coeffs_12b = [1, 1, 3, 2, -1, -3, -11, -8, -12, -4, -4]
coeffs_12c = [1, 1j, -1, -1j, 1]

# Znajdywanie pierwiastków wielomianów metodą Laguerre'a
def laguerre(poly, z0, tol=1e-12, max_iter=100):
    n = len(poly) - 1  
    
    for _ in range(max_iter):
        P = sum(c * (z0 ** (n - i)) for i, c in enumerate(poly))  
        
        if abs(P) < tol:
            return z0 
        
        # Obliczanie pochodnych
        dP = sum(c * (n - i) * (z0 ** (n - i - 1)) for i, c in enumerate(poly[:-1])) 
        d2P = sum(c * (n - i) * (n - i - 1) * (z0 ** (n - i - 2)) for i, c in enumerate(poly[:-2]))
        G = dP / P
        H = G ** 2 - d2P / P
       
        denom1 = G + cmath.sqrt((n - 1) * (n * H - G ** 2))
        denom2 = G - cmath.sqrt((n - 1) * (n * H - G ** 2))
       
        if abs(denom1) > abs(denom2):
            dz = n / denom1
        else:
            dz = n / denom2
        z0 -= dz
        
        if abs(dz) < tol:
            return z0  
    return z0 

# Deflacja wielomianu względem znalezionego pierwiastka
def deflate_polynomial(poly, root):
    n = len(poly) - 1
    new_poly = [0] * n
    new_poly[0] = poly[0]
    for i in range(1, n):
        new_poly[i] = poly[i] + new_poly[i - 1] * root
    return new_poly

# Znajdowanie wszystkich pierwiastków wielomianu
def find_all_roots(poly, tol=1e-12, max_iter=100):
    roots = []
    while len(poly) > 1:
        z0 = complex(1, 1)
        root = laguerre(poly, z0, tol, max_iter)
        roots.append(root)
        poly = deflate_polynomial(poly, root)
    return roots

# Znajdowanie pierwiastków podanych wielomianów
roots_12a = find_all_roots(coeffs_12a)
roots_12b = find_all_roots(coeffs_12b)
roots_12c = find_all_roots(coeffs_12c)

# Wyświetlanie wyniku
def format_roots(roots):
    return ', '.join(f'{root:.6g}' for root in roots)

print("Skompresowane pierwiastki równania 12a:", format_roots(roots_12a))
print("Skompresowane pierwiastki równania 12b:", format_roots(roots_12b))
print("Skompresowane pierwiastki równania 12c:", format_roots(roots_12c))