import numpy as np

def laguerre_method(coeffs, tol=1e-9, max_iter=100):
    def poly_eval(coeffs, x):
        return np.polyval(coeffs, x)
    
    def poly_deriv(coeffs):
        return np.polyder(coeffs)
    
    n = len(coeffs) - 1
    roots = []
    
    while n > 0:
        x = np.random.randn()
        for _ in range(max_iter):
            f = poly_eval(coeffs, x)
            f1 = poly_eval(poly_deriv(coeffs), x)
            f2 = poly_eval(poly_deriv(poly_deriv(coeffs)), x)
            
            if abs(f) < tol:
                break
            
            G = f1 / f
            H = G**2 - f2 / f
            denom1 = G + np.sqrt((n - 1) * (n * H - G**2))
            denom2 = G - np.sqrt((n - 1) * (n * H - G**2))
            a = n / (denom1 if abs(denom1) > abs(denom2) else denom2)
            x -= a
        
        roots.append(x)
        coeffs = np.polydiv(coeffs, [1, -x])[0]
        n -= 1
    
    return np.array(roots)

# Example usage: Solve polynomial equation from (12a)
coeffs_12a = [243, -486, 783, -990, 558, -28, -72, 16]
roots_12a = laguerre_method(coeffs_12a)
print("Roots of equation (12a):", roots_12a)
