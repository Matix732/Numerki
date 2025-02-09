import math

# Formowanie funkcji
def f(x):
    return (1/4) * x**4 - (1/2) * x**2 - (1/16) * x

# Metoda złotego podziału do znajdowania minimum funkcji
def golden_section_search(f, a, b, tol=1e-6):
    gr = (math.sqrt(5) - 1) / 2 
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = f(c)
    fd = f(d)
    nfev = 2  # Licznik ewaluacji funkcji
    while abs(c - d) > tol:
        if fc < fd:
            b = d
            d = c
            fd = fc
            c = b - gr * (b - a)
            fc = f(c)
        else:
            a = c
            c = d
            fc = fd
            d = a + gr * (b - a)
            fd = f(d)
        nfev += 1
    return (a + b) / 2, nfev

# Metoda Brenta do znajdowania minimum funkcji
def brent_method(f, a, b, tol=1e-6, max_iter=100):
    x = w = v = a + 0.5 * (b - a)
    fx = fw = fv = f(x)
    d = e = b - a
    nfev = 1  # Licznik ewaluacji funkcji
    for i in range(max_iter):
        m = 0.5 * (a + b)
        tol1 = tol * abs(x) + 1e-10
        if abs(x - m) <= 2 * tol1 - 0.5 * (b - a):
            break
        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2 * (q - r)
            if q > 0:
                p = -p
            q = abs(q)
            etemp = e
            e = d
            if abs(p) >= abs(0.5 * q * etemp) or p <= q * (a - x) or p >= q * (b - x):
                e = b - x if x >= m else a - x
                d = 0.5 * e
            else:
                d = p / q
                u = x + d
                if u - a < 2 * tol1 or b - u < 2 * tol1:
                    d = tol1 if x < m else -tol1
        else:
            e = b - x if x >= m else a - x
            d = 0.5 * e
        u = x + d if abs(d) >= tol1 else x + (tol1 if d > 0 else -tol1)
        fu = f(u)
        nfev += 1
        if fu <= fx:
            if u >= x:
                a = x
            else:
                b = x
            v, w, x = w, x, u
            fv, fw, fx = fw, fx, fu
        else:
            if u >= x:
                b = u
            else:
                a = u
            if fu <= fw or w == x:
                v, w = w, u
                fv, fw = fw, fu
            elif fu <= fv or v == x or v == w:
                v = u
                fv = fu
    return x, nfev

# Zdefiniowanie przedziału poszukiwań dla uzyskania tych samych wyników
a = -2
b = 2

# Znalazienie minimum funkcji za pomocą metody złotego podziału
minimum_golden, nfev_golden = golden_section_search(f, a, b)
value_golden = f(minimum_golden)

# Znalezienie minimum funkcji za pomocą metody Brenta
minimum_brent, nfev_brent = brent_method(f, a, b)
value_brent = f(minimum_brent)

# Wyświetlenie wyników
print(f"Golden Section Search - Minimum at x = {minimum_golden:.6f}, y = {value_golden:.6f}")
print(f"Brent's Method - Minimum at x = {minimum_brent:.6f}, y = {value_brent:.6f}")

# Porównanie liczby obliczeń funkcji
print(f"Golden Section Search evaluations: {nfev_golden}")
print(f"Brent's Method evaluations: {nfev_brent}")