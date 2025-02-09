import numpy as np

# Definicja układu równań
def equations(vars):
    x, y = vars
    eq1 = 2*x**2 + y**2 - 2
    eq2 = (x - 0.5)**2 + (y - 1)**2 - 0.25
    return np.array([eq1, eq2])

# Definicja macierzy Jacobiego
def jacobian(vars):
    x, y = vars
    df1_dx = 4*x
    df1_dy = 2*y
    df2_dx = 2*(x - 0.5)
    df2_dy = 2*(y - 1)
    return np.array([[df1_dx, df1_dy], [df2_dx, df2_dy]])

# Metoda Newtona-Raphsona do rozwiązywania układów nieliniowych
def newton_raphson(initial_guess, tol=1e-6, max_iter=100):
    vars = initial_guess
    for i in range(max_iter):
        f = equations(vars)
        J = jacobian(vars)
        delta = np.linalg.solve(J, -f)
        vars += delta
        if np.linalg.norm(delta) < tol:
            break
    return vars

# Początkowe przybliżenie
initial_guess = np.array([1.0, 1.0])

# Rozwiązanie układu równań
solution = newton_raphson(initial_guess)

# Wyświetlenie wyników
print(f"x = {solution[0]}, y = {solution[1]}")