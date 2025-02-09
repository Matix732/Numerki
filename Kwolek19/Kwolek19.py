import numpy as np

# Definicja funkcji Rosenbrocka
def rosenbrock(x, y):
    return (1 - x)**2 + 100 * (y - x**2)**2

# Obliczenie gradientu
def gradient(x, y):
    df_dx = -2 * (1 - x) - 400 * x * (y - x**2)
    df_dy = 200 * (y - x**2)
    return np.array([df_dx, df_dy])

# Metoda gradientu prostego do znalezienia minimum
def gradient_descent(starting_point, learning_rate, tol, max_iterations):
    point = np.array(starting_point)
    for i in range(max_iterations):
        grad = gradient(point[0], point[1])
        new_point = point - learning_rate * grad
        
        if np.linalg.norm(new_point - point) < tol:
            break
        
        point = new_point
    
    return point

# Parametry
starting_point = [0.5, 1.5]
learning_rate = 0.001
tol = 1e-7
max_iterations = 10000

# Znajdź minimum
minimum_point = gradient_descent(starting_point, learning_rate, tol, max_iterations)
print("Minimum found at:", minimum_point)
print("Function value at minimum:", rosenbrock(minimum_point[0], minimum_point[1]))