import numpy as np
import matplotlib.pyplot as plt

# Zadane punkty
nodes = np.array([-7/8, -5/8, -3/8, -1/8, 1/8, 3/8, 5/8, 7/8])

# Formowanie funkcji
def function(x):
    return 1 / (1 + 5 * x**2)

# Obliczanie wag dla interpolacji Hormanna-Floatera
def compute_weights(nodes, d_parameter):
    n = nodes.shape[0]
    weights = np.zeros(n)
    for k in range(n):
        for i in range(max(0, k-d_parameter), min(k, n-d_parameter) + 1):
            weight = (-1)**i
            for j in range(i, i + d_parameter):
                if j != k:
                    weight /= (nodes[j] - nodes[k])
            weights[k] += weight
    return weights

# Funkcja interpolująca
def interpolate_value(x, weights, nodes, values):
    if x in nodes:
        return function(x)
    numerator = sum(weights[i] * values[i] / (x - nodes[i]) for i in range(len(nodes)))
    denominator = sum(weights[i] / (x - nodes[i]) for i in range(len(nodes)))
    return numerator / denominator

# Parametry interpolacji
d_parameter = 3
values = function(nodes)
weights = compute_weights(nodes, d_parameter)

# Generowanie punktów do wykresu
x_plot = np.arange(nodes.min(), nodes.max() + 0.01, 0.01)
y_plot = [interpolate_value(x, weights, nodes, values) for x in x_plot]

# Rysowanie wykresu
plt.plot(x_plot, y_plot, label='Interpolacja Floatera-Hormanna')
plt.scatter(nodes, values, color='red', label='Punkty danych')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Interpolacja Floatera-Hormanna')
plt.grid(True)
plt.show()