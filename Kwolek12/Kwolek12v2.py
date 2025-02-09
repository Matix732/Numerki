import numpy as np
import matplotlib.pyplot as plt

# Zadane punkty 
nodes = np.array([-7/8, -5/8, -3/8, -1/8, 1/8, 3/8, 5/8, 7/8])

# Formowanie funkcji
def function(x):
    return 1.0 / (1 + 5 * x**2)

# Obliczanie wag dla interpolacji Hormanna-Floatera
def calc_weights(nodes, parameter):
    n = nodes.shape[0]
    weights = np.zeros(n)
    k = 0
    while k < n:
        for i in range(max(0, k-parameter), min(k, n-parameter) + 1):
            acc = (-1.0) ** i
            for j in range(i, i + parameter):
                if j != k:
                    acc /= (nodes[j] - nodes[k])
            weights[k] += acc
        k += 1
    return weights

# Interpolacja 
def interpoluj_wartosc(x, weights, nodes, node_values):
    i = 0
    while i < len(nodes):
        if x == nodes[i]:
            return function(x)
        i += 1
    numerator = sum(weights[i] / (x - nodes[i]) for i in range(len(nodes)))
    denominator = sum((weights[i] / (x - nodes[i])) * node_values[i] for i in range(len(nodes)))
    return denominator / numerator

# Parametry interpolacji
parameter = 3
node_values = function(nodes)
weights = calc_weights(nodes, parameter)

# Generowanie punktów do wykresu
x_values = np.arange(nodes.min(), nodes.max() + 0.01, 0.01)
y_values = [interpoluj_wartosc(x, weights, nodes, node_values) for x in x_values]

# Rysowanie wykresu
plt.scatter(nodes, node_values, color="red", label="Points")
plt.plot(x_values, y_values, color="blue", label="Interpolation polynomial")
plt.xticks(np.arange(nodes.min(), nodes.max() + abs(nodes[0] - nodes[1]), abs(nodes[0] - nodes[1])))
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.legend()
plt.title("Hormann-Floater Interpolation")
plt.savefig("wykres.png")
plt.show()