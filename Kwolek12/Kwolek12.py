import matplotlib.pyplot as plt
import numpy as np

# Formowanie funkcji
def function(x):
    return 1 / (1 + 5 * x * x)

# Zadane punkty funkcji
x = np.array([-7/8, -5/8, -3/8, -1/8, 1/8, 3/8, 5/8, 7/8])
y = np.array([function(xi) for xi in x])

# Interpolacja Floatera-Hormanna
def floater_hormann(x, y, d, x_plot):
    n = len(x)
    c = np.zeros((n, n))
    for i in range(n):
        c[i, 0] = y[i]
    for j in range(1, d + 1):
        for i in range(n - j):
            c[i, j] = ((x_plot - x[i + j]) * c[i, j - 1] - (x_plot - x[i]) * c[i + 1, j - 1]) / (x[i] - x[i + j])
    return c[0, d]

# Tworzenie punktów do wykresu
x_plot = np.linspace(-1, 1, 400)
y_plot = np.array([floater_hormann(x, y, 3, xi) for xi in x_plot])

# Rysowanie wykresu
plt.plot(x_plot, y_plot, label='Interpolacja Floatera-Hormanna')
plt.scatter(x, y, color='red', label='Punkty danych')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Interpolacja Floatera-Hormanna')
plt.grid(True)
plt.show()