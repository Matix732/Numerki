import numpy as np
import matplotlib.pyplot as plt

# Definiujemy funkcję podcałkową
def integrand(t):
    return np.cos((1 + t) / (t**2 + 0.04)) * np.exp(-t**2)

# Własna funkcja do całkowania metodą trapezów
# n - liczba podprzedziałów (większa -> większa dokładność)
def integrate_trapezoid(f, a, b, n=10000):
    x = np.linspace(a, b, n)
    dx = x[1] - x[0]
    y = f(x)
    return np.sum((y[:-1] + y[1:]) / 2 * dx)

# Definiujemy funkcję F(x) jako całkę od -∞ do x
def F(x):
    lower_limit = -100 # Przybliżenie -nieskończoności
    return integrate_trapezoid(integrand, lower_limit, x)

# Generowanie punktów do wykresu
x_values = np.linspace(-10, 10, 400)
y_values = [F(x) for x in x_values]

# Rysowanie wykresy
plt.plot(x_values, y_values, label='y')
plt.xlabel('x')
plt.ylabel('y')
plt.title('F(x) graph')
plt.legend()
plt.grid(True)
plt.show()

# Obliczamy granicę lim x->inf F(x)
# Używamy dużego x jako przybliżenie nieskończoności
upper_limit = 100 
lower_limit = -100  
limit = integrate_trapezoid(integrand, lower_limit, upper_limit)
print(f"Granica lim x->inf F(x) wynosi: {limit:.8f}")