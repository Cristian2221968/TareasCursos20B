import sympy as sp
import math
import numpy as np
import matplotlib.pyplot as plt

x = sp.Symbol('x')
a, b = -1, 1  

def baseVectores(n):
    base = []
    for i in range(n):
        base.append(x**i)
    return base 

def productoNormal(f, g):
    return sp.integrate(f*g, (x, a, b))

n = int(input("cantidad de vectores de la base"))

hFuncion = sp.sin(3*x)*(1-x**2)

PolyLegendre = [sp.legendre(k, x) for k in range(n)]

coeficientesLegendre = []

for i, pI in enumerate(PolyLegendre):
    c_i = productoNormal(hFuncion, pI) / productoNormal(pI, pI)
    coeficientesLegendre.append(c_i)
    print(f"Coeficiente c_{i} para el polinomio de Legendre: {c_i}")

hFuncionLengre = 0

for i in range(n): 
    hFuncionLengre += coeficientesLegendre[i] * PolyLegendre[i]
print("\nExpansión en Legendre:", sp.simplify(hFuncionLengre))

PolyChebysheu = [sp.chebyshevu(l, x) for l in range(n)]

coeficientesChebysheu = []

for j, pJ in enumerate(PolyChebysheu):
    c_j = productoNormal(hFuncion, pJ) / productoNormal(pJ, pJ)
    coeficientesChebysheu.append(c_j)
    print(f"Coeficiente c_{j} para el polinomio de Chebyshev: {c_j}")

hFuncionChebysheu = 0

for i in range(n): 
    hFuncionChebysheu += coeficientesChebysheu [i] * PolyChebysheu[i]
print("\nExpansión en Chebyshev:", sp.simplify(hFuncionChebysheu))

baseMono = baseVectores(n)

ecuacionGram = [] 
coeficientesMono = sp.symbols('k_0:n-1')

for j, bk in enumerate(baseMono):
    inc = 0
    for i, bi in enumerate(baseMono):
        inc+= coeficientesMono[i] * productoNormal(bi, bk)
        columIgual = productoNormal(hFuncion, bk)
    ecuacionGram.append(sp.Eq(inc, columIgual))

solucionMono = sp.solve(ecuacionGram, coeficientesMono)
print("Solución del sistema de ecuaciones de Gram:", solucionMono)

hFuncionMono = 0

for i in range(n): 
    hFuncionMono += solucionMono[coeficientesMono[i]] * baseMono[i]
print("\nExpansión en monomios:", sp.simplify(hFuncionMono))

f_exacta = sp.lambdify(x, hFuncion, 'numpy')
f_leg = sp.lambdify(x, sp.simplify(hFuncionLengre), 'numpy')
f_cheb = sp.lambdify(x, sp.simplify(hFuncionChebysheu), 'numpy')
f_mono = sp.lambdify(x, sp.simplify(hFuncionMono), 'numpy')
xx = np.linspace(-1,1,400)

plt.figure(figsize=(8,5))
plt.plot(xx, f_exacta(xx), 'k', linewidth=2, label='h(x) exacta')
plt.plot(xx, f_mono(xx), 'g--', linewidth=2, label='Expansión Monomios')
plt.plot(xx, f_cheb(xx), 'b-.', linewidth=2, label='Expansión Chebyshev')

ax = plt.gca()
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.title("Comparación: h(x) vs Monomios vs Chebyshev")
plt.xlabel("x")
plt.ylabel("h(x)")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()

plt.figure(figsize=(8,5))
plt.plot(xx, f_exacta(xx), 'k', linewidth=2, label='h(x) exacta')
plt.plot(xx, f_mono(xx), 'g--', linewidth=2, label='Expansión Monomios')
plt.plot(xx, f_leg(xx), 'r--', linewidth=2, label='Expansión Legendre')

ax = plt.gca()
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.title("Comparación: h(x) vs Monomios y Legendre")
plt.xlabel("x")
plt.ylabel("h(x)")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()

plt.figure(figsize=(8,5))
plt.plot(xx, f_exacta(xx), 'k', linewidth=2, label='h(x) exacta')
plt.plot(xx, f_leg(xx), 'r--', linewidth=2, label='Expansión Legendre')
plt.plot(xx, f_cheb(xx), 'b-.', linewidth=2, label='Expansión Chebyshev')

ax = plt.gca()
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.title("Comparación: h(x) vs Legendre y Chebyshev")
plt.xlabel("x")
plt.ylabel("h(x)")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.show()