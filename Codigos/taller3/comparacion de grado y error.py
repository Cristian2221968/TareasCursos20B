import sympy as sp
import math
import numpy as np

x = sp.Symbol('x')
a, b = -1, 1  

def baseVectores(n):
    base = []
    for i in range(n):
        base.append(x**i)
    return base 

def productoNormal(f, g):
    return sp.integrate(f*g, (x, a, b))

def errorDeAprox(f,g):
    error = sp.integrate((f - g)**2, (x, a, b))
    return math.sqrt(error)

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

for grado in range(1, n):
    # Truncar cada expansión al grado actual
    hcheb_trunc  = sum(coeficientesChebysheu[i]*PolyChebysheu[i] for i in range(grado))
    hLeg_trunc  = sum(coeficientesLegendre[i]*PolyLegendre[i] for i in range(grado))
    
    # Comparar diferencia
    diff = sp.simplify(hcheb_trunc - hLeg_trunc)
    
    if diff != 0:
        print(f"Las expansiones empiezan a diferir en grado {grado-1}")
        print("Diferencia =", diff)
        break
else:
    print("Coinciden hasta el grado", n)

for grado in range(1, n):
    
    hMono_trunc = sum(solucionMono[coeficientesMono[i]]*baseMono[i] for i in range(grado))
    hLeg_trunc  = sum(coeficientesLegendre[i]*PolyLegendre[i] for i in range(grado))
    
   
    diff = sp.simplify(hMono_trunc - hLeg_trunc)
    
    if diff != 0:
        print(f"Las expansiones empiezan a diferir en grado {grado-1}")
        print("Diferencia =", diff)
        break
else:
    print("Coinciden hasta el grado", n-1)

for grado in range(1, n):

    hMono_trunc = sum(solucionMono[coeficientesMono[i]]*baseMono[i] for i in range(grado))
    hcheb_trunc  = sum(coeficientesChebysheu[i]*PolyChebysheu[i] for i in range(grado))
    
   
    diff = sp.simplify(hMono_trunc - hcheb_trunc)
    
    if diff != 0:
        print(f"Las expansiones empiezan a diferir en grado {grado-1}")
        print("Diferencia =", diff)
        break
else:
    print("Coinciden hasta el grado", n-1)

errorLegendre = errorDeAprox(hFuncion, hFuncionLengre)
errorChebysheu = errorDeAprox(hFuncion, hFuncionChebysheu)
errorMono = errorDeAprox(hFuncion, hFuncionMono)
print(f"\nError de la aproximación con Legendre: {errorLegendre}")
print(f"Error de la aproximación con Chebysheu: {errorChebysheu}")
print(f"Error de la aproximación con Monomios: {errorMono}")