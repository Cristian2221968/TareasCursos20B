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

def productoCheb(f, g):
    return sp.integrate(f*g*sp.sqrt(1 - x**2), (x, a, b))

n = int(input("cantidad de vectores de la base"))

baseComprobar = baseVectores(n)
for i in range(n):
    for j in range(i+1, n):
            comprobacion = productoNormal(baseComprobar[i],baseComprobar[j])
            if comprobacion == 0:
                print(f"Los vectores son ortogonales{i}: {baseVectores(n)[i]} y {baseVectores(n)[i+1]}") 
            else:
                print(f"Los vectores no son ortogonales{i}: {baseVectores(n)[i]} y {baseVectores(n)[i+1]}")
print("Base inicial:", baseComprobar)
baseAOrtogonalizar = baseVectores(n) 

print("Base inicial:", baseAOrtogonalizar)

ortogonalizados = []
for v in baseAOrtogonalizar:
    u = v
    for o in ortogonalizados:
        u -= (productoNormal(v, o) / productoNormal(o, o)) * o
    ortogonalizados.append(sp.simplify(u))
    
print("Base ortogonalizada:")
for i, poly in enumerate (ortogonalizados, 1):
    print(f"Polinomio de Legendre {i}: {poly}")

ortogonalizadosCheb = []
for v in baseAOrtogonalizar:
    u = v
    for o in ortogonalizadosCheb:
        num = productoCheb(v, o)
        den = productoCheb(o, o)
        print(f"Integrando: <{v},{o}> = {num}, <{o},{o}> = {den}") 
        if den != 0:
            u -= (num/den) * o
    ortogonalizadosCheb.append(sp.simplify(u))
    
print("Base ortogonalizada por Cheb:")
for i, poly in enumerate (ortogonalizadosCheb, 1):
    print(f"Polinomio de Chebyshey {i}: {poly}")

