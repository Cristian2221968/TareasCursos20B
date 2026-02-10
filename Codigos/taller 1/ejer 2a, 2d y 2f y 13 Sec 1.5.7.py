import sympy as sp
import math
from sympy.vector import CoordSys3D, gradient, divergence, curl

x, y, z = sp.symbols('x y z')
N = CoordSys3D('N')

phi = sp.Function('phi')(x, y, z)
psi = sp.Function('psi')(x, y, z)

a1 = sp.Function('a1')(x, y, z)
a2 = sp.Function('a2')(x, y, z)
a3 = sp.Function('a3')(x, y, z)
a = a1*N.i + a2*N.j + a3*N.k

###Ejercicio 2a ###
lhs_a = gradient(phi*psi, N)
rhs_a = phi*gradient(psi, N) + psi*gradient(phi, N)
assert (lhs_a - rhs_a).to_matrix(N) == sp.Matrix([0,0,0])

###Ejercicio 2d ###

lhs_d = divergence(curl(a, N), N)
assert sp.simplify(lhs_d) == 0

###Ejercicio 2f ###

d = sp.diff

c1 = d(a3, y) - d(a2, z)
c2 = d(a1, z) - d(a3, x)
c3 = d(a2, x) - d(a1, y)

L1 = d(c3, y) - d(c2, z)
L2 = d(c1, z) - d(c3, x)
L3 = d(c2, x) - d(c1, y)

div_a = d(a1, x) + d(a2, y) + d(a3, z)
lap = lambda s: d(s, x, 2) + d(s, y, 2) + d(s, z, 2)
R1 = d(div_a, x) - lap(a1)
R2 = d(div_a, y) - lap(a2)
R3 = d(div_a, z) - lap(a3)

assert sp.simplify(L1-R1) == 0 and sp.simplify(L2-R2) == 0 and sp.simplify(L3-R3) == 0

###Ejercicio 13 ###

t = sp.symbols('t', real=True)
x, y = sp.cos(t), sp.sin(t)

Fx = -y/(x**2 + y**2)
Fy =  x/(x**2 + y**2)

dxdt, dydt = sp.diff(x,t), sp.diff(y,t)
integrand = sp.simplify(Fx*dxdt + Fy*dydt)  

Wcampo_a = sp.integrate(integrand, (t, 0, sp.pi))     
Wcontra_a = -Wcampo_a                                 

Wcampo_b = sp.integrate(integrand, (t, 0, -sp.pi)) 
Wcontra_b = -Wcampo_b                                  

print(Wcampo_a, Wcontra_a, Wcampo_b, Wcontra_b)