import sympy as sp

# Definimos símbolos
a = sp.Symbol('a', positive=True)
m = sp.Symbol('m', positive=True)

# -------------------------
# (a) Cuadrado en el plano xy
# -------------------------

# Posiciones de los vértices
r1 = sp.Matrix([0, 0])
r2 = sp.Matrix([a, 0])
r3 = sp.Matrix([a, a])
r4 = sp.Matrix([0, a])

# Lista de posiciones
posiciones_cuadrado = [r1, r2, r3, r4]

# Centro de masa
r_cm_a = sum((m * r for r in posiciones_cuadrado),
             sp.Matrix([0, 0])) / (4 * m)


# Sustituimos a = 2
r_cm_a = r_cm_a.subs(a, 2)

print("Centro de masa (a):", r_cm_a)

# -------------------------
# (b) Cubo
# -------------------------

# Vértices inferiores (z = 0)
v1 = sp.Matrix([0, 0, 0])
v2 = sp.Matrix([a, 0, 0])
v3 = sp.Matrix([a, a, 0])
v4 = sp.Matrix([0, a, 0])

# Vértices superiores (z = a)
v5 = sp.Matrix([0, 0, a])
v6 = sp.Matrix([a, 0, a])
v7 = sp.Matrix([a, a, a])
v8 = sp.Matrix([0, a, a])

# Lista de posiciones
posiciones_cubo = [v1, v2, v3, v4, v5, v6, v7, v8]

# Centro de masa
r_cm_b = sum((m * r for r in posiciones_cuadrado),
             sp.Matrix([0, 0])) / (8 * m)


# Sustituimos a = 2
r_cm_b = r_cm_b.subs(a, 2)

print("Centro de masa (b):", r_cm_b)
