import sympy as sp
import math

e1, e2, e3 = sp.Matrix([1,0,0]), sp.Matrix([0,1,0]), sp.Matrix([0,0,1])
# vectores
a = sp.Matrix([1,2,3])
b = sp.Matrix([4,5,6])
c = sp.Matrix([3,2,1])
d = sp.Matrix([6,5,4])

a_ei = sp.Matrix([ a.dot(e1), a.dot(e2), a.dot(e3)])
b_ei = sp.Matrix([ b.dot(e1), b.dot(e2), b.dot(e3)])
c_ei = sp.Matrix([ c.dot(e1), c.dot(e2), c.dot(e3)])
d_ei = sp.Matrix([ d.dot(e1), d.dot(e2), d.dot(e3)])

print (a_ei)
print (b_ei)
print (c_ei)
print (d_ei)

sum_1 = a_ei + b_ei + c_ei + d_ei
sum_2 = a_ei + b_ei - c_ei - d_ei
sum_3 = a_ei - b_ei + c_ei - d_ei
sum_4 = - a_ei + b_ei - c_ei + d_ei

print (sum_1)
print (sum_2)
print (sum_3)
print (sum_4)

vectors = {'a': a, 'b': b, 'c': c, 'd': d}

def angle(u,v):
    nu = sp.sqrt(u.dot(u))
    nv = sp.sqrt(v.dot(v))
    cosang = sp.simplify(u.dot(v)/(nu*nv))
    cosang_numeric = float(cosang.evalf())
    cosang_numeric = max(-1.0, min(1.0, cosang_numeric))
    deg = math.degrees(math.acos(cosang_numeric))
    return sp.simplify(cosang), deg

angles_with_basis = {}
for name, vec in vectors.items():
    angles_with_basis[name] = {}
    for bname, basis in [('e1',e1),('e2',e2),('e3',e3)]:
        cosv, deg = angle(vec,basis)
        angles_with_basis[name][bname] = {'cos(theta)': cosv, 'theta_deg': round(deg,6)}


magnitudes = {name: sp.sqrt(vec.dot(vec)) for name, vec in vectors.items()}

cos_ab, deg_ab = angle(a,b)
cos_cd, deg_cd = angle(c,d)

proj_a_on_b = (a.dot(b) / b.dot(b)) * b

M = sp.Matrix.hstack(a,b,c,d)
rank_vectors = M.rank()  # If <= 2 -> coplanar; if <=3 -> they lie in 3D space of course
are_coplanar = rank_vectors <= 2

dot_expr = (a+b).dot(c+d)

cross_ab = a.cross(b)
cross_bc = b.cross(c)
cross_cd = c.cross(d)

def angle_with_d(u):
    cosud, degud = angle(u,d)
    return {'vector': u, 'cos(theta)': sp.simplify(u.dot(d)/(sp.sqrt(u.dot(u))*sp.sqrt(d.dot(d)))), 'theta_deg': round(degud,6)}

angles_cross_with_d = {
    'a×b': angle_with_d(cross_ab),
    'b×c': angle_with_d(cross_bc),
    'c×d': angle_with_d(cross_cd)
}

scalar_triple = c.dot(a.cross(b))