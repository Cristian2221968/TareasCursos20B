import sympy as sp
import math

x, y = sp.symbols('x y')
a, b = -1 , 1

def innprod(f,g, x):
    return sp.integrate(f * g , (x,a,b))
    
#bases iniciales del programa
baseEspa1 = [x**i for i in range(0,3)]
baseEspa2 = [y**i for i in range(0,3)]

print("las bases de los dos espacios vectoriales son ")
sp.pprint(baseEspa1)
print()
sp.pprint(baseEspa2)
print()

#base del tensor( producto tensorial)
baseTensor = []
for i in range(len(baseEspa1)):
    for j in range(len(baseEspa2)):

        baseTensor.append(baseEspa1[i] * baseEspa2[j])

print("base del tensor")
sp.pprint(baseTensor)
print()

#polinomios del inciso b
p_x = x**2 + x + 3
p_y = y + 1

#coeficientes del espacio 1 
coefEsp1 = []
for i in range(len(baseEspa1)):
    coeficiente = innprod(p_x,baseEspa1[i],x)/innprod(baseEspa1[i],baseEspa1[i],x)
    coefEsp1.append(coeficiente)

#coefientes del espacio 2
coefEsp2 = []
for j in range(len(baseEspa2)):
    coeficiente = innprod(p_y,baseEspa2[j],y)/innprod(baseEspa2[i],baseEspa2[i],y)
    coefEsp2.append(coeficiente)

#coeficientes del de tensor
coefTensor = []
for i in range(len(coefEsp1)):
    for j in range(len(coefEsp2)):
        coefTensor.append(coefEsp1[i] * coefEsp2[j]) 

print("coeficientes del tensor")
sp.pprint(coefTensor)
print()

#tensor expandido
ten_p_xy = 0
for i in range(len(baseTensor)):
    ten_p_xy += coefTensor[i]*baseTensor[i]

sp.pprint(f" tensor con los p(x) ={p_x} y g(y) ={p_y} ")
sp.pprint(ten_p_xy)
print()

#polinomios de legendre en P 
polyLegende_X = []
for v in baseEspa1:
    u = v
    for o in polyLegende_X:
        u -= (innprod(v, o, x ) / innprod(o, o, x )) * o
    polyLegende_X.append(sp.simplify(u))

print("polinomios de legendre ")
sp.pprint(polyLegende_X)
print()

#coeficentes del polinomio en base de legendre
coefP_xL = []
for l in range(len(polyLegende_X)):
    coeficienteL = innprod(p_x, polyLegende_X[l], x) / innprod(polyLegende_X[l], polyLegende_X[l], x)
    coefP_xL.append(coeficienteL)

#polinomio p(x) en base de legendre
p_xy_Le = 0
for i in range(len(polyLegende_X)):
    p_xy_Le += coefP_xL[i]*polyLegende_X[i]

print(f" polinomio p(x) expresado en la base de Legendre ")
sp.pprint(p_xy_Le)
coefP_yL = []
print()

#polinomios de legendre en g
polyLegende_y = []
for v in baseEspa2:
    u = v
    for o in polyLegende_y:
        u -= (innprod(v, o, y ) / innprod(o, o, y )) * o
    polyLegende_y.append(sp.simplify(u))

#polinomio g(y) en base de legendre
for l in range(len(polyLegende_y)):
    coeficienteL = innprod(p_y, polyLegende_y[l], y) / innprod(polyLegende_y[l], polyLegende_y[l], y)
    coefP_yL.append(coeficienteL)

#coeficientes de tensor en la base de los polinomios de legendre
coefTensorPo = []
for i in range(len(coefP_xL)):
    for j in range(len(coefP_yL)):
        coefTensorPo.append(coefP_xL[i] * coefP_yL[j]) 

print("componentes de tensor en base de legendre ")
sp.pprint(coefTensorPo)
print()

#base del tensor en la base de los polinomios de legrendre
baseTensor_Po = []
for i in range(len(polyLegende_X)):
    for j in range(len(polyLegende_y)):

        baseTensor_Po.append(polyLegende_X[i] * polyLegende_y[j])

print("base de tensor en base de legendre ")
sp.pprint(baseTensor_Po)
print()

#tensor expandido en la base de los polinomios de legrendre
ten_P_xy_Le = 0
for i in range(len(baseTensor_Po)):
    ten_P_xy_Le += coefTensorPo[i]*baseTensor_Po[i]

print("tensor expresado en la base de los polinomios de Lengendre")
sp.pprint(ten_P_xy_Le)