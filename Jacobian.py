import matplotlib.pyplot as plt
import numpy as np
import sympy
from sympy import *



x,y=symbols('x y')
vars=[x,y]
exp1=sympy.sympify("x**2+y**2-1")
exp2=sympy.sympify("x**2 - y")
f=[exp1,exp2]
J=sympy.zeros(len(f),len(vars))
for i,fi in enumerate(f):
    for j,s in enumerate(vars):
        J[i,j]=sympy.diff(fi,s)
J_inv=J.inv()

value_x=0.5
value_y=0.5
eps=0.0001
matrix_corn=Matrix([value_x,value_y])
matrix_func=Matrix([exp1.subs({x:value_x,y:value_y}),exp2.subs({x:value_x,y:value_y})])
#print(matrix_func)
count=0
while(abs(matrix_func[0])>eps and abs(matrix_func[1])>eps):

    pr=J_inv*matrix_func
    matrix_corn[0]=matrix_corn[0]-(pr[0].subs({x:value_x,y:value_y}))
    matrix_corn[1]=matrix_corn[1]-(pr[1].subs({x:value_x, y:value_y}))
    matrix_func[0]=exp1.subs({x:matrix_corn[0],y:matrix_corn[1]})
    matrix_func[1] = exp2.subs({x: matrix_corn[0], y: matrix_corn[1]})
    count+=1
    print("count ",count)
    print("корни ",matrix_corn)
    print("значение  ",matrix_func)


