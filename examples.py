import draft as d
import matplotlib.pyplot as plt
import numpy as np
import flint  as ft
import sympy as sp
from copy import copy, deepcopy
import sys
import time
from sympy import *
import inspect
import math
from pprint import pprint

x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
t1= sp.Symbol('t1')
t2= sp.Symbol('t2')
q1= sp.Symbol('q1')
q2= sp.Symbol('q2')
q3= sp.Symbol('q3')
X=[x1,x2]

##############################################
# Example of solver method 
############################################
P1=x1**2-x2**2
P2=x1**2+x2**2-1


P1=d.poly_normal_to_list(P1,X)
P2=d.poly_normal_to_list(P2,X)
P=[P1,P2]

B=[ft.arb(0.1,2),ft.arb(0.1,2)]
B=[ft.arb(-0.7,0.2),ft.arb(-0.7,0.3)]
x_teld=[-0.7,-0.7]


jac=d.jacobian_of_function_list(P)
jac=d.matrixval(jac,B)
b=[d.evaluation_poly_list(Pi,B) for Pi in P]


T=d.hansen_hengupta(x_teld,jac,b,B,B)
d.ftprint(T)
#d.ftprint(d.hansen_hengupta(x_teld,jac,b,B,B))




#hansen_hengupta and then come back to solver
	
a=ft.arb(1,1)
b=ft.arb(0,1)
x=ft.arb(0,5)

#d.ftprint([d.gauss_seidel_dim1(a,b,x)])






















"""
# Example of the  paper 
a1=8
a2=5
a3=8
a4=5  
b=9

P1=sp.expand((a1**2+2*a1*x1-a2**2+x1**2+x2**2)*t1**2-4*a1*t1*x2+a1**2-2*a1*x1-a2**2+x1**2+x2**2)
P2=sp.expand((a4**2+(-2*b+2*x1)*a4-a3**3+b**2-2*b*x1+x1**2+x2**2)*t2**2-4*a4*t2*x2+a4**2+a4*(2*b-2*x1)-a3**3+b**2-2*b*x1+x1**2+x2**2)
func=Matrix([P1,P2])
var=Matrix(X)
M=func.jacobian(var)

P3= sp.expand((M[:,2:]).det()* (M[:,0:2]).det()) #P3 is the pruduct of minors 



P=[draft.poly_normal_to_list(P1,X),draft.poly_normal_to_list(P2,X),draft.poly_normal_to_list(P3,X)]

B=[ft.arb(0,20),ft.arb(0,20),ft.arb(0,1),ft.arb(0,1)]
jac=draft.jacobian_of_function_list(P)


T= draft.curve_tracer(P,B,jac,wth=0.01)
print(len(T[0]))
print(len(T[1]))

TY=[[[ round(float(Tij.lower()),5),round(float(Tij.upper()) ,5) ]  for Tij in Ti]  for Ti in T[0]]


projection=[ Ti[:2] for Ti in TY ]
#projection=[[Ti[0],Ti[1]] for Ti in T[0]]
  

fig, ax = plt.subplots()
plt.grid(True)
ax.set_ylim(-50, 50)
ax.set_xlim(-50, 50)
ax.set_xlabel('x')
ax.set_ylabel('y')
c=0

for box in projection:
    c+=1
    rectangle= plt.Rectangle((box[0][0], box[1][0] ), box[0][1]-box[0][0],box[1][1]-box[1][0], fc='g')
    plt.gca().add_patch(rectangle)
    #rectangle= plt.Rectangle((round(float(box[0].lower()),3), round(float(box[1].lower())),3), 2*float(box[0].rad()),2*float(box[1].rad()), fc='g')
    #plt.gca().add_patch(rectangle)

plt.show()
"""








#Example 2
"""
P1=x1**2+x2**2-q1**2
P2=sp.expand(simplify(sp.expand((x1 + 17*((-t**2+1)*(t**2+1)**(-1))-15.9)**2 + (x2 + 17*((2*t)*(t**2+1)**(-1)) )**2-q2**2)*(t**2+1)))
P3=sp.expand(simplify((t**2+1)**3* sp.expand((x1 + 20.8*((-t**2+1)*cos(.8822)*(t**2+1)**(-1)-2*t*sin(.8822)*(t**2+1)**(-1)) )**2 + (x2 + 20.8*(2*t*cos(0.8822)*(t**2+1)**(-1)+(-t**2+1)*sin(0.8822)*(t**2+1)**(-1))- 10)**2 -q3**2  )))
P5=t-1
#For finding P4(the aspects boundary equation, we do the following)

func=Matrix([P1,P2,P3])
var=Matrix(X)
M=func.jacobian(var)
P31= (M[:,:3]).det()
P32=(M[:,3:]).det()
P4=sp.expand(P32*P31)
"""





#Example 2 case x3=0
"""
P1=x1**2+x2**2-q1**2
P2=sp.expand((x1 + 1.1)**2 + x2**2 - q2**2)
P3=sp.expand((x1 + 20.8*cos(0.8822))**2 + (x2 + 20.8*sin(0.8822) - 10)**2 - q3**2)


#For finding P4(the aspects boundary equation, we do the following)

func=Matrix([P1,P2,P3])
var=Matrix(X)
M=func.jacobian(var)
minor02=(M[:,[1,3,4]]).det() #equal to 8*q2*q3*x2, but q1, q2 are never zeros
minor03= (M[:,[1,2,4]]).det()   #equal to -8*q1*q3*x2    q1, q3 are never zeros
minor04= (M[:,[1,2,3]]).det()   #equal to 4*q1*q2*(2*x2 + 12.1209717960827)  

minor12=(M[:,[0,3,4]]).det() #equal to 8*q2*q3*x1, but q1, q2 are never zeros
minor13=(M[:,[0,2,4]]).det() #equal to -4*q1*q3*(2*x1 + 2.2)
minor14=(M[:,[0,2,4]]).det() #equal to -4*q1*q3*(2*x1 + 2.2)
"""
