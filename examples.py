import draft as d
import function_version as fv
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
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')
X=[[x1,x2,x3,x4],[r3,r4],t]
#Defining the curve:
#P1=x2**3*x1**2+2*x2**2*x1**4+3*x3**2+ 4*x4**2+1
P1=x1+2*x1*x2+3*x1*x3+5*x1*x4
P2=x2-x4**3
P3=x3-x4
Ball1=d.inBall(P1,X)
Ball2=d.inBall(P2,X)
Ball3=d.inBall(P3,X)

U=[ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1)]
	  		
P=[d.poly_normal_to_list(Pi,X[0]) for Pi in [P1,P2,P3] ]
Jet_P_list=d.Jet_poly_list(P)

U_pluse=fv.F_Ballplus(U)
U_minus=fv.F_Ballminus(U)


Jet_P_func=[ {T: fv.poly_list_tofunc(Jet_P_list[i][T]) for T in Jet_P_list[i] } for i in range(3) ]

DP1=d.poly_normal_to_list(Ball1[1],X[0]+X[1]+[X[2]])


Ball=d.Ball_interval(P)



d.ftprint([ d.evaluation_poly_list(Pi,U) for Pi in Ball ])
d.ftprint(fv.Ball_system(Jet_P_func,U))





##########################################
#Comparing SP for func and the classical input##
##########################################
#classical input
"""x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')
X=[[x1,x2,x3,x4],[r3,r4],t]
#Defining the curve:
#P1=x2**3*x1**2+2*x2**2*x1**4+3*x3**2+ 4*x4**2+1
P1=x1-x4**2
P2=x2-x4**3
P3=x3-x4
Ball1=d.inBall(P1,X)
Ball2=d.inBall(P2,X)
Ball3=d.inBall(P3,X)

U=[ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0.5,0.5)]
	  		
P=[d.poly_normal_to_list(Pi,X[0]) for Pi in [P1,P2,P3] ]
Jet_P_list=d.Jet_poly_list(P)

U_pluse=fv.F_Ballplus(U)
U_minus=fv.F_Ballminus(U)


Jet_P_func={T: fv.poly_list_tofunc(Jet_P_list[0][T]) for T in Jet_P_list[0] }
SP1=d.poly_normal_to_list(Ball1[0],X[0]+X[1]+[X[2]])

print(Ball1)
d.ftprint([d.evaluation_poly_list(SP1,U)] )
d.ftprint([ 1/2 *( fv.poly_list_tofunc(Jet_P_list[0][(0,0,0,0)])(U_pluse) \
	   +fv.poly_list_tofunc(Jet_P_list[0][(0,0,0,0)])(U_minus) ) ])"""

#########################################
# general  test for  polynomials ... not ready yet 
"""x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')

X=[[x1,x2,x4],[r4],t]

#Defining the curve:
P1=x1-x4**2+1
P2=x2-x4**3+x4

U=[ft.arb(0.03,0.2),ft.arb(0.03,0.2),ft.arb(0.03,0.2),ft.arb(1.03,0.2),ft.arb(1.53,1)]
T1=d.inBall(P1,X)
T2=d.inBall(P2,X)



#changing the polynomials data to lists
P1=d.poly_normal_to_list(P1,X[0])
P2=d.poly_normal_to_list(P2,X[0])
P=[P1,P2]



func1=fv.poly_list_tofunc(P1)
func2=fv.poly_list_tofunc(P2)




func=lambda B : [func1(B),func2(B)]


jac_poly=d.jacobian_of_function_list(P)
jac_func=lambda B: [[fv.poly_list_tofunc(Pij)(B) for Pij in Pi ]  for Pi in jac_poly ]



H_P1=d.jacobian_of_function_list(jac_poly[0])
H_P1_func=[ [fv.poly_list_tofunc(partial_Pij)  for \
             partial_Pij  in H_P1[i]]  for i in range(len(H_P1)) ]  


H_P2=d.jacobian_of_function_list(jac_poly[1])
H_P2_func=[ [fv.poly_list_tofunc(partial_Pij)  for \
             partial_Pij  in H_P2[i]]  for i in range(len(H_P2)) ]  

H_func= lambda B: [ [[Pij(B) for Pij in Pi  ] for Pi in Hi  ] for Hi in [H_P1_func,H_P2_func] ]




#ball system
ball_func=fv.Ball_func(func,jac_func,U)

print(fv.jac_Ball(func,jac_func,H_func,U))
input()                  

jac_ball_func=lambda U: fv.complete_jac_ball(func,jac_func,H_func,U)




print(fv.func_solver(ball_func,jac_ball_func1,U))

"""
##################################
#simple example solver_func#######
##################################
"""x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
X=[x1,x2,x3]

P1=x1**2-x2**2
P2=x1**2+x2**2+x3**2-1
P3=x3-0.5

P1=d.poly_normal_to_list(P1,X)
P2=d.poly_normal_to_list(P2,X)
P3=d.poly_normal_to_list(P3,X)
P=[P1,P2,P3]
func= [fv.poly_list_tofunc(Pi)  for Pi in P]

jac_poly=d.jacobian_of_function_list(P)
jac_func=[[fv.poly_list_tofunc(Pij)  for Pij in Pi ]  for Pi in jac_poly ]



B=[ft.arb(0,3),ft.arb(0,3),ft.arb(0,3)]

for Ti in  fv.func_solver(func,jac_func,B):
	d.ftprint(Ti)

"""
##################################
#simple example solver_func#######
##################################
"""
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
X=[x1,x2]

P1=x1**2-x2**2
P2=x1**2+x2**2-1

P1=d.poly_normal_to_list(P1,X)
P2=d.poly_normal_to_list(P2,X)
P=[P1,P2]
func= [fv.poly_list_tofunc(Pi)  for Pi in P]

jac_poly=d.jacobian_of_function_list(P)
jac_func=[[fv.poly_list_tofunc(Pij)  for Pij in Pi ]  for Pi in jac_poly ]



B=[ft.arb(0,3),ft.arb(0,3)]

for Ti in  fv.func_solver(func,jac_func,B):
	d.ftprint(Ti)
"""
###########################################################################################
# Comparing the func version with the poly_list (not ready)
###########################################################################################

"""
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')

X=[[x1,x2,x3,x4],[r3,r4],t]

#Defining the curve:
P1=x1-x4**2+1
P2=x2-x4**3+x4
P3=x4-x3



#changing the polynomials data to lists
P1=d.poly_normal_to_list(P1,X[0])
P2=d.poly_normal_to_list(P2,X[0])
P3=d.poly_normal_to_list(P3,X[0])
P=[P1,P2,P3]

func1=fv.poly_list_tofunc(P1)
func2=fv.poly_list_tofunc(P2)
func3=fv.poly_list_tofunc(P3)
func=[func1,func2,func3]


jac_func=[ [fv.poly_list_tofunc(Pij) for Pij in d.jacobian_of_function_list(P)[i] ] for i in range(len(d.jacobian_of_function_list(P)))  ]

partial_P1=d.jacobian_of_function_list(P)[0]

Jac_funcP1=[ fv.poly_list_tofunc(Pi) for Pi in partial_P1 ]

H_P1=d.jacobian_of_function_list(partial_P1)
H_P1_func=[ [fv.poly_list_tofunc(partial_Pij)  for partial_Pij  in H_P1[i]]  for i in range(len(H_P1))      ]   

partial_P2=d.jacobian_of_function_list(P)[1]
H_P2=d.jacobian_of_function_list(partial_P2)
H_P2_func=[ [fv.poly_list_tofunc(partial_Pij)  for partial_Pij  in H_P2[i]]  for i in range(len(H_P2))  ]

partial_P3=d.jacobian_of_function_list(P)[2]
H_P3=d.jacobian_of_function_list(partial_P3)
H_P3_func=[ [fv.poly_list_tofunc(partial_Pij)  for partial_Pij  in H_P3[i]]  for i in range(len(H_P3))]  

B_Ball=[ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(2,1)]


B=[ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1)]

d.ftprint(fv.jac_Ball_of_onefunction(func1,Jac_funcP1,H_P1_func,B_Ball))
ballsystem=d.Ball_for_interval_poly(P[0],X)[0]
jacballsystem=d.jacobian_of_function_list([ballsystem])[0]

input()
d.ftprint([ d.evaluation_poly_list(P_i,B_Ball) for P_i in   jacballsystem ] )
input()

print(jacballsystem)

#for Ti in H_P1_func:
#	print(Ti[0](B))

#changing the polynomials to funcs """
###########################################################################################
# Example of finding a node using Ball system I manually computed B_Ball ##################
#because with big boxes, the algorithm does not stop in a realistic time (337.74 seconds) #
###########################################################################################
"""x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
X=[x1,x2,x3,x4]

#Defining the curve:
P1=x1-x4**2+1
P2=x2-x4**3+x4
P3=x4-x3

#changing the polynomials data to lists
P1=d.poly_normal_to_list(P1,X)
P2=d.poly_normal_to_list(P2,X)
P3=d.poly_normal_to_list(P3,X)
P=[P1,P2,P3]
B=[ft.arb(0.03,0.2),ft.arb(0.03,0.2),ft.arb(0.03,1.1),ft.arb(0.03,1.1)]

import math
# computing  the Ball system and B_Ball
Ball=d.Ball_interval(P)
B_Ball=d.B_Ball_calculator(B)

B_Ball=[ft.arb(0.003,1.01),ft.arb(0.003,1.01),ft.arb(0.003,1.01),ft.arb(0.003,1.01), ft.arb(1/(math.sqrt(2))+0.03,1.01),ft.arb(1/(math.sqrt(2)+0.03),1.01),ft.arb(2.03,1.1)]
jac_Ball= d.jacobian_of_function_list(Ball)   #the Jacobian of Ball
answer=d.solver2(Ball,jac_Ball,B_Ball)
for Ti in answer:
	d.ftprint(Ti) """
##############################################################################
# Example of finding a cusp using Ball system (27.3 seconds for this example)#
##############################################################################
"""x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
X=[x1,x2,x3,x4]

#Defining the curve:
P1=x1-x4**2
P2=x2-x4**3
P3=x4-x3

#changing the polynomials data to lists
P1=d.poly_normal_to_list(P1,X)
P2=d.poly_normal_to_list(P2,X)
P3=d.poly_normal_to_list(P3,X)
P=[P1,P2,P3]
B=[ft.arb(0.01,0.2),ft.arb(0.01,0.2),ft.arb(0.1,0.2),ft.arb(0.01,0.2)]


# computing  the Ball system and B_Ball
Ball=d.Ball_interval(P)
B_Ball=d.B_Ball_calculator(B)
jac_Ball= d.jacobian_of_function_list(Ball)   #the Jacobian of Ball

answer=d.solver2(Ball,jac_Ball,B_Ball)


for Ti in answer:
	d.ftprint(Ti)"""
##############################################
# Example of solver method ###################
############################################
"""x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
X=[x1,x2,x3]
P1=x1-x3**2
P2=x2-x3**3
P3=x1-x3

P1=d.poly_normal_to_list(P1,X)
P2=d.poly_normal_to_list(P2,X)
P3=d.poly_normal_to_list(P3,X)
P=[P1,P2,P3]
B=[ft.arb(0,3),ft.arb(0,3),ft.arb(0,3)]
x_teld=[0,0]
jac=d.jacobian_of_function_list(P)
answer=d.solver2(P,jac,B)
for Ti in answer:
	d.ftprint(Ti)  """                          
                                             
#############################################
##############################################
# Example of solver method (0.065 seconds)################
############################################
"""x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
X=[x1,x2]
P1=x1**2-x2**2
P2=x1**2+x2**2-1


P1=d.poly_normal_to_list(P1,X)
P2=d.poly_normal_to_list(P2,X)
P=[P1,P2]
B=[ft.arb(0,3),ft.arb(0,3)]
x_teld=[0,0]
jac=d.jacobian_of_function_list(P)
answer=d.solver2(P,jac,B)
for Ti in answer:
	d.ftprint(Ti) """ 


















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
