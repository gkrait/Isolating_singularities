
import math
from sympy import *
from pprint import pprint
import flint  as ft



import draft as d
import function_version as fv
import matplotlib.pyplot as plt
import numpy as np
import sys

import sympy as sp
from copy import copy, deepcopy
import time

import inspect





x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')
X=[[x1,x2,x4],[r4],t]

P1=lambda B: d.power_interval(B[0]-8*ft.arb.cos(B[2]),2)+ d.power_interval(B[1]-8*ft.arb.sin(B[2]),2)-25
P2=lambda B:  d.power_interval(B[0]-9-5*ft.arb.cos(B[2]),2)+ d.power_interval(B[1]-5*ft.arb.sin(B[3]),2)-64
P3=lambda B: (B[0]-8*ft.arb.cos(B[2]))*ft.arb.sin(B[2])-(16*(B[1]-8*ft.arb.sin(B[2])))*ft.arb.cos(B[2])*  \
(B[0]-9-5*ft.arb.cos(B[3]))*ft.arb.sin(B[3])-(10*(B[1]-5*ft.arb.sin(B[3])))*ft.arb.cos(B[3])


P=[P1,P2,P3]

B=[ft.arb(0.2),ft.arb(0.2),ft.arb(0,1),ft.arb(0,1)]

P1_x1=lambda B: 2*B[0]-16*ft.arb.cos(B[2])
P1_x2=lambda B: 2*B[1]-16*ft.arb.sin(B[2])
P1_q1=lambda B: 16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) - 16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2])
P1_q2= lambda B :ft.arb(0)

P2_x1=lambda B:2*B[0]-18-10*ft.arb.cos(B[3])
P2_x2=lambda B: 2*B[1]-10*ft.arb.sin(B[3])
P2_q1=lambda B: ft.arb(0)
P2_q2=lambda B: 10*(B[0]-9-5*ft.arb.cos(B[3]))*ft.arb.sin(B[3])-(10*(B[1]-5*ft.arb.sin(B[3])))*ft.arb.cos(B[3])

B=[ft.arb(0.2),ft.arb(0.2),ft.arb(0,2),ft.arb(0,2)]


P3_x1=lambda B: 10*(16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) -\
 16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*((2*B[0] - 16*ft.arb.cos(B[2]))*\
 (2*B[1] - 10*ft.arb.sin(B[3]))-(2*B[1] - 16*ft.arb.sin(B[2]))*(2*B[0] - 10*ft.arb.cos(B[3]) - \
 18))*ft.arb.sin(B[3])+(16*(B[0]-8*ft.arb.cos(B[2]))*ft.arb.sin(B[2])-16*(B[1] -\
 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*(-10*(B[1]-5*ft.arb.sin(B[3]))*ft.arb.cos(B[3]) +\
 10*(B[0]-5*ft.arb.cos(B[3])-9)*ft.arb.sin(B[3]))*(32*ft.arb.sin(B[2])-\
 20*ft.arb.sin(B[3]))+16*((2*B[0]-16*ft.arb.cos(B[2]))*(2*B[1]-10*ft.arb.sin(B[3])) -\
 (2*B[1]-16*ft.arb.sin(B[2]))*(2*B[0]-10*ft.arb.cos(B[3])-18))*(-10*(B[1] - \
 	5*ft.arb.sin(B[3]))*ft.arb.cos(B[3])+10*(B[0]-5*ft.arb.cos(B[3])-9)*\
 ft.arb.sin(B[3]))*ft.arb.sin(B[2])


P3_x2=lambda B: -10*(16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) -\
16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*((2*B[0] - 16*ft.arb.cos(B[2]))\
*(2*B[1] - 10*ft.arb.sin(B[3])) - (2*B[1] - 16*ft.arb.sin(B[2]))*(2*B[0] - \
	10*ft.arb.cos(B[3]) - 18))*ft.arb.cos(B[3]) + (16*(B[0] - 8*ft.arb.cos(B[2]))*\
ft.arb.sin(B[2]) - 16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*(-10*(B[1] - \
	5*ft.arb.sin(B[3]))*ft.arb.cos(B[3]) + 10*(B[0] - 5*ft.arb.cos(B[3]) - 9)*\
ft.arb.sin(B[3]))*(-32*ft.arb.cos(B[2]) + 20*ft.arb.cos(B[3]) + 36) - 16*((2*B[0] -\
 16*ft.arb.cos(B[2]))*(2*B[1] - 10*ft.arb.sin(B[3])) - (2*B[1] - 16*ft.arb.sin(B[2]))*\
(2*B[0] - 10*ft.arb.cos(B[3]) - 18))*(-10*(B[1] - 5*ft.arb.sin(B[3]))*\
ft.arb.cos(B[3]) + 10*(B[0] - 5*ft.arb.cos(B[3]) - 9)*ft.arb.sin(B[3]))*ft.arb.cos(B[2])

P3_q1=lambda B:  (16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) - 16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*(-10*(B[1] - 5*ft.arb.sin(B[3]))*ft.arb.cos(B[3]) + 10*(B[0] - 5*ft.arb.cos(B[3]) - 9)*ft.arb.sin(B[3]))*(16*(2*B[1] - 10*ft.arb.sin(B[3]))*ft.arb.sin(B[2]) + 16*(2*B[0] - 10*ft.arb.cos(B[3]) - 18)*ft.arb.cos(B[2])) + ((2*B[0] - 16*ft.arb.cos(B[2]))*(2*B[1] - 10*ft.arb.sin(B[3])) - (2*B[1] - 16*ft.arb.sin(B[2]))*(2*B[0] - 10*ft.arb.cos(B[3]) - 18))*(-10*(B[1] - 5*ft.arb.sin(B[3]))*ft.arb.cos(B[3]) + 10*(B[0] - 5*ft.arb.cos(B[3]) - 9)*ft.arb.sin(B[3]))*((16*B[0] - 128*ft.arb.cos(B[2]))*ft.arb.cos(B[2]) - (-16*B[1] + 128*ft.arb.sin(B[2]))*ft.arb.sin(B[2]) + 128*d.power_interval(ft.arb.sin(B[2]),2) + 128*d.power_interval(ft.arb.cos(B[2]),2))


P3_q2=lambda B: (16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) - 16*(B[1] - \
 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*((2*B[0] - 16*ft.arb.cos(B[2]))*(2*B[1] -\
  10*ft.arb.sin(B[3])) - (2*B[1] - 16*ft.arb.sin(B[2]))*(2*B[0] - 10*ft.arb.cos(B[3])\
   - 18))*(-(-10*B[1] + 50*ft.arb.sin(B[3]))*ft.arb.sin(B[3]) + (10*B[0] - \
   50*ft.arb.cos(B[3]) - 90)*ft.arb.cos(B[3]) + 50*d.power_interval(ft.arb.sin(B[3]),2) +\
    50*ft.arb.cos(B[3])**2) + (16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) - \
    16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*(-10*(2*B[0] - \
    	16*ft.arb.cos(B[2]))*ft.arb.cos(B[3]) + 10*(-2*B[1] + 16*ft.arb.sin(B[2]))*\
    ft.arb.sin(B[3]))*(-10*(B[1] - 5*ft.arb.sin(B[3]))*ft.arb.cos(B[3]) + 10*(B[0] -\
     5*ft.arb.cos(B[3]) - 9)*ft.arb.sin(B[3]))

jac=[[P1_x1,P1_x2,P1_q1,P1_q2],[P2_x1,P2_x2,P2_q1,P2_q2],[P3_x1,P3_x2,P3_q1,P3_q2]]


print(fv.curve_tracer(P,B,jac,wth=0.0001,wth2=0.001))


############################################################
## Solver with analytic functions.. it works but slowly ####
############################################################
"""
P1=lambda B: B[0]- d.intervals_multi(ft.arb.cos(B[2]),(3+d.power_interval(ft.arb.sin(B[2]),4)))+3
P2=lambda B: B[1] - d.intervals_multi(d.power_interval(ft.arb.sin(B[2]),2),(3 + ft.arb.sin(8*B[2])))

P=[P1,P2]
B=[ft.arb(0),ft.arb(0),ft.arb(0)]

P11=lambda B: ft.arb(1)
P111=lambda B: ft.arb(0)




P14=lambda B:  4* d.intervals_multi( d.power_interval(ft.arb.sin(B[2]),3 ), \
 d.power_interval (ft.arb.cos(B[2]), 2) )\
  -d.intervals_multi( (3 + d.power_interval(ft.arb.sin(B[2]),4)), ft.arb.sin(B[2]))

P144=lambda B:  d.intervals_multi( cos(B[2]), (-13* d.power_interval( ft.arb.sin(B[2]),4) + \
	12* d.intervals_multi(d.power_interval(ft.arb.sin(B[2]),2), d.power_interval(ft.arb.cos(B[2]),2)) - 3) )

P1444=lambda B: d.intervals_multi( ft.arb.sin(B[2]), \
 (13* d.power_interval(ft.arb.sin(B[2]),4) + 24* d.power_interval(ft.arb.cos(B[2]),4)\
 -88*d.intervals_multi(d.power_interval(ft.arb.sin(B[2]),2), d.power_interval(ft.arb.cos(B[2]),2) )+ 3) )


P24=lambda B:2* d.intervals_multi( ft.arb.sin(B[2]), \
 (d.intervals_multi( (3 + ft.arb.sin(8*B[2])), ft.arb.cos(B[2]) ) \
 	+ 4* d.intervals_multi(ft.arb.sin(B[2]),  ft.arb.cos(8*B[2])))  )

P244= lambda B: -6* d.intervals_multi(d.power_interval(ft.arb.sin(B[2]),2), (11*ft.arb.sin(8*B[2]) + 1) )+\
2* d.intervals_multi(d.power_interval(ft.arb.cos(B[2]),2),(3+ft.arb.sin(8*B[2])) ) +\
32*d.intervals_multi(d.intervals_multi(ft.arb.cos(B[2]), \
	ft.arb.sin(B[2])), ft.arb.cos(8*B[2])   )

P2444= lambda B:-8*(-6*d.intervals_multi((d.power_interval(ft.arb.cos(B[2]),2)),ft.arb.cos(8*B[2])) +\
70*d.intervals_multi((d.power_interval(ft.arb.sin(B[2]),2)),ft.arb.cos(8*B[2])) +\
d.intervals_multi(d.intervals_multi(ft.arb.cos(B[2]),ft.arb.sin(B[2])), \
48*ft.arb.sin(B[2]) +3 ) )





JetP=[{(0,0,0):P1, (1,0,0): P11,(0,0,1): P14, (0,0,2):P144,(0,0,3):P1444  },\
{(0,0,0):P2, (0,1,0): P11,(0,0,1): P24, (0,0,2):P244,(0,0,3):P2444  } ]
 


U=[ft.arb(-1,1.1),ft.arb(1,1.1),ft.arb(0,1),ft.arb(1,0.000005),ft.arb(1,1.01)]
#U=[ft.arb(-1.5,1.51),ft.arb(1.5,1.51),ft.arb(0,2),\
#ft.arb(1.00000003,0.01),ft.arb(0,1.51)]


Ball=fv.Ball_system(JetP,U)
Jac= fv.Jacobian_of_Ball(JetP,U)


func_Ball=[lambda U1,i=i:fv.Ball_system(JetP,U1)[i] for i in range(len(Ball))]
func_Ball_alt=[lambda U1,i=i:fv.Alternative_Ball_system(JetP,U1)[i] for i in range(len(Ball))]
func_Jac=lambda U1:fv.Jacobian_of_Ball(JetP,U1)


T=fv.func_solver(func_Ball,func_Ball_alt,func_Jac,U,3)



print(T)

"""



########################################
### path tracer example by analytic functions &bout 7 minutes
##############################################
"""
P1=lambda B: B[0]- d.intervals_multi(ft.arb.cos(B[2]),(3+d.power_interval(ft.arb.sin(B[2]),4)))+3
P2=lambda B: B[1] - d.intervals_multi(d.power_interval(ft.arb.sin(B[2]),2),(3 + ft.arb.sin(8*B[2])))
P=[P1,P2]
#defining the Jacobian
P1_x=lambda B: ft.arb(1)
P1_xx=lambda B: ft.arb(0)
P1_z=lambda B:  4* d.intervals_multi( d.power_interval(ft.arb.sin(B[2]),3 ), \
 d.power_interval (ft.arb.cos(B[2]), 2) )\
  -d.intervals_multi( (3 + d.power_interval(ft.arb.sin(B[2]),4)), ft.arb.sin(B[2]))
P2_z=lambda B:2* d.intervals_multi( ft.arb.sin(B[2]), \
 (d.intervals_multi( (3 + ft.arb.sin(8*B[2])), ft.arb.cos(B[2]) ) \
 	+ 4* d.intervals_multi(ft.arb.sin(B[2]),  ft.arb.cos(8*B[2])))  )

B=[ft.arb(-2.5,2.5),ft.arb(2,2),ft.arb(0,5)]
jac_mat=[[P1_x,P1_xx,P1_z],[P1_xx,P1_x,P2_z]]
T=fv.curve_tracer(P,B,jac_mat)

projection=[Ti[:2] for Ti in T[0]]


fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlim(-5, 0)
ax.set_ylim(0,4)
ax.set_xlabel('x')
ax.set_ylabel('y')
c=0


for box in projection:
    rectangle= plt.Rectangle((float(box[0].lower()),float(box[1].lower()) ), \
    	float(box[0].upper())-float(box[0].lower()),float(box[1].upper())-float(box[1].lower()), fc='g')
    plt.gca().add_patch(rectangle)

plt.show()





"""
###########
#solver##########
"""
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')
X=[[x1,x2,x4],[r4],t]

P1=lambda B: B[0]- d.intervals_multi(ft.arb.cos(B[3]),(3+d.power_interval(ft.arb.sin(B[3]),4)))+3
P2=lambda B: B[1] - d.intervals_multi(d.power_interval(ft.arb.sin(B[3]),2),(3 + ft.arb.sin(8*B[3])))
P3=lambda B: B[2]- B[3]
P=[P1,P2,P3]
B=[ft.arb(0),ft.arb(0),ft.arb(0),ft.arb(0)]

P11=lambda B: ft.arb(1)
P111=lambda B: ft.arb(0)

P14=lambda B:  4* d.intervals_multi( d.power_interval(ft.arb.sin(B[3]),3 ), \
 d.power_interval (ft.arb.cos(B[3]), 2) )\
  -d.intervals_multi( (3 + d.power_interval(ft.arb.sin(B[3]),4)), ft.arb.sin(B[3]))

P144=lambda B:  d.intervals_multi( cos(B[3]), (-13* d.power_interval( ft.arb.sin(B[3]),4) + \
	12* d.intervals_multi(d.power_interval(ft.arb.sin(B[3]),2), d.power_interval(ft.arb.cos(B[3]),2)) - 3) )

P1444=lambda B: d.intervals_multi( ft.arb.sin(B[3]), \
 (13* d.power_interval(ft.arb.sin(B[3]),4) + 24* d.power_interval(ft.arb.cos(B[3]),4)\
 -88*d.intervals_multi(d.power_interval(ft.arb.sin(B[3]),2), d.power_interval(ft.arb.cos(B[3]),2) )+ 3) )


P24=lambda B:2* d.intervals_multi( ft.arb.sin(B[3]), \
 (d.intervals_multi( (3 + ft.arb.sin(8*B[3])), ft.arb.cos(B[3]) ) \
 	+ 4* d.intervals_multi(ft.arb.sin(B[3]),  ft.arb.cos(8*B[3])))  )

P244= lambda B: -6* d.intervals_multi(d.power_interval(ft.arb.sin(B[3]),2), (11*ft.arb.sin(8*B[3]) + 1) )+\
2* d.intervals_multi(d.power_interval(ft.arb.cos(B[3]),2),(3+ft.arb.sin(8*B[3])) ) +\
32*d.intervals_multi(d.intervals_multi(ft.arb.cos(B[3]), \
	ft.arb.sin(B[3])), ft.arb.cos(8*B[3])   )

P2444= lambda B:-8*(-6*d.intervals_multi((d.power_interval(ft.arb.cos(B[3]),2)),ft.arb.cos(8*B[3])) +\
70*d.intervals_multi((d.power_interval(ft.arb.sin(B[3]),2)),ft.arb.cos(8*B[3])) +\
d.intervals_multi(d.intervals_multi(ft.arb.cos(B[3]),ft.arb.sin(B[3])), \
48*ft.arb.sin(B[3]) +3 ) )


P34=lambda B: ft.arb(-1)


JetP=[{(0,0,0,0):P1, (1,0,0,0): P11,(0,0,0,1): P14, (0,0,0,2):P144,(0,0,0,3):P1444  },\
{(0,0,0,0):P2, (0,1,0,0): P11,(0,0,0,1): P24, (0,0,0,2):P244,(0,0,0,3):P2444  },\
{(0,0,0,0):P3, (0,0,1,0):P11, (0,0,0,1):P34} ]
 


U=[ft.arb(0.03,2),ft.arb(0.03,2),ft.arb(0.03,2),ft.arb(0.03,2),ft.arb(0.03,1),ft.arb(0.03,1),ft.arb(0.03,3)]
#U=[ft.arb(0.00003,0.001),ft.arb(0.00003,0.001),ft.arb(0.00003,0.001),\
#ft.arb(0.00003,0.001),ft.arb(0.707,0.001),ft.arb(0.707,0.001),ft.arb(0.000003,0.001)]




Ball=fv.Ball_system(JetP,U)
Jac=  fv.Jacobian_of_Ball(JetP,U)


func_Ball=[lambda U1,i=i:fv.Ball_system(JetP,U1)[i] for i in range(len(Ball))]

func_Jac=lambda U1:fv.Jacobian_of_Ball(JetP,U1)


T=fv.func_solver(func_Ball,func_Jac,U)
d.ftprint(T) 
"""
############################################
#draft#######################################
#######################################

"""
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')


P1=x3**3-x3*x1-x2
P2=3*x3**2-x1
P3=x4-x3

X=[[x1,x2,x3,x4],[r3,r4],t]
P=[d.poly_normal_to_list(Pi,X[0]) for Pi in [P1,P2,P3] ]

Ball= d.Ball_interval(P)
Ball_jacob= d.jacobian_of_function_list(Ball)	
sol=[0,0,0,0,0.7,0.7,0]
M=Matrix(d.matrixval(Ball_jacob,sol))
pprint(M)
#pprint(M[[0,2,3,4],0:4])
input()


B=[ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1)]


P1=lambda B: d.power_interval(B[0] -8* ft.arb.cos(B[2]),2) + d.power_interval(B[1] -8* ft.arb.sin(B[2]),2) -25
P2=lambda B: d.power_interval(B[0]-9 -5* ft.arb.cos(B[3]),2) + d.power_interval(B[1] -5* ft.arb.sin(B[3]),2) -64
P3=lambda B: 16*(B[0]-8*ft.arb.cos(B[2]))*ft.arb.sin(B[2])-\
(16*B[1]-128*ft.arb.sin(B[2]))*ft.arb.cos(B[2])*10*(B[0]-9-5*ft.arb.cos(B[3]))*ft.arb.sin(B[3])-(10*B[1]-50*ft.arb.sin(B[3]))*ft.arb.cos(B[3])
P=[P1,P2,P3]

partielP13= lambda B:16*(B[0]-8*ft.arb.cos(B[2]))*ft.arb.sin(B[2])-\
(16*B[1]-128*ft.arb.sin(B[2]))*ft.arb.cos(B[2])
partielP24= lambda B:10*(B[0]-9-5*ft.arb.cos(B[3]))*ft.arb.sin(B[3])-(10*B[1]-50*ft.arb.sin(B[3]))*ft.arb.cos(B[3])
partielP14=lambda B:ft.arb(0)

partielP11=2*B[0]-16*ft.arb.cos(B[2])

partielP12=2*B[1]-16*ft.arb.sin(B[2])

partielP21=2*B[0]-18-10*ft.arb.cos(B[3])

jac=[[partielP11,partielP12],[partielP11,partielP22]]

print(fv.curve_tracer(P,B,jac))



P=[]
"""
##############################
#curve tracer function example
##############################
"""
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')
X=[[x1,x2,x4],[r4],t]
P1=x1-x4**2
P2=x2-x4**3

P=[d.poly_normal_to_list(Pi,X[0]) for Pi in [P1,P2] ]
func_P=[ fv.poly_list_tofunc(Pi) for Pi in P  ]
Jet_P_list=d.Jet_poly_list(P)
Jet_P_list[1].update({(0,0,3):[ [[0,0,0],-6] ]})

B=[ft.arb(0,1),ft.arb(0,1),ft.arb(0,1)]
def jet_to_jac(Jet_P):
 jac=[]
 Id_n=list(np.eye(len(Jet_P)+1, dtype=int))
 for j in range(len(Jet_P)):
   jac.append([])
   for i in range(len(Jet_P)+1):
     jac[j].append(fv.poly_list_tofunc( Jet_P[j][tuple(Id_n[i])]) )
 return jac 

T=jet_to_jac(Jet_P_list) 

func_jac=lambda U: [ [Pij(U)  for Pij in Pi] for Pi in T   ]


d.ftprint(fv.curve_tracer(func_P,B,func_jac)[0])


"""


##############################################
#proving that the evaluation of  jac of ball #
#converges to the eval at the solution########
##############################################
"""
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')
X=[[x1,x2,x4],[r4],t]
P1=x1-x4**2
P2=x2-x4**3


print(d.inBall(P1,X))

n=1

U=[ft.arb(0.001,1/n),ft.arb(0.001,1/n),ft.arb(0.001,1/n),ft.arb(1.001,1/n),ft.arb(0.001,1/n)]	  		
P=[d.poly_normal_to_list(Pi,X[0]) for Pi in [P1,P2] ]
func_P=[ fv.poly_list_tofunc(Pi) for Pi in P  ]
Jet_P_list=d.Jet_poly_list(P)
Jet_P_list[1].update({(0,0,3):[ [[0,0,0],-6] ]})
Ball= d.Ball_interval(P)
Ball_jacob= d.jacobian_of_function_list(Ball)



Jet_P_func=[ {T: fv.poly_list_tofunc(Jet_P_list[i][T]) for T in Jet_P_list[i] } for i in range(2) ]

func_Ball_jacob= fv.Jacobian_of_Ball(Jet_P_func,U)

func_jac=lambda U1: fv.Jacobian_of_Ball(Jet_P_func,U1)
Ball_func=[]
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[0])
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[1])
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[2])
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[3])
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[4])



sol=[0,0,0,1,0]
eval_jac_ball=d.matrixval(Ball_jacob,sol)

m=0
for i in range(5):
	for j in range(5):
		#print(eval_jac_ball[i])
		#print( func_Ball_jacob[i])
		if m < ft.arb(func_Ball_jacob[i][j]).rad():
			m=ft.arb(func_Ball_jacob[i][j]).rad()
		if eval_jac_ball[i][j] not in ft.arb( func_Ball_jacob[i][j]) :
			print('false')
        		

d.ftprint(fv.func_solver(Ball_func,func_jac,U))

"""

################################################
# example of the solver where the input is the jet of P 
#################################################### 
"""
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')
r4= sp.Symbol('r4')
t= sp.Symbol('t')
X=[[x1,x2,x4],[r4],t]
P1=x1-x4**2
P2=x2-x4**3





U=[ft.arb(0.01,0.1),ft.arb(0.01,0.1),ft.arb(0.01,0.1),ft.arb(1.01,0.1),ft.arb(1.01,0.1)]	  		
P=[d.poly_normal_to_list(Pi,X[0]) for Pi in [P1,P2] ]
func_P=[ fv.poly_list_tofunc(Pi) for Pi in P  ]
Jet_P_list=d.Jet_poly_list(P)
Jet_P_list[1].update({(0,0,3):[ [[0,0,0],-6] ]})
Ball= d.Ball_interval(P)
Ball_jacob= d.jacobian_of_function_list(Ball)
eval_jac_ball=d.matrixval(Ball_jacob,U)


Jet_P_func=[ {T: fv.poly_list_tofunc(Jet_P_list[i][T]) for T in Jet_P_list[i] } for i in range(2) ]

func_Ball_jacob= fv.Jacobian_of_Ball(Jet_P_func,U)

func_jac=lambda U1: fv.Jacobian_of_Ball(Jet_P_func,U1)
Ball_func=[]
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[0])
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[1])
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[2])
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[3])
Ball_func.append(lambda U1: fv.Ball_system(Jet_P_func,U1)[4])




#d.ftprint(fv.func_solver(Ball_func,func_jac,U))
x=[ft.arb(0),ft.arb(0),ft.arb(0),ft.arb(1),ft.arb(0)]
b=[Balli(x) for Balli in Ball_func]

A=func_Ball_jacob= fv.Jacobian_of_Ball(Jet_P_func,x)


print(d.hansen_hengupta(x,A,b,x,x))
"""
##############################
#the solver with analytic maps
#################################

"""
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')

X=[x1,x2]
#Defining the curve:

P1=lambda B: ft.arb.sin(B[0]+B[1])
P2= lambda B: d.intervals_multi(B[0],ft.arb.exp(B[0]))
P=[P1,P2]
P11=lambda B:ft.arb.cos(B[0]+B[1])
P21=lambda B: d.intervals_multi(B[0],ft.arb.exp(B[0]))+ft.arb.exp(B[0])
P22=lambda B:0
jac=lambda U: [[P11(U),P11(U)],[P21(U),P22(U)]]
B=[ft.arb(0.1,1),ft.arb(0.1,1)]

for Ti in fv.func_solver(P,jac,B):
 d.ftprint(Ti)

"""
################################################
#the jacobian of Ball(P)   func & poly versions#
################################################
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
#P1=x2**3*x1**2+2*x2**2*x1**4+3*x3**2+ 4*x4**2+1
P1=x1-x4**2+1
P2=x2-x4**3+x4
P3=x3-x4

U=[ft.arb(0.1,1),ft.arb(0.1,1),ft.arb(0.1,1),ft.arb(0.1,1),ft.arb(0.1,1),ft.arb(0.1,1),ft.arb(0.5,0.5)]	  		
P=[d.poly_normal_to_list(Pi,X[0]) for Pi in [P1,P2,P3] ]
Jet_P_list=d.Jet_poly_list(P)
Jet_P_list[0].update({(0,0,0,3):[ [[0,0,0,0],-6] ]})
Jet_P_list[1].update({(0,0,0,3):[ [[0,0,0,0],-6] ]})
Ball= d.Ball_interval(P)
Ball_jacob= d.jacobian_of_function_list(Ball)
eval_jac_ball=d.matrixval(Ball_jacob,U)



Jet_P_func=[ {T: fv.poly_list_tofunc(Jet_P_list[i][T]) for T in Jet_P_list[i] } for i in range(3) ]
#d.ftprint(fv.derivatives_of_SDPi(Jet_P_func[0],U)[0])
#d.ftprint(fv.derivatives_of_SDPi(Jet_P_func[0],U)[1])

print('the func version:')
d.ftprint(fv.Jacobian_of_Ball(Jet_P_func,U))
print('the polynomial version:')
d.ftprint(eval_jac_ball)
"""
#########################################################
#Comparing the ball system when the in the input is a function
###############################################################
"""x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
x3= sp.Symbol('x3')
x4= sp.Symbol('x4')
r3= sp.Symbol('r3')""
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
"""
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
X=[x1,x2]

P1=x1**2-x2**2
P2=x1**2+x2**2-1


P1=d.poly_normal_to_list(P1,X)
P2=d.poly_normal_to_list(P2,X)

P=[P1,P2]
func= [fv.poly_list_tofunc(Pi)  for Pi in P]
B=[ft.arb(0,3),ft.arb(0,3)]
jac_poly=d.jacobian_of_function_list(P)
jac_func=lambda U:  [[d.evaluation_poly_list(jac_poly[0][0],U),\
 d.evaluation_poly_list(jac_poly[0][1],U)],\
 [d.evaluation_poly_list(jac_poly[1][0],U),d.evaluation_poly_list(jac_poly[1][1],U)]]

print(jac_func)
input()
d.ftprint(fv.func_solver(func,jac_func,B))


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
