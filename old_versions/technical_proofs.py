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


a11= sp.Symbol('a11')
a12= sp.Symbol('a12')
a13= sp.Symbol('a13')
a14= sp.Symbol('a14')

a21= sp.Symbol('a21')
a22= sp.Symbol('a22')
a23= sp.Symbol('a23')
a24= sp.Symbol('a24')

a31= sp.Symbol('a31')
a32= sp.Symbol('a32')
a33= sp.Symbol('a33')
a34= sp.Symbol('a34')

P1=sp.expand(x1-(a11*x4+a12*x4**2+a13*x4**3+a14*x4**4))
P2=sp.expand(x2-(a21*x4+a22*x4**2+a23*x4**3+a24*x4**4))
P3=sp.expand(x3-(a31*x4+a32*x4**2+a33*x4**3+a34*x4**4))
P=[P1,P2,P3]
Ball=[]
X=d.genvars(4)
S=[]
D=[]
for i in range(3):
	S.append(d.inBall(P[i],X)[0])
	D.append(d.inBall(P[i],X)[1])
Ball=S+D+[r3**2+r4**2-1]
Jac_ball=Matrix(Ball).jacobian(Matrix([*X[0],*X[1],X[2]]))	

R=Jac_ball.subs([(x1,0),(x2,0),(x3,0),(x4,0),(t,0),(a11,0),(a21,0)])
determin=R[[0,1,3,4,5,6],[0,1,3,4,5,6]].det()

pprint(factor(determin))