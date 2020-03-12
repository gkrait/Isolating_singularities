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



def compose(f, g):
    return lambda x: f(g(x))
def sqrt_t(I):
  if 0 in I:
    return d.ftconstructor(0,math.sqrt(float(I.upper())))
  else:
       return ftconstructor(math.sqrt(float(I.lower())),math.sqrt(float(I.upper())))
def F_Ballplus(U):    # len(U) is odd
  n=int((len(U)+1)/2)
  Y=U[:n]
  r_times_sqrtt=[yi *sqrt_t(U[2*n-2]) for yi in U[n:2*n-2] ]
  F1=U[:2]
  for i in range(2,n):
      F1.append(Y[i]+r_times_sqrtt[i-2])
  return F1
def F_Ballminus(U):    # len(U) is odd
  n=int((len(U)+1)/2)
  Y=U[:n]
  r_times_sqrtt=[yi *sqrt_t(U[2*n-2]) for yi in U[n:2*n-2] ]
  F2=U[:2]
  for i in range(2,n):
      F2.append(Y[i]-r_times_sqrtt[i-2])
  return F2          

def Ball_func(func,Jac_func): # func is a function that sends a list of intervals to an internal ....   Jac_func(i) is the pratial dervative of func wrt the i-variable 
   S_func=lambda U: 1/2*(compose(func,F_Ballplus)(U)+compose(func,F_Ballplus)(U))
   D_func= lambda U: 1/(2*sqrt_t(U[len(U)-1]))*(compose(func,F_Ballplus)(U)+compose(func,F_Ballplus)(U)) if 0 not \
    in U[len(U)-1] else sum([nabla_funci(U)*ri for Yi,ri in zip(Jac_func,[0,0]+U[int((len(U)+1)/2):len(U)-1])  ])
   return [S_func,D_func] 


def poly_list_tofunc(P):
  return lambda B: d.evaluation_poly_list(P,B)
