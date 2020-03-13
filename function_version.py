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
def eval_func(func,U):
  return func(U)
def func_matrixval(Mat,X): #Mat  is matrix and X is a list of ft.arb
         evaluation_jac_X=jac[:]
         for i in range(len(jac)):
             for j in range(len(jac[0])):
                 evaluation_jac_X[i][j]=eval_func(jac[i][j],X)
                 if evaluation_jac_X[i][j]==0:       #fixing a bug in flint since Python cannot create arb from type <class 'sympy.core.numbers.Zero'>
                     evaluation_jac_X[i][j]=ft.arb(0)
         return evaluation_jac_X

def invertability_jac_Ball(func,Jac_func,U): #still not ready 
  if 0 not in U[len(U)-1]:
    F_Ball=[F_Ballplus(U),F_Ballminus(U)]
    # checking whether a point is in L_c
    first_minor=d.i_minor(eval_jac,0)
    second_minor=d.i_minor(eval_jac,1)
    first_minor_q1=d.invertibility_of_a_matrix(func_matrixval(first_minor,F_Ball[0]))
    second_minor_q1= d.invertibility_of_a_matrix(func_matrixval(second_minor,F_Ball[0]))
    if  first_minor_q1!=1 or first_minor_q1 !=1:
      return -1
    else:
      first_minor_q2=d.invertibility_of_a_matrix(func_matrixval(first_minor,F_Ball[1]))
      second_minor_q2= d.invertibility_of_a_matrix(func_matrixval(second_minor,F_Ball[1]))
      if  first_minor_q2!=1 or first_minor_q2 !=1:
        return -1

      
def jac_Ball(func,Jac_func,H_func,U): #func is from R^n to R
  n=int((len(U)+1)/2)
  if 0 not in U[2*n-2]:
    answer=[]
    for j in range(n):
      answer.append(Ball_func(Jac_func[j],H_func[j])[0](U))
    for j in range(2,n):
      answer.append(Ball_func(Jac_func[j],H_func[j])[1](U)* U[2*n-2])
    D_P=[Ball_func(Jac_func[j],H_func[j])[1](U) for j in range(2,n)]
    sum=ft.arb(0)
    for i in range(2,n):
      sum += D_P[i]*U[i+n-2]
    answer.appen(sum)
  return answer    






