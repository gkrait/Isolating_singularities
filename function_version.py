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
def composenbox(func,g):  # input of f is a list of n ft.arb and g is a list of 2n-1 ft.arb
  return lambda U:[funci(g(U)) for funci in func ]


def sqrt_t(I):
  if 0 in I:
    return d.ftconstructor(0,math.sqrt(float(I.upper())))
  else:
       return d.ftconstructor(math.sqrt(float(I.lower())),math.sqrt(float(I.upper())))
def F_Ballplus(U):    # len(U) is odd
  n=int((len(U)+1)/2)
  Y=U[:n]
  r_times_sqrtt=[d.intervals_multi(yi ,sqrt_t(U[2*n-2])) for yi in U[n:2*n-2] ]
  F1=U[:2]
  for i in range(2,n):
      F1.append(Y[i]+r_times_sqrtt[i-2])
  return F1
def F_Ballminus(U):    # len(U) is odd
  n=int((len(U)+1)/2)
  Y=U[:n]
  r_times_sqrtt=[d.intervals_multi(yi ,sqrt_t(U[2*n-2])) for yi in U[n:2*n-2] ]
  F2=U[:2]
  for i in range(2,n):
      F2.append(Y[i]-r_times_sqrtt[i-2])
  return F2          


def poly_list_tofunc(P):
  return lambda B: d.evaluation_poly_list(P,B)

def SD_Pi(JetPi,U): # JetPi (dictionary) is the jet of one function Pi
  n=int((len(U)+1)/2)
  Pi=JetPi[(0,)*n]
  S_Pi= 1/2*(Pi(F_Ballplus(U))+Pi(F_Ballminus(U)))
  if 0 not in U[2*n-2]:
    D_Pi=(Pi(F_Ballplus(U))-Pi(F_Ballminus(U)))/(2*sqrt_t(U[2*n-2]))
  else:
    n=int((len(U)+1)/2)
    Id_n=list(np.eye(n, dtype=int))
    nablaP=[  JetPi[tuple(Ii)] for Ii in  Id_n ]  
    D_Pi=ft.arb(0)
    for i in range(len(Id_n)-2):
      D_Pi += d.intervals_multi(nablaP[i+2](U),U[n+i])

  return [S_Pi, D_Pi]
def Ball_system(JetP,U):
  S_func=[]
  D_func=[]
  n=int((len(U)+1)/2)
  for i in range(n-1):
     SDPi=SD_Pi(JetP[i],U)
     S_func.append(SDPi[0])
     D_func.append(SDPi[1])
  Ball=S_func+ D_func
  last_eq=sum( [d.power_interval(Ui,2) for Ui in U[n:2*n-2] ] )-ft.arb(1) 
  Ball.append(last_eq)
  return Ball
    


"""
def Ball_func(func,Jac_func,U): # func is a function that sends a list of intervals to an internal ....   Jac_func(i) is the pratial dervative of func wrt the i-variable 
 S_func=[]
 D_func=[]
 n=int((len(U)+1)/2 )
 for i in range(len(F_Ballplus(U))-1):

   S_func.append (1/2*(func(F_Ballplus(U))[i]+func(F_Ballminus(U))[i]) )
   if 0 not in U[len(U)-1] :
    D_func.append(1/(2*sqrt_t(U[len(U)-1]))*(func(F_Ballplus(U))[i]+func(F_Ballminus(U))[i]))
   else:
    D_func.append( sum([d.intervals_multi(nabla_funci,ri) for \
      nabla_funci,ri in zip(Jac_func(U)[i],[0,0]+U[int((len(U)+1)/2):len(U)-1])  ]) )
 last_eq=sum( [d.intervals_multi(Ui,Ui) for Ui in U[n:2*n-2] ] )-1
 return S_func+D_func+[last_eq]





def complete_jac_ball(func,jac_func,H_func,U): 
  jac_ball=[]
  SP=[]
  DP=[]
  for i in range(len(func)):
     SP.append(jac_Ball_of_onefunction(func[i],jac_func[i],H_func[i],U)[0] )
     DP.append(jac_Ball_of_onefunction(func[i],jac_func[i],H_func[i],U)[1] )
  
  jac_ball=SP+DP
  last_eq= lambda U:  [0]*int((len(U)+1)/2) + U[ int((len(U)+1)/2): len(U)-1  ] +[0] 
  jac_ball.append(last_eq(U))
  return jac_ball

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

      #the following method is not ready
def jac_Ball(func,Jac_func,H_func,U): #func is from R^n to R.. Jac_func is the Jacobian of func 
  write again this 
  n=int((len(U)+1)/2)
  answer=[]
  #Computing the partial derivatives of S.func (see the proof of Lemma 55)
  for j in range(n-1):
    nabla_SD_Pi=Ball_func(lambda U1: Jac_func(U1)[j],lambda U1: H_func(U1)[j],U )[:2*n-2]
    D_nblaPi=nabla_SD_Pi[n:]
    for P in nabla_SD_Pi[n:]:
      P=d.intervals_multi(P,U[2*n-2])
    sum1=ft.arb(0) 
    for i in range(2,n):
        sum1 += d.intervals_multi(D_nblaPi[i-2],U[i+n-2])
    nabla_SD_Pi.append(sum1) 
    answer.append(nabla_SD_Pi)
  ############################################
  #Computing the partial derivatives of D.func 
  for j in range(n-1):
     nabla_SD_Pi=Ball_func(lambda U1: Jac_func(U1)[j],lambda U1: H_func(U1)[j],U )[:2*n-2]
     D_nblaPj=nabla_SD_Pi[n:]+nabla_SD_Pi[]
     sum1=ft.arb(0) 
     if 0 not in U[2*n-2]:
       T=sum([d.intervals_multi(xi,yi) for xi, yi  in zip(nabla_SD_Pi[2:n], U[n:2*n-2] )  ]   )
       T=d.intervals_multi(T,1/(2*U[2*n-2]))
       nbla_DPi.append(T)
     answer.append(nbla_DPi)  


  for j in range(2,n):
    answer.append(d.intervals_multi(Ball_func(Jac_func[j],H_func[j])[1](U), U[2*n-2]))
    D_P=[Ball_func(Jac_func[j],H_func[j])[1](U) for j in range(2,n)]
  #for the derivative wrt t: 
  ############################################
  #Computing the partial derivatives of D.func    



def func_matrixval(jac,X): #jac as i_minor and X is a list of ft.arb
         #for i in range(len(jac)):
          #   for j in range(len(jac[0])):

          #       evaluation_jac_X[i][j]=jac[i][j](X)
          #       if evaluation_jac_X[i][j]==0:       #fixing a bug in flint since Python cannot create arb from type <class 'sympy.core.numbers.Zero'>
          #           evaluation_jac_X[i][j]=ft.arb(0)
         
          T=[]
          for i in range(len(jac)):
             T.append([])
             for j in range(len(jac[i])):
              T[i].append(jac[i][j](X))
          return T

def func_solver(P,jac,B,k=2): #k is the number of parts in which every interval is divided  

    it=0
    Solutions=[]
    L=[B]
    while len(L) !=0:
        it=it+1
        current_box=L[0]  #evaluating P at the currecnt_box:
        value_of_P_current_box= [Pi(current_box) for Pi in P]
        solution_in_current_box=1
        for value in value_of_P_current_box:  #checking whether the box has no solution
            try:                    # notice that if value ==0 then a TypeError occurs
               if 0 not in (ft.arb(1)*(value)):
                 solution_in_current_box=0
                 break
            except TypeError:
                 if value !=0:
                      solution_in_current_box=0
                      break            
        if solution_in_current_box==0:
             L.remove(current_box)
        else:
                   jac_eval_current_box=func_matrixval(jac,current_box)
                   mid_box=[ft.arb(float(interval.mid()))  for interval in current_box]
                   b=[Pi(mid_box) for Pi in P]
                   print(jac[0][0](current_box))
                   d.ftprint(current_box)
                   d.ftprint(value_of_P_current_box)
                   pprint(jac_eval_current_box)
                   input()
                   if d.invertibility_of_a_matrix(jac_eval_current_box)==1:
                    Image_of_current_box=d.hansen_hengupta(mid_box,jac_eval_current_box,b,current_box,current_box)
                    if Image_of_current_box !='empty':
                                    currecnt_box_contains_its_image=1
                                    for i  in range(len(Image_of_current_box)):
                                          if  (Image_of_current_box[i]).lower() <= (current_box[i]).lower() or (Image_of_current_box[i]).upper() >= (current_box[i]).upper():
                                                currecnt_box_contains_its_image=0

                                    if currecnt_box_contains_its_image==1 :
                                          #print('Hi',Solutions)
                                          Solutions.append(Image_of_current_box)    
                                          L.remove(current_box)
                                    else:
                                           try:    # to check whether the intersection of current_box with its image (Image_of_current_box) is non empty
                                               for i in range(len(current_box)):
                                                   Intersection = current_box[i].intersection(Image_of_current_box[i])
                                               new_children=d.k_subdivide(current_box,k)

                                               L.remove(current_box)
                                               L =   L +new_children
                                           except:
                                               L.remove(current_box)
                    else:
                          L.remove(L[0])

            
                   else:
                    new_children=d.k_subdivide(current_box,k)
                    L.remove(current_box)
                    L =   L +new_children

                                               

    #print(Solutions)
    return Solutions    

"""

