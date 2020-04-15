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
import operator
import random

def compose(f, g):
    return lambda x: f(g(x))
def composenbox(func,g):  # input of f is a list of n ft.arb and g is a list of 2n-1 ft.arb
  return lambda U:[funci(g(U)) for funci in func ]
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
              T[i].append(ft.arb(jac[i][j](X))) 
          return T

def sqrt_t(I):
  if 0 in I:
    return d.ftconstructor(0,math.sqrt(float(I.upper())))
  elif 0 < float(I.lower()):  
       sqrt1=math.sqrt(float(I.lower()))
       sqrt2=math.sqrt(float(I.upper()))
       return ft.arb(0.5*sqrt1+0.5*sqrt2,0.5*sqrt2-0.5*sqrt1)  #d.ftconstructor(math.sqrt(float(I.lower())),math.sqrt(float(I.upper())))
  else: 

    return -sqrt_t(-I)

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

  S_Pi= 0.5*Pi(F_Ballplus(U))+0.5*Pi(F_Ballminus(U))
  if 0 not in U[2*n-2]:

    D_Pi=(Pi(F_Ballplus(U))-Pi(F_Ballminus(U)))/(2*sqrt_t(U[2*n-2]))
  else:
    n=int((len(U)+1)/2)
    Id_n=list(np.eye(n, dtype=int))
    zero_function=lambda U:0
    nablaP=[JetPi[tuple(Ii)] if tuple(Ii) in JetPi else zero_function  for Ii in  Id_n ]  
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
    
def derivatives_of_SDPi(JetPi,U):

  n=int((len(U)+1)/2)
  Id_n=list(np.eye(n, dtype=int))
  Jet_nabla_Pi=[{}]*n
  for dervative_index in JetPi:
    for i in range(n):
      if dervative_index[i]>0:
        Jet_nabla_Pi[i]={**Jet_nabla_Pi[i], **{tuple(map(operator.sub,\
         dervative_index, Id_n[i])):JetPi[dervative_index] } }      
  #computing S_nablaPi and D_nablaPi
  S_nablaPi=[]
  D_nablaPi=[]

  for Jet_Pi_xj in Jet_nabla_Pi:
    if Jet_Pi_xj != {}:
      SD_Pi_xj=SD_Pi(Jet_Pi_xj,U)
    else:
       SD_Pi_xj=[ft.arb(0),ft.arb(0)]
    S_nablaPi.append(SD_Pi_xj[0])
    D_nablaPi.append(SD_Pi_xj[1])
  
 
  #computing the rows of the Jacobian of the Ball system 
  first_row=S_nablaPi[:]
  r=U[n:2*n-2]
  sum1=ft.arb(0)
  for i in range(2,n):
    first_row.append(d.intervals_multi(D_nablaPi[i],U[2*n-2]))
    sum1 +=d.intervals_multi(D_nablaPi[i],r[i-2])
  first_row.append(0.5*sum1)

  second_row=D_nablaPi+S_nablaPi[2:]
  if 0 not in U[2*n-2]:
    sum1=ft.arb(0)
    for i in range(2,n):
      sum1 +=d.intervals_multi(S_nablaPi[i],r[i-2])
    sum1 -= SD_Pi(JetPi,U)[1]
    second_row.append(d.intervals_multi(sum1, 0.5/(sqrt_t(U[2*n-2])) ))
  else:
    sum1=ft.arb(0)
    for i in range(2,n):
      for j in range(2,n):
        ri_rj=d.intervals_multi(r[i-2],r[j-2])
        for k in range(2,n):
          ri_rj_rk=d.intervals_multi(ri_rj,r[k-2])
          s=tuple([sum(x) for x in zip(tuple(Id_n[i]),tuple(Id_n[j]),tuple(Id_n[k]))])
          if s in JetPi:
            t=U[2*n-2]
            sqrt_of_t=extension_t(t)  #computing f'''(-sqrt(t),sqrt(t))
            U_prime=U[:2]
            for i in range(2,n):
              U_prime.append(U[i]+d.intervals_multi(r[i-2],sqrt_of_t) )
            sum1 +=d.intervals_multi(JetPi[s](U_prime),ri_rj_rk) 
    
    second_row.append(sum1/6)        
  return [first_row,second_row]

def extension_t(t):
  if 0 in t:
    sqrt_t=max( -float(t.lower()), float(t.upper()))
    return d.ftconstructor(-sqrt_t,sqrt_t )



def Jacobian_of_Ball(JetP,U):
  nabla_SP=[]
  nabla_DP=[]
  n=int((len(U)+1)/2)
  for i in range(n-1):

     partial_SDPi=derivatives_of_SDPi(JetP[i],U)

     
     nabla_SP.append(partial_SDPi[0])
     nabla_DP.append(partial_SDPi[1])
  Ball=nabla_SP+nabla_DP
  last_eq=[0]*n+[2*ri for ri in U[n:2*n-2]]+[0]
  Ball.append(last_eq)
  return Ball

def func_solver(P,jac,B,k=3): #k is the number of parts in which every interval is divided  
    it=0
    Solutions=[]
    L=[B]
    while len(L) !=0:
        
        it=it+1

        current_box=L[0]  #evaluating P at the currecnt_box:
        value_of_P_current_box=[ Pi(current_box) for Pi in P ]
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
             #L.remove(current_box)
             L=L[1:]
        
        else:
                   
                   jac_eval_current_box=jac(current_box)
                   mid_box=[ft.arb(float(interval.mid()))  for interval in current_box]
                   b=[Pi(mid_box) for Pi in P] 
                   if d.invertibility_of_a_matrix(jac_eval_current_box)==1:

                    Image_of_current_box=d.hansen_hengupta(mid_box,jac_eval_current_box,b,current_box,current_box)
                    #d.ftprint(current_box)
                    #print(Image_of_current_box)
                    #input()
                    
  
                    if Image_of_current_box !='empty':
                                    currecnt_box_contains_its_image=1
                                    
                                    for i  in range(len(Image_of_current_box)):
                                          if  (Image_of_current_box[i]).lower() <= (current_box[i]).lower() or (Image_of_current_box[i]).upper() >= (current_box[i]).upper():
                                                currecnt_box_contains_its_image=0

                                    if currecnt_box_contains_its_image==1 :
                                          Solutions.append(Image_of_current_box)    

                                          L=L[1:]
                          
                                    else:
                                           try:    # to check whether the intersection of current_box with its image (Image_of_current_box) is non empty
                                               for i in range(len(current_box)):
                                                   Intersection = current_box[i].intersection(Image_of_current_box[i])
                                               new_children=d.k_subdivide(current_box,k)

                                               #L.remove(current_box)
                                               L=L[1:]
                                               L =   L +new_children
                                           except:
                                               
                                               L=L[1:]
                    else:
                          L.remove(L[0])
                         

            
                   else:
                    new_children=d.k_subdivide(current_box,k)
                    #L.remove(current_box)
                    L=L[1:]
                    L =   L +new_children

                                               


    return Solutions  


def curve_tracer(P,B,jac,wth=0.001,wth2=0.01):  #returns all list of boxes that contains smooth parts of a curve and it stops if the curve is smooth
    list_of_boxes=[B]
    smoothness=1
    regular_boxes=[]
    smoothness=1
    while len(list_of_boxes)!=0 and d.width(list_of_boxes[0])> wth:
              membership=1
              eval_P=[Pi(list_of_boxes[0]) for Pi in P]           #checking whether the box contains a point of the curve
              for eval_Pi_at_B in eval_P:
                  if 0 not in ft.arb(1)*(eval_Pi_at_B):
                      membership=0
              if membership==0:
                list_of_boxes.remove(list_of_boxes[0])
                  #print("empty")
              else:
                   eval_jac= [[jacij(list_of_boxes[0]  ) for jacij in jaci] for \
                   jaci in jac] 
                   full_rankness= d.checking_full_rank(eval_jac)
                   if   full_rankness ==1 and d.width(list_of_boxes[0])<wth2 :

                        regular_boxes.append(list_of_boxes[0])
                        list_of_boxes.remove(list_of_boxes[0])
                   else:

                        new_children=d.subdivide(list_of_boxes[0])
                        list_of_boxes.remove(list_of_boxes[0])
                        list_of_boxes = list_of_boxes +new_children
    if len(list_of_boxes)!=0:
      smoothness=-1                    
      #print(list_of_boxes[0])
    return [regular_boxes,list_of_boxes]

def box_membership(B1,B2):
    comparing=8
    B1_contains_B2=1
    B2_contains_B1=1
    for i in range(len(B1)):

        if B2[i] not in B1[i]:
                B1_contains_B2=0
        if B1[i] not in B2[i]:
            B2_contains_B1=0
    if  B1_contains_B2==1 and   B2_contains_B1==1:
        comparing=3
    elif   B1_contains_B2==1:
        comparing=1
    elif   B2_contains_B1==1:
        comparing=2
    else:
        comparing=0
    return comparing


def jet_to_jac(Jet_P):
 jac=[]
 Id_n=list(np.eye(len(Jet_P)+1, dtype=int))
 for j in range(len(Jet_P)):
   jac.append([])
   for i in range(len(Jet_P)+1):
     jac[j].append( Jet_P[j][tuple(Id_n[i])] )
 return jac 
"""
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


  

"""

