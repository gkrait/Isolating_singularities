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
import functools


def sevr_mul(*arg):
  answer=ft.arb(1)
  for a in arg:
   answer *= a
  return answer

def sevr_add(*arg):
  answer=ft.arb(0)
  for a in arg:
   answer += a
  return answer

class interv(ft.arb):
  def __pow__(self,n): 
   
   answer=interv(1)
   for i in range(n):
    answer *= self
    return interv(answer)
  def __str__(self):
    return  str([float(self.lower()),float(self.upper())])
  def __rmul__(self,other):
    intervals_multi(self ,other)

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno
def abs_value(interval):
    a=abs((interval.lower()))
    b=abs((interval.upper()))
    c=max(a,b)
    return ft.arb(c/2, c/2)
def gershgorin_circles(M):   #returns Gershgorin circles of an intrival matrix M
    radius= [0]*(len(M))
    for i in range(len(M)):
        for j in range(len(M)):
            if i !=j:
              radius[i]=radius[i]+ abs_value(M[i][j])
    T=[]   #Now we have that radius is the list of the radiuses of  The gershgorin circles
    i=0

    for i in range(len(radius)):
         m_lower=float((M[i][i]).lower())
         m_upper=float((M[i][i]).upper())
         r_upper=float((radius[i]).upper())
         gershgorin_lower=m_upper+r_upper
         gershgorin_upper =m_lower-r_upper
         T.append(ft.arb(0.5*(gershgorin_upper+gershgorin_lower),0.5*(gershgorin_upper-gershgorin_lower) ))
    return T



def gauss_seidel_dim1p(a,b,x):
    return b/a
def eliminate_the_last_variable(M):
    N=M
    if type(M[0])!= list:
        N.remove(M[len(N)-1])
    else:
        N=[eliminate_the_last_variable(Mi) for Mi in M]
    return N # No use so far

    #Evaluating polynomials \not used yet

def empty_list_generator(m,n):   #generates a n-dim list of m zeros
    M=[0]*(m+1)
    N=M
    if n>1:
      for j in range(len(M)):
               N[j]=empty_list_generator(m,n-1)

      M=N

    return M
def replacing_an_elements_of_list(the_index,the_list,x):
     if len(the_index)==1:
         the_list[the_index[0]]=x

     else:
          the_new_index=the_index[1:len(the_index)]
          the_new_list=the_list[the_index[0]]
          replacing_an_elements_of_list(the_new_index,the_new_list,x)
def reversing_list(M):
    N=M
    if type(N[0])==type([]):
        for T in N:
            T=(reversing_list(T))


    return(N)
def coefficient_matrix(f):
       non_zero_terms= f.terms()
       max_degree=0
       for T in non_zero_terms:
           if max_degree < max(list(T[0])):
               max_degree=max(list(T[0]))
       M=empty_list_generator(max_degree, len(list(T[0])))
       for T in non_zero_terms:
           N=list(T[0])
           N.reverse()
           replacing_an_elements_of_list(N,M,T[1])
           #L=reversing_list(M)
       return M
def coefficient_matrix_list(list_poly):
       max_degree=0
       for T in list_poly:
           if max_degree < max(list(T[0])):
               max_degree=max(list(T[0]))
       M=empty_list_generator(max_degree, len(list(T[0])))
       for T in list_poly:
           N=list(T[0])
           N.reverse()
           replacing_an_elements_of_list(N,M,T[1])
       return M

       #Ball Sysyem
def genvars(n):     #already exists in Ball lib
  x=[]
  r=[]
  i = 1
  while i < int(n)+1:
    x.insert(i, sp.Symbol('x'+str(i)) )
    if i>2 :
        r.insert(i, sp.Symbol('r'+str(i)) )
    i= i+1
  t=sp.Symbol('t')
  return [x,r,t]
def inBall(P,X):     #already exists in Ball lib
 i=2
 x=X[0]
 r=X[1]
 t=X[2]
 SPp=P
 DPm=P
 while i< len(x):
    SPp=SPp.subs(x[i],x[i]+r[i-2]*sp.sqrt(t))
    DPm=DPm.subs(x[i],x[i]-r[i-2]*sp.sqrt(t))
    i=i+1

 SPp=sp.expand(SPp)
 DPm=sp.expand(DPm)
 SP=sp.expand(0.5*(SPp+DPm))
 DP=sp.expand(0.5*(SPp-DPm)/sp.sqrt(t))
 return [SP,DP]
def removing_duplicated_terms(list_poly):
    list_poly_reduced=[]
    duplicated_index=[]
    k1=0
    k2=0
    k3=0
    while k1< len(list_poly):   #searching for equal monomials in  poly_SP and poly_DP
             if k1 not in duplicated_index:
                  list_poly_reduced.append(list_poly[k1])
                  k3+=1
                  k2=k1+1
                  while k2< len(list_poly) :
                     if list_poly[k3-1][0]==list_poly[k2][0]:
                       duplicated_index.append(k2)
                       list_poly_reduced[k3-1][1]+=list_poly[k2][1]
                       k2 =k2+1
                     else:
                        k2=k2+1

                  k1=k1+1

             else :
                k1+=1
    return(list_poly_reduced)
def Ball_for_terms(L,X):  #L=[list of integers,a coefficient] X is list of variables
    P=1
    if type(L[0])==list:
      for i in range(len(L[0])) :
        P=P*X[0][i]**(L[0][i])
      M= inBall(P,X)
      S=[sp.poly(T,*X[0],*X[1],X[2]).terms() for T in  M]
      S[0][0]=list(S[0][0])
      for i in  range(len(S)):
        for j in range(len(S[i])):
              S[i][j]=list(S[i][j])
              S[i][j][0]=list(S[i][j][0])
              S[i][j][1]=L[1]*S[i][j][1]
              
    
    return S
def Ball_for_interval_poly(list_of_terms,X):
    poly_SP=[]
    poly_DP=[]
    for T in list_of_terms:
        SDP=Ball_for_terms(T,X)
        poly_SP=poly_SP+SDP[0]
        poly_DP=poly_DP+SDP[1]

        
    

    return [removing_duplicated_terms(poly_SP),removing_duplicated_terms(poly_DP)]
def gauss_seidel_dimnp(A,b,x,z):
            x_prime=[xi for xi in x]
            for i in range(len(x)):
                x_i_prime=gauss_seidel_dim1p(A[i][i],b[i]-np.dot(A[i],x_prime)+A[i][i]*x_prime[i],z[i])
                x_prime[i]=x_i_prime
            return(x_prime)


def subdivide(B): #B is a list of ft.arb
     if len(B)==1:
         children=[[ft.arb(B[0].mid()+(B[0].rad()/2), B[0].rad()/2)],[ft.arb(B[0].mid()-(0.5*B[0].rad()),0.5* B[0].rad())]]
     else:
         B_prime=list(B)
         B_prime.remove(B_prime[len(B_prime)-1])
         B_pri_sub= subdivide(B_prime)
         B_one_sub=subdivide([B[len(B)-1]])
         children=cartesian_product(B_pri_sub,B_one_sub)
     return children
def ftconstructor(a,b): #if an interval is [a,b] it returns ft.arb((a+b)/2,(b-a)/2)
     return ft.arb((a+b)/2,(b-a)/2)

def k_subdivide(B,k=2): #B is a list of ft.arb
     children=[]
     if len(B)==1:
      for i in range(k):
        children.append([ftconstructor(B[0].lower()+2*i*(B[0].rad())/k,B[0].lower()+2*(i+1)*(B[0].rad())/k )])
     else:
         B_prime=B[:]
         B_prime.remove(B_prime[len(B_prime)-1])
         B_pri_sub= k_subdivide(B_prime,k)
         B_one_sub=k_subdivide([B[len(B)-1]],k)
         children=cartesian_product(B_pri_sub,B_one_sub)
     return children     
def i_minor(jac,i): #jac is a list of lists of  list_poly. i_minor(jac,i) is the i-th minor of jac
         jac_array=copy(np.array(jac))
         minor_i_array=np.delete(jac_array, i, axis=1)
         minor_i=minor_i_array.tolist()
         return minor_i
def matrixval(jac,X): #jac as i_minor and X is a list of ft.arb
         #coeffs_matrix_matrix=[copy(jaci) for jaci in jac]
         evaluation_jac_X=[copy(jaci) for jaci in jac]
         for i in range(len(jac)):
             for j in range(len(jac[0])):
                 evaluation_jac_X[i][j]=evaluation_poly_list(jac[i][j],X)
                 if evaluation_jac_X[i][j]==0:       #fixing a bug in flint since Python cannot create arb from type <class 'sympy.core.numbers.Zero'>
                     evaluation_jac_X[i][j]=ft.arb(0)
         return evaluation_jac_X
def width(B):
    if type(B[0])==type(ft.arb(1)):
     maximal=max([float(Bi.rad()) for Bi in B])
     return maximal
    if type(B[0])==type([]):
     B_ft=[ftconstructor(Bi[0],Bi[1]) for Bi in B]
     return width(B_ft) 
def checking_full_rank(M): #checks whether a list of lists of arb is full rank
            minors=[i_minor(M,t) for t in range(len(M[0]))]
            full_rank=0
            invertibility_of_eval_minors=[invertibility_of_a_matrix(Mi) for Mi in minors]
            #print((minors[len(M[0])-1][0]).lower(),(minors[len(M[0])-1][0]).upper())

            if  1 in invertibility_of_eval_minors:
                full_rank=1
            return full_rank

def derivative_poly_list(poly_list,i): #computes the derivative of poly_list wrt i-th variable
    derivative=[]
    for term in poly_list:
       if term[0][i-1]!=0:
            derivative_mon=copy(term[0])
            derivative_mon[i-1]=term[0][i-1]-1
            derivative_coeff=term[1]*(term[0][i-1])
            derivative_term=[derivative_mon,derivative_coeff]
            derivative.append(derivative_term)
    if len(derivative)==0:
           derivative=[[[0]*len(poly_list[0][0]),0]]
    return derivative
def jacobian_of_function_list(P):
     jacobian=[]
     n=len(P[0][0][0])
     for i in range(len(P)):
         jacobian.append([])
         for j in range(n):
              jacobian[i].insert(j,derivative_poly_list(P[i],j+1))
     return jacobian
def Ball_interval(list_polys):
    X=genvars(len(list_polys[0][0][0]))
    Ball_system1=[]
    Ball_system2=[]
    for list_poly in list_polys:
        SD_list_poly=Ball_for_interval_poly(list_poly,X)
        Ball_system1.append(SD_list_poly[0])  #SP's
        Ball_system2.append(SD_list_poly[1])   #DP's
    Ball_system= [*Ball_system1,*Ball_system2]   #removing unnessary terms
    for Pi in Ball_system:
        for Pij in Pi:
            if Pij[1]==0:
                Pi.remove(Pij)
    last_eq=[]
    for r in X[1]:
          Ter=sp.poly(r**2,*X[0],*X[1],X[2]).terms()
          Ter[0]=list(Ter[0])
          Ter[0][0]=list(Ter[0][0])
          last_eq.append(Ter[0])

    last_eq.append([[0]*(2*len(X[0])-1),-1])


    return [*Ball_system,last_eq]


def two_intervals_compare(a,b):
    Answer=0
    if a.mid()==b.mid() and a.rad()==b.rad():
        Answer =1
    return Answer
def interval_difference_enclosure(a,b):
     if a in b :
         Difference= 'empty'
     else:
           try :
              Intersection=a.intersection(b)
              if Intersection.upper() >= a.upper():
                  Difference=ft.arb(0.5*a.lower()+0.5*Intersection.lower(),0.5*Intersection.lower()-0.5*a.lower())
              elif   Intersection.lower() <= a.lower():
                  Difference=ft.arb(0.5*a.upper()+0.5*Intersection.upper(),0.5*a.upper()-0.5*Intersection.upper())
              else:
                  Difference=a
           except:
              Difference=a
     return  Difference
def boxes_intersection(B1,B2):
  if type(B1[0])==type(ft.arb(1)):
    inters=[]
    for i in range(len(B1)):
       try:
         inters.append(B1[i].intersection(B2[i]))
       except:
         inters=[]
         break
    return  inters     
  elif type(B1[0])==type([]):
    B1_ft=[ftconstructor(Bi[0],Bi[1]) for Bi in B1]
    B2_ft=[ftconstructor(Bi[0],Bi[1]) for Bi in B2]
    intersec= boxes_intersection(B1_ft,B2_ft)
    return [[float(Bi.lower()),float(Bi.upper())] for Bi in intersec]
 
def components_intersection(c1,c2,com=2): 
  for b1 in c1:
    for b2 in c2:
      if boxes_intersection(b1[:com],b2[:com]) !=[]:
        return True
  return False      


def intersection_of_bounded_with_unbounded_interval(x,boundary,infty): #returns the intersection of x with [boundary, infty] where infty= 1 or -1  is corresponding to + infty or -infty respectevly
                 Intersection=x
                 if infty==1:
                     if x.upper()< boundary:
                            Intersection='empty'
                     elif  x.upper() >= boundary and x.lower() < boundary:
                         Intersection=ft.arb(0.5*x.upper()+ 0.5*boundary,0.5*x.upper()-0.5*boundary)
                     elif x.lower() >= boundary:
                            Intersection=x
                 else:
                         if x.lower() >  boundary:
                                Intersection='empty'
                         elif  x.lower() <= boundary and x.upper() > boundary:
                            Intersection=ft.arb((0.5*x.lower())+ (0.5*boundary),(0.5*boundary)-(0.5*x.lower()))
                         elif x.upper() <= boundary:
                                 Intersection=x
                 return Intersection
def special_cases_gauss_seidel(a,b,x,precision=0.00000001): # it computes Gamma operator for the case 0 in a (see Neumaier Proposition 4.3.1, eq (5))
     Gamma =x
     if a*ft.arb(1) in ft.arb(0):
         Gamma='empty'
     elif b.lower()> 0 :
        if  (-a.lower() > precision  and a.upper()> precision):  # that is if 0 is not in the boundary of a
            interval = ft.arb((a.upper()* b.lower()+ (a.lower()* b.lower()) )/(2*a.upper()*a.lower()), (a.upper()* b.lower()- (a.lower()* b.lower()) )/(2*a.upper()*a.lower()))
            Gamma=interval_difference_enclosure(x,interval)
        elif a.upper()> precision :     #that is a.lower()=0 then the interval is [b.lower()/a.upper(), infty ]

              Gamma=intersection_of_bounded_with_unbounded_interval(x,b.lower()/a.upper(),1)
        elif -a.lower()> precision :     #that is a.upper()=0 then the interval is [-infty ,b.lower()/a.lower()  ]
              Gamma=intersection_of_bounded_with_unbounded_interval(x,b.lower()/a.lower(),-1)
     elif b.upper()< 0 :
         if  -a.lower() > precision  and a.upper()> precision :
            interval =ft.arb(0.5*(b.upper()/a.upper())+0.5*(b.upper()/a.lower()), 0.5*(b.upper()/a.lower())-0.5*(b.upper()/a.upper()) )
            Gamma=interval_difference_enclosure(x,interval)
         elif  a.upper()> precision :
               Gamma=intersection_of_bounded_with_unbounded_interval(x,b.upper()/a.upper(),-1)
         elif -a.lower()> precision :
                Gamma=intersection_of_bounded_with_unbounded_interval(x,b.upper()/a.lower(),1)
     elif  -b.lower() > precision  and b.upper()> precision: # if 0 is in interor of b
         Gamma=x
     elif b.upper()> precision:
                if  -a.lower() > precision  and a.upper()> precision :
                    Gamma=x
                elif  a.upper()> precision :
                      Gamma=intersection_of_bounded_with_unbounded_interval(x,0,1)
                elif -a.lower()> precision :
                       Gamma=intersection_of_bounded_with_unbounded_interval(x,0,-1)
     elif -b.lower() > precision:
        if  -a.lower() > precision  and a.upper()> precision :
            Gamma=x
        elif  a.upper()> precision :
              Gamma=intersection_of_bounded_with_unbounded_interval(x,0,-1)
        elif -a.lower()> precision :
               Gamma=intersection_of_bounded_with_unbounded_interval(x,0,1)
     return Gamma
def cartesian_product(B1,B2):
    the_product=[]
    for i in range(len(B1)):
        for j in range(len(B2)):
            """print(B1[i])
            print(B2[j])
            input()"""
            the_product.append(B1[i]+B2[j])
    return the_product
def polyvalnd(M,X):
    evaluation=ft.arb(0)

    if len(X)==1:
        N=[]
        for Mi in M:     #to avoid the problem: sp.Integer * arb = not correct answer
            if type(Mi)!= ft.arb:
                N.append(float(Mi))
            else:
                N.append(Mi)
        evaluation= np.polynomial.polynomial.polyval(X[0],N)
    else:
        Y=copy(X[0:len(X)-1])
        N=copy([polyvalnd(Mi,Y) for Mi in M])

        evaluation=np.polynomial.polynomial.polyval(X[len(X)-1],N)
    return evaluation

def gauss_seidel_dim1(a,b,x):
    Answer=x
    if 0 not in (a * ft.arb(1)):
          try:
            Answer=x.intersection((b * ft.arb(1))/a)
          except:
             Answer='empty'
    else:

        Answer=special_cases_gauss_seidel(a,b,x)
    return Answer
def gauss_seidel_dimnorigin(A,b,x,z): #A,b,x,z=lists (matrix) of ft.arb
            #ftprint(x)
            #ftprint(b)
            #L=[ftprint(Ai) for Ai in A]
            #input()
            x_prime=[xi for xi in x]
            if invertibility_of_a_matrix(A)==1:
              A_arb_mat=ft.arb_mat(A)
              inv_m_A=(A_arb_mat.mid()).inv()
              A_precon=inv_m_A*ft.arb_mat(A)
              b_arb_mat=ft.arb_mat([[bi] for bi in b])
              b_precon=inv_m_A*b_arb_mat
              for i in range(len(x)):
                      sum=0
                      for j in range(len(x)):
                          if i !=j:
                              sum+=A_precon[i,j]*x_prime[j]
                      x_prime[i]=gauss_seidel_dim1(A_precon[i,i],b_precon[i,0]-sum,z[i])
                      if x_prime[i]=='empty':
                                x_prime='empty'
                                break
              #ftprint(x_prime)
              #input()                  
            return x_prime

def gauss_seidel_dimn(A,b,x,z): #A,b,x,z=lists (matrix) of ft.arb
            x_prime=[xi for xi in x]
            if invertibility_of_a_matrix(A)==1:
              A_arb_mat=ft.arb_mat(A)
              inv_m_A=(A_arb_mat.mid()).inv()
              A_precon=inv_m_A*ft.arb_mat(A)
              b_arb_mat=ft.arb_mat([[bi] for bi in b])
              b_precon=inv_m_A*b_arb_mat
              """X=A_precon.solve(b_precon)
              sol=[]
              for i in range(len(x)):
                sol.append(X[i,0])
              #print(type(sol[0]))
              #print(ft.arb.solve)
              #input()  


              """
              for i in range(len(x)):
                      sum=0
                      for j in range(len(x)):
                          if i !=j:
                              sum += intervals_multi(A_precon[i,j],x_prime[j])
       
                      x_prime[i]=gauss_seidel_dim1(A_precon[i,i],b_precon[i,0]-sum,z[i])
                      if x_prime[i]=='empty':
                                x_prime='empty'
                                break
              #ftprint(x_prime)
              #input()                  
            return x_prime

def poly_normal_to_list(f,X):  # f is a sympy expression, X is  the  set of variables
      S=sp.poly(f,X).terms()
      M=[]
      #S[0][0]=list(S[0][0])
      for i in  range(len(S)):
                M.append([0,0])
                M[i][0]=list(S[i][0])
                M[i][1]=S[i][1]
      return M





def checking_regularity_of_sysyem1(P,B,jac,wth=0.2):
    M=[coefficient_matrix_list(Pi) for Pi in P]
    list_of_boxes=[B]
    smoothness=1
    regular_boxes=[]
    regular_volume=0
    empty_volume=0
    unknownvolume=n_volume(B)
    it=0
    while len(list_of_boxes)!=0 and width(list_of_boxes[0])> wth:
              #print(width(list_of_boxes[0]))
              membership=1
              eval_P=[polyvalnd(Mi,list_of_boxes[0]) for Mi in M]

              for eval_Pi_at_B in eval_P:
                  if 0 not in ft.arb(1)*(eval_Pi_at_B):
                      membership=0
                      break
              if membership==0:
                 empty_volume +=n_volume(list_of_boxes[0])
                 unknownvolume=unknownvolume-n_volume(list_of_boxes[0])
                 regular_boxes.append(list_of_boxes[0])
                 list_of_boxes.remove(list_of_boxes[0])


              else :
                     eval_jac=matrixval(jac,list_of_boxes[0])
                     full_rankness= invertibility_of_a_matrix(eval_jac)
                     if   full_rankness ==1:
                         regular_volume +=n_volume(list_of_boxes[0])
                         unknownvolume=unknownvolume-n_volume(list_of_boxes[0])
                         regular_boxes.append(list_of_boxes[0])
                         list_of_boxes.remove(list_of_boxes[0])

                     else:
                        new_children=subdivide(list_of_boxes[0])
                        list_of_boxes.remove(list_of_boxes[0])
                        list_of_boxes = list_of_boxes +new_children
              print("regular volume=",regular_volume)
              print("empty volume=",empty_volume)
              print("unknown volume=",unknownvolume)
                        #print('Unknown box detected')
              #if len(list_of_boxes)!=0:
                #smoothness=-1

    return regular_boxes

def checking_regularity_of_sysyem(P,B,jac,wth=0.2):
    list_of_boxes=[B]
    smoothness=1
    regular_boxes=[]
    empty_boxes=[]
    regular_volume=0
    empty_volume=0
    unknownvolume=n_volume(B)
    while len(list_of_boxes)!=0 and width(list_of_boxes[0])> wth:
              current_box=list_of_boxes[0]
              membership=1
              eval_P=[evaluation_poly_list(Pi,current_box) for Pi in P]
              for eval_Pi_at_B in eval_P:
                  if type(eval_Pi_at_B) != type(ft.arb(1)):
                       eval_Pi_at_B=float(eval_Pi_at_B)
                  if 0 not in ft.arb(1)*(eval_Pi_at_B):
                      membership=0
                      break
              if membership==0:
                 empty_boxes.append(list_of_boxes[0])
                 empty_volume +=n_volume(list_of_boxes[0])
                 unknownvolume=unknownvolume-n_volume(list_of_boxes[0])
                 list_of_boxes.remove(list_of_boxes[0])
              else :
                     eval_jac=[[ft.arb(evaluation_poly_list(jac_ij,list_of_boxes[0])) for jac_ij in jac_i ] for jac_i in jac ]
                     if  invertibility_of_a_matrix(eval_jac)==1: # 0 not in (ft.arb_mat(eval_jac)).det():
                         #print(ft.arb_mat(eval_jac))
                         regular_boxes.append(list_of_boxes[0])
                         regular_volume=regular_volume+n_volume(list_of_boxes[0])
                         unknownvolume=unknownvolume-n_volume(list_of_boxes[0])
                         list_of_boxes.remove(list_of_boxes[0])
                     else:
                        new_children=subdivide(list_of_boxes[0])
                        list_of_boxes.remove(list_of_boxes[0])
                        list_of_boxes = list_of_boxes +new_children
    
    if len(list_of_boxes)!=0:
                smoothness=-1

    return smoothness
def n_volume(B):
    volume=1
    for Bi in B :
        volume=volume*(2*Bi.rad())
    return volume
def power_interval(a,n):
  the_power=ft.arb(1)
  if type(a) != type(ft.arb(1)):
    return float(a)**n
  elif n!=0:
    if a.lower()>=0:
        the_power=ft.arb(0.5*(a.lower())**n+0.5*(a.upper())**n,0.5*(a.upper())**n-0.5*(a.lower())**n)
    elif a.upper()<= 0:
        if n %2 ==0:
            
            the_power=ft.arb(0.5*(a.lower())**n+0.5*(a.upper())**n,0.5*(a.lower())**n-0.5*(a.upper())**n)
            
        else:
            the_power=ft.arb(0.5*(a.lower())**n+0.5*(a.upper())**n,0.5*(a.upper())**n-0.5*(a.lower())**n)
    else:
        a1= ft.arb(0.5*a.lower(),-0.5*a.lower())
        a2= ft.arb(0.5*a.upper(),0.5*a.upper())
        the_power2=ft.arb(0.5*(float(a2.upper()))**n,0.5*(float(a2.upper()))**n)
        if n %2 ==0:
            the_power1=ft.arb(0.5*(a1.lower())**n,0.5*(a1.lower())**n)
        else:
            the_power1=ft.arb(0.5*(a1.lower())**n,-0.5*(a1.lower())**n)
        the_power=the_power1.union(the_power2)
  return the_power
def intervals_multi(B1,B2): #does interval multiplication for B1,B2 if they are of different sign to avoid a bug in flint
    if type(B1)!=type(ft.arb(1)) :
        B1=ft.arb(float(B1))
    if type(B2)!=type(ft.arb(1)) :
            B2=ft.arb(float(B2))
    end_points=[B1.lower()*B2.lower(),B1.lower()*B2.upper(),B1.upper()*B2.lower(),B1.upper()*B2.upper() ]
    if 0 not in B1 and 0 not in B2 :
        lower_point=min(end_points)
        upper_point=max(end_points)
        the_product=ft.arb(0.5*(upper_point+lower_point),0.5*(upper_point-lower_point))
    else:
        the_product=B1*B2 
    return the_product

def evaluation_poly_list(poly_list, B):
    evaluation=ft.arb(0)
    for term in poly_list:
        value_of_term=evaluation_term(term, B)
        if type(value_of_term)==type(ft.arb(1)):
           evaluation=evaluation+value_of_term
        else:
                evaluation=evaluation+ft.arb(float(value_of_term))
       
    return evaluation

def evaluation_term(term, B):    #evaluate a term at a box B... term
       if [str(T) for T in term[0]]==["0"]*len(term[0]) :
           
           if type(term[1]) ==type(ft.arb(1)): #to avoid an error by sympy  zeros
              the_evaluation=term[1]
           else:
                the_evaluation =ft.arb(float(term[1]))
      
       else:
           the_evaluation=ft.arb(1)
           for i in range(len(term[0])):
                if term[0][i]!= 0:
                   the_evaluation=intervals_multi(the_evaluation,(power_interval(B[i],term[0][i])))
                   #print(B[i], term[0],i)
                   #print(the_evaluation.lower())
                   #ftprint([intervals_multi(the_evaluation,(power_interval(B[i],term[0][i])))])
                   #input()
                   #the_evaluation=(the_evaluation)*(power_interval(B[i],term[0][i]))
           if type(term[1]) ==type(ft.arb(1)): #to avoid an error by sympy  zeros
                      the_evaluation=the_evaluation*term[1]
           else:
               the_evaluation =intervals_multi(the_evaluation,ft.arb(float(term[1])) )

       return the_evaluation
def matrix_multi(M1,M2):
    the_product=[[]]*(len(M1))

    for i in range(len(M1)):
        for j in range(len(M2[0])):
            sum=0
            for k in range(len(M1[0])):
                sum=sum+intervals_multi(M1[i][k],M2[k][j])
            the_product[i]=the_product[i]+[sum]
    return the_product

def invertibility_of_a_matrix(M):    # M is interval matrix .. invertibility_of_a_matrix(M) returns 0 if we are sure that M is singular, 1 if we are sure that M invertable and -1 if we do not  know.. The smaller width M has the more sure we are
 Answer =3

 L=[[Mij.mid() if type(Mij)==type(ft.arb(1)) else ft.arb(float(Mij)) for Mij in Mi] for Mi in M ]
 T=[]
 try:
          T=(Matrix(L)).inv()
 except:
        Answer=0
 if Answer==3:
        T_list=[[float(T[j,i]) for i in range(len(M))] for j  in range(len(M))]
        Pre_cond=matrix_multi(T_list,M)
        E=ft.arb_mat(T_list)*ft.arb_mat(M)
        zero_eigenvalue_detected=0
        Gersh=gershgorin_circles(Pre_cond)
        
        for R in Gersh:
             if 0 in R:
                 zero_eigenvalue_detected=1
                 break        
        if zero_eigenvalue_detected==0:
            Answer =1
        if zero_eigenvalue_detected==1:
            Answer=-1
        #print('the inver: ', Answer)
        #print([[[Mij.lower(),Mij.upper()] for Mij in Mi] for Mi in Pre_cond  ])
        #print('Gersh=',[[R.lower(),R.upper()] for R in Gersh ])
          

        #if 0 not in E.det():
        #    Answer=1
        #else:
        #    Answer = -1

 return Answer    #here is a problem

def checking_smoothness1(P,B,jac,wth=0.1):
    M=[coefficient_matrix_list(Pi) for Pi in P]
    list_of_boxes=[B]
    smoothness=1
    regular_boxes=[]
    while len(list_of_boxes)!=0 and width(list_of_boxes[0])> wth:
              membership=1
              eval_P=[polyvalnd(Mi,list_of_boxes[0]) for Mi in M]
              for eval_Pi_at_B in eval_P:
                  if 0 not in ft.arb(1)*(eval_Pi_at_B):
                      membership=0
              eval_jac=matrixval(jac,list_of_boxes[0])
              full_rankness= checking_full_rank(eval_jac)

              if   membership==0  or full_rankness ==1:
                  regular_boxes.append(list_of_boxes[0])
                  list_of_boxes.remove(list_of_boxes[0])
              else:
                 new_children=subdivide(list_of_boxes[0])
                 list_of_boxes.remove(list_of_boxes[0])
                 list_of_boxes = list_of_boxes +new_children
    if len(list_of_boxes)!=0:
        smoothness=-1

    return smoothness



def curve_tracer(P,B,jac,wth=0.001,wth2=1):  #returns all list of boxes that contains smooth parts of a curve and it stops if the curve is smooth
    list_of_boxes=[B]
    smoothness=1
    regular_boxes=[]
    smoothness=1
    while len(list_of_boxes)!=0 and width(list_of_boxes[0])> wth:
              membership=1
              eval_P=[evaluation_poly_list(Pi,list_of_boxes[0]) for Pi in P]           #checking whether the box contains a point of the curve
              for eval_Pi_at_B in eval_P:
                  if 0 not in ft.arb(1)*(eval_Pi_at_B):
                      membership=0
              if membership==0:
                list_of_boxes.remove(list_of_boxes[0])
                  #print("empty")
              else:
                   eval_jac=[[evaluation_poly_list(jac_ij,list_of_boxes[0]) for jac_ij in Jac_i ]  for Jac_i in jac  ]
                   full_rankness= checking_full_rank(eval_jac)
                   if   full_rankness ==1 and width(list_of_boxes[0])<wth2 :
                        #print([ [round(float(jac_ij.lower()),5),round(float(jac_ij.upper()),5)] for jac_ij in list_of_boxes[0]  ])
                        #input()
                        #print("regular")
                        #input()
                        regular_boxes.append(list_of_boxes[0])
                        list_of_boxes.remove(list_of_boxes[0])
                   else:
                        #print('width=',width(list_of_boxes[0]))
                        #print('the box',[ [round(float(jac_ij.lower()),5),round(float(jac_ij.upper()),5)] for jac_ij in list_of_boxes[0]  ])
                        #pprint(Matrix([ [[round(float(jac_ij.lower()),5),round(float(jac_ij.upper()),5)] for jac_ij in jac_i] for jac_i in eval_jac  ]))
                        #input()
                        #print("we do not know")
                        new_children=subdivide(list_of_boxes[0])
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
        """print(B1[i])
        print(B2[i])
        input()"""
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
def connected_components(list_of_boxes):  #Returns the connected components of list_of_boxes
  connected_components={}
  list_removable_ind=[]
  flaged_boxes=[]
  for box_index in range(len(list_of_boxes)):
    if box_index not in flaged_boxes:
      connected_components[box_index]=[box_index]
      for diff_box in range(box_index+1,len(list_of_boxes)):
        intersection_certificate=1
        intersection_box=[]
        for j in range(len(list_of_boxes[0])):
          try:
            intersection_j=list_of_boxes[box_index][j].intersection(list_of_boxes[diff_box][j])
            intersection_box.append(intersection_j)
          except:
             intersection_certificate=0
             break
        if intersection_certificate ==1:
            connected_components[box_index].append(diff_box)
            flaged_boxes.append(diff_box)              
  return connected_components

def distance(B1,B2):
  if len(B1)==0 or len(B2)==0:
    return []
  elif type(B1[0])==type(ft.arb()):
    norm=ft.arb(0)
    for i in range(len(B1)):
      norm += power_interval(B1[i]-B2[i],2)
    if norm.lower()>=0:
        return ftconstructor(math.sqrt(float(norm.lower())),math.sqrt(float(norm.upper())))
    else :
        return ftconstructor(0,math.sqrt(float(norm.upper()))) 
  elif  type(B1[0])==type([]): 
      B1_ft=[ftconstructor(Bi[0],Bi[1]) for Bi in B1]
      B2_ft=[ftconstructor(Bi[0],Bi[1]) for Bi in B2]
      dis= distance(B1_ft,B2_ft)
      return [float(dis.lower()),float(dis.upper())] 
        
    

                                             
    

def B_Ball_calculator(B):   # returns B_Ball
    norm_diff=ft.arb(0)                                            #computing \Xi (Remark 5.2.10) 
    for i in range(2,len(B)):
    	norm_diff+= power_interval(B[i]-B[i],2)
    
    max_t= float(norm_diff.upper())
    B_Ball= B[:2]  
    for i in range(2,len(B)):
    	B_Ball.append(B[i]+ft.arb(0,max_t))
    B_Ball+=[ft.arb(0,1)]*(len(B)-2)
    B_Ball+= [ft.arb(0.5*max_t, 0.5*max_t)]
    return  B_Ball

def box_union(B1,B2):
    the_union=[]
    
    if B1==[]:
      return B2
    if B2==[]:
     return B1
    #print(B1);input() 
    if type(B1[0])==type(ft.arb(1)):
      for i in range(len(B1)) :
            the_union.append(B1[i].union(B2[i]))
      return the_union
    if  type(B1[0])==type([1]): 
      B1_ft=[ftconstructor(Bi[0],Bi[1]) for Bi in B1]
      B2_ft=[ftconstructor(Bi[0],Bi[1]) for Bi in B2]
      uni= box_union(B1_ft,B2_ft)
      return [[float(Bi.lower()),float(Bi.upper())] for Bi in uni]

def checking_smoothness(P,B,jac,wth=0.1):
     M=[coefficient_matrix_list(Pi) for Pi in P]
     list_of_boxes=[B]
     smoothness=1
     regular_boxes=[]
     while len(list_of_boxes)!=0 and width(list_of_boxes[0])> wth:
               membership=1
               eval_P=[polyvalnd(Mi,list_of_boxes[0]) for Mi in M]
               for eval_Pi_at_B in eval_P:
                   if 0 not in ft.arb(1)*(eval_Pi_at_B):
                       membership=0
               eval_jac=matrixval(jac,list_of_boxes[0])
               full_rankness= checking_full_rank(eval_jac)
               if   membership==0  or full_rankness ==1:
                   regular_boxes.append(list_of_boxes[0])
                   list_of_boxes.remove(list_of_boxes[0])
               else:
                  new_children=subdivide(list_of_boxes[0])
                  list_of_boxes.remove(list_of_boxes[0])
                  list_of_boxes = list_of_boxes +new_children
     if len(list_of_boxes)!=0:
         smoothness=-1
     print("Empty or regular boxes:",regular_boxes )
     print("Unknown-behavior boxes",list_of_boxes)
     return smoothness


def ftprint(B,k=3):
 answer=[]
 if type(B[0])==type(ft.arb(1)):   
  for Bi in B:
    if type(Bi)==type(ft.arb(0,1)):
      answer.append([round(float(Bi.lower()),k),round(float(Bi.upper()),k) ] )
    else:
      answer.append([Bi,Bi])
  print(answer)
 elif  type(B[0])==type([]):
  pprint(Matrix([[ [round(float(Bij.lower()),k),round(float(Bij.upper()),k) ] if \
    type(Bij)==type(ft.arb(1)) else [Bij,Bij  ] for Bij in Bi ] for Bi in B ] ))
def hansen_hengupta(x_teld,A,b,x,z):   
    """
    It returns the output of Hansen Hengupta operator
    A=[[ft.arb(2,1),ft.arb(0,1),ft.arb(0,1)],[ft.arb(0,1),ft.arb(3,1),ft.arb(0,1)],[ft.arb(0,1),ft.arb(0,1),ft.arb(4,1)]]
    b=[ft.arb(4,1),ft.arb(6,1),ft.arb(5,1)]
    x=[ft.arb(0,1),ft.arb(0,1.5),ft.arb(0,1),ft.arb(1,1),ft.arb(0,1)]
    z=[ft.arb(1,0.1),ft.arb(1,0.5),ft.arb(0,1),ft.arb(1,1),ft.arb(0,1)]
    b=[ft.arb(0,0.1),ft.arb(0,0.1),ft.arb(0,0.1),ft.arb(0,0.1),ft.arb(0,0.1)]
    x_teld=[0.01,0.01,0.01,1.01,0.01]
    >>> hansen_hengupta(x_teld,A,b,x,x)
    [ft.arb(0.0100000000000000 +/- 2.10e-19,0.100000000558794 +/- 4.56e-16), ft.arb(0.0100000000000000 +/- 2.10e-19,0.715100008994341 +/- 1.04e-16),ft.arb(0.0100000000000000 +/- 2.10e-19,0.0500000002793968 +/- 2.77e-17),ft.arb(1.01000000000000 +/- 9.0e-18,0.0500000002793968 +/- 2.77e-17),ft.arb(0.0100000000000000 +/- 2.10e-19,0.126582279801369 +/- 2.87e-16)]
    >>> hansen_hengupta(x_teld,A,b,z,z)
    'empty'
    """
    x_prime=[]
    z_prime=[]
    b_prime=[-bi for bi in b]
    for i in range(len(x)):
      x_prime.append(x[i]-x_teld[i])
      z_prime.append(z[i]-x_teld[i])
    
    Gamma=[]
    Gamma1=[]
    Gamma=gauss_seidel_dimn(A,b_prime,x_prime,z_prime)
    """Sol=ft.acb_mat(A).solve(ft.acb_mat(b_prime))

  
    for i in range(len(b_prime)):
      Gamma1.append(Sol[1,0].real)

    for i in range(len(b_prime)):
      try:
        Gamma.append(Gamma1[i].intersection(x_prime[i]))  
      except:
        Gamma="empty"
        break"""
        

    
    S=[]
    if Gamma!= 'empty':
      #ftprint(x)

      """print("Gamma:")
      #pprint(x_teld)
      #print('x')
      ftprint(x,5)
      #print('x_prime') 
      ftprint(x_prime,5)
      pprint(Matrix(A))
      #ftprint(b,5)
      ftprint(b_prime,5)
    
      ftprint(Gamma,5)
      input()"""
      for i in range(len(Gamma)):
        S.append(Gamma[i]+x_teld[i])
      Answer=S[:]  
      #ftprint(S)
      #input()  
        
    else:
        Answer ='empty'

    
    return Answer
def solver(P,jac,B,k=2): #k is the number of parts in which every interval is divided  
    it=0
    Solutions=[]
    L=[B]
    while len(L) !=0:
        it=it+1
        current_box=L[0]  #evaluating P at the currecnt_box:
        value_of_P_current_box= [evaluation_poly_list(Pi,current_box) for Pi in P]
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
                   jac_eval_current_box=matrixval(jac,current_box)
                   mid_box=[ft.arb(float(interval.mid()))  for interval in current_box]
                   b=[evaluation_poly_list(Pi,mid_box) for Pi in P]
                   """print('mid',mid_box)
                   print('the invertability',invertibility_of_a_matrix(jac_eval_current_box))
                   try:
                    pprint(jac_eval_current_box)
                   except:
                    pass
                   print('b',b)
                   ftprint(current_box)
                   input()"""
                   if invertibility_of_a_matrix(jac_eval_current_box)==1:
                    Image_of_current_box=hansen_hengupta(mid_box,jac_eval_current_box,b,current_box,current_box)
                    #ftprint(Image_of_current_box)
                    #input()
                    if Image_of_current_box !='empty':
                                    """ftprint(current_box)
                                    ftprint(Image_of_current_box)
                                    input()"""
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
                                               new_children=k_subdivide(current_box,k)

                                               L.remove(current_box)
                                               L =   L +new_children
                                           except:
                                               L.remove(current_box)
                    else:
                          L.remove(L[0])

            
                   else:
                    new_children=k_subdivide(current_box,k)
                    L.remove(current_box)
                    L =   L +new_children

                                               

    #print(Solutions)
    return Solutions     

def solver2(P,jac,B):
   t1=time.time()
   S=solver(P,jac,B)
   S2=[]
   k=2
   connected_components1=connected_components(S)
   while len(connected_components1)!= len(S) :
    for component in connected_components1:
      if len(connected_components1[component])==1:
        S2.append(S[connected_components1[component][0]])
      else:
        union_box=S[connected_components1[component][0]]
        for i in range(1,len(connected_components1[component])):
          union_box=box_union(union_box,S[connected_components1[component][i]])
        union_box_inflated= [x + y for x, y in zip(union_box, [ft.arb(0,0.001)]*len(union_box))]      

        S3=solver(P,jac,union_box_inflated,k)
        S2=S2+S3
    k+=1    
    S=S2[:]
    S2=[]
    connected_components1= connected_components(S) 

   t2=time.time()
   print(t2-t1)
   return S    




    

def compose(f, g):
    return lambda x: f(g(x))
def sqrt_t(I):
  if 0 in I:
    return ftconstructor(0,math.sqrt(float(I.upper())))
  else:
       return ftconstructor(math.sqrt(float(I.lower())),math.sqrt(float(I.upper())))
def F_Ballplus(U):    # len(U) is odd
 if type(U[0])==type(ft.arb(1)): 
  n=int((len(U)+1)/2)
  Y=U[:n]
  r_times_sqrtt=[yi *sqrt_t(U[2*n-2]) for yi in U[n:2*n-2] ]
  F1=U[:2]
  for i in range(2,n):
      F1.append(Y[i]+r_times_sqrtt[i-2])
  return F1
 if type(U[0])==type([]):
   U_ft=[ftconstructor(Ui[0],Ui[1]) for Ui in U]
   F1=F_Ballplus(U_ft)
   return [[ float(F1i.lower()),float(F1i.upper()) ] for F1i in  F1]
def F_Ballminus(U):    # len(U) is odd
  if type(U[0])==type(ft.arb(1)): 
    n=int((len(U)+1)/2)
    Y=U[:n]
    r_times_sqrtt=[yi *sqrt_t(U[2*n-2]) for yi in U[n:2*n-2] ]
    F1=U[:2]
    for i in range(2,n):
      F1.append(Y[i]-r_times_sqrtt[i-2])
    return F1
  if type(U[0])==type([]):
   U_ft=[ftconstructor(Ui[0],Ui[1]) for Ui in U]
   F1=F_Ballplus(U_ft)
   return [[ float(F1i.lower()),float(F1i.upper()) ] for F1i in  F1]        

def Ball_func(func,Jac_func): # func is a function that sends a list of intervals to an internal ....   Jac_func(i) is the pratial dervative of func wrt the i-variable 
   S_func=lambda U: 1/2*(compose(func,F_Ballplus)(U)+compose(func,F_Ballplus)(U))
   D_func= lambda U: 1/(2*sqrt_t(U[len(U)-1]))*(compose(func,F_Ballplus)(U)+compose(func,F_Ballplus)(U)) if 0 not \
    in U[len(U)-1] else sum([nabla_funci(U)*ri for Yi,ri in zip(Jac_func,[0,0]+U[int((len(U)+1)/2):len(U)-1])  ])
   return [S_func,D_func] 


def poly_list_tofunc(P):
  return lambda B: evaluation_poly_list(P,B)


def Jet_poly_list(P):  # P is map from R^n to R^{n-1} and returns JetP
  jac_P=jacobian_of_function_list(P)
  Id_n=np.eye(len(P)+1, dtype=int) #the identity matrix of size n
  JetP=[]
  for j in range(len(P)):
    JetP.append({ tuple(Id_n[i]):jac_P[j][i] for i in range(len(P)+1)  })
  for k in range(len(P)):
    J2=jacobian_of_function_list(jac_P[k])
    Jet2={}
    for i in range(len(J2)):
      for j in range(i,len(J2[0])):
        Jet2.update({tuple(map(operator.add, Id_n[i], Id_n[j])):J2[i][j] })
    JetP[k].update(Jet2)
    JetP[k].update({(0,)*(len(P)+1):P[k] })
  return JetP 

def decimal_str(x: float, decimals: int = 10) -> str:
    return format(x, f".{decimals}f").lstrip().rstrip('0')

def ft_normal(B):
   return [ [round(float(Bi.lower()),6),round(float(Bi.upper()),6)] for Bi in B  ]
     



"""
  
#example to transfer a sympy expression to a poly_list
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
t1= sp.Symbol('t1')
t2= sp.Symbol('t2')
X=[x1,x2,t1,t2]


P1=x1
P2=x2
P3=t1
P4=t2
P=[P1,P2,P3,P4]
P=[poly_normal_to_list(Pi,X) for Pi in P] 
B=[ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1)]
jac=jacobian_of_function_list(P)
#print(solver(P,jac,B))
#input()








#Example 1 R_RRRR_ robot
a1=8.1
a2=5.1
a3=8.2
a4=5.4  
b=9

P1=sp.expand((a1**2+2*a1*x1-a2**2+x1**2+x2**2)*t1**2-4*a1*t1*x2+a1**2-2*a1*x1-a2**2+x1**2+x2**2)
P2=sp.expand((a4**2+(-2*b+2*x1)*a4-a3**3+b**2-2*b*x1+x1**2+x2**2)*t2**2-4*a4*t2*x2+a4**2+a4*(2*b-2*x1)-a3**3+b**2-2*b*x1+x1**2+x2**2)
func=Matrix([P1,P2])
var=Matrix(X)
M=func.jacobian(var)

P3= (M[:,2:]).det()* (M[:,0:2]).det() #P3 is the pruduct of minors 

#P3=sp.expand(a1*a4*(2*t1*x1+t1*x2-x2)*(-t2**2*x2+2*b*t2-2*t2*x1+x2))
B=[ft.arb(5,2),ft.arb(5,2),ft.arb(0,1),ft.arb(0,1)]   #the defining domain
P=[P1,P2,P3]      
P=[poly_normal_to_list(Pi,X) for Pi in P]   #representing P as a list  
jac=jacobian_of_function_list(P) # the Jacobian matrix of P



Ball=(Ball_interval(P))
B_Ball=B_Ball_calculator(B)
"""






"""
x1= sp.Symbol('x1')
x2= sp.Symbol('x2')
q1= sp.Symbol('q1')
q2= sp.Symbol('q2')
X=[x1,x2,q1,q2]
P1=sp.expand(x1**2+x2**2-q1**2)
P2=sp.expand((x1-9)**2+x2**2-q2**2)
P3=x2
P=[P1,P2,P3]      
P=[poly_normal_to_list(Pi,X) for Pi in P]   #representing P as a list  
jac=jacobian_of_function_list(P) # the Jacobian matrix of P
B=[ft.arb(0,20),ft.arb(0,20),ft.arb(4,2),ft.arb(6.5,2.5)]   #the defining domain
"""



#print(len(no_rep))

#Ball=Ball_interval(L)   # computing the ball system
#jac_Ball=jacobian_of_function_list(Ball)
#B_Ball=B_Ball_calculator(B) 


#M=checking_regularity_of_sysyem(Ball,B_Ball,jac_Ball,0.3)
#print(len(M[0]),len(M[1]))

#checking the smoothness of the curve:







"""
#Example2
P1=sp.expand(x1-t2**2)
P2=sp.expand(x2-t2**3)
P3=sp.expand(t1-t2)
U=[ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1)]
B=[ft.arb(0,1),ft.arb(0,1),ft.arb(0,1),ft.arb(0,1)]
"""
import pickle


"""

recall=open("list_midpoints.txt","rb")
T_mid=pickle.load(recall)

recall=open("list_raduces.txt","rb")
T_rad=pickle.load(recall)

regular=[[ft.arb(T_mid[i][j],T_rad[i][j]) for j in range(len(T_mid[0]))]
                      for i in range(len(T_mid)) ]
proj=[[regulari[0],regulari[1]] for regulari in regular  ]


recall=open("overlapped_boxes_2D_mid.txt","rb")
T_mid=pickle.load(recall)
recall=open("overlapped_boxes_2D_rad.txt","rb")
T_rad=pickle.load(recall)
union_overlopping_boxes= [[ft.arb(T_mid[i][j],T_rad[i][j]) for j in range(len(T_mid[0]))]
                      for i in range(len(T_mid)) ]
out_list=[]
for Ui in union_overlopping_boxes:
    U=B_Ball_calculator(Ui)
    P=[poly_normal_to_list(Pi,X) for Pi in [P1,P2,P3]]
    Ball=Ball_interval(P)
    Jac_Ball=jacobian_of_function_list(Ball)
    T=checking_regularity_of_sysyem(Ball,U,Jac_Ball,wth=0.1)
    if T != 1 :
      out_list=1
print(len(out_list))
"""

"""

Bi_pickled= open("overlapped_boxes_2D_rad.txt","wb")
pickle.dump(T_rad,Bi_pickled)
Bi_pickled.closed()


Bi_pickled= open("list_midpoints.txt","wb")
pickle.dump(T_mid,Bi_pickled)
Bi_pickled.closed()

Bi_pickled= open("list_raduces.txt","wb")
pickle.dump(T_rad,Bi_pickled)
Bi_pickled.closed()

T=overlapping_boxes_2d(regular)
T_mid=[[float(Tij.mid()) for Tij in Ti ] for Ti in T]
T_rad=[[float(Tij.mid()) for Tij in Ti ] for Ti in T]
overlapped_mid_p= open("overlapped_mid.txt","wb")
pickle.dump(T_mid,overlapped_mid_p)
Bi_pickled.closed()
#T=curve_tracer(P,B,jac,0.00005)
#T_mid=[[Tij.mid() for Tij in Ti] for Ti in T]
#T_rad=[[Tij.rad() for Tij in Ti] for Ti in T]
"""

#print([round(float(Ti.upper()),3) for Ti in T])


#B=[ft.arb(0.03,3.15),ft.arb(0.01,3.15),ft.arb(0.02,1),ft.arb(0.02,1)]
#Checking Assumption 1 :
#print(checking_smoothness(P,B,Jac,0.01))












"""
projection0= [[T[0],T[1]] for T in unknownboxes]
projection=[]
bad_ind=[]
for i in range(len(projection0)):
    if i not in bad_ind:
      for j in range(i+1,len(projection0)):
          membership1=1
     
      for k in range(len(projection0[0])):
            if not projection0[i][k].contains(projection0[j][k]) :
                membership1=0
                break
          if membership1==1:
              projection.append(projection0[i])
              bad_ind.append(j)
          membership2=1
          for k in range(len(projection0[0])):
            if not projection0[j][k].contains(projection0[i][k]) :
                membership2=0
                break
          if membership1==0 and membership2==1:
              projection.append(projection0[j])
"Ho")

print(len(projection))




"""



#print(max([ T.rad() for T in  subdivide(U)[0]]))

#print(max([ T.rad() for T in  subdivide(subdivide(U)[0])[0]]))



#B=[ft.arb(0,5),ft.arb(0,5),ft.arb(0,5),ft.arb(0,5),ft.arb(0,5)]
#P=[[[[1,0],1],[[0,2],-1]],[[[0,1],1],[[3,0],-1]]]
#P=[[[[2,0],1],[[0,1],-1]],[[[0,1],1],[[0,0],-2]]]
#P=[[[[3,0],1],[[2,0],-6],[[1,0],11],[[0,0],-6],[[0,1],-1]],[[[0,1],1]]]
#P=[[[[1,0,0],1],[[0,0,3],-1]],[[[0,1,0],1],[[0,0,2],-1]]]
#Ball=Ball_interval(P)
#M=coefficient_matrix_list(Ball[4]) #a problem in polyvalnd for this example
#print(polyvalnd(M,B))
#jac=jacobian_of_function_list(Ball)
#T= solver(Ball,jac,B)
#print(T)

#print([[Bij.mid() for Bij in Bi]  for Bi in subdivide(B)])
#print([Ti.lower() for Ti in  T[0]])
#print([Ti.upper() for Ti in  T[0]])
#M_list=[coefficient_matrix_list(Pi) for Pi in P]
#print([polyvalnd(Mi, T) for Mi in M_list])
#print([T[0].upper() for T in solver(P,jac,B)])
#print(jac_B)
#print(hansen_hengupta(x_teld,jac_B,b,B,B)[0].lower(), hansen_hengupta(x_teld,jac_B,b,B,B)[0].upper() )
#print(B[0].lower(),B[0].upper())
#print([ T.lower() for T in hansen_hengupta(x_teld,jac_B,b,B,B)])
#print([ T.upper() for T in hansen_hengupta(x_teld,jac_B,b,B,B)])
#print(hansen_hengupta(x_teld,jac_B,b,B,B)[0] in B[0])
#print(hansen_hengupta(x_teld,jac_B,b,B,B))
#print(B[0] in hansen_hengupta(x_teld,jac_B,b,B,B )[0])
#example checking Assumption 1
#Example 1
"""
P1= [[[[1,0,0],1],[[0,0,2],-1],[[0,0,0],1] ],   [[[0,1,0],1],[[0,0,3],-1],[[0,0,1],1] ]  ] # the curve x-z^3=y-z^2 =0
Ball=Ball_interval(P1)
b=[0,0,0,0,0]
x=[ft.arb(5,0.1),ft.arb(6,0.1),ft.arb(-7.1,0.3),ft.arb(9,0.1),ft.arb(8,1.5)]
z=[ft.arb(5,0.1),ft.arb(6,0.1),ft.arb(-7.1,0.3),ft.arb(9,0.1),ft.arb(8,1.5)]
jac_Ball=jacobian_of_function_list(Ball)
eval_jac_Ball_x=matrixval(jac_Ball,x)
x_teld=[0,0,-1,1.01,0]
#print(hansen_hengupta(x_teld,eval_jac_Ball_x,b,x,z)[0] in x[0] )


#Example 2
P2=[[[[1,0,0],1],[[0,0,3],-1] ],   [[[0,1,0],1],[[0,0,2],-1] ]  ]
Ball=Ball_interval(P2)
b=[0,0,0,0,0]
x=[ft.arb(0,0.001),ft.arb(0,0.001),ft.arb(0,0.0003),ft.arb(1,0.1),ft.arb(0,0.005)]
z=[ft.arb(0,0.001),ft.arb(0,0.001),ft.arb(0,0.0003),ft.arb(1,0.1),ft.arb(0,0.005)]
jac_Ball=jacobian_of_function_list(Ball)
eval_jac_Ball_x=matrixval(jac_Ball,x)
x_teld=[0,0,0,1.01,0]
#print(hansen_hengupta(x_teld,eval_jac_Ball_x,b,x,z)[1] in x[1])"""


"""







#eva_jac_Ball=matrixval(jac_Ball,B)

#print(jac_Ball[0][0])
#print(subdivide(B)[0][0].mid())
#example Hansen-Sengupta operator
#b=[3,3,3]
#x=[ft.arb(1,1),ft.arb(1,1),ft.arb(1,1)]
#x_teld=[0,0,0]
#z=[ft.arb(1,1),ft.arb(1,1),ft.arb(1,1)]
#A=np.array([[16,3],[7,-11]])
#A=[[ft.arb(3,1),ft.arb(1,1),ft.arb(1,1)],[ft.arb(1,1),ft.arb(3,1),ft.arb(1,1)],[ft.arb(1,1),ft.arb(1,1),ft.arb(3,1)]]
#print(hansen_hengupta(x_teld,A,b,x,z))
#
#t4=gauss_seidel_dimn(A,b,t3,t3)
#print(hansen_hengupta(x_teld,A,b,x,z))
#example Ball_interval(L)
#L=[[[[1,0,0],1],[[0,0,3],-1]],[[[0,1,0],1],[[0,0,2],-1]]] # L is the curve x1-x^3, x2-x^2
#B=Ball_interval(L)  #the Ball system of that curve
#print(B)
#print(polyvalnd(coefficient_matrix_list(B[0]),[0,0,0,1,0]))  checking the evaluation



#print(F)
#print(Ball_for_monmial(L[1],X)[0][0][0])
#Example of invertibility_of_a_matrix
#M=ft.arb_mat([[ft.arb(0,1),ft.arb(-2,3),ft.arb(-2,5)],[ft.arb(0,1),ft.arb(-7,-5),ft.arb(-1,8)],[ft.arb(0,1),ft.arb(-2,3),ft.arb(0,1)]])
#N=ft.arb_mat([[ft.arb(1,0.01),ft.arb(0,0.01)],[ft.arb(0,0.01),ft.arb(1,0.01)]])
#K=ft.arb_mat([[ft.arb(1,2),ft.arb(0,2)],[ft.arb(0,2),ft.arb(1,3)]])
#R=invertibility_of_a_matrix(K)
#print(R,invertibility_of_a_matrix(K))
#print(R)
#L=[[[0,0,0,3],2],[[0,1,0,3],2]]
#print(Ball_for_interval_poly(L))
#M=coefficient_matrix(f)
#M=[[1, 2, 3], [4, 5, 6], [6, 7, 8], [9, 10, 11]]
#M=polyvalnd(M,[ft.arb(1,1),ft.arb(1,1),ft.arb(1,1)])
#M=empty_list_generator(1,3)
#reaching_the_elements_of_list([1,1,1],M,ft.arb(2,3))
#print(M)"""
