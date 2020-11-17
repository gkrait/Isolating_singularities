import math
import matplotlib.pyplot as plt
import os
import pickle 
import interval_arithmetic as d

from pprint import pprint
from sympy.parsing.sympy_parser import parse_expr
import sympy as sp 
import os 
from tcusp import cusp_Ball_solver, evaluation_exp

import matplotlib.patches as mpatches
import csv
from scipy import spatial
import flint as ft
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import itertools
import timeit
import time
import numpy as np
from evaluation_file1 import evlist
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components




def ploting_boxes(boxes,uncer_boxes, var=[0,1], B=[[-20,20],[-20,20]],x=0.1,nodes=[], cusps=[],uncer_Solutions=[],Legend=False,color="green",variabel_name="x" ):
   fig, ax = plt.subplots()
   #plt.grid(True)
   ax.set_xlim(B[0][0], B[0][1])
   ax.set_ylim(B[1][0], B[1][1])
   ax.set_xlabel(variabel_name+str(1))
   ax.set_ylabel(variabel_name+str(2))
   """try:
    ax.title(open("system.txt","r").read())
   except:
    pass"""
   
   #textstr = open("system.txt","r").read()
   #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
   #ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
   #     verticalalignment='top', bbox=props)
   c=0
   green_patch = mpatches.Patch(color=color, label='smooth part') 
   red_patch = mpatches.Patch(color='red', label='unknown part')
   node_patch = mpatches.Patch(color='black', label='Certified nodes',fill=None)
   cusp_patch = mpatches.Patch(color='blue', label='Projection of certified solution with t=0 ',fill=None)
   if Legend==True:
     plt.legend(handles=[green_patch,red_patch,node_patch,cusp_patch])
   for box in boxes:
     rectangle= plt.Rectangle((box[var[0]][0],box[var[1]][0]) , \
      (box[var[0]][1]-box[var[0]][0]),(box[var[1]][1]-box[var[1]][0]),color=color)
     plt.gca().add_patch(rectangle)
   for box in uncer_boxes:

      rectangle= plt.Rectangle((box[var[0]][0],box[var[1]][0]) , \
      (box[var[0]][1]-box[var[0]][0]),(box[var[1]][1]-box[var[1]][0]), fc='r')
      plt.gca().add_patch(rectangle)
   for box in nodes:
     rectangle= plt.Rectangle((box[0][0]-x,box[1][0]-x) ,\
      2*x+box[0][1]-box[0][0],2*x+box[1][1]-box[1][0], fc='y',fill=None)
     plt.gca().add_patch(rectangle) 
   for box in cusps:
     rectangle= plt.Rectangle((box[0][0]-x,box[1][0]-x) ,\
      2*x+box[0][1]-box[0][0],2*x+box[1][1]-box[1][0], fc='y',color="blue",fill=None)
     plt.gca().add_patch(rectangle) 
   for box in uncer_Solutions:
     rectangle= plt.Rectangle((box[0][0]-x,box[1][0]-x) ,\
      2*x+box[0][1]-box[0][0],2*x+box[1][1]-box[1][0], fc='y',color="red",fill=None)
     plt.gca().add_patch(rectangle)   
   plt.savefig("fig.jpg",dpi=1000) 
   plt.show()
def Ball_node_gen(equations,B_Ball,X):
    P=open(equations,"r").readlines()
    P=[Pi.replace('\n','') for Pi in P]
    n=len(X)
    V=""" Variables \n """
    for i in range(n):
        V += "x" +str(i+1) + " in " + str(B_Ball[i]) +" ; \n"
    for i in range(n,2*n-2):
        V += "r" +str(i-n+3) + " in " + str(B_Ball[i]) +" ; \n" 
    V += "t" + " in " + str(B_Ball[2*n-2]) +" ; \n"       
    V +="Constraints \n"   
    for Pi in P:
        t1=time.time();H=SDP_str(Pi,X);print(time.time()-t1)
        V += H[0]
        V += H[1]
    last_eq=""
    for i in range(3,n):
        last_eq += "r"+str(i)+"^2+"
    last_eq += "r" +str(n)+"^2 -1=0;"    
    V += last_eq +"\n"
    f= open("eq.txt","w+")
    f.write(V) 
    f.write("end")
    f.close() 
def Ball_solver(equations,B_Ball,X):     #the width condition needs to be added  Do not suse this one 
  L=[B_Ball]
  certified_boxes=[]
  uncertified_boxes=[]
  n=len(X)
  while len(L) !=0: 
    solvability=1
    if B_Ball[2*n-2][0] <= 0 <=   B_Ball[2*n-2][1] and \
    d.width([ d.ftconstructor(Bi[0],Bi[1]) for Bi in L[0] ] ) <0.1 :
      Ball_cusp_gen(equations,B_Ball,X,eps_min=eps_min,eps_max=eps_max)
    elif (B_Ball[2*n-2][0] > 0 or 0 >   B_Ball[2*n-2][1] ) \
    and d.width([ d.ftconstructor(Bi[0],Bi[1]) for Bi in L[0] ] ) <0.1:
      Ball_node_gen(equations,B_Ball,X)
    else:
      children=cb.plane_subdivision(L[0])
      L.remove(L[0])
      L += children
      solvability=0
    if solvability==1:
      ibex_output=cb.solving_with_ibex()
      if ibex_output[0]== "Empty":
        
        L.remove(L[0])
      elif len(ibex_output[0]) !=0:  
        
         certified_boxes +=cb.computing_boxes(ibex_output[0])
         L.remove(L[0])
      elif len(ibex_output[1])!=0:   
        
        uncertified_boxes +=cb.computing_boxes(ibex_output[1])
        L.remove(L[0])
      else:  
        children=cb.plane_subdivision(L[0])
        L.remove(L[0])
        L += children
    
  return [certified_boxes,uncertified_boxes]      
def SDP_str(P,X):
    n=len(X)
    P_pluse=P
    P_minus=P
    for i in range(2,n):
        P_pluse=P_pluse.replace("x"+str(i+1),"(x"+str(i+1) + "+ r"+str(i+1) +"*sqrt(t))")
        P_minus=P_minus.replace("x"+str(i+1),"(x"+str(i+1) + "- r"+str(i+1) +"*sqrt(t))")
    
    if parse_expr(P.replace("^","**")).is_polynomial():
     SP= str(sp.simplify(sp.expand("0.5*(" + P_pluse + "+" +P_minus+") \n"))); SP=(SP+"=0; \n").replace("**","^")
     DP=str(sp.simplify(sp.expand("0.5*(" + P_pluse + "- (" +P_minus+") )/(sqrt(t)) \n"))); DP=(DP+"=0; \n").replace("**","^")
    else:
     SP= ("0.5*(" + P_pluse + "+" +P_minus+")  =0; \n").replace("**","^")
     DP= ("(" + P_pluse + "- (" +P_minus+")) =0;  \n").replace("**","^")
    return [SP,DP]
def Ball_generating_system(P,B_Ball,X,eps_min=0.001):
    n=len(X)
    V=""" Variables \n """
    for i in range(n):
        if B_Ball[i][0] != B_Ball[i][1]:
         V += "x" +str(i+1) + " in " + str(B_Ball[i]) +" ; \n"
        else: 
         V += "x" +str(i+1) + " in " + str([B_Ball[i][0]-eps_min, B_Ball[i][1]+eps_min]) +" ; \n"
    for i in range(n,2*n-2):
        V += "r" +str(i-n+3) + " in " + str(B_Ball[i]) +" ; \n" 
    V += "t" + " in " + str(B_Ball[2*n-2]) +" ; \n"       
    V +="Constraints \n"   
    for Pi in P:
        V +=Pi

    f= open("eq.txt","w+")
    f.write(V) 
    f.write("end")
    f.close()
def intersting_boxes1(f,b):
    pickle_in=open(f,"rb")
    curve=pickle.load(pickle_in)
    pickle_in.close()
    intersting_boxes=[]
    uncer_boxes=[]
    for box in curve[0]:
        if  b[0][0] <= box[0][0] <= box[0][1] <=b[0][1] and \
          b[1][0] <= box[1][0] <= box[1][1] <=b[1][1]:
             intersting_boxes.append(box)
    for box in curve[1]:
        if  b[0][0] <= box[0][0] <= box[0][1] <=b[0][1] and \
          b[1][0] <= box[1][0] <= box[1][1] <=b[1][1]:
             uncer_boxes.append(box)
    return [intersting_boxes,uncer_boxes]  
def intersting_boxes(curve,b):
    cer_intersting_boxes=[]
    uncer_intersting_boxes=[]
    for box in curve[0]:
        if  b[0][0] <= box[0][0] <= box[0][1] <=b[0][1] and \
          b[1][0] <= box[1][0] <= box[1][1] <=b[1][1]:
             cer_intersting_boxes.append(box)
    for box in curve[1]:
        if  b[0][0] <= box[0][0] <= box[0][1] <=b[0][1] and \
          b[1][0] <= box[1][0] <= box[1][1] <=b[1][1]:
             uncer_intersting_boxes.append(box)
    return [cer_intersting_boxes,uncer_intersting_boxes]  
def ibex_output(P,B,X):
    os.system("ibexsolve   --eps-max=0.1 -s  eq.txt  > output.txt")
    g=open('output.txt','r')
    result=g.readlines()
    T=computing_boxes(result)

    return T    
def estimating_t1(components,upper_bound=200000):  #it works only if len(components)
  t1=upper_bound
  t2=0
  for box1 in components[0]:
    for box2 in components[1]:
        a=d.distance(box1,box2).lower()
        b=d.distance(box1,box2).upper()
    if t1 > a:
     t1=a 
    if t2<b:
     t2=b 
  t=d.ftconstructor(t1,t2)
  t=0.25*d.power_interval(t,2)     

  return [float(t.lower()),float(t.upper())]     
def estimating_t(components):  #it works only if len(components)==2
  uni=d.distance(components[0][0],components[1][0])
  
  for box1 in components[0]:
    for box2 in components[1]:
        #print(uni);input()
        if d.width(box1[2:])< 1e-6:
          box1[2:]=[[Bi[0] - 1e-6, Bi[1]+1e-6 ]  for Bi in   box1[2:]]
        if d.width(box2[2:])< 1e-6:
          box2[2:]=[[Bi[0] - 1e-6, Bi[1]+1e-6 ]  for Bi in   box2[2:]]  
        a=d.distance(box1[2:],box2[2:])
        #print(a);input()
        uni=d.box_union([uni],[a])[0]
        """if t1 > a[0]:
          t1=a[0]
        if t2<a[1]:
           t2=a[1]"""  
  t=d.ftconstructor(uni[0],uni[1])
  t=0.25*d.power_interval(t,2) 
  return [float(t.lower()),float(t.upper())] 

def boxes_compare(box1,box2):
    flage=0
    for i in range(len(box1)-1,-1,-1):

        if box1[i][0] > box2[i][0]: 
            return 1
        if  box1[i][0] < box2[i][0]: 
          return -1
    return 0      
def boxes_sort(boxes):
    sorted_boxes=boxes[:]
    for i in range(len(boxes)-1):
        for j in range(i+1,len(boxes)):
            if boxes_compare(sorted_boxes[i],sorted_boxes[j]) ==1:
                sorted_boxes[i], sorted_boxes[j] =sorted_boxes[j], sorted_boxes[i]
    return sorted_boxes            
def connected_compnants(boxes):
    #ftboxes=[ [d.ftconstructor(boxi[0],boxi[1]) for boxi in box ] for box in boxes ]
    ftboxes=boxes[:]
    components=[[ftboxes[0]]]
    for i in range(1,len(ftboxes)):
        boxi_isused=0
        for j in  range(len(components)):
            membership=0
            for k in range(len(components[j])):   
                if d.boxes_intersection(ftboxes[i],components[j][k]) !=[] :
                    components[j].append(ftboxes[i])
                    membership=1
                    boxi_isused=1
                    break
            if membership==1:
              break 
        if boxi_isused==0:
          components.append([ftboxes[i]])
    unused=list(range(len(components)))
    components1=components[:]
    components2=[]
    while len(components1) != len(components2) :  
        for i in unused:
            for j in   [j for j in list(range(i+1,len(components))) if j in unused ]:
                intersection_exists=False
                is_looping=True
                for boxi in components[i]:
                    for boxj in components[j]:
                        if d.boxes_intersection(boxi,boxj)!=[]:
                            
                            is_looping = False
                            intersection_exists=True
                            break
                    if is_looping==False:
                       break
                if intersection_exists== True:
                    components[i] += components[j]
                    unused.remove(j)

        components2=components1[:]
        components1=[components[k] for k in unused ]
                            
    return components1                        

def planner_connected_compnants(boxes): 
    if len(boxes)==0:
      return []
    ftboxes=boxes[:]
    #ftboxes=[ [d.ftconstructor(boxi[0],boxi[1]) for boxi in box ] for box in boxes ]
    components=[[ftboxes[0]] ]
    for i in range(1,len(ftboxes)):
        boxi_isused=0
        for j in  range(len(components)):
            membership=0
            for k in range(len(components[j])):   
                if d.boxes_intersection(ftboxes[i][:2],components[j][k][:2]) !=[]: # and \
                #d.boxes_intersection(ftboxes[i],components[j][k]) ==[]:
                    components[j].append(ftboxes[i])
                    membership=1
                    boxi_isused=1
                    break  
            if membership==1:
              break 
        if boxi_isused==0:
          components.append([ftboxes[i]])
        
    unused=list(range(len(components)))
    components1=components[:]
    components2=[]
    while len(components1) != len(components2) :
        for i in unused:
            for j in   [j for j in list(range(i+1,len(components))) if j in unused ]:
                intersection_exists=False
                is_looping=True
                for boxi in components[i]:
                    for boxj in components[j]:
                        if d.boxes_intersection(boxi[:2],boxj[:2])!=[] :#and \
                        #d.boxes_intersection(boxi[:2],boxj[:2]) != []  :
                            is_looping = False
                            intersection_exists=True
                            break
                    if is_looping==False:
                       break
                if intersection_exists== True:
                    components[i] += components[j]
                    unused.remove(j)
            components2=components1[:]
            components1=[components[k] for k in unused ]                    
    
    return components1                        
def estimating_yandr(components,upper_bound=100000):
  r_bounds=[[upper_bound,0]]*(len(components[0][0])-2)
  r_list=[]
  y_list=[]
  for box1 in components[0]:
    for box2 in components[1]:
        ft_box1= [d.ftconstructor(Bi[0],Bi[1]) for Bi in box1  ]
        ft_box2= [d.ftconstructor(Bi[0],Bi[1]) for Bi in box2  ]
        
        y_list.append([0.5*(q1+q2) for q1,q2 in zip(ft_box1[2:],ft_box2[2:])])
        norm_q1q2=d.distance(box1[2:],box2[2:])
        norm_q1q2=d.ftconstructor(norm_q1q2[0],norm_q1q2[1])
        q1q2=[ft_box1[i]-ft_box2[i] for i in range(2,len(box1)) ]
        
        r=[ ri/norm_q1q2 for ri in q1q2 ]
        r_list.append(r)
  r=[]
  y=[]
  for i in range(len(y_list[0])):
    yi1=min([float(y[i].lower()) for y in y_list  ])
    yi2=max([float(y[i].upper()) for y in y_list  ])
    y.append([yi1,yi2])
  for i in range(len(r_list[0])):
        ri1=min([float(r[i].lower()) for r in r_list  ])
        ri2=max([float(r[i].upper()) for r in r_list  ])
        r.append([ri1,ri2])
  r=d.boxes_intersection(r,[[-1,1]]*(len(components[0][0])-2))          

  return y+r       
def detecting_nodes(boxes,B,f,X,eps): #boxes are list of cer and uncer curve
    mixes_boxes= [[1,box ] for box in boxes[0] ] +[[0,box ] for box in boxes[1]] #putting flaggs for cer and uncer boxes
    ftboxes=[ [box[0], [d.ftconstructor(boxi[0],boxi[1])  for boxi in box[1]] ] for box in mixes_boxes ]    
    nodes_lifting=[]
    used=[]
    P=[ Pi.replace("\n","") for Pi in   open(f,"r").readlines()  ]
    for i in range(len(ftboxes)):
        for j in range(i+1,len(ftboxes)):
            Mariam_ft=d.boxes_intersection(ftboxes[i][1],ftboxes[j][1])
            Mariam=[[float(Bi.lower()),float(Bi.upper()) ] for Bi in Mariam_ft]
            if (Mariam ==[] and \
             d.boxes_intersection(ftboxes[i][1][:2],ftboxes[j][1][:2])) or\
               (Mariam != [] and enclosing_curve(f,Mariam,X,eps_max=0.1) ==[[],[]] ): #needs to work more
                if i not in used:
                    used.append(i)
                    nodes_lifting.append(ftboxes[i])
                if j not in used:
                    used.append(j)
                    nodes_lifting.append(ftboxes[j])

    components= planner_connected_compnants(nodes_lifting)
    cer_components=[]
    uncer_components=[]
    component_normal=[]
    for component in components:
        boxes_component=[box[1] for box in component]
        component_normal =[ [[ float(Bi.lower()),  float(Bi.upper()) ] for Bi in box[1] ] for box in component ]
        if 0  not in [ box[0] for box in  component]  and eval_file_gen(f,component_normal,X) =="[]\n" :
            cer_components.append(boxes_component)
        else: 
            uncer_components.append(boxes_component)
    return [cer_components,uncer_components]         
def intersect_in_2D(class1,class2,monotonicity=1):
  pl_intesected_pairs=[]
  if monotonicity==1:
    for i in range(len(class1)):
      for j in range(len(class2)):
        if d.boxes_intersection(class1[i][:2],class2[j][:2]) !=[] and   d.boxes_intersection(class1[i],class2[j]) ==[] :
          if [class2[j],class1[i]] not in pl_intesected_pairs:
            pl_intesected_pairs.append([class1[i],class2[j]])
  elif  monotonicity==0:       
       for i in range(len(class1)):
         for j in range(len(class2)):
           if d.boxes_intersection(class1[i][:2],class2[j][:2]) !=[]:
             if [class2[j],class1[i]] not in pl_intesected_pairs:
                 pl_intesected_pairs.append([class1[i],class2[j]])
  elif  monotonicity==2:       
       inters_indic=[]
       for i in range(len(class1)):
         inters_indic.append([])
         for j in range(len(class2)):
           if d.boxes_intersection(class1[i][:2],class2[j][:2]) !=[]:
             inters_indic[i]=  inters_indic[i] +[j] 
       for k in range(len(class1)):
         if len(inters_indic[k])> 3:
          for j in range(len(inters_indic[k])):
            if  [class2[j],class1[k]] not in pl_intesected_pairs:
              pl_intesected_pairs.append([class1[k], class2[j]])
       
                      
  return pl_intesected_pairs     
def solving_fornodes(equations,boxes,B,X,eps=0.1):
    plane_components=detecting_nodes(boxes,B,equations,X,eps)#[0]
    g=open(equations,'r')
    P=[ Pi.replace("\n","") for Pi in   g.readlines()  ]
    Ball_solutions=[]
    for plane_component in plane_components:
        x1=float(min([ai[0].lower() for ai in plane_component]))
        x2=float(max([ai[0].upper() for ai in plane_component]))
        y1=float(min([ai[1].lower() for ai in plane_component]))
        y2=float(max([ai[1].upper() for ai in plane_component]))
        components=connected_compnants(plane_component)
        r=[ [float(ri[0]),float(ri[1])] for ri in    estimating_r(components)   ]
        t=estimating_t(components)
        t=[float(t[0]),float(t[1])]
        B_Ball=[[x1,x2],[y1,y2]]+r +[t]
        Ball_generating_system(P,B_Ball,X)
        solutionsi=ibex_output(P,B_Ball,X)
        Ball_solutions +=solutionsi
    return Ball_solutions
def normal_subdivision(B):
  ft_B=d.subdivide([d.ftconstructor(Bi[0],Bi[1]) for Bi in B[:]])
  return [d.ft_normal(Bi)  for Bi in ft_B]
def plane_subdivision(B):
  
  ft_B2=d.subdivide([d.ftconstructor(Bi[0],Bi[1]) for Bi in B[:2]])
  normal_B2=[d.ft_normal(Bi)  for Bi in ft_B2]
  return d.cartesian_product(normal_B2,[B[2:]])
def system_generator(f,B,X):
    g = open(f, "r")
    L = g.readlines()
    g.close()
    f = open("eq.txt", "w+")
    f.write("Variables \n")
    for i in range(len(X)):
        f.write(str(X[i]) + " in " + str(B[i]) + " ; \n")
    f.write("Constraints \n")
    for Li in L:
        f.write(Li.replace("\n", "") + "=0; \n")
    f.write("end ")
    f.close()

    return f
def solving_with_ibex(eps=0.1):
  uncer_content=[]
  cer_content=[]
  os.system("ibexsolve   --eps-max="+ str(eps) +" -s  eq.txt  > output.txt")
  g=open('output.txt','r')
  result=g.read()
  with open('output.txt') as f:
    if "successful" in result:
      cer_content = f.readlines()
    elif  "infeasible" not in result and "done! but some boxes" in result:
      uncer_content = f.readlines()
    elif "infeasible problem" in result:
      uncer_content="Empty"
      cer_content="Empty"
  return [cer_content,uncer_content]      
def computing_boxes():
  content=open("output.txt","r").readlines()
  for line in content[:]:
    if "number of cells" in line:
      
      tree_size= 2*int(line.split()[3])+1 #to compute tree size
      
      break
  if "infeasible" in open("output.txt","r").read():
    return ["Empty","Empty",tree_size]
  cer=[]; uncer=[]
  i=0
  Answer=[]
  for fi in content:
    try:
      a=fi.index('(')
      b=fi.index(')')
      T=(fi[a:b+1]).replace('(','[')
      T=(fi[a:b+1]).replace('(','[')
      T=T.replace(')',']')
      T=T.split(";")
      E=[]
      i=0
      for Ti in T:
        Ti= Ti.replace('[',"")
        Ti= Ti.replace(']',"")
        Ti=Ti.replace('<','')
        Ti=Ti.replace('>','')
        x=Ti.index(",")
        a=float(Ti[:x])
        b=float(Ti[x+1:])
        E.append([])
        E[i]=[a,b]
        i+=1
      if "solution n" in fi or "boundary n" in fi:
        cer.append(E)
      elif "unknown n" in fi:
        uncer.append(E)
    except ValueError:
          pass 
  if len(cer)==0:
    return [[],uncer,0] 
  return [cer,uncer,tree_size]        
def enclosing_curve(system,B,X,eps_min=0.1,eps_max=0.1, info=False): 
  start_enc=time.time()
  tree_size=0
  L=[B]
  certified_boxes=[]
  uncertified_boxes=[]
  while len(L) !=0: 
    system_generator(system,L[0],X)
    os.system("ibexsolve   --eps-max="+ str(eps_max)+"  --eps-min="+ str(eps_min) + " -s  eq.txt  > output.txt")
    content=open("output.txt","r").readlines()
    res=computing_boxes()

    ibex_output=res[:2]
    tree_size += res[2]
    if ibex_output ==[[],[]] and max([Bi[1]-Bi[0] for Bi in L[0]  ]) < eps_min :  
      uncertified_boxes.append(L[0])
      L.remove(L[0]);
    elif ibex_output ==[[],[]] :
      children=normal_subdivision(L[0]) 
      #children=plane_subdivision(L[0])
      L.remove(L[0]);
      L += children  # print warning ################################################################""
      print("warning!")
      tree_size +=1

    elif ibex_output[1]== "Empty":
      L.remove(L[0])
    else:

      if len(ibex_output[0]) !=0:
       certified_boxes += ibex_output[0]
      if len(ibex_output[1])!=0: 
       uncertified_boxes += ibex_output[1]
      L.remove(L[0])  
  if len(certified_boxes) ==0:
    return [[],uncertified_boxes,tree_size]   
  end_enc=time.time()
  if info==True:
    print("Timing of enclosing curve ", end_enc-start_enc )      
  return [certified_boxes,uncertified_boxes,tree_size]                        
def loopsfree_checker(f,certified_boxes,uncer_boxes,P): #Assumption: no cusps
  L=eval_file_gen(f,certified_boxes,X)
  while L.replace('\n',"") != "[]":
    L=L.replace('[','')
    L=L.replace(']','')
    L=L.replace('\n','')
    L=L.split(",")
    for i in L:
      children=normal_subdivision(certified_boxes[int(i)])
      certified_boxes.remove(certified_boxes[int(i)])
      for child in children:
        cer_children, uncer_children= enclosing_curve(f,child,X)
        certified_boxes +=cer_children
        uncer_boxes +=uncer_children
    L =  eval_file_gen(f,certified_boxes,X)
  return L     
def eval_file_gen(f,boxes,X,special_function=[]): #condition: len(boxes[0]) is even
  functions=["sin","cos","tan","exp"]+special_function
  if len(boxes[0])==0:
    return []
  n=len(boxes[0])
  m=len(boxes)
  g=open(f,'r')
  P_str=g.readlines()
  P_str= [Pi.replace('\n','') for Pi in P_str]
  P_str= [Pi.replace('^','**') for Pi in P_str]
  P_exp= [parse_expr(Pi) for Pi in P_str]
  #computing jac and the minors
  jac=sp.Matrix(P_str).jacobian(sp.Matrix(X))
  minor1=jac[:,1:].det()
  minor2=jac[:,[i for i in range(n) if i != 1]  ].det()
  fil=open("evaluation_file1.py","w")
  fil.write("import flint as ft \n")
  fil.write("import sympy as sp \n")
  fil.write("import interval_arithmetic as d \n")
  fil.write("def  evlist(boxes): \n")
  #fil.write("boxes="+str(boxes)+"\n")
  fil.write(" ftboxes=[ [d.ftconstructor(Bi[0],Bi[1]) for Bi in B ] for B in boxes ] \n"    )
  fil.write(" n=len(boxes[0])\n")
  fil.write(" m=len(boxes)\n")
  fil.write(" m1=[]\n")
  fil.write(" m2=[]\n")
  minor1_str=str(minor1)
  minor2_str=str(minor2)
  for i in range(n):
    minor1_str= minor1_str.replace("x"+str(i+1),"B["+str(i)+"]" )
    minor2_str= minor2_str.replace("x"+str(i+1),"B["+str(i)+"]" )
  for func in functions:
    minor1_str=minor1_str.replace(func,"ft.arb."+func)
    minor2_str=minor2_str.replace(func,"ft.arb."+func)
  fil.write(" for B in ftboxes: \n")
  fil.write("   m1.append(ft.arb("+ minor1_str + ")) \n")
  fil.write("   m2.append( ft.arb("+ minor2_str + ")) \n")  
  fil.write(" innrer_loops=[i for i in range(m) if 0 in m1[i] and 0 in m2[i] ]\n")
  fil.write(" x1mon=[i for i in range(m) if 0  not in m1[i] ]\n")
  fil.write(" x2mon=[i for i in range(m) if 0 in m1[i] and 0 not in m2[i] ]\n")
  fil.write(" return [x1mon,x2mon , innrer_loops] \n")
  fil.close()
  
  #t=os.popen("python3 evaluation_file1.py ").read()
  return evlist(boxes)
def boxes_classifier(system,boxes,X,special_function=[]):
  if len(boxes[0])==0:
    return [[],[],boxes[1]]
  certified_boxes ,uncer_boxes, tree_size =boxes
  L=eval_file_gen(system,certified_boxes,X)
  if L==[[],[],[]]:
    return [[],[],[],uncer_boxes]
  else:
      return  [ [certified_boxes[i]  for i in range(len(certified_boxes)) if i in L[0]] ,\
      [certified_boxes[i]  for i in range(len(certified_boxes)) if i in L[1]],\
      [certified_boxes[i]  for i in range(len(certified_boxes)) if i in L[2]], uncer_boxes ] #can be enhanced
def projection_checker(solutions):
  if len(solutions)==0:
    return [[],[]] 
  m=len(solutions[0])
  n=int((m+1)/2)
  rep=[]
  unaccepted=[]
  intersect_in2d=[[] for i in range(len(solutions))]
  for i in range(len(solutions)-1):
    for j in range(i+1,len(solutions)):
      if d.boxes_intersection(solutions[i],solutions[j]) != [] or d.boxes_intersection(solutions[i],solutions[j][:n]+[ [-Bi[1],-Bi[0]] for Bi in  solutions[j][n:m-1]]+[solutions[j][m-1]]) != [] :
        rep.append(j)
      elif   d.boxes_intersection(solutions[i][:2],solutions[j][:2]) != []:
        if i not in unaccepted:
          unaccepted.append(i)
        if j not in unaccepted:
          unaccepted.append(j)     
  #intersect_in2d[i] = intersect_in2d[i]+[ j]
  """elif d.boxes_intersection(solutions[i][:2],solutions[j][:2]) !=[] and (\
      (d.boxes_intersection(solutions[i][n:2*n-2],[[-Bi[1],-Bi[0]] for Bi in solutions[j][n:2*n-2]]) ==[] and \
      d.boxes_intersection(solutions[i][n:2*n-2],[[Bi[0],Bi[1]] for Bi in solutions[j][n:2*n-2]]) ==[] ) \
          or \
       d.boxes_intersection(solutions[i][2:n]+[solutions[i][2*n-2]], solutions[j][2:n]+[solutions[j][2*n-2]]) ==[]) : 
        intersect_in2d[i] = intersect_in2d[i]+[ j]"""

  
  """accepted=[]
  acc_ind=[]
  unaccepted=[]
  unacc_ind=[]
  for i in range(len(solutions)):
    if len(intersect_in2d[i]) ==0 and i not in unacc_ind+acc_ind+rep:
      accepted.append(solutions[i])
      acc_ind.append(i)
      continue
    elif  i not in unacc_ind+acc_ind+rep:
      unaccepted.append(solutions[i])
      unacc_ind.append(i)
    for k in intersect_in2d[i]:
     if k not in unacc_ind:  
       unaccepted.append(solutions[k]) 
       unacc_ind.append(k)   """        
  return [[ solutions[i]  for i in range(len(solutions)) if i not in rep + unaccepted ], [ solutions[i]  for i in range(len(solutions)) if i  in unaccepted ]]     
def Ball_given_2nboxes(system,X, B1,B2, monotonicity_B1=1,monotonicity_B2=1):
  B1_ft=[d.ftconstructor(Bi[0],Bi[1]) for Bi in B1]
  B2_ft=[d.ftconstructor(Bi[0],Bi[1]) for Bi in B2]
  P=[Pi.replace("\n","") for Pi in  open(system,"r").readlines()]
  sol="Empty"
  if d.boxes_intersection(B1_ft, B2_ft) ==[] and monotonicity_B1== monotonicity_B2==1:
    t=estimating_t([[B1_ft], [B2_ft]])
    y_and_r=estimating_yandr([[B1_ft], [B2_ft]])
    intersec_B1B2_in2d=d.boxes_intersection(B1_ft[:2],B2_ft[:2])
    intersec_B1B2_in2d=[ [float(Bi.lower()),float(Bi.upper())]  for Bi in intersec_B1B2_in2d ]
    B_Ball=intersec_B1B2_in2d +y_and_r +[t]
    Ball_node_gen(system,B_Ball,X)
    os.system("ibexsolve   --eps-max=0.1 -s  eq.txt  > output.txt")
    sol=computing_boxes()
  #if d.boxes_intersection(B1_ft, B2_ft) ==[]:
  #  pass
  return sol  
def all_pairs_oflist(L):
  pairs=[]
  for i in range(len(L)-1):
    for j in range(i+1,len(L)):
      pairs.append([L[i],L[j]])
  return pairs    
def checking_assumptions(curve_data): #the input of this function is the output of Ball_solver
  if len(curve_data[0][1]) !=0 :
    return 0
  Ball_sols_ft=[[d.ftconstructor(Bi[0],Bi[1]) for Bi in B] for B in  curve_data[1][0]]+[[d.ftconstructor(Bi[0],Bi[1]) for Bi in B] for B in  curve_data[1][1]]
  alph3=assum_alph3_checker(Ball_sols_ft)
  if alph3==1 :
    return 1
  else:
    return 0
def csv_saver(L,type_L="Ball"):
  dic=[]
  if type_L== "Ball" :
    n=int((len(L[0])+1)/2)
    for j in range(len(L)):
       dic.append({})
       for i in range(n):
         dic[j]["x"+str(i+1)]=L[j][i]
       for i in range(n,2*n-2):
          dic[j]["r"+str(i+3-n)]=L[j][i]
       dic[j]["t"]= L[j][2*n-2]
  return dic     
def dict2csv(dictlist, csvfile):
    """
    Takes a list of dictionaries as inpu and outputs a CSV file.
    """
    f = open(csvfile, 'wb')

    fieldnames = dictlist[0].keys()

    csvwriter = csv.DictWriter(f, delimiter=',', fieldnames=fieldnames)
    csvwriter.writerow(dict((fn, fn) for fn in fieldnames))
    for row in dictlist:
        csvwriter.writerow(row)
    fn.close()          
def assum_alph3_checker(solutions):
  comparing_list=[[]]*len(solutions)
  for i in range(len(solutions)-1):
    for j in range(i+1,len(solutions)):
      if d.boxes_intersection(solutions[i][:2],solutions[j][:2]) !=[]:
        comparing_list[i].append(j)
        comparing_list[j].append(i)
  matching=[len(T) for T in comparing_list]
  if max(matching) <=2:
    return 1
  else:
    return 0

def plotting_3D(boxes,Box,var=[0,1,2]):
  ax = plt.figure().add_subplot(111, projection='3d')
  ax.set_xlim(Box[0][0], Box[0][1])
  ax.set_ylim(Box[1][0], Box[1][1])
  ax.set_zlim(Box[2][0], Box[2][1])
  ax.set_xlabel("x"+str(var[0]+1))
  ax.set_ylabel("x"+str(var[1]+1))
  ax.set_zlabel("x"+str(var[2]+1))
  for box in boxes : 
    V=[[box[j][0] for j in range(3)] ,  [box[j][1] for j in range(3)]]
    #ax.scatter3D(box[0], box[1], box[2])
    points =list(itertools.product(*box))
    faces=[[points[0],points[2],points[6],points[4]],
    [points[0],points[2],points[3],points[1]],
 [points[0],points[1],points[5],points[4]], 
 [points[2],points[3],points[7],points[6]], 
 [points[1],points[3],points[7],points[5]]]
    ax.add_collection3d(Poly3DCollection(faces, 
 facecolors='green', linewidths=1,edgecolors='green', alpha=.25))
  plt.show()  
def check_for_pre(B,eps=1e-6):
  for Bi in B:
    
    if d.width([Bi])< eps:
      Bi[0]=Bi[0]-eps
      Bi[1]=Bi[1]+eps
  return B    


def robot(a,b,c,d,e):
  
  P1= "(x1 -" +str(a) +  "*cos(x3))^2 + (x2 -"+ str(a)+ " *sin(x3))^2 - " +str(b**2) 
  P2="(x1 -(" +str(e)+ ") - (" +str(-d) +  ")*cos(x4))^2 + (x2 - ("+ str(d)+ " )*sin(x4))^2 - " +str(c**2) 
  jac=sp.Matrix([P1,P2]).jacobian(sp.Matrix(["x1","x2","x3","x4"]))
  P_ser=sp.simplify( sp.expand(jac[:,2:].det()))
  P_par=sp.simplify( sp.expand(jac[:,:2].det()))
  f=open("serial.txt","w")
  f.write(P1+"\n")
  f.write(P2+"\n")
  f.write(str(P_ser).replace("**","^"))
  f.close()


  f=open("parallel.txt","w")
  f.write(P1+"\n")
  f.write(P2+"\n")
  f.write(str(P_par).replace("**","^"))
  f.close()



def enclosing_singularities(system,boxes,B,X,eps_max=0.1,eps_min=0.01, threshold=0.01,info=False): #there still computing Ball  On the case where tow monotonic boxes intersect
  #Defining variables######################################
  ##########################################################
  t1=time.time()
  
  tree_size_ball=0
  n=len(B);
  P=[Pi.replace("\n","") for Pi in  open(system,"r").readlines()]
  is_poly= False not in [ parse_expr(Pi.replace("^","**")).is_polynomial()  for Pi in P ]
  ##### Computing Ball for nodes 
  Ball_P=[]
  for Pi in P:
    H=SDP_str(Pi,X)
    Ball_P.append(H[0])
    Ball_P.append(H[1]) 
  last_eq=""
  for i in range(3,n):
        last_eq += "r"+str(i)+"^2+"
  last_eq += "r" +str(n)+"^2 -1=0; \n"
  Ball_P.append(last_eq) 
  certified_boxes, uncertified_boxes,tree_size_curve= boxes
  time_ball=[] 
  time_combinatorics=[] 
  ## Classifying boxes into x1-monotonous, x2-but-not-x1 monotonous, non-monotonous and non-smooth boxes 
  #####################################################################################################
  com_start=time.time()
  classes= boxes_classifier(system,boxes,X)
  
  Lcc=[[ci] for ci in  classes[2]]  #Lcc is tended to be the set of connected components
  
  mon_boxes=classes[0]+classes[1]
  mon_mid=[[0.5*(Bij[1]+Bij[0]) for Bij in Bi[:2] ] for Bi in mon_boxes ]
  mon_rad=[ 2*d.width(Bi[:2]) for Bi in  mon_boxes ]
  tree = spatial.KDTree(mon_mid)
  Kdtree_boxes=[tree.query_ball_point(m,r=(math.sqrt(2))*r ) for m,r in zip(mon_mid,mon_rad)] 
  L=[]
  for i in range(len(Kdtree_boxes)):
    if len(Kdtree_boxes[i]) >3:
      for j in Kdtree_boxes[i]:
        if j not in L:
          L.append(j)


  L0=[i for i in L if i < len(classes[0])]
  L1=[i for i in L if i >= len(classes[0])]

  graph0=np.zeros((len(L0), len(L0)))
  for i in L0:
    
    for j in Kdtree_boxes[i]:
      if   j in L0 and  d.boxes_intersection( mon_boxes[i], mon_boxes[j]) !=[] :
        graph0[L0.index(i),L0.index(j)]=1

        
  graph0 = csr_matrix(graph0)
  n_components0, labels0 = connected_components(csgraph=graph0, directed=False, return_labels=True)
  comp_class_0=empty_lists = [ [] for i in range(n_components0) ]
  for i in range(len(L0)):
    comp_class_0[labels0[i]].append(classes[0][L0[i]]) 
  Lcc=Lcc+ comp_class_0 


  if len(classes[1])>0:
   graph1=np.zeros((len(L1), len(L1)))
   for i in L1:
    
    for j in Kdtree_boxes[i]:
      if   j in L1 and  d.boxes_intersection( mon_boxes[i], mon_boxes[j]) !=[] :
        graph1[L1.index(i),L1.index(j)]=1     
   graph1 = csr_matrix(graph1)
   n_components1, labels1 = connected_components(csgraph=graph1, directed=False, return_labels=True)
   comp_class_1=empty_lists = [ [] for i in range(n_components1) ]
   for i in range(len(L1)):
    comp_class_1[labels1[i]].append(mon_boxes[L1[i]])


   Lcc=Lcc+ comp_class_1 
  
  

  t2=time.time() 
  #############################################################################
  #Computing intersting_boxes= {box in boxes[0] with box has at least three neighbors  in 2D}#
  #intersting_boxes is equivalent to L in Marc's Algorithm
  #############################################################################

  ###############################################################################################################################
  #For every two distinct connected components of Lcc, we check if they intersect in 2D, if so we compute the ball system ans solve it  #
  ################################################################################################################################ 

  Solutions=[[],[]] # list of two lists, certified and uncertified solutions
  for c1_index in range(len(Lcc)-1):
    for c2_index in range(c1_index+1,len(Lcc)):
      ball_start=time.time()
      
      if d.components_intersection(Lcc[c1_index],Lcc[c2_index] ) == True:
        
        ####################################################################
        ##Now we discuss several cases depending on how far  0 from t is#####
        #################################################################### 
        uni=[]
        for box in Lcc[c1_index]+Lcc[c2_index] :
            uni = d.box_union(uni,box)
        t=estimating_t([Lcc[c1_index],Lcc[c2_index]]); t1p = d.ftconstructor(t[0],t[1]); t=[float(t1p.lower()),float(t1p.upper())];
        # t contains 0  
        if  t[0] <= 0<= t[1]   :
          t1po=[t[0],min(threshold,t[1])]
          B_Ball=uni[:] +[[-1-eps_min,+1+eps_min]]*(n-2) +[t1po]
      
          if n==3:
            B_Ball[n]=[0.999,1.001]
          B_Ball= check_for_pre(B_Ball) 
          if is_poly:
            Ball_generating_system(Ball_P,B_Ball,X,eps_min)
            os.system("ibexsolve   --eps-max="+ str(eps_max)+"  --eps-min="+ str(eps_min) + " -s  eq.txt  > output.txt")
            res=computing_boxes()
          else:
            res=cusp_Ball_solver(P,B_Ball,X,eps_min=eps_min,eps_max=eps_max)
          Sol=res[:2]
          tree_size_ball += res[len(res)-1]
          if Sol[0] != "Empty" and Sol !=[[],[]]:
            Solutions[0]=Solutions[0]+Sol[0] 
            Solutions[1]=Solutions[1]+Sol[1] 
          elif Sol == [[],[]]   and d.width(uni[:2]) > eps_min:
                
                new_B=B_Ball[:n]
                res=enclosing_curve(system,new_B,X,eps_max=0.1*eps_max, eps_min=eps_min)
                resul=enclosing_singularities(system,res,new_B,X,eps_max=eps_max,eps_min=eps_min)
                Solutions[0] += resul[0]+resul[1] 
                Solutions[1] += resul[2]
                boxes[1] += res[1]
                tree_size_ball += resul[3]
                tree_size_curve +=resul[4]
          elif   Sol == [[],[]]   and d.width(uni) <= eps_min:   
                Solutions[1].append(B_Ball)
        
          if t[1] > threshold:
           
           t2pos=[threshold,t[1]]
           r=[ [float(ri[0]),float(ri[1])] for ri in  estimating_yandr([Lcc[c1_index],Lcc[c2_index]])]
           B_Ball=uni[:2] +r[:n-2]+[[-1-eps_min,+1+eps_min]]*(n-2) +[t2pos] 
           if n==3:
            B_Ball[n]=[0.999,1.001]
         
           Ball_generating_system(Ball_P,B_Ball,X,eps_min)
           os.system("ibexsolve   --eps-max="+ str(eps_max)+"  --eps-min="+ str(eps_min) + " -s  eq.txt  > output.txt")
           res=computing_boxes()
           Sol=res[:2]
           tree_size_ball += res[len(res)-1]
           if Sol[0] != "Empty" and Sol !=[[],[]]:
            Solutions[0]=Solutions[0]+Sol[0] 
            Solutions[1]=Solutions[1]+Sol[1]
           elif Sol == [[],[]]   and d.width(uni[:2]) > eps_min:
                new_B=uni
                res=enclosing_curve(system,new_B,X,eps_max=0.1*eps_max, eps_min=eps_min)
                resul=enclosing_singularities(system,res,new_B,X,eps_max=eps_max,eps_min=eps_min)
                Solutions[0] += resul[0]+resul[1] 
                Solutions[1] += resul[2]
                boxes[1] += res[1]
                tree_size_ball += resul[3]
                tree_size_curve +=resul[4]
           elif   Sol == [[],[]]   and d.width(uni[:2]) <= eps_min:   
                Solutions[1].append(B_Ball)         
          
        # t does not contain 0   but close enough       
        elif t[1]< threshold :
          
          r=[ [float(ri[0]),float(ri[1])] for ri in  estimating_yandr([Lcc[c1_index],Lcc[c2_index]])]
          B_Ball=uni[:2] +r[:] +[t] 
          if n==3:
            B_Ball[n]=[0.999,1.001]
          B_Ball= check_for_pre(B_Ball)
          if is_poly:
            Ball_generating_system(Ball_P,B_Ball,X,eps_min)
            os.system("ibexsolve   --eps-max="+ str(eps_max)+"  --eps-min="+ str(eps_min) + " -s  eq.txt  > output.txt")
            res=computing_boxes()
          else:
            res=cusp_Ball_solver(P,B_Ball,X,eps_min=eps_min,eps_max=eps_max)
          Sol=res[:2]
          tree_size_ball += res[len(res)-1]
          if Sol[0] != "Empty" and Sol !=[[],[]]:
            Solutions[0]=Solutions[0]+Sol[0] 
            Solutions[1]=Solutions[1]+Sol[1] 
          elif Sol == [[],[]]   and d.width(uni[:2]) > eps_min:
                
                new_B=uni
                res=enclosing_curve(system,new_B,X,eps_max=0.1*eps_max, eps_min=eps_min)
                resul=enclosing_singularities(system,res,new_B,X,eps_max=eps_max,eps_min=eps_min)
                Solutions[0] += resul[0]+resul[1] 
                Solutions[1] += resul[2]
                boxes[1] += res[1]
                tree_size_ball += resul[3]
                tree_size_curve +=resul[4]
          elif   Sol == [[],[]]   and d.width(uni[:2]) <= eps_min:   
                Solutions[1].append(B_Ball)  
         # t is a away from 0
        elif t[0]> threshold or  (t[0] < threshold and t[1] > threshold ):
           r=[ [float(ri[0]),float(ri[1])] for ri in  estimating_yandr([Lcc[c1_index],Lcc[c2_index]])]
           B_Ball=uni[:2] +r[:] +[t] 
           if n==3:
            B_Ball[n]=[0.999,1.001]
           B_Ball= check_for_pre(B_Ball)
           Ball_generating_system(Ball_P,B_Ball,X,eps_min)
           os.system("ibexsolve   --eps-max="+ str(eps_max)+"  --eps-min="+ str(eps_min) + " -s  eq.txt  > output.txt")
           ibex_output=computing_boxes()
           Sol=ibex_output[:2]
           tree_size_ball += ibex_output[len(ibex_output)-1]
           if Sol[0] != "Empty" and Sol != [[],[]] :
            Solutions[0]=Solutions[0]+Sol[0] 
            Solutions[1]=Solutions[1]+Sol[1]
           elif Sol == [[],[]]   and d.width(uni[:2]) > eps_min:
                new_B=uni
                res=enclosing_curve(system,new_B,X,eps_max=0.1*eps_max, eps_min=eps_min)
                resul=enclosing_singularities(system,res,new_B,X,eps_max=eps_max,eps_min=eps_min)
                Solutions[0] += resul[0]+resul[1] 
                Solutions[1] += resul[2]
                boxes[1] += res[1]
                tree_size_ball += resul[3]
                tree_size_curve +=resul[4]
           elif   Sol == [[],[]]   and d.width(uni[:2]) <= eps_min:   
                Solutions[1].append(B_Ball)        
         
      ball_end=time.time()
      time_ball.append(ball_end-ball_start)   

  ####################################################################################################
  #Now we have all ball solutions. We check whether all of them satisfy Conditions alpha2 and alpha3 
  ####################################################################################################
  com_start=time.time()
  nodes=[]
  Sols_with_0in_t=[]
  checker=projection_checker(Solutions[0])
  Solutions[1]= Solutions[1] +checker[1]
  Solutions[0]=[Bi for Bi in checker[0] if Bi[2*n-2][1] >= 0   ] 
  for solution in Solutions[0] :
    if 0 >= solution[2*n-2][0] and 0 <= solution[2*n-2][1]:
      Sols_with_0in_t.append(solution)
    else:
      nodes.append(solution) 
  com_end=time.time()
  time_combinatorics.append(com_end-com_start)  


  if info ==True:
   print("tree_size_curve: ",tree_size_curve )
   print("tree_size_ball: ",tree_size_ball )


   print("Timing of solving Ball system  : ", sum(time_ball))
   print("Timing of combinatorial steps: ", t2-t1)
  return [nodes,Sols_with_0in_t, Solutions[1], tree_size_ball, tree_size_curve]    
"""
#robot(9,3,5.1,7.1,5)   robot 1 
#robot(8,5,8,5,15) robot 2     
robot(8,5,8,5,5) # robot 3   
robot(5,5,5,5,5) #robot 4 
robot(10,8,7,8,2) # robot 5 
robot(10,8,7,8,6) # robot 6 

robot(10,6,2,4,8)"""


