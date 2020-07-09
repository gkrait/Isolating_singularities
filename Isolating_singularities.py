
import matplotlib.pyplot as plt
import os
import pickle 
import interval_arithmetic as d
from pprint import pprint
from sympy.parsing.sympy_parser import parse_expr
import sympy as sp 
import os 
from cusp import cusp_Ball_solver, evaluation_exp
import matplotlib.patches as mpatches


def ploting_boxes(boxes,uncer_boxes, var=[0,1], B=[[-20,20],[-20,20]],a=1,b=10,nodes=[], cusps=[],uncer_Solutions=[],color="g"):
   fig, ax = plt.subplots()
   #plt.grid(True)
   ax.set_xlim(B[0][0], B[0][1])
   ax.set_ylim(B[1][0], B[1][1])
   ax.set_xlabel('x'+str(var[0]+1))
   ax.set_ylabel('x'+str(var[1]+1))
   try:
    ax.set_title(open("system.txt","r").read())
   except:
    pass
   

   c=0
   green_patch = mpatches.Patch(color='green', label='smooth part')
   red_patch = mpatches.Patch(color='red', label='unknown part')
   node_patch = mpatches.Patch(color='black', label='Certified nodes',fill=None)
   cusp_patch = mpatches.Patch(color='blue', label='cusps or small nodes',fill=None)
   plt.legend(handles=[green_patch,red_patch,node_patch,cusp_patch])
   for box in boxes:
     rectangle= plt.Rectangle((box[var[0]][0],box[var[1]][0]) , \
      a*(box[var[0]][1]-box[var[0]][0]),a*(box[var[1]][1]-box[var[1]][0]), fc=color)
     plt.gca().add_patch(rectangle)
   for box in uncer_boxes:
      rectangle= plt.Rectangle((box[var[0]][0],box[var[1]][0]) , \
      a*(box[var[0]][1]-box[var[0]][0]),a*(box[var[1]][1]-box[var[1]][0]), fc='r')
      plt.gca().add_patch(rectangle)
   for box in nodes:
     rectangle= plt.Rectangle((box[0][0]-0.05,box[1][0]-0.05) ,\
      (0.09),(0.09), fc='y',fill=None)
     plt.gca().add_patch(rectangle) 
   for box in cusps:
     rectangle= plt.Rectangle((box[0][0]-0.05,box[1][0]-0.05) ,\
      (0.09),(0.09),color="blue", fc='y',fill=None)
     plt.gca().add_patch(rectangle)
   for box in uncer_Solutions:
     rectangle= plt.Rectangle((box[0][0]-0.05,box[1][0]-0.05) ,\
      (0.09),(0.09),color="red", fc='y',fill=None)
     plt.gca().add_patch(rectangle)  
   
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
        V += SDP_str(Pi,X)[0]
        V += SDP_str(Pi,X)[1]
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
			Ball_cusp_gen(equations,B_Ball,X)
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
    P_pluse=P[:]
    P_minus=P[:]
    for i in range(2,n):
        P_pluse=P_pluse.replace("x"+str(i+1),"(x"+str(i+1) + "+ r"+str(i+1) +"*sqrt(t))")
        P_minus=P_minus.replace("x"+str(i+1),"(x"+str(i+1) + "- r"+str(i+1) +"*sqrt(t))")
    SP= "0.5*(" + P_pluse + "+" +P_minus+")=0; \n"
    DP= "0.5*(" + P_pluse + "- (" +P_minus+") )/(sqrt(t))=0; \n"
    return [SP,DP]
def Ball_generating_system(P,B_Ball,X):
    n=len(X)
    V=""" Variables \n """
    for i in range(n):
        V += "x" +str(i+1) + " in " + str(B_Ball[i]) +" ; \n"
    for i in range(n,2*n-2):
        V += "r" +str(i-n+3) + " in " + str(B_Ball[i]) +" ; \n" 
    V += "t" + " in " + str(B_Ball[2*n-2]) +" ; \n"       
    V +="Constraints \n"   
    for Pi in P:
        V += SDP_str(Pi,X)[0]
        V += SDP_str(Pi,X)[1]

    last_eq=""
    for i in range(3,n):
        last_eq += "r"+str(i)+"^2+"
    last_eq += "r" +str(n)+"^2 -1=0;"    

    V += last_eq +"\n"

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
def estimating_t(components,upper_bound=19000.8):  #it works only if len(components)==2
  t1=upper_bound
  t2=0
  for box1 in components[0]:
    for box2 in components[1]:

        a=d.distance(box1[2:],box2[2:])

        if t1 > a[0]:
          t1=a[0]
        if t2<a[1]:
           t2=a[1]  
  t=d.ftconstructor(t1,t2)
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
                if d.boxes_intersection(ftboxes[i],components[j][k]) !=[]:
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

                if d.boxes_intersection(ftboxes[i][:2],components[j][k][:2]) !=[] and \
                d.boxes_intersection(ftboxes[i],components[j][k]) ==[]:
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
                        if d.boxes_intersection(boxi,boxj)==[] and \
                       d.boxes_intersection(boxi[:2],boxj[:2]) != []  :
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




    return components      
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
               (Mariam != [] and enclosing_curve(f,Mariam,X,eps*0.1) ==[[],[]] ): #needs to work more
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
  if "infeasible" in open("output.txt","r").read():
    return "Empty"
  content=open("output.txt","r").readlines()
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
  return [cer,uncer]        
def enclosing_curve(system,B,X,eps=0.1): 
  L=[B]
  certified_boxes=[]
  uncertified_boxes=[]
  while len(L) !=0: 
    system_generator(system,L[0],X)
    content=open("output.txt","r").readlines()
    os.system("ibexsolve   --eps-max=" + str(eps) +" -s  eq.txt  > output.txt")
    ibex_output=computing_boxes()
    #ibex_output=solving_with_ibex(eps)
    if ibex_output ==[[],[]] and max([Bi[1]-Bi[0] for Bi in L[0]  ]) < eps*0.01 :  
      uncertified_boxes.append(L[0])
      L.remove(L[0]);
    elif ibex_output ==[[],[]] :
      children=plane_subdivision(L[0])
      L.remove(L[0]);
      L += children
    elif ibex_output== "Empty":
      L.remove(L[0])
    else:
      if len(ibex_output[0]) !=0:
       certified_boxes += ibex_output[0]
      if len(ibex_output[1])!=0: 
       uncertified_boxes += ibex_output[0]
      L.remove(L[0]) 
    
  return [certified_boxes,uncertified_boxes]                        
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
  fil=open("evaluation_file.py","w")
  fil.write("import flint as ft \n")
  fil.write("import sympy as sp \n")
  fil.write("import interval_arithmetic as d \n")
  fil.write("def eval_func():\n  pass \n")
  fil.write("boxes="+str(boxes)+"\n")
  fil.write("ftboxes=[ [d.ftconstructor(Bi[0],Bi[1]) for Bi in B ] for B in boxes ] \n"    )
  fil.write("n=len(boxes[0])\n")
  fil.write("m=len(boxes)\n")
  fil.write("m1=[]\n")
  fil.write("m2=[]\n")
  minor1_str=str(minor1)
  minor2_str=str(minor2)
  for i in range(n):
    minor1_str= minor1_str.replace("x"+str(i+1),"B["+str(i)+"]" )
    minor2_str= minor2_str.replace("x"+str(i+1),"B["+str(i)+"]" )
  for func in functions:
    minor1_str=minor1_str.replace(func,"ft.arb."+func)
    minor2_str=minor2_str.replace(func,"ft.arb."+func)
  fil.write("for B in ftboxes: \n")
  fil.write("   m1.append("+ minor1_str + ") \n")
  fil.write("   m2.append("+ minor2_str + ") \n")  
  fil.write("innrer_loops=[i for i in range(m) if 0 in m1[i] and 0 in m2[i] ]\n")
  fil.write("print(innrer_loops)\n")
  fil.close()
  t=os.popen("python3 evaluation_file.py ").read()
  return t
def boxes_classifier(system,boxes,X,special_function=[]):
  if len(boxes[0])==0:
    return [[],[],boxes[1]]
  certified_boxes ,uncer_boxes =boxes
  L=eval_file_gen(system,certified_boxes,X)
  if L==[]:
    return [[],[],uncer_boxes]
  it=0
  L=L.replace('[','')
  L=L.replace(']','')
  L=L.replace('\n','')
  L=L.split(",")
  if L !=[""]:
    L=[int(li) for li in L]
    return   [ [certified_boxes[i]  for i in range(len(certified_boxes)) if i not in L] ,\
    [certified_boxes[i]  for i in L ], \
    uncer_boxes ]
  else:
      return  [ [certified_boxes[i]  for i in range(len(certified_boxes)) if i not in L] ,[], uncer_boxes ] #can be enhanced
def projection_checker(solutions):
  if len(solutions)==0:
    return [[],[]]
  m=len(solutions[0])
  n=int((m+1)/2)
  intersect_in2d=[[]]*len(solutions)
  for i in range(len(solutions)-1):
    for j in range(i+1,len(solutions)):
      """pprint(sp.Matrix(solutions[i]))
      pprint(sp.Matrix(solutions[j]));"""
      if d.boxes_intersection(solutions[i][:2],solutions[j][:2]) !=[] and \
      d.boxes_intersection(solutions[i][n:2*n-2],[[-Bi[1],-Bi[0]] for Bi in solutions[j][n:2*n-2]]) ==[] and \
       d.boxes_intersection(solutions[i][2:n], solutions[j][2:n]) ==[] : 
        intersect_in2d[i] = intersect_in2d[i]+[ j]
      """print(d.boxes_intersection(solutions[i][:2],solutions[j][:2]) !=[] and \
      d.boxes_intersection(solutions[i][n:2*n-2],[[-Bi[1],-Bi[0]] for Bi in solutions[j][n:2*n-2]]) ==[] and \
       d.boxes_intersection(solutions[i][2:n], solutions[j][2:n]) ==[] );input() """
  accepted=[]
  acc_ind=[]
  unaccepted=[]
  unacc_ind=[]
  for i in range(len(solutions)):

    if len(intersect_in2d[i]) ==0 and i not in unacc_ind:
      accepted.append(solutions[i])

      continue
    unaccepted.append(solutions[i])
    unacc_ind.append(i)
    for k in intersect_in2d[i]:
     if k not in unacc_ind:
       
       unaccepted.append(solutions[k]) 
       unacc_ind.append(k)  
  return [accepted, unaccepted] 		
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
def enclosing_singularities(system,boxes,B,X,eps=0.1): #there still computing Ball  On the case where tow monotonic boxes intersect
  n=len(B);
  P=[Pi.replace("\n","") for Pi in  open(system,"r").readlines()]
  certified_boxes, uncertified_boxes= boxes
  classes= boxes_classifier(system,boxes,X,special_function=[])
  cer_Solutions=[]
  uncer_Solutions=[]
  #############################################################################
  #Solving Ball for B1 and B2 in R^n such that C is monotonic in B1 and B2
  #######################################################################
  monotonic_pairs=intersect_in_2D(classes[0],classes[0])
  monotonic_componants=[ Bi[0] for Bi in  monotonic_pairs ] +[ Bi[1] for Bi in  monotonic_pairs ]
  """checked_boxes=[]
  for potential_node in monotonic_pairs:
    if potential_node[0] in checked_boxes and potential_node[1] in checked_boxes:
      continue"""
  plane_components= planner_connected_compnants(monotonic_componants)
  H=[]
  for plane_component in plane_components:  
        x1=float(min([ai[0][0] for ai in plane_component]))
        x2=float(max([ai[0][1] for ai in plane_component]))
        y1=float(min([ai[1][0] for ai in plane_component]))
        y2=float(max([ai[1][1] for ai in plane_component]))
        components=connected_compnants(plane_component)
        pairs_of_branches=all_pairs_oflist(components)

        for pair_branches in  pairs_of_branches:
          all_boxes=pair_branches[0]+pair_branches[1]
          uni=[]
          for box in all_boxes:
            uni = d.box_union(uni,box)

          t=estimating_t(pair_branches)
          r=[ [float(ri[0]),float(ri[1])] for ri in  estimating_yandr(pair_branches)]
          B_Ball=uni[:2] +r +[t] 
          H.append(B_Ball)    
          Ball_generating_system(P,B_Ball,X)
          os.system("ibexsolve   --eps-max="+ str(eps)+" -s  eq.txt  > output.txt")
          Solutions=computing_boxes()
          if Solutions != "Empty":
            cer_Solutions += Solutions[0]
            uncer_Solutions += Solutions[1]  
          

            
    #There still the case B1B2[0],B1B2[1] are not disjoint 
  ########################################################################################################
  #Solving Ball for potential_cusp, a box in  R^n such that C is not monotonic 
  ########################################################################################################
  checked_boxes=[]
  for potential_cusp in classes[1]:
    ###finding cusps (or small loops) in potential_cusp####
    plane_intersecting_boxes= intersect_in_2D([potential_cusp],classes[0]+classes[1]+classes[2],monotonicity=0)
    
    intersecting_boxes= [pair_i[1] for pair_i in plane_intersecting_boxes \
     if  d.boxes_intersection(pair_i[1], potential_cusp)!=[] ] 
    
    ##########
    uni= potential_cusp[:]
    checked_boxes.append(potential_cusp)
    for box in intersecting_boxes:
     if box in checked_boxes:
      continue
     uni = d.box_union(uni,box)
     checked_boxes.append(box)
    max_q1q2=d.distance(uni[2:],uni[2:])
    max_q1q2=d.ftconstructor(max_q1q2[0],max_q1q2[1])
    t=d.power_interval(max_q1q2,2)/4
    t=[float(t.lower()),float(t.upper())]
    B_Ball=uni +[[-1.01,1.01]]*(n-2)+[t]
    sol=cusp_Ball_solver(P,B_Ball,X)
    if sol != "Empty":
         cer_Solutions += sol[0]
         uncer_Solutions += sol[1]
    ####finding nodes that have the same projection with potential_cusp
    non_intersecting_boxes= [pair_i[1] for pair_i in plane_intersecting_boxes \
     if  d.boxes_intersection(pair_i[1], potential_cusp)==[] ] 
    
    for aligned in non_intersecting_boxes:
      if aligned  in checked_boxes:
        continue
      boxes_intersect_aligned=[B  for B in non_intersecting_boxes if d.boxes_intersection(aligned,B) != []  ]
      uni=aligned[:]
      for boxi  in boxes_intersect_aligned:
          if boxi in checked_boxes:
            continue
          uni=d.box_union(uni,boxi)
          checked_boxes.append(boxi)
      t=estimating_t([[potential_cusp],[uni]])
      r=[ [float(ri[0]),float(ri[1])] for ri in  estimating_yandr([[potential_cusp],[uni]])]
      B_Ball=potential_cusp[:2]+r +[t]      
      Ball_generating_system(P,B_Ball,X)
      os.system("ibexsolve   --eps-max="+ str(eps)+" -s  eq.txt  > output.txt")
      Solutions=computing_boxes()
      if Solutions != "Empty":
          cer_Solutions += Solutions[0]
          uncer_Solutions += Solutions[1]               
  nodes=[]
  cups_or_smallnodes=[]
  checker=projection_checker(cer_Solutions)
  uncer_Solutions= uncer_Solutions +checker[1]
  cer_Solutions=[Bi for Bi in checker[0] if Bi[2*n-2][1] >= 0   ] 
  for solution in cer_Solutions :
    if 0 >= solution[2*n-2][0] and 0 <= solution[2*n-2][1]:
      cups_or_smallnodes.append(solution)
    else:
      nodes.append(solution) 
  return [nodes,cups_or_smallnodes, uncer_Solutions ]    
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
def checking_assumptions(curve_data): #the input of this function is the output of Ball_solver
	if len(curve_data[0][1]) !=0 :
		return 0
	Ball_sols_ft=[[d.ftconstructor(Bi[0],Bi[1]) for Bi in B] for B in  curve_data[1][0]]+[[d.ftconstructor(Bi[0],Bi[1]) for Bi in B] for B in  curve_data[1][1]]
	alph3=assum_alph3_checker(Ball_sols_ft)
	if alph3==1 :
		return 1
	else:
		return 0
System="system.txt" 
#Box=[[-5, 15], [-15, 15],[-3.14,3.14],[-3.14,3.14]]
Box=[[-2.01,3.03],[-1.03,5.03],[-1.6,1.6]]
X=[sp.Symbol("x"+str(i)) for i in range(1,4)]
boxes =enclosing_curve(System,Box,X,eps=0.1)

nodes, cups_or_smallnodes,uncer_Solutions=enclosing_singularities(System,boxes,Box,X)
#plotting the singularities
ploting_boxes(boxes[0],boxes[1] , nodes = nodes, cusps= cups_or_smallnodes,uncer_Solutions=uncer_Solutions )
##################################
#Declaring parameters #######
##################################
"""System="system.txt" 
Box=[[-5,15],[-15,15],[-3.14,3.14],[-3.14,3.14]]
X=[sp.Symbol("x"+str(i)) for i in range(1,5)]

##################################
#Applying the function #######
##################################
boxes =enclosing_curve(System,Box,X)

"""