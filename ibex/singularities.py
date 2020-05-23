import matplotlib.pyplot as plt
import os
import pickle 
from pprint import pprint
import sympy
import computing_boxes as cb
import draft as d
import flint as ft
import math
from sympy.parsing.sympy_parser import parse_expr 


def SDP_str(P):
    n=len(P)+1
    P_pluse=P[:]
    P_minus=P[:]
    for i in range(2,n):
        P_pluse=P_pluse.replace("x"+str(i+1),"(x"+str(i+1) + "+ r"+str(i+1) +"*sqrt(t))")
        P_minus=P_minus.replace("x"+str(i+1),"(x"+str(i+1) + "- r"+str(i+1) +"*sqrt(t))")
    SP= "0.5*(" + P_pluse + "+" +P_minus+")=0; \n"
    DP= "0.5*(" + P_pluse + "- (" +P_minus+") )/(sqrt(t))=0; \n"
    return [SP,DP]
def generating_system(P,B_Ball):
    n=len(P)+1
    V=""" Variables \n """
    for i in range(n):
        V += "x" +str(i+1) + " in " + str(B_Ball[i]) +" ; \n"
    for i in range(n,2*n-2):
        V += "r" +str(i-n+3) + " in " + str(B_Ball[i]) +" ; \n" 
    V += "t" + " in " + str(B_Ball[2*n-2]) +" ; \n"       
    V +="Constraints \n"   
    
    for Pi in P:
        V += SDP_str(Pi)[0]
        V += SDP_str(Pi)[1]

    last_eq=""
    for i in range(3,n):
        last_eq += "r"+str(i)+"^2+"
    last_eq += "r" +str(n)+"^2 -1=0;"    

    V += last_eq +"\n"

    f= open("eq.txt","w+")
    f.write(V) 
    
    
    f.write("end")
    f.close()
def intersting_boxes(f,b):
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
def finding_nodes(P,B):
    generating_system(P,B)
    os.system("ibexsolve   --eps-max=0.1 -s  eq.txt  > output.txt")
    g=open('output.txt','r')
    result=g.readlines()

    T=cb.computing_boxes(result)

    return T    
def estimating_t1(components,upper_bound=200):  #it works only if len(components)
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

def estimating_t(components,upper_bound=19.8):  #it works only if len(components)
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
    ftboxes=boxes[:]
    #ftboxes=[ [d.ftconstructor(boxi[0],boxi[1]) for boxi in box ] for box in boxes ]
    components=[[ftboxes[0]]]
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
def estimating_r(components,upper_bound=1000):
  r_bounds=[[upper_bound,0]]*(len(components[0][0])-2)
  r_list=[]
  y_list=[]
  for box1 in components[0]:
    for box2 in components[1]:
        y_list.append([0.5*(q1+q2) for q1,q2 in zip(box1[2:],box2[2:])])
        norm_q1q2=d.distance(box1[2:],box2[2:])
        q1q2=[box1[i]-box2[i] for i in range(2,len(box1)) ]
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
  
def detecting_nodes(f,B):
    boxes=intersting_boxes(f,B)[0]

    boxes=[ [d.ftconstructor(boxi[0],boxi[1]) for boxi in box ] for box in boxes ]
    nodes_lifting=[]
    used=[]

    for i in range(len(boxes)):
        for j in range(i+1,len(boxes)):
            if d.boxes_intersection(boxes[i],boxes[j]) ==[] and \
             d.boxes_intersection(boxes[i][:2],boxes[j][:2]) :
                if i not in used:
                    used.append(i)
                    nodes_lifting.append(boxes[i])
                if j not in used:
                    used.append(j)
                    nodes_lifting.append(boxes[j])
    components= planner_connected_compnants(nodes_lifting)

    return components           

def solving_fornodes(equations,f,B):
    plane_components=detecting_nodes(f,B)
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

        generating_system(P,B_Ball)
        solutionsi=finding_nodes(P,B_Ball)

        Ball_solutions +=solutionsi
        #pprint(B_Ball[:])
        #print(solutionsi)

        


    return Ball_solutions




        


"""
equations="equations.txt"
f="boxes_second_branch"
B=[[-20,20],[-20,20],[-3.14,3.14],[-3.14,3.14]]
T=solving_fornodes(equations,f,B)
print(len(T),T)
"""
#T=[d.ft_normal(Ti) for Ti in T ]





#cb.ploting_boxes(T,[],a=10)


"""
x1=sympy.Symbol("x1")
x2=sympy.Symbol("x2")
q1=sympy.Symbol("q1")
q2=sympy.Symbol("q2")
r3=sympy.Symbol("r3")
r4=sympy.Symbol("r4")
t=sympy.Symbol("t")
X=[x1,x2,q1,q2,r3,r4,t]
S=[]
for Pi in P:
    S.append(SDP_str(Pi)[0].replace('=0; \n',""))
    S.append(SDP_str(Pi)[1].replace('=0; \n',""))
    S.append("r3^2+r4^2-1")
S1=[]
for Pi in S:
        Pi=Pi.replace("^","**")
        S1.append(parse_expr(Pi))    

jac=sympy.Matrix(S1).jacobian(sympy.Matrix(X)) 
eval_jac=jac.subs([(x1,0.9207659841484748),(x2,6.176745907220626),(q1,1.422816280529649),\
    (q2,1.60413583605439),(r3,1),(r4,0)])
input()"""
















"""
sorted=boxes_sort(T)

components=connected_compnants(sorted)
print(len(components))
t=estimating_t(components)
r=estimating_r(components)




S=finding_nodes(P,b[0],b[1],[t[0],t[1]])
"""







"""
for box in components[0]:
    rectangle= plt.Rectangle((float(box[2].lower()),float(box[3].lower()) ), \
        float(box[2].upper())-float(box[2].lower()),float(box[3].upper())-float(box[3].lower()), fc='r')
    plt.gca().add_patch(rectangle)
for box in components[1]:
    rectangle= plt.Rectangle((float(box[2].lower()),float(box[3].lower()) ), \
        float(box[2].upper())-float(box[2].lower()),float(box[3].upper())-float(box[3].lower()), fc='g')
    plt.gca().add_patch(rectangle)        
for box in components[2]:
    rectangle= plt.Rectangle((float(box[2].lower()),float(box[3].lower()) ), \
        float(box[2].upper())-float(box[2].lower()),float(box[3].upper())-float(box[3].lower()), fc='b')
    plt.gca().add_patch(rectangle)  
"""





"""for box in T:
    rectangle= plt.Rectangle((float(box[2].lower()),float(box[3].lower()) ), \
        float(box[2].upper())-float(box[2].lower()),float(box[3].upper())-float(box[3].lower()), fc='b')
    plt.gca().add_patch(rectangle)   
plt.show()"""




"""
pickle_in=open("boxes_second_branch","rb")
branch1=pickle.load(pickle_in)
pickle_in.close()



pprint(sympy.Matrix(intersting_boxes))


     

#print(finding_nodes(P))
"""




"""







X=[]
Y=[]
for box in  T:
    X.append(box[4][0])
    Y.append(box[5][0])

plt.plot(X,Y,"ro")
plt.show()    



"""





            

    


      