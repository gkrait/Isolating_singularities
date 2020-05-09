import matplotlib.pyplot as plt
import os
import pickle 
from pprint import pprint
import sympy as sp
import computing_boxes as cb
import draft as d
import flint as ft
import math
#import singularities as s
from sympy.parsing.sympy_parser import parse_expr

def inter_intersect(I1,I2):
    T=[]
    for i in range(len(I1)):
      T.append([max(I1[i][0],I2[i][0]),min(I1[i][1],I2[i][1])] )
    return T  
def DP_str(P):
    D=[]
    x1=sp.Symbol("x1")
    x2=sp.Symbol("x2")
    q1=sp.Symbol("q1")
    q2=sp.Symbol("q2")
    r3=sp.Symbol("r3")
    r4=sp.Symbol("r4")
    X=[x1,x2,q1,q2]
    for Pi in P:
        Pi=Pi.replace("^","**")
        D.append(parse_expr(Pi))
    jacobian_P= sp.Matrix(D).jacobian(sp.Matrix(X)) 
    minor_P=jacobian_P[:,2:]
    T=minor_P *sp.Matrix([r3,r4])
    answer=[]
    for i in range(len(P)):
        t=str(T[i,0])
        t=t.replace("**","^")
        answer.append(t)
    return answer

def generating_system(P,x1,x2,t):
    V=""" Variables \n
    x1 in """ + str(x1) + """; 
    x2 in """ + str(x2) +"""; 
    q1 in [-3.14,3.14];
    q2 in [-3.14,3.14];
    r3 in [-1,1];
    r4 in [-1,1];
    t in """ +str(t)+ """ ; 
    Constraints \n """
    f= open("eq.txt","w+")
    f.write(V)
    for Pi in P:
        f.write(Pi+"=0;\n")
    S=DP_str(P)
    for Si in S:
          f.write(Si+"=0; \n")
    f.write("r3^2+r4^2-1=0;\n")
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
def finding_nodes(P,x1,x2,t):
    generating_system(P,x1,x2,t)
    os.system("ibexsolve   --eps-max=0.1 -s  eq.txt  > output.txt")
    g=open('output.txt','r')
    result=g.readlines()
    T=cb.computing_boxes(result)
    return T    
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
    for i in range(len(box1)):

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
    ftboxes=[ [d.ftconstructor(boxi[0],boxi[1]) for boxi in box ] for box in boxes ]
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
    return components      
                         

"""
P1="(x1 - 8*cos(q1))^2 + (x2 - 8*sin(q1) )^2 - 25"
P2="(x1 - 9 - 5* cos(q2) )^2 + (x2 - 5* sin(q2))^2 - 64"
P3="(2*x1 - 16*cos(q1))*(2*x2 - 10*sin(q2)) - (2*x2 - 16*sin(q1))*(2*x1 - 10*cos(q2) - 18)"
"""
P1="x1-q2^3"
P2="x2-q2^2"
P3="q1-q2"
P=[P1,P2,P3]
b=[[-0.9,0.2],[-0.3,0.6],[-0.9,0.2],[-0.9,0.2],[-1,1],[-1,1],[-0.1,0.92] ]
generating_system(P,b[0],b[1],b[6])

    
#f="boxes_first_branch_silhouette"



"""fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlim(0.75,1)
ax.set_ylim(6,6.25)
#plt.xticks(np.arange(-20, 20, 2.0))
#plt.yticks(np.arange(-20, 20, 2.0))
ax.set_xlabel('q1')
ax.set_ylabel('q2')



for box in T:
     rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
        box[0][1]-box[0][0],box[1][1]-box[1][0], fc='g')
     plt.gca().add_patch(rectangle)
plt.show()
input() 
"""



"""
components=connected_compnants(sorted)
t=estimating_t(components)

S=finding_nodes(P,b[0],b[1],t)
"""





"""
fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlim(3,6)
ax.set_ylim(10,14)
#plt.xticks(np.arange(-20, 20, 2.0))
#plt.yticks(np.arange(-20, 20, 2.0))
ax.set_xlabel('q1')
ax.set_ylabel('q2')



for box in T[:]:
    rectangle= plt.Rectangle((float(box[0].lower()),float(box[1].lower()) ), \
        float(box[0].upper())-float(box[0].lower()),float(box[1].upper())-float(box[1].lower()), fc='r')
    plt.gca().add_patch(rectangle)"""
"""for box in T[:]:
    rectangle= plt.Rectangle((float(box[2].lower()),float(box[3].lower()) ), \
        float(box[2].upper())-float(box[2].lower()),float(box[3].upper())-float(box[3].lower()), fc='r')
    plt.gca().add_patch(rectangle) """   


"""
for box in C[1]:
    rectangle= plt.Rectangle((float(box[2].lower()),float(box[3].lower()) ), \
        float(box[2].upper())-float(box[2].lower()),float(box[3].upper())-float(box[3].lower()), fc='g')
    plt.gca().add_patch(rectangle)
for box in C[2]:
    rectangle= plt.Rectangle((float(box[0].lower()),float(box[1].lower()) ), \
        float(box[0].upper())-float(box[0].lower()),float(box[1].upper())-float(box[1].lower()), fc='b')
    plt.gca().add_patch(rectangle)"""    
#plt.show()




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


"""
def connected_components(boxes): #not working
    index=[set()]*len(boxes)
    for i in range(len(boxes)-1):
       for j in range(i+1,len(boxes)):
        if d.boxes_intersection(boxes[i],boxes[j]):
            index[i].add(j)
            index[j].add(i)
    

    components=[index[0]]
    for i in index[0]:
        index[0] +=index[i]
    


    for i in range(1,len(boxes)):
        for component in components:
            if component.intersection(index[i]) !=[]:
                component += index[i]"""





            

    


      