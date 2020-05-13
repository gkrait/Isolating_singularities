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
    P_pluse=P.replace("q1","(q1+r3*sqrt(t))")
    P_pluse=P_pluse.replace("q2","(q2+r4*sqrt(t))")
    P_minus=P.replace("q1","(q1-r3*sqrt(t))")
    P_minus=P_minus.replace("q2","(q2-r4*sqrt(t))")
    SP= "0.5*(" + P_pluse + "+" +P_minus+")=0; \n"
    DP= "0.5*(" + P_pluse + "- (" +P_minus+") )/(sqrt(t))=0; \n"
    return [SP,DP]
def generating_system(P,x1,x2,t):
    V=""" Variables \n
    x1 in """ + str(x1) + """; 
    x2 in """ + str(x2) +"""; 
    q1 in [-3.14,3.14];
    q2 in [-3.14,3.14];
    r3 in [-1.1, 1.1] ;
    r4 in [-1.1,1.1] ;
    t in """ +str(t)+ """ ; 
    Constraints \n """
    f= open("eq.txt","w+")
    f.write(V)
    for Pi in P:
        f.write(SDP_str(Pi)[0])
        f.write(SDP_str(Pi)[1])
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
def estimating_r(components):
  r31=5
  r32=0
  r41=5
  r42=0
  for box1 in components[0]:
    for box2 in components[1]:
        norm_q1q2=d.distance(box1[2:],box2[2:])
        q1q2=[box1[i]-box2[i] for i in range(2,len(box1)) ]
        r3,r4=[ ri/norm_q1q2 for ri in q1q2 ]
        if r31 > r3.lower():
            r31=float(r3.lower())
        if r32 < r3.upper():
            r32=float(r3.upper()) 
        if   r41 > r4.lower(): 
          r41=float(r4.lower())
        if r42 < r4.upper():
           r42= r4.upper()
   
  return [[r31,r32],[r41,r42]]
P1="(x1 - 8*cos(q1))^2 + (x2 - 8*sin(q1) )^2 - 25"
P2="(x1 - 9 - 5* cos(q2) )^2 + (x2 - 5* sin(q2))^2 - 64"
P3="(2*x1 - 16*cos(q1))*(2*x2 - 10*sin(q2)) - (2*x2 - 16*sin(q1))*(2*x1 - 10*cos(q2) - 18)"
P=[P1,P2,P3]




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








f="boxes_first_branch_silhouette"
b=[[4.5,5.1],[11,13.15],[-3.15,3.15],[-3.15,3.15],[-1,1],[-1,1],[0.1,9.2] ]

"""
with open('output.txt') as f:
            cer_content = f.readlines()  
T=cb.computing_boxes(cer_content)            

pickle_out=open("additional1","wb")
pickle.dump(T,pickle_out)
pickle_out.close()

"""




pickle_in=open("gabe_firstbranch","rb")
T=pickle.load(pickle_in)
pickle_in.close()

pickle_in=open("additional","rb")
T1=pickle.load(pickle_in)
pickle_in.close()

pickle_in=open("additional1","rb")
T2=pickle.load(pickle_in)
pickle_in.close()






S=T+T1+T2

pickle_out=open("gape1","wb")
pickle.dump(S,pickle_out)
pickle_out.close()




"""
sorted=boxes_sort(T)

components=connected_compnants(sorted)
print(len(components))
t=estimating_t(components)
r=estimating_r(components)




S=finding_nodes(P,b[0],b[1],[t[0],t[1]])
"""






fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlim(1, 1.44)
ax.set_ylim(1.9, 2)
#plt.xticks(np.arange(-20, 20, 2.0))
#plt.yticks(np.arange(-20, 20, 2.0))
ax.set_xlabel('q1')
ax.set_ylabel('q2')

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





            

    


      