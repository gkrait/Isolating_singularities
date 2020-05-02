import matplotlib.pyplot as plt
import os
import pickle 
from pprint import pprint
import sympy
import computing_boxes as cb
import draft as d


def inter_intersect(I1,I2):
    T=[]
    for i in range(len(I1)):
      T.append([max(I1[i][0],I2[i][0]),min(I1[i][1],I2[i][1])] )
    return T  
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
    r3 in [-1,1];
    r4 in [-1,1];
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
    os.system("ibexsolve   --eps-max=0.1  -s eq.txt  > output.txt")
    g=open('output.txt','r')
    result=g.readlines()
    T=cb.computing_boxes(result)
    return T    

def estimating_t(boxes):
  if len(boxes)> 1:  
    t1=d.distance(boxes[0],boxes[1]).lower()
    t2=d.distance(boxes[0],boxes[1]).upper()
    for i in range(len(boxes)-1):
     for j in range(i+1,len(boxes)):
        length_ij=d.distance(boxes[i][2:],boxes[j][2:])
        if t1 > length_ij.lower():
            t1=length_ij.lower()
        if t2 < length_ij.upper():
         t1=length_ij.upper()
    t=d.ftconstructor(t1,t2)
    t=0.25*d.power_interval(t,2)     
  return t         




                          





P1="(x1 - 8*cos(q1))^2 + (x2 - 8*sin(q1) )^2 - 25"
P2="(x1 - 9 - 5* cos(q2) )^2 + (x2 - 5* sin(q2))^2 - 64"
P3="(2*x1 - 16*cos(q1))*(2*x2 - 10*sin(q2)) - (2*x2 - 16*sin(q1))*(2*x1 - 10*cos(q2) - 18)"
P=[P1,P2,P3]

f="boxes_second_branch"



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
b=[[0.8,1],[6.125,6.25],[-3.15,3.15],[-3.15,3.15],[-1,1],[-1,1],[0.4,2] ]

T=intersting_boxes(f,b)[0]
T=[[d.ftconstructor(Tij[0],Tij[1])  for Tij in Ti  ]  for Ti in T]
print(d.connected_components(T))
input()
t=[float(estimating_t(T).lower()),float(estimating_t(T).upper())]
print(t)
input()
t=[6,6.1]
print(finding_nodes(P,b[0],b[1],t))

"""
fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlim(-3.15,3.15)
ax.set_ylim(-3.15,3.15)
#plt.xticks(np.arange(-20, 20, 2.0))
#plt.yticks(np.arange(-20, 20, 2.0))
ax.set_xlabel('q1')
ax.set_ylabel('q2')
for box in C[0]:
    rectangle= plt.Rectangle((float(box[0].lower()),float(box[1].lower()) ), \
        float(box[0].upper())-float(box[0].lower()),float(box[1].upper())-float(box[1].lower()), fc='r')
    plt.gca().add_patch(rectangle)
for box in C[1]:
    rectangle= plt.Rectangle((float(box[2].lower()),float(box[3].lower()) ), \
        float(box[2].upper())-float(box[2].lower()),float(box[3].upper())-float(box[3].lower()), fc='g')
    plt.gca().add_patch(rectangle)
for box in C[2]:
    rectangle= plt.Rectangle((float(box[0].lower()),float(box[1].lower()) ), \
        float(box[0].upper())-float(box[0].lower()),float(box[1].upper())-float(box[1].lower()), fc='b')
    plt.gca().add_patch(rectangle)    
plt.show()
"""



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




def connected_components(boxes): #not working
    index=[set()*len(boxes)]
    for i in range(len(boxes)-1):
       for j in range(i+1,len(boxes)):
        if d.boxes_intersection(boxes[i],boxes[j]):
            index[i].add(j)
            index[j].add(i)
    
    used=[]
    components=[index[0]]
    for i in index[0]:
        for j in range(len(boxes)):
            if i in index[j]:
                component += index[j]
                used.append(j)
    


    """for i in range(1,len(boxes)):
        for component in components:
            if component.intersection(index[i]) !=[]:
                component += index[i]"""





            

    



    """ftboxes=[ [d.ftconstructor(boxi[0],boxi[1]) for boxi in box ] for box in boxes ]
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
          #print(i)
          #d.ftprint(ftboxes[i])
          #input()
          components.append([ftboxes[i]])
    S=[]
    index=[[]]*len(components)
    print(index)
    for i in range(len(components)-1):
      for j in range(i,len(components)):
         flage_k1=0
         for k1 in range(len(components[i])):
           flage_k2=0
           for k2 in range(len(components[j])):
              if d.boxes_intersection(components[i][k1],components[j][k2])!=[]:
                S.append(components[i]+components[j])
                index[i].append(j)
                #flage_k2=1
                #flage_k1=1
                #break
           if flage_k1 ==1:
             break

    return index """        