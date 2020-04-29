import matplotlib.pyplot as plt
import os
import pickle 
from pprint import pprint
import sympy

a1=25
b1=-25
a2=25
b2=-25

pickle_in=open("boxes_second_branch","rb")
branch1=pickle.load(pickle_in)
pickle_in.close()

intersting_boxes=[]
b=[[3.2,3.6],[-0.1,0.15]]

for box in branch1[0]:
	if    b[0][0] <= box[0][0] <= box[0][1] <=b[0][1] and \
	b[1][0] <= box[1][0] <= box[1][1] <=b[1][1]:
	  intersting_boxes.append(box)
	  a1=min(a1,box[2][0])
	  b1=max(b1,box[2][1])
	  a2=min(a2,box[3][0])
	  b2=max(b2,box[3][1])


q1=[a1,b1]
q2=[a2,b2]

q1_minusq2=[[a1-b1,b1-a1] ,[a2-b2,b2-a2] ]


P1="(x1 - 8*cos(q1))^2 + (x2 - 8*sin(q1) )^2 - 25"
P2="(x1 - 9 - 5* cos(q2) )^2 + (x2 - 5* sin(q2))^2 - 64"
P3="(16*(x1 - 8*cos(q1))*sin(q1) - 16*(x2 - 8*sin(q1))*cos(q1))*(-10*(x2 - 5*sin(q2))*cos(q2) + 10*(x1 - 5*cos(q2) - 9)*sin(q2))"


The solver does not work


def SDP_str(P):
	P_pluse=P.replace("q1","q1+r3*sqrt(t)")
	P_pluse=P_pluse.replace("q2","q2+r4*sqrt(t)")
	P_minus=P.replace("q1","q1+r3*sqrt(t)")
	P_minus=P_minus.replace("q2","q2+r4*sqrt(t)")
	SP= "0.5*(" + P_pluse + "+" +P_minus+")=0;"
	DP= "0.5*(" + P_pluse + "-" +P_minus+")/(sqrt(t))=0;"
	return [SP,DP]

def generating_system():
    V=""" Variables 
    x1 in [3.2,3.6]; 
    x2 in [-0.1,0.15]; 
    q1 in [-0.5,0.5];
    q2 in [-1.5,1.5];
    r3 in [-1,1];
    r4 in [-1,1];
    t in [0.0000001,4];
    Constraints """
    f= open("eq.txt","w+")
    f.write(V)
    f.write(SDP_str(P1)[0])
    f.write(SDP_str(P1)[1])
    f.write(SDP_str(P2)[0])
    f.write(SDP_str(P2)[1])
    f.write(SDP_str(P3)[0])
    f.write(SDP_str(P3)[1])
    f.write("r3^2+r4^2-1=0;")
    f.write("end")

generating_system()
os.system("ibexsolve   --eps-max=0.1 -s  eq.txt  > output.txt")
g=open('output.txt','r')
result=g.read()
print(result)

"""
fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlim(-3.15,3.15)
ax.set_ylim(-3.15,3.15)
#plt.xticks(np.arange(-20, 20, 2.0))
#plt.yticks(np.arange(-20, 20, 2.0))
ax.set_xlabel('x')
ax.set_ylabel('y')
c=0

for box in intersting_boxes:
     rectangle= plt.Rectangle((a1,a2) , \
    	b1-a1,b2-a2, fc='g')
     plt.gca().add_patch(rectangle)
plt.show()     """
