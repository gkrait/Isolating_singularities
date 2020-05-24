
import matplotlib.pyplot as plt
import os
import pickle 
import draft as d
from pprint import pprint
from sympy.parsing.sympy_parser import parse_expr
import sympy as sp 
import os 

def normal_subdivision(B):
	ft_B=d.subdivide([d.ftconstructor(Bi[0],Bi[1]) for Bi in B[:]])
	return [d.ft_normal(Bi)  for Bi in ft_B]

def plane_subdivision(B):
	ft_B2=d.subdivide([d.ftconstructor(Bi[0],Bi[1]) for Bi in B[:2]])
	normal_B2=[d.ft_normal(Bi)  for Bi in ft_B2]
	return d.cartesian_product(normal_B2,[B[2:]])
def system_generator(f,B):
	g=open(f,"r")
	L=g.readlines()
	g.close()
	f= open("eq.txt","w+")
	f.write("Variables \n")
	f.write("x1 in " + str(B[0]) +" ; \n")
	f.write("x2 in " + str(B[1])+" ; \n")
	for i in range(2,len(B)):
		f.write("x" +str(i+1) + " in " + str(B[i]) +" ; \n")
	f.write("Constraints \n")
	for Li in L:
		f.write(Li.replace("\n","") +"=0; \n")
	f.write("end ")
	f.close()
	return f 
def solving_with_ibex():
	uncer_content=[]
	cer_content=[]
	os.system("ibexsolve   --eps-max=0.1 -s  eq.txt  > output.txt")
	g=open('output.txt','r')
	result=g.read()
	with open('output.txt') as f:
		if "successful" in result:
			cer_content = f.readlines()
		elif  "infeasible problem" not in result and "done! but some boxes" in result:
			uncer_content = f.readlines()
		elif "infeasible problem" in result:
			uncer_content="Empty"
			cer_content="Empty"




	return [cer_content,uncer_content]			
def computing_boxes(content):
 i=0
 Answer=[]

 for fi in content:
 	try:
 	 a=fi.index('(')
 	 b=fi.index(')')
 	 T=(fi[a:b+1]).replace('(','[')
 	 T=T.replace(')',']')
 	 T=T.split(";")

 	 E=[]
 	 i=0
 	 for Ti in T:
 	 	#print(Ti)
 	 	Ti= Ti.replace('[',"")
 	 	Ti= Ti.replace(']',"")
 	 	Ti=Ti.replace('<','')
 	 	Ti=Ti.replace('>','')
 	 	x=Ti.index(",")
 	 	#print(Ti)
 	 	#input()
 	 	a=float(Ti[:x])
 	 	b=float(Ti[x+1:])
 	 	E.append([])
 	 	E[i]=[a,b]
 	 	i+=1
 	 Answer.append(E)
 	except ValueError:
 		#print("Hi")
 		k=1
 	  
 return Answer
def ploting_boxes(boxes,uncer_boxes, var=[0,1], B=[[-20,20],[-20,20]],a=1,b=10,nodes=[],color="g"):
   fig, ax = plt.subplots()
   plt.grid(True)
   ax.set_xlim(B[0][0], B[0][1])
   ax.set_ylim(B[1][0], B[1][1])
   ax.set_xlabel(r'$x_1$')
   ax.set_ylabel(r'$x_2$')
   c=0
   for box in boxes:
     rectangle= plt.Rectangle((box[var[0]][0],box[var[1]][0]) , \
    	a*(box[var[0]][1]-box[var[0]][0]),a*(box[var[1]][1]-box[var[1]][0]), fc=color)
     plt.gca().add_patch(rectangle)
   for box in uncer_boxes:
    	rectangle= plt.Rectangle((box[var[0]][0],box[var[1]][0]) , \
    	a*(box[var[0]][1]-box[var[0]][0]),a*(box[var[1]][1]-box[var[1]][0]), fc='r')
    	plt.gca().add_patch(rectangle)

   for box in nodes:
     rectangle= plt.Rectangle((box[0][0]-0.25,box[1][0]-0.25) , \
    	(box[0][1]-box[0][0]+0.5),(box[1][1]-box[1][0]+0.5), fc='y',fill=None)
     plt.gca().add_patch(rectangle)   



   plt.show()
def solver(f,B): #Assumption: no cusps
	L=[B]
	certified_boxes=[]
	uncertified_boxes=[]
	while len(L) !=0:
		system=system_generator(f,L[0])
		ibex_output=solving_with_ibex()
		if ibex_output[0]== "Empty":
			L.remove(L[0])
		elif len(ibex_output[0]) !=0:
			certified_boxes +=computing_boxes(ibex_output[0])
			L.remove(L[0])
		elif len(ibex_output[1])!=0:
			uncertified_boxes +=computing_boxes(ibex_output[1])
			L.remove(L[0])
		else:
			children=plane_subdivision(L[0])
			L.remove(L[0])
			L += children
	L=eval_file_gen(f,certified_boxes)
	
	
	print([Li for Li in L])
	while L.replace('\n',"") != "[]":
		L=L.replace('[','')
		L=L.replace(']','')
		L=L.replace('\n','')
		L=L.split(",")
		print(L)
		input()
		for i in L:
			#print(i)
			#input()
			children=normal_subdivision(certified_boxes[int(i)])
			certified_boxes.remove(certified_boxes[int(i)])
			certified_boxes +=children
		L =  eval_file_gen(f,certified_boxes)
		#print(L)
		#input()
 	    

	return [certified_boxes,uncertified_boxes]		

def eval_file_gen(f,boxes,special_function=[]): #condition: len(boxes[0]) is even
	functions=["sin","cos","tan","exp"]+special_function
	#computing P as sympy exprision
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
	#m2=sp.lambdify(X,minor2)
	fil=open("evaluation_file.py","w")
	fil.write("import flint as ft \n")
	fil.write("import sympy as sp \n")
	fil.write("import draft as d \n")
	#fil.write("from computing_boxes import ploting_boxes \n")
	fil.write("boxes="+str(boxes)+"\n")
	fil.write("ftboxes=[ [d.ftconstructor(Bi[0],Bi[1]) for Bi in B ] for B in boxes ] \n"    )
	fil.write("n=len(boxes[0])\n")
	fil.write("m=len(boxes)\n")
	fil.write("m1=[]\n")
	fil.write("m2=[]\n")
	#fil.write("X=[] \n")
	#fil.write("for i in range(n): \n")
	#fil.write("  X.append(sp.Symbol(\"x\"  +str(i+1)))\n")
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
	#cmd=os.system("python3 evaluation_file.py ")
	t=os.popen("python3 evaluation_file.py ").read()
	return t


X=[]
for i in range(4):
	X.append(sp.Symbol("x"+str(i+1)))

f="equations.txt" 
B=[[-5,15],[-15,15],[-3.14,3.14],[-3.14,3.14]]
#B=[d.ftconstructor(Bi[0],Bi[1]) for Bi in B ]
T=solver(f,B)

#print(len(T[0]),len(T[1]))



"""
pickle_in=open("boxes_first_branch_silhouette","rb")
boxes=pickle.load(pickle_in)[0]
pickle_in.close()


(eval_file_gen(f,boxes))
"""

"""T=solver(f,B)
pickle_out=open("boxes_first_branch_silhouette","wb")
pickle.dump(T,pickle_out)
pickle_out.close()

ploting_boxes(T[0],T[1])"""

"""
for i in range(-20,20):
	for j in range(-20,20):
		f=system_generator(i,j)
		T=solving_with_ibex(f)
		content=T[0]
		uncer_content=T[1]
		if len(content) != 0:
			Answer=computing_boxes(content)
			projection += [Ti for Ti in Answer ]
		if len(uncer_content) !=0:
		   Answer=computing_boxes(uncer_content)
		   uncer_projection += Answer







ploting_boxes(projection,uncer_projection)


"""




#os.system("cd Documents")

"""
f= open("eq.txt","w+")
f.write("Variables \n")
f.write("x1 in [" + str(4)+ "," + str(5) +" ]; \n")
f.write("x2 in [" + str(11.9)+ "," + str(12.2) +" ]; \n")
L=q1 in [0.5,1.2];
q2 in [1.2,2.4];
Constraints
(x1 - 8*cos(q1))^2 + (x2 - 8*sin(q1) )^2 - 25=0;
(x1 - 9 - 5* cos(q2) )^2 + (x2 - 5* sin(q2))^2 - 64 = 0;
(16*(x1 - 8*cos(q1))*sin(q1) - 16*(x2 - 8*sin(q1))*cos(q1))*(-10*(x2 - 5*sin(q2))*cos(q2) + 10*(x1 - 5*cos(q2) - 9)*sin(q2))=0;
end  
f.write(L)
f.close()
content=solving_with_ibex(f)[0]
T=computing_boxes(content)
print(len(T))
input()

fig, ax = plt.subplots()
plt.grid(True)

ax.set_xlim(4,5)
ax.set_ylim(11.2,12.4)
#plt.xticks(np.arange(-6,14, 2.0))
#plt.yticks(np.arange(-13,13,2))
ax.set_xlabel('x1')
ax.set_ylabel('x2')

for box in T:
     rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='g',label='branch 1')
     plt.gca().add_patch(rectangle)
plt.show()     """