
import matplotlib.pyplot as plt
import os

def system_generator(i,j):
	f= open("eq.txt","w+")
	f.write("Variables \n")
	f.write("x1 in [" + str(i)+ "," + str(i+1) +" ]; \n")
	f.write("x2 in [" + str(j)+ "," + str(j+1) +" ]; \n")
	L=""" q1 in [-3.15,3.15];
	q2 in [-3.15,3.15];
	Constraints
	//you can use C++ comments
	(x1 - 8*cos(q1))^2 + (x2 - 8*sin(q1) )^2 - 25=0;
	(x1 - 9 - 5* cos(q2) )^2 + (x2 - 5* sin(q2))^2 - 64 = 0;
	(2*x1 - 16*cos(q1))*(2*x2 - 10*sin(q2)) - (2*x2 - 16*sin(q1))*(2*x1 - 10*cos(q2) - 18)=0;
	end """ 
	f.write(L)
	f.close()
	return f  
def solving_with_ibex(f):
	os.system("ibexsolve   --eps-max=0.1 -s  eq.txt  > output.txt")
	g=open('output.txt','r')
	result=g.read()
	if "successful" in result:
		with open('output.txt') as f:
			content = f.readlines()
		return content
	else:
		return 1			
def computing_boxes(content):
 i=0
 Answer=[]

 for fi in content:
 	try:
 	 a=fi.index('[')
 	 b=fi.index(')')
 	 T=(fi[a:b+1]).replace('(','[')
 	 T=T.replace(')',']')
 	 T=T.split(";")
 	 E=[[]]*(content[len(content)-4].count('['))
 	 i=0
 	 for Ti in T:
 	 	Ti= Ti.replace('[',"")
 	 	Ti= Ti.replace(']',"")
 	 	x=Ti.index(",")
 	 	a=float(Ti[:x])
 	 	b=float(Ti[x+1:])
 	 	E[i]=[a,b]
 	 	i+=1
 	 Answer.append(E)
 	except ValueError:
 		k=1
 	  
 return Answer
def ploting_boxes(boxes):
   fig, ax = plt.subplots()
   plt.grid(True)
   ax.set_xlim(-20, 20)
   ax.set_ylim(-20,20)
   ax.set_xlabel('x')
   ax.set_ylabel('y')
   c=0

   for box in projection:
     rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='g')
     plt.gca().add_patch(rectangle)


   plt.show()
        
projection=[]

for i in range(-20,20):
	for j in range(-20,20):
		f=system_generator(i,j)
		content=solving_with_ibex(f)
		if content != 1:
			Answer=computing_boxes(content)
			projection += [Ti[:2] for Ti in Answer ]

ploting_boxes(projection)







#os.system("cd Documents")




		







