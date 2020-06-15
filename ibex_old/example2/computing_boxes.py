
import matplotlib.pyplot as plt
import os
import pickle 

def system_generator(i,j,B):
	g=open("equations.txt","r")
	L=g.read()
	g.close()
	f= open("eq.txt","w+")
	f.write("Variables \n")
	f.write("x1 in [" + str(i)+ "," + str(i+1) +" ]; \n")
	f.write("x2 in [" + str(j)+ "," + str(j+1) +" ]; \n")
	f.write("q1 in [" + str(B[2][0])+ "," + str(str(B[2][1])) +" ]; \n")
	f.write("q2 in [" + str(B[2][0])+ "," + str(str(B[2][1])) +" ]; \n")
	f.write("q3 in [" + str(B[2][0])+ "," + str(str(B[2][1])) +" ]; \n")
	f.write("Constraints \n")
	f.write(L)
	f.write("end")


	f.close()
	return f  
def solving_with_ibex(f):
	uncer_content=[]
	cer_content=[]
	os.system("ibexsolve   --eps-max=0.1 -s  eq.txt  > output.txt")
	g=open('output.txt','r')
	result=g.read()
	if "successful" in result:
		with open('output.txt') as f:
			cer_content = f.readlines()
		
	elif  "infeasible problem" not in result and "done! but some boxes have 'unknown' status" in result:
		with open('output.txt') as f:
			uncer_content = f.readlines()

	


	
	return [cer_content,uncer_content]			
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
def ploting_boxes(boxes,uncer_boxes):
   fig, ax = plt.subplots()
   plt.grid(True)
   ax.set_xlim(-20, 20)
   ax.set_ylim(-20,20)
   ax.set_xlabel('x')
   ax.set_ylabel('y')
   c=0

   for box in boxes:
     rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='g')
     plt.gca().add_patch(rectangle)

   for box in uncer_boxes:
    	rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='r')
    	plt.gca().add_patch(rectangle) 



   plt.show()
        
projection=[]
uncer_projection=[]

B=[[1,2],[1,2],[10,32],[10,32],[10,32] ]

for i in range(-20,20):
	for j in range(-20,20):
		f=system_generator(i,j,B)
		T=solving_with_ibex(f)
		content=T[0]
		uncer_content=T[1]
		if len(content) != 0:
			Answer=computing_boxes(content)
			projection += [Ti for Ti in Answer ]
		if len(uncer_content) !=0:
		   Answer=computing_boxes(uncer_content)
		   uncer_projection += Answer
pprint(Matrix(projection))		   


fig, ax = plt.subplots()
plt.grid(True)

ax.set_xlim(-20,20)
ax.set_ylim(-20,20)
#plt.xticks(np.arange(-6,14, 2.0))
#plt.yticks(np.arange(-13,13,2))
ax.set_xlabel('x1')
ax.set_ylabel('x2')

for box in projection:
     rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='g')
     plt.gca().add_patch(rectangle)
plt.show() 




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





pickle_out=open("boxes_second_branch","wb")
pickle.dump([projection,uncer_projection],pickle_out)
pickle_out.close()


ploting_boxes(projection,uncer_projection)


"""




#os.system("cd Documents")

