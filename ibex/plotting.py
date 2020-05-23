import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
#import singularities as s
from pprint import pprint
import computing_boxes as cb
pickle_in=open("boxes_first_branch_silhouette","rb")
branch1=pickle.load(pickle_in)
pickle_in.close()
B=[[-20,20],[-20,20],[-3.14,3.14],[-3.14,3.14]]
pickle_in=open("boxes_second_branch","rb")
branch2=pickle.load(pickle_in)
pickle_in.close()

cb.ploting_boxes(branch1[0],branch2[0],B=B[2:],var=[2,3],b=10000000)

input()
equations="equations.txt"
f="boxes_second_branch"




nodes=s.solving_fornodes(equations,f,B)

fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlim(-20, 20)
ax.set_ylim(-20, 20)
#plt.xticks(np.arange(-6,14, 2.0))
#plt.yticks(np.arange(-13,13,2))
ax.set_xlabel('x1')
ax.set_ylabel('x2')





green_patch = mpatches.Patch(color='green', label='Branch 1 (Silhouette)')
blue_patch = mpatches.Patch(color='blue', label='Branch 2')
black_patch = mpatches.Patch(color='black', label='Certified nodes',fill=None)
plt.legend(handles=[green_patch,blue_patch,black_patch])

for box in branch1[0]:
     rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='g',label='branch 1')
     plt.gca().add_patch(rectangle)
    
for box in branch2[0]:
     rectangle1= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='b')
     plt.gca().add_patch(rectangle1) 
    

for box in branch2[1]:
     rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='r')
     plt.gca().add_patch(rectangle)  

for box in branch1[1]:
     rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='r')
     plt.gca().add_patch(rectangle)            


 
for box in nodes:
     rectangle= plt.Rectangle((box[0][0]-0.25,box[1][0]-0.25) , \
    	(box[0][1]-box[0][0]+0.5),(box[1][1]-box[1][0]+0.5), fc='y',fill=None)
     plt.gca().add_patch(rectangle)  
   
plt.show()     
