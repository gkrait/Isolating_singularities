import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches

pickle_in=open("boxes_first_branch_silhouette","rb")
branch1=pickle.load(pickle_in)
pickle_in.close()

pickle_in=open("boxes_second_branch","rb")
branch2=pickle.load(pickle_in)
pickle_in.close()



fig, ax = plt.subplots()
plt.grid(True)
ax.set_xlim(3,4)
ax.set_ylim(-0.2,0.2)
"""plt.xticks(np.arange(-20, 20, 2.0))
plt.yticks(np.arange(-20, 20, 2.0))"""
ax.set_xlabel('x')
ax.set_ylabel('y')
c=0

blue_patch = mpatches.Patch(color='blue', label='Branch 1 (Silhouette)')
plt.legend(handles=[blue_patch])
green_patch = mpatches.Patch(color='green', label='Branch 2')
red_patch = mpatches.Patch(color='red', label='proj of unknown')
plt.legend(handles=[blue_patch,green_patch,red_patch])

for box in branch1[0]:
     rectangle= plt.Rectangle((box[0][0],box[1][0]) , \
    	box[0][1]-box[0][0],box[1][1]-box[1][0], fc='g',label='branch 1')
     plt.gca().add_patch(rectangle)
rectangle.set_label("first branche")     
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
plt.show()     
