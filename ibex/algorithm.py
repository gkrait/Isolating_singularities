import matplotlib.pyplot as plt
import os
import pickle 
import draft as d
from pprint import pprint
import computing_boxes as cb
import singularities as s

equations="equations1.txt" 
B=[[-3,3],[-1,6.1],[-3.14,3.14]]


"""
f=open("output.txt","r")
g=f.readlines()
print(cb.computing_boxes(g))
"""


T=cb.solver(equations,B)
pickle_out=open("example2","wb")
pickle.dump(T,pickle_out)
pickle_out.close()

f="example2"
nodes_boxes=s.solving_fornodes(equations,f,B)
#print(T)
#input()

cb.ploting_boxes(T[0],[],B=[[-3.1,3.1],[-1,6.1]],nodes=nodes_boxes)
