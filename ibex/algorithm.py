import matplotlib.pyplot as plt
import os
import pickle 
import draft as d
from pprint import pprint
import computing_boxes as cb
import singularities as s

equations="equations.txt" 
B=[[-20,20],[-20,20],[-3.14,3.14],[-3.14,3.14]]

#B=[[-3.1,3.1],[-1,6.1],[-3.14,3.14]]





T=cb.solver(equations,B)
pickle_out=open("enclosing boxes","wb")
pickle.dump(T,pickle_out)
pickle_out.close()

f="enclosing boxes"
pickle_in=open(f,"rb")
curve=pickle.load(pickle_in)
pickle_in.close()


nodes_boxes=s.solving_fornodes(equations,f,B)

#input()

cb.ploting_boxes(curve[0],[],B=B[2:],var=[2,3], nodes=nodes_boxes,b=10000000)
