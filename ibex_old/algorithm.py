import matplotlib.pyplot as plt
import os
import pickle 
import draft as d
from pprint import pprint
import computing_boxes as cb
import sympy as sp 
from sympy.parsing.sympy_parser import parse_expr

def kinematic_equations(equations,X):
	manifold_str=[ Pi.replace("\n","") for Pi in   open(equations,"r").readlines()  ]
	n=len(X)
	manifold_exp= [Pi.replace('^','**') for Pi in manifold_str]
	manifold_exp= [parse_expr(Pi) for Pi in manifold_exp]
	jac_manifold=sp.Matrix(manifold_exp).jacobian(sp.Matrix(X))
	serial=open("serial.txt","w")
	serial.writelines([Pi+"\n" for Pi in  manifold_str])
	serial.write(str(jac_manifold[:,:int(n/2)].det()))
	serial.close()
	parallel=open("parallel.txt","w")
	parallel.writelines([Pi+"\n" for Pi in  manifold_str])
	parallel.write(str(jac_manifold[:,int(n/2):].det()))
	parallel.close()


def isolating_singularities(equations,B,X):
	kinematic_equations(equations,X) 


X=[]
for i in range(4):
	X.append(sp.Symbol("x"+str(i+1)))



equations="equations.txt" 
kinematic_equations(equations,X)
B=[[-5,15],[-15,15],[-3.14,3.14],[-3.14,3.14]]
#B=[[0.3,1.1],[-0.7,0.7],[-5,5]]


"""
sing=(cb.branche_sing(f,B,X))
cb.ploting_boxes(sing[0][0],sing[0][1], nodes=sing[1][0], cusps=sing[1][1])
#print(len(T[0]),len(T[1]))"""