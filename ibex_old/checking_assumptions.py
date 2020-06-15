import sympy as sp
import computing_boxes as cb
from sympy.parsing.sympy_parser import parse_expr
from pprint import pprint
import flint as ft
import draft as d
import os 
import numpy as np
from evaluation_file import eval_func

def evaluation_exp(expr,B,X):
	expr_int=expr
	n=len(B_Ball)
	n=int((n+1)/2)
	for i in range(n):
		expr_int=expr_int.replace("x"+str(i+1),"B_f["+str(i)+"]")
	for i in range(n,len(X)-1):
		expr_int=expr_int.replace("r"+str(i+3-n),"B_f["+str(i)+"]")
	expr_int=expr_int.replace("t","B_f["+str(len(X)-1)+"]")	
	f=open("evaluation_file.py","w")
	f.write("import flint as ft\n")
	f.write("import draft as d \n")
	f.write("from sympy.parsing.sympy_parser import parse_expr \n")
	f.write("from sympy import * \n")
	f.write("from numpy import * \n")
	f.write("import operator \n")
	f.write("def eval_func(): \n")
	for i in range(n):
	 f.write(" x"+str(i+1)+"=Symbol(\""+ str(X[i])+"\"  )\n")
	for i in range(n,len(X)-1):
		f.write(" r"+str(i-n+3)+"=Symbol(\""+ str(X[i])+"\"  )\n")
	f.write(" t"+"=Symbol(\""+ str(X[len(X)-1])+"\"  )\n")
	f.write(" B="+str(B)+"\n")
	f.write(" f=srepr( parse_expr ( \"  " +str(expr)+ "  \") ) \n")
	f.write(" B_f=[ d.ftconstructor(Bi[0],Bi[1]) for Bi in B ]  \n"    )
	for i in range(len(X)):
		f.write(" f=f.replace(\"Symbol(\'"+str(X[i])+\
			"\')\", \" B_f[ "+str(i) +"] \")  \n")
	f.write(" f=f.replace(\"Add\",\"d.sevr_add\")\n")	
	f.write(" f=f.replace(\"Mul\",\"d.sevr_mul\")\n")
	f.write(" f=f.replace(\"Pow\",\"d.power_interval\")\n")
	f.write(" f=f.replace(\"Integer\",\"int\")\n")
	f.write(" f=f.replace(\"Float\",\"float\")\n")
	f.write(" f=f.replace(\", precision=53\", \"\")\n")
	#f.write("B_f=[ d.interv(d.ftconstructor(Bi[0],Bi[1])) for Bi in B ]  \n"    )
	f.write(" return eval(f) \n")
	f.close() 
	#from  evaluation_file import  eval_func
	answer=ft.arb(eval_func())

	return [float(answer.lower()),float(answer.upper())]
def Cauchy_form_poly(equations,X):
	#####Computing the Taylor polynomial#######################
	n=len(X)
	P_str=[Pi.replace('\n','') for Pi in  open(equations,'r').readlines()]
	P_pluse=[P_stri for P_stri in P_str]
	P_minus=[P_stri for P_stri in P_str]
	for i in range(2,n):
		P_pluse=[P_plusei.replace("x"+str(i+1),"(x"+str(i+1) + "+ r"+str(i+1) +"*t)")  for P_plusei in P_pluse]
		P_minus=[P_minusi.replace("x"+str(i+1),"(x"+str(i+1) + "- r"+str(i+1) +"*t)") for P_minusi in P_minus]
	DP= [("0.5*(" + P_plusei + "- (" +P_minusi+") )").replace("^","**") for P_plusei,P_minusi in zip(P_pluse,P_minus) ]
	#X_Ball=list(parse_expr(DP[0]).free_symbols)
	#t_in= X_Ball.index(sp.Symbol('t'))
	#t=X_Ball[t_in]
	t=sp.Symbol('t')
	D_ser=[sp.series(parse_expr(DPi),t).removeO() for DPi in DP]
	D_ser=[sp.expand(seri/t)  for seri in  D_ser ]
	D_ser=[seri.subs(t,sp.sqrt(t)) for seri in D_ser]
	return D_ser
def Ball_cusp_gen(equations,B_Ball,X):
	n=len(X)
	DP=Cauchy_form_poly(equations,X)
	t=sp.Symbol("t")
	coeff_t2= [ DPi.coeff(t**2)  for DPi in DP]
	eval_coefft2=[]
	D_list=[]
	X_Ball=X+[sp.Symbol("r"+str(i+1)) for i in range(2,n)]+[t]
	###Adding the remainder ##################
	for i  in range(len(DP)):

		coeff_t2=DP[i].coeff(t**2)

		ev=evaluation_exp(str(coeff_t2),B_Ball,X_Ball)
		eval_coefft2.append(ev)
		DP[i] -= coeff_t2 *t**2
		ci=sp.Symbol("c"+str(i+1))
		DP[i]=sp.simplify(DP[i])+0.5*ci *t**2

	V=""" Constants \n """
	for i in range(len(DP)):
		V+="c"+str(i+1)+" in " + str(eval_coefft2[i])+" ; \n"
	V +=""" Variables \n """
	for i in range(n):
		V += "x" +str(i+1) + " in " + str(B_Ball[i]) +" ; \n"
	for i in range(n,2*n-2):	
		V += "r" +str(i-n+3) + " in " + str(B_Ball[i]) +" ; \n"
	V += "t" + " in " + str(B_Ball[2*n-2]) +" ; \n" 
	V +="Constraints \n"   
	P=open(equations,'r').readlines()
	for Pi in P:
		V += Pi.replace('\n','') +"=0; \n"
	for Di in DP:
		D= str(Di).replace('\n','') +"=0; \n"
		D=D.replace("**","^")
		V +=D
	last_eq=""
	for i in range(3,n):
		last_eq += "r"+str(i)+"^2+"
	last_eq += "r" +str(n)+"^2 -1=0;" 
	V += last_eq +"\n"  
	f= open("eq.txt","w+") 
	f.write(V) 
	f.write("end")
	f.close()  	
def Ball_node_gen(equations,B_Ball,X):
    P=open(equations,"r").readlines()
    P=[Pi.replace('\n','') for Pi in P]
    n=len(X)
    V=""" Variables \n """
    for i in range(n):
        V += "x" +str(i+1) + " in " + str(B_Ball[i]) +" ; \n"
    for i in range(n,2*n-2):
        V += "r" +str(i-n+3) + " in " + str(B_Ball[i]) +" ; \n" 
    V += "t" + " in " + str(B_Ball[2*n-2]) +" ; \n"       
    V +="Constraints \n"   
    for Pi in P:
        V += cb.SDP_str(Pi,X)[0]
        V += cb.SDP_str(Pi,X)[1]
    last_eq=""
    for i in range(3,n):
        last_eq += "r"+str(i)+"^2+"
    last_eq += "r" +str(n)+"^2 -1=0;"    
    V += last_eq +"\n"
    f= open("eq.txt","w+")
    f.write(V) 
    f.write("end")
    f.close() 
def Ball_solver(equations,B_Ball,X):     #the width condition needs to be added
	L=[B_Ball]
	certified_boxes=[]
	uncertified_boxes=[]
	n=len(X)
	while len(L) !=0: 
	
		solvability=1
		if B_Ball[2*n-2][0] <= 0 <=   B_Ball[2*n-2][1] and \
		d.width([ d.ftconstructor(Bi[0],Bi[1]) for Bi in L[0] ] ) <0.1 :
			Ball_cusp_gen(equations,B_Ball,X)
		elif (B_Ball[2*n-2][0] > 0 or 0 >   B_Ball[2*n-2][1] ) \
		and d.width([ d.ftconstructor(Bi[0],Bi[1]) for Bi in L[0] ] ) <0.1:
			Ball_node_gen(equations,B_Ball,X)
		else:
			children=cb.plane_subdivision(L[0])
			L.remove(L[0])
			L += children
			solvability=0
		if solvability==1:
			ibex_output=cb.solving_with_ibex()
			if ibex_output[0]== "Empty":
		    
			  L.remove(L[0])
			elif len(ibex_output[0]) !=0:  
		    
			   certified_boxes +=cb.computing_boxes(ibex_output[0])
			   L.remove(L[0])
			elif len(ibex_output[1])!=0:   
		    
			  uncertified_boxes +=cb.computing_boxes(ibex_output[1])
			  L.remove(L[0])
			else:  
			  children=cb.plane_subdivision(L[0])
			  L.remove(L[0])
			  L += children
		
	return [certified_boxes,uncertified_boxes]		
		
		
    


X=[]

for i in range(3):
	X.append(sp.Symbol("x"+str(i+1)))
equations="equations3.txt"
B=[[-0.5,1],[-0.5,1],[-0.5,1]]
B_Ball=[[-0.1,3.1],[-0.1,3.1],[-1.6,1.6],[-1.1,1.1],[-0.1,6.2]]
#B=[ d.ftconstructor(Bi[0],Bi[1]) for Bi in B ]
"""T=Ball_solver(equations,B_Ball,X)

cb.ploting_boxes(T[0],T[1],B=[[-0.1,3.1],[-0.1,3.1]],a=1000)
"""