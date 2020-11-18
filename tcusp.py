 
import matplotlib.pyplot as plt
import os
import pickle 
from pprint import pprint
import sympy as sp
#import computing_boxes as cb
import interval_arithmetic as d
import flint as ft
import math
from sympy.parsing.sympy_parser import parse_expr
from evaluation_file import eval_func
#from computing_boxes import  computing_boxes
#import computing_boxes as cb
def inter_intersect(I1,I2):
    T=[]
    for i in range(len(I1)):
      T.append([max(I1[i][0],I2[i][0]),min(I1[i][1],I2[i][1])] )
    return T  
def DP_str(P):
    D=[]
    X=[ (sp.Symbol("x"+str(i))) for i in range(1,len(P)+2) ]
    r=[(sp.Symbol("r"+str(i))) for i in range(3,len(P)+2)]
    for Pi in P:
        Pi=Pi.replace("^","**")
        D.append(parse_expr(Pi))
    jacobian_P= sp.Matrix(D).jacobian(sp.Matrix(X)) 
    minor_P=jacobian_P[:,2:]
    T=minor_P *sp.Matrix([[ri] for ri in r  ])
    answer=[]
    for i in range(len(P)):
        t=str(T[i,0])
        t=t.replace("**","^")
        answer.append(t)
    return answer
def generating_system(P,B_Ball):
    n=len(P)+1
    V=""" Variables \n """
    for i in range(n):
        V += "x" +str(i+1) + " in " + str(B_Ball[i]) +" ; \n"
    for i in range(n,2*n-2):
        V += "r" +str(i-n+3) + " in " + str(B_Ball[i]) +" ; \n"

    #V += "t" + " in " + str(B_Ball[2*n-2]) +" ; \n"       
    V +="Constraints \n" 
    f= open("eq.txt","w+")
    f.write(V)
    for Pi in P:
        f.write(Pi+"=0;\n")
    S=DP_str(P)

    for Si in S:
          f.write(Si+"=0; \n")
    last_eq=""
    for i in range(3,len(P)+1):
        last_eq += "r"+str(i)+"^2+"
    last_eq += "r" +str(n)+"^2 -1=0;" 
    f.write(last_eq+"\n")
    f.write("end")
    f.close()
def intersting_boxes(f,b):
    pickle_in=open(f,"rb")
    curve=pickle.load(pickle_in)
    pickle_in.close()
    intersting_boxes=[]
    uncer_boxes=[]
    for box in curve[0]:
        if  b[0][0] <= box[0][0] <= box[0][1] <=b[0][1] and \
          b[1][0] <= box[1][0] <= box[1][1] <=b[1][1]:
             intersting_boxes.append(box)
    for box in curve[1]:
        if  b[0][0] <= box[0][0] <= box[0][1] <=b[0][1] and \
          b[1][0] <= box[1][0] <= box[1][1] <=b[1][1]:
             uncer_boxes.append(box)
    return [intersting_boxes,uncer_boxes]  
def finding_nodes(P,x1,x2,t):
    generating_system(P,x1,x2,t)
    os.system("ibexsolve   --eps-max=0.1 -s  eq.txt  > output.txt")
    g=open('output.txt','r')
    result=g.readlines()
    T=computing_boxes(result)
    return T    
def estimating_t(components,upper_bound=19.8):  #it works only if len(components)
  t1=upper_bound
  t2=0
  for box1 in components[0]:
    for box2 in components[1]:
        a=d.distance(box1,box2).lower()
        b=d.distance(box1,box2).upper()
    if t1 > a:
     t1=a 
    if t2<b:
     t2=b 
  t=d.ftconstructor(t1,t2)
  t=0.25*d.power_interval(t,2)     

  return [float(t.lower()),float(t.upper())]     
def boxes_compare(box1,box2):
    flage=0
    for i in range(len(box1)):

        if box1[i][0] > box2[i][0]: 
            return 1
        if  box1[i][0] < box2[i][0]: 
          return -1
    return 0      
def boxes_sort(boxes):
    sorted_boxes=boxes[:]
    for i in range(len(boxes)-1):
        for j in range(i+1,len(boxes)):
            if boxes_compare(sorted_boxes[i],sorted_boxes[j]) ==1:
                sorted_boxes[i], sorted_boxes[j] =sorted_boxes[j], sorted_boxes[i]
    return sorted_boxes            
def connected_compnants(boxes):
    ftboxes=[ [d.ftconstructor(boxi[0],boxi[1]) for boxi in box ] for box in boxes ]
    components=[[ftboxes[0]]]
    for i in range(1,len(ftboxes)):
        boxi_isused=0
        for j in  range(len(components)):
            membership=0
            for k in range(len(components[j])):   
                if d.boxes_intersection(ftboxes[i],components[j][k]) !=[]:
                    components[j].append(ftboxes[i])
                    membership=1
                    boxi_isused=1
                    break
                
            if membership==1:
              break 
        if boxi_isused==0:
          components.append([ftboxes[i]])
    return components      
def decimal_str(x: float, decimals: int = 10) -> str:
    return format(x, f".{decimals}f").lstrip().rstrip('0')                        


def evaluation_exp(expr,B,X):
    expr_int=expr
    n=len(B)
    n=int((n+1)/2)
    for i in range(n):
        expr_int=expr_int.replace("x"+str(i+1),"B_f["+str(i)+"]")
    for i in range(n,len(X)-1):
        expr_int=expr_int.replace("r"+str(i+3-n),"B_f["+str(i)+"]")
    expr_int=expr_int.replace("t","B_f["+str(len(X)-1)+"]") 
    f=open("evaluation_file.py","w")
    f.write("import flint as ft\n")
    f.write("import interval_arithmetic as d \n")
    f.write("from sympy.parsing.sympy_parser import parse_expr \n")
    f.write("from sympy import * \n")
    f.write("from numpy import * \n")
    f.write("import operator \n")
    f.write("def eval_func(expr,B): \n")
    for i in range(n):
     f.write(" x"+str(i+1)+"=Symbol(\""+ str(X[i])+"\"  )\n")
    for i in range(n,len(X)-1):
        f.write(" r"+str(i-n+3)+"=Symbol(\""+ str(X[i])+"\"  )\n")
    f.write(" t"+"=Symbol(\""+ str(X[len(X)-1])+"\"  )\n")
    #f.write(" B="+str(B)+"\n")
    f.write(" f=srepr( parse_expr ("  +"expr"+  ") ) \n")
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
    #f.write(" print(f) \n")
    f.write(" return eval(f) \n")
    f.close() 
    #from  evaluation_file import  eval_func
    #print(X);print(expr); print(eval_func(expr,B));input()
    answer=ft.arb(eval_func(expr,B))

    return [float(answer.lower()),float(answer.upper())]
def Cauchy_form_poly(P,X,var=[0,1]):
    #####Computing the Taylor polynomial#######################

    n=len(X)
    P_str=P[:]
    P_pluse=[P_stri for P_stri in P_str]
    P_minus=[P_stri for P_stri in P_str]
    for i in range(n):
        if i not in var:
         P_pluse=[P_plusei.replace("x"+str(i+1),"(x"+str(i+1) + "+ r"+str(i+1) +"*t)")  for P_plusei in P_pluse]
         P_minus=[P_minusi.replace("x"+str(i+1),"(x"+str(i+1) + "- r"+str(i+1) +"*t)") for P_minusi in P_minus]
    DP= [("0.5*(" + P_plusei + "- (" +P_minusi+") )").replace("^","**") for P_plusei,P_minusi in zip(P_pluse,P_minus) ]
    SP= [("0.5*(" + P_plusei + "+ (" +P_minusi+") )").replace("^","**") for P_plusei,P_minusi in zip(P_pluse,P_minus) ]
    #X_Ball=list(parse_expr(DP[0]).free_symbols)
    #t_in= X_Ball.index(sp.Symbol('t'))
    #t=X_Ball[t_in]
    t=sp.Symbol('t')
    D_ser=[sp.series(parse_expr(DPi),t).removeO() for DPi in DP]
    D_ser=[sp.expand(seri/t)  for seri in  D_ser ]
    D_ser=[seri.subs(t,sp.sqrt(t)) for seri in D_ser]
    return P_pluse
def Ball_cusp_gen(equations,B_Ball,X,var=[0,1]):
    
    n=len(X)
    P_pluse=Cauchy_form_poly(equations,X,var)
    t=sp.Symbol("t")
    X_Ball=X+[sp.Symbol("r"+str(i+1)) for i in range(2,n)]+[t]
    coeffs=[]
    coeffd=[]
    SP=[]
    DP=[]
    for i in range(len(P_pluse)):
     a0=sp.expand(P_pluse[i].replace("t","0"))
     a1=sp.diff(P_pluse[i],t)
     
     a2 =sp.diff(a1,t)
     a3=sp.diff(a2,t)
     a4=sp.diff(a3,t)
     a5=sp.diff(a4,t)
     a6=sp.diff(a5,t)
     a1=sp.expand(a1.subs(t,0))
     a2=sp.expand(a2.subs(t,0))
     a3=sp.expand(a3.subs(t,0))
     a4=sp.expand(a4.subs(t,0))

     cs=evaluation_exp(str(a6),B_Ball,X_Ball)
     cd=evaluation_exp(str(a5),B_Ball,X_Ball)
     coeffs.append(cs)
     coeffd.append(cd)
     SP.append(str(a0+0.5*(a2)*t + (1/24)*(a4)*t**2   ) +"+(1/(2*3*4*5*6))*t^3*cs"+ str(i+1) )
     DP.append(str(a1+(1/6)*(a3)*t   ) +"+(1/(2*3*4*5))*t^2*cd"+ str(i+1) )
    """ci1=sp.Symbol("cs"+str(1))
    ci2=sp.Symbol("cs"+str(2))
    ci1d=sp.Symbol("cd"+str(1))
    ci2d=sp.Symbol("cd"+str(2))
    X_Balls=X_Ball+[ci1,ci2]
    B_Balls=B_Ball +coeffs
    X_Balld=X_Ball+[ci1d,ci2d]
    B_Balld=B_Ball +coeffd
    for Pi in SP:   
        print(B_Ball)
        print(Pi);input()
        print(evaluation_exp(Pi.replace("^","**"),B_Balls,X_Balls));input()
    for Pi in DP:   
        print(B_Ball)
        print(Pi);input()
        print(evaluation_exp(Pi.replace("^","**"),B_Balld,X_Balld));input()    
    #print([evaluation_exp(Pi.replace("^","**"),B_Ball,X_Ball) for Pi in SP+DP ]);input() """
    

    """
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
        DP[i]=sp.simplify(DP[i])+0.5*ci *t**2 """ 
    V=""" Constants \n """
    for i in range(len(DP)):
        V+="cs"+str(i+1)+" in " + str(coeffs[i])+" ; \n" 
    for i in range(len(DP)):
        V+="cd"+str(i+1)+" in " + str(coeffd[i])+" ; \n"    

    V +=""" Variables \n """
    for i in range(n):
        V += "x" +str(i+1) + " in " + str(B_Ball[i]) +" ; \n"
    k=0
    for i in range(n):    
        if i not in var:
          V += "r" +str(i+1) + " in " + str(B_Ball[k+n]) +" ; \n"
          k+=1
    V += "t" + " in " +str(B_Ball[2*n-2]) +" ; \n" 
    #V += "t" + " in [" +str(0) +","+ str(B_Ball[2*n-2][1]) +"] ; \n" 
    V +="Constraints \n"   
    P=equations[:]
    for i in range(len(SP)):
        #V += P[i] +"=0 ; \n "
        V += SP[i].replace("**","^") +"=0; \n"
        V += DP[i].replace("**","^") +"=0; \n"

        """Pi= Pi.replace('\n','') +"=0; \n"
        Pi= Pi.replace('t','sqrt(t)') 
        V += Pi.replace("**","^")
    for Di in DP:
        D= str(Di).replace('\n','') +"=0; \n"
        D=D.replace("**","^")
        V +=D"""
    last_eq=""
    for i in range(n):
        if i not in var:
         last_eq += "r"+str(i+1)+"^2+"
    last_eq=last_eq[:-1]
    last_eq += " -1=0;" 
    V += last_eq +"\n"  
    f= open("eq.txt","w+") 
    f.write(V) 
    f.write("end")
    f.close() 
    


def computing_boxes():
  content=open("output.txt","r").readlines()
  for line in content[:]:
    if "number of cells" in line:
      tree_size= 2*int(line.split()[3])+1 #to compute tree size
      break
  if "infeasible" in open("output.txt","r").read():
    return ["Empty","Empty",tree_size]
  
  cer=[]; uncer=[]
  i=0
  Answer=[]
  for fi in content:
    try:
      a=fi.index('(')
      b=fi.index(')')
      T=(fi[a:b+1]).replace('(','[')
      T=(fi[a:b+1]).replace('(','[')
      T=T.replace(')',']')
      T=T.split(";")
      E=[]
      i=0
      for Ti in T:
        Ti= Ti.replace('[',"")
        Ti= Ti.replace(']',"")
        Ti=Ti.replace('<','')
        Ti=Ti.replace('>','')
        x=Ti.index(",")
        a=float(Ti[:x])
        b=float(Ti[x+1:])
        E.append([])
        E[i]=[a,b]
        i+=1
      if "solution n" in fi or "boundary n" in fi:
        cer.append(E)
      elif "unknown n" in fi:
        uncer.append(E)
    except ValueError:
          pass 
  return [cer,uncer,tree_size] 
def cusp_Ball_solver(P,B,X,var=[0,1],eps_min=1e-7,eps_max=0.1):
    Ball_cusp_gen(P,B,X,var) #check the evaluation_exp function
    os.system("ibexsolve   --eps-max="+ str(eps_max)+"  --eps-min="+ str(eps_min) + " -s  eq.txt  > output.txt")
    
    return computing_boxes()

#generating_system(P,B)
#print(ibex_output(P,B))


