import flint as ft 
import sympy as sp 
import interval_arithmetic as d 
def  evlist(boxes): 
 ftboxes=[ [d.ftconstructor(Bi[0],Bi[1]) for Bi in B ] for B in boxes ] 
 n=len(boxes[0])
 m=len(boxes)
 m1=[]
 m2=[]
 for B in ftboxes: 
   m1.append(ft.arb(0.0002*ft.arb.sin(B[2])*ft.arb.cos(B[2]))) 
   m2.append( ft.arb(-2*ft.arb.sin(B[2])**2*ft.arb.cos(B[2]))) 
 innrer_loops=[i for i in range(m) if 0 in m1[i] and 0 in m2[i] ]
 x1mon=[i for i in range(m) if 0  not in m1[i] ]
 x2mon=[i for i in range(m) if 0 in m1[i] and 0 not in m2[i] ]
 return [x1mon,x2mon , innrer_loops] 
