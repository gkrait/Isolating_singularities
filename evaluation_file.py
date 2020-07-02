import flint as ft
import interval_arithmetic as d 
from sympy.parsing.sympy_parser import parse_expr 
from sympy import * 
from numpy import * 
import operator 
def eval_func(): 
 x1=Symbol("x1"  )
 x2=Symbol("x2"  )
 x3=Symbol("x3"  )
 r3=Symbol("r3"  )
 t=Symbol("t"  )
 B=[[2.9897822524640874, 3.0000162132915507], [-3.96301012839384e-11, 0.024572495447888507], [-0.06349793290964928, 0.08255727118566965], [-1.01, 1.01], [-1.018555118720288e-11, 0.0053330307317090665]]
 f=srepr( parse_expr ( "  -361.066666666667*r3**5*sin(x3)**2*cos(8*x3) - 384.266666666667*r3**5*sin(x3)*sin(8*x3)*cos(x3) - 0.8*r3**5*sin(x3)*cos(x3) + 88.0*r3**5*cos(x3)**2*cos(8*x3)  ") ) 
 B_f=[ d.ftconstructor(Bi[0],Bi[1]) for Bi in B ]  
 f=f.replace("Symbol('x1')", " B_f[ 0] ")  
 f=f.replace("Symbol('x2')", " B_f[ 1] ")  
 f=f.replace("Symbol('x3')", " B_f[ 2] ")  
 f=f.replace("Symbol('r3')", " B_f[ 3] ")  
 f=f.replace("Symbol('t')", " B_f[ 4] ")  
 f=f.replace("Add","d.sevr_add")
 f=f.replace("Mul","d.sevr_mul")
 f=f.replace("Pow","d.power_interval")
 f=f.replace("Integer","int")
 f=f.replace("Float","float")
 f=f.replace(", precision=53", "")
 return eval(f) 
