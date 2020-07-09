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
 B=[[2.9794870263337376, 3.000001332189351], [-1.3773451823517746e-10, 0.033837973759145185], [-0.11701202015186296, 0.09107253984754339], [-1.01, 1.01], [-3.3704916461985146e-11, 0.010824796171565189]]
 f=srepr( parse_expr ( "  -9426.13333333333*r3**5*sin(x3)**2*cos(16*x3) - 5632.26666666667*r3**5*sin(x3)*sin(16*x3)*cos(x3) - 0.8*r3**5*sin(x3)*cos(x3) + 688.0*r3**5*cos(x3)**2*cos(16*x3)  ") ) 
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
