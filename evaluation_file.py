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
 B=[[-1.000000000144927, -0.9566492756997367], [-0.03801741441704251, 0.03801741448704236], [1.4850635712606615, 1.656529082461611], [-1.01, 1.01], [-1.5891889800639936e-11, 0.007350105470801903]]
 f=srepr( parse_expr ( "  -1561.6*r3**5*sin(4*x3)**2*cos(4*x3) + 512.0*r3**5*cos(4*x3)**3  ") ) 
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
