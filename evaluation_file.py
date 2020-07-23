import flint as ft
import interval_arithmetic as d 
from sympy.parsing.sympy_parser import parse_expr 
from sympy import * 
from numpy import * 
import operator 
def eval_func(expr,B): 
 x1=Symbol("x1"  )
 x2=Symbol("x2"  )
 x3=Symbol("x3"  )
 r3=Symbol("r3"  )
 t=Symbol("t"  )
 f=srepr( parse_expr (expr) ) 
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
