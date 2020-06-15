import draft as d
import matplotlib.pyplot as plt
import numpy as np
import flint  as ft
import sympy as sp
from copy import copy, deepcopy
import sys
import time
from sympy import *
import inspect
import math
from pprint import pprint

"""B=[ft.arb(0.5,0.5),ft.arb(2.5,0.5),ft.arb(3.5,0.5)]
P=[[[1,1,1],2],[[0,2,1],3]]
P=['sin',P]
X=d.genvars(4)
print(d.Ball_for_interval_poly([[[1,0,0,0],1],[[0,0,0,2],-1],[[0,0,0],1]],X) )"""


A = ft.acb_mat([[1,1],[1,-1]])
B=ft.acb_mat([[0],[0]])
print(A.solve(B))
