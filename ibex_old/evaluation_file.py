import flint as ft 
import sympy as sp 
import draft as d 
def eval_func():
  pass 
boxes=[[[-0.009913982108861342, -0.0025158218616127435], [3.00503940841096, 3.019805054981555], [1.5714252840296008, 1.5732748255538804]], [[-0.0042373384016543255, 0.005435537798402309], [2.9915357528894244, 3.0108605107175688], [-1.5718556596522972, -1.5694374440838288]], [[-0.002524768433437398, 0.008492953862509044], [2.9830193374427036, 3.005041025468847], [1.5686730898987435, 1.5714275167801817]]]
ftboxes=[ [d.ftconstructor(Bi[0],Bi[1]) for Bi in B ] for B in boxes ] 
n=len(boxes[0])
m=len(boxes)
m1=[]
m2=[]
for B in ftboxes: 
   m1.append((-ft.arb.sin(B[2])**4 - 3)*ft.arb.sin(B[2]) + 4*ft.arb.sin(B[2])**3*ft.arb.cos(B[2])**2) 
   m2.append(-2*(ft.arb.sin(8*B[2]) + 3)*ft.arb.sin(B[2])*ft.arb.cos(B[2]) - 8*ft.arb.sin(B[2])**2*ft.arb.cos(8*B[2])) 
innrer_loops=[i for i in range(m) if 0 in m1[i] and 0 in m2[i] ]
print(innrer_loops)
