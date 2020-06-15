import flint as ft 
import sympy as sp 
import draft as d 
def eval_func():
  pass 
boxes=[[[1.0812616307283927, 1.1193018077518966], [-6.295928488219059, -6.231079219690426], [-0.7565329157807386, -0.7516742890513373], [-1.627207229064572, -1.6179016586153576]], [[1.058480302626916, 1.1382109395327948], [-6.3407838156576855, -6.285324006051722], [-2.0475167411168167, -2.0352115821114802], [-1.633594976810668, -1.61570744046166]], [[1.0761512038281142, 1.0993217435008336], [-6.3357881324778305, -6.295928488165885], [-0.7588748405845178, -0.7559319465375081], [-1.6297229951369392, -1.6240352838421988]], [[1.0585715175600674, 1.1016244818218004], [-6.384158740805242, -6.313115216729011], [-0.762197533623846, -0.7567454377655007], [-1.6348797169353544, -1.624324363104786]]]
ftboxes=[ [d.ftconstructor(Bi[0],Bi[1]) for Bi in B ] for B in boxes ] 
n=len(boxes[0])
m=len(boxes)
m1=[]
m2=[]
for B in ftboxes: 
   m1.append(-(2*B[1] - 16*ft.arb.sin(B[2]))*(-10*(B[1] - 5*ft.arb.sin(B[3]))*ft.arb.cos(B[3]) + 10*(B[0] - 5*ft.arb.cos(B[3]) - 9)*ft.arb.sin(B[3]))*(16*(2*B[1] - 10*ft.arb.sin(B[3]))*ft.arb.sin(B[2]) + 16*(2*B[0] - 10*ft.arb.cos(B[3]) - 18)*ft.arb.cos(B[2])) - (2*B[1] - 10*ft.arb.sin(B[3]))*(16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) - 16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*(-10*(2*B[0] - 16*ft.arb.cos(B[2]))*ft.arb.cos(B[3]) + 10*(-2*B[1] + 16*ft.arb.sin(B[2]))*ft.arb.sin(B[3])) + (16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) - 16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*(-10*(B[1] - 5*ft.arb.sin(B[3]))*ft.arb.cos(B[3]) + 10*(B[0] - 5*ft.arb.cos(B[3]) - 9)*ft.arb.sin(B[3]))*(-32*ft.arb.cos(B[2]) + 20*ft.arb.cos(B[3]) + 36)) 
   m2.append(-(2*B[0] - 16*ft.arb.cos(B[2]))*(-10*(B[1] - 5*ft.arb.sin(B[3]))*ft.arb.cos(B[3]) + 10*(B[0] - 5*ft.arb.cos(B[3]) - 9)*ft.arb.sin(B[3]))*(16*(2*B[1] - 10*ft.arb.sin(B[3]))*ft.arb.sin(B[2]) + 16*(2*B[0] - 10*ft.arb.cos(B[3]) - 18)*ft.arb.cos(B[2])) - (16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) - 16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*(-10*(2*B[0] - 16*ft.arb.cos(B[2]))*ft.arb.cos(B[3]) + 10*(-2*B[1] + 16*ft.arb.sin(B[2]))*ft.arb.sin(B[3]))*(2*B[0] - 10*ft.arb.cos(B[3]) - 18) + (16*(B[0] - 8*ft.arb.cos(B[2]))*ft.arb.sin(B[2]) - 16*(B[1] - 8*ft.arb.sin(B[2]))*ft.arb.cos(B[2]))*(-10*(B[1] - 5*ft.arb.sin(B[3]))*ft.arb.cos(B[3]) + 10*(B[0] - 5*ft.arb.cos(B[3]) - 9)*ft.arb.sin(B[3]))*(32*ft.arb.sin(B[2]) - 20*ft.arb.sin(B[3]))) 
innrer_loops=[i for i in range(m) if 0 in m1[i] and 0 in m2[i] ]
print(innrer_loops)
