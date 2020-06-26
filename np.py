import numpy as np
from numpy import linalg as la
import scipy.linalg
from numpy.linalg import inv
import math

def smoth_f(N_x, j):
 if j<N_x/2:
  d=j+1
 if j>=N_x/2:
  d=N_x-j
  
 x=1.0-d/3.0
 
 smoth_f_val=0
 if x<=0:
  smoth_f_val=1.0
 elif x>=1:
  smoth_f_val=0.0
 elif x>0 and x<1:
  smoth_f_val=0.5*(1.0-math.tanh((x-0.5)/(x-x**2)))

 return smoth_f_val




def Smo_1d(N_x, i, Model):
 value_val=0
 if Model[0]=="off":
  value_val=1.0
 if Model[0]=="Sin":
  value_val=math.sin( (math.pi/N_x)*i)**2
 if Model[0]=="Tanh":
   value_val=smoth_f(N_x, i)

 return value_val



def make_obc_h_1d(L, N, N_e,Model):
 N=int(N)
 E_total=[]
 Epsilon=L/(N+1.)
 h = np.zeros([N, N])
 for x1 in range(N):
  for x2 in range(N):
   if x1==x2+1:
      h[x1,x2] = -1.0*Smo_1d(N, x1, Model)
   if x1==x2-1:
      h[x1,x2] = -1.0*Smo_1d(N, x2, Model)
   if x1==x2:
      h[x1,x2] = 2.0*Smo_1d(N, x1, Model)
      #print x1, Smo_1d(N, x1, Model)

 #h = h #+ np.eye(N)*2
 #print h
 #h = (Epsilon**-2) * h

 w=np.linalg.eigvalsh(h)
 w_sorted=np.sort(w)
 #print w_sorted
 sum_val=0.
 for i in xrange(N_e):
  sum_val=sum_val+w_sorted[i]

 return sum_val 


def make_pbc_h_1d(L, N, N_e):
 N=int(N)
 E_total=[]
 Epsilon=L/N
 h = np.zeros([N, N])
 for x1 in range(N):
  for x2 in range(N):
   if x1==x2+1 or x1==x2-1:
      h[x1,x2] = -1.0
 h = h + np.eye(N)*2
 h[0,-1]=-1.0
 h[-1,0]=-1.0
 h = (Epsilon**-2) * h

 w=np.linalg.eigvalsh(h)
 w_sorted=np.sort(w)
 sum_val=0.
 for i in xrange(N_e):
   sum_val=sum_val+w_sorted[i]

 return sum_val 




def make_obc_h_2d(L, N, N_e):
 N=int(N)
 E_total=[]
 Epsilon=L/(N+1.)
 h = np.zeros([N*N, N*N])
 for x1 in range(N):
  for x2 in range(N):
   for y1 in range(N):
    for y2 in range(N):
     if ((abs(x1-x2)==1) & (abs(y1-y2)==0)) or ((abs(x1-x2)==0) & (abs(y1-y2)==1)):
      h[x1*N+y1,x2*N+y2] = -1.0


 h = h + np.eye(N*N)*4
 #print h
 h = (Epsilon**-2) * h

 w=np.linalg.eigvalsh(h)
 w_sorted=np.sort(w)
 sum_val=0.
 for i in xrange(N_e):
   sum_val=sum_val+w_sorted[i]

 return sum_val 




def make_pbc_h_2d(L, N, N_e):
 N=int(N)
 E_total=[]
 Epsilon=L/N
 h = np.zeros([N*N, N*N])
 for x1 in range(N):
  for x2 in range(N):
   for y1 in range(N):
    for y2 in range(N):
     if  ( ((abs(x1-x2)==1) or (abs(x1-x2+N)==1) or (abs(x1-x2-N)==1)) & (abs(y1-y2)==0) ) or  ( ( (abs(y1-y2)==1) or (abs(y1-y2+N)==1) or (abs(y1-y2-N) == 1)) & (abs(x1-x2)==0) ):
      h[x1*N+y1,x2*N+y2] = -1.0


 h = h + np.eye(N*N)*4
 #print h
 h = (Epsilon**-2) * h

 w=np.linalg.eigvalsh(h)
 w_sorted=np.sort(w)
 sum_val=0.
 for i in xrange(N_e):
   sum_val=sum_val+w_sorted[i]

 return sum_val 


L=6.
N=40
N_e=10
Model=["off"]
E_t=make_obc_h_1d(L, N, N_e, Model)
#print E_t*2.0



L=6.
N=16
N_e=10
#E_t=make_pbc_h_2d(L, N, N_e)
#print E_t



L=6.
N=16
N_e=14
E_t=make_obc_h_2d(L, N, N_e)
print E_t



# 
# 
# L=1.
# N=20
# N_e=3
# E_t=make_pbc_h_1d(L, N, N_e)
# print E_t









L=6.0
Start=4
End=10
E_list=[]
N_list=[]
for N in xrange(Start,End):
 template = np.zeros([N, N])
 for idx in np.ndindex(template.shape):
    if idx[0]==idx[1]+1 or idx[0]==idx[1]-1:
      template[idx] = (-16.*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2))
    if idx[0]==idx[1]+2 or idx[0]==idx[1]-2:
      template[idx] = (+1.*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2))
    if idx[0]==idx[1] and idx[1]==0 :
      template[idx] = (((30.0-1)*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2)))+((((idx[0]+1)*(L/(N+1)))-(L/2.0))**2)*0.0
    elif idx[0]==idx[1] and idx[1]==N-1 :
      template[idx] = (((30.0-1)*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2)))+((((idx[0]+1)*(L/(N+1)))-(L/2.0))**2)*0.0
    elif idx[0]==idx[1]:
      template[idx] = ((30.0*(N+1)*(N+1)*(1.0/12.0))*(1/(L**2)))+((((idx[0]+1)*(L/(N+1)))-(L/2.0))**2)*0.0

 w, v = la.eig(template)
 E_list.append(min(w))
 #E_list.append(sorted(w)[1])
 #E_list.append(sorted(w)[2])
# print  N,"2", min(w)*3+sorted(w)[1]
 #print  N, min(w)*3+sorted(w)[1]


# print  N,"8", min(w)*6+sorted(w)[1]*6+sorted(w)[2]*4
# print  N,"6", min(w)*6+sorted(w)[1]*4+sorted(w)[2]*2
# print  N,"14", min(w)*8+sorted(w)[1]*8+sorted(w)[2]*7+sorted(w)[3]*5

# print  N, min(w)*6+sorted(w)[1]*4+sorted(w)[2]*2


# print  N,"4", min(w)*4+sorted(w)[1]*4
# print  N,"10", min(w)*8+sorted(w)[1]*6+sorted(w)[2]*4+sorted(w)[3]*2

# print  N, min(w)*8+sorted(w)[1]*6+sorted(w)[2]*4+sorted(w)[3]*2

# print  N, min(w)*6+sorted(w)[1]*2














