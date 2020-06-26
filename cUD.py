import pyUni10 as uni10
import copy
import numpy as np
import scipy as sp
from numpy import linalg as npLA
from scipy import linalg as LA
import MPSclass 
import math
import TruncateU as TU
import time



def Short_TrotterSteps(N_iterF):
 List_delN=[]

 Delta_N=(0.02, N_iterF)
 List_delN.append(Delta_N)


# Delta_N=(0.005, N_iterF)
# List_delN.append(Delta_N)


 #Delta_N=(0.05, N_iterF)
 #List_delN.append(Delta_N)


 return List_delN








# def fermionicOPT(Sys,bdi,bdi1):
# 
#  bdo = uni10.Bond(uni10.BD_OUT, bdi.Qlist())
#  bdo1 = uni10.Bond(uni10.BD_OUT, bdi1.Qlist())
#  bdii = uni10.Bond(uni10.BD_IN, bdi.Qlist())
#  bdii1 = uni10.Bond(uni10.BD_IN, bdi1.Qlist())
# 
# 
#  T = uni10.UniTensor([bdii,bdii1,bdo,bdo1])
#  T.setLabel([1,2,3,4])
#  if Sys[0]=="Fer":
#   template = np.zeros([bdi.dim(), bdi1.dim(), bdo.dim(), bdo1.dim()])
# 
#   for idx in np.ndindex(template.shape):
#    if idx[0]==idx[2] and idx[1]==idx[3]:
#     if bdi1.Qlist()[idx[1]].prt() == uni10.PRT_ODD and bdi.Qlist()[idx[0]].prt() == uni10.PRT_ODD:
#      template[idx] = -1.
#     else:
#      template[idx] = +1.
# 
#   T.setRawElem(template.reshape(-1))
#  elif Sys[0]=="Bos":
#   T.identity()
# 
#  T.permute([1,2,3,4],2)
#  #print T.printDiagram(), "T"
#  return T
# 





def fermionicOPT(Sys,bdi, bdi1):


 bdo = uni10.Bond(uni10.BD_OUT, bdi.Qlist())
 bdo1 = uni10.Bond(uni10.BD_OUT, bdi1.Qlist())
 bdii = uni10.Bond(uni10.BD_IN, bdi.Qlist())
 bdii1 = uni10.Bond(uni10.BD_IN, bdi1.Qlist())


 T = uni10.UniTensor([bdii,bdii1,bdo,bdo1])
 T.setLabel([1,2,3,4])

 if Sys[0]=="Fer":
  template = np.zeros([bdi.dim(), bdi1.dim(), bdo.dim(), bdo1.dim()])

  for idx in np.ndindex(template.shape):
   if idx[0]==idx[2] and idx[1]==idx[3]:
    if Sys[9]=="Z2":
     if bdi1.Qlist()[idx[1]].prt() == uni10.PRT_ODD and bdi.Qlist()[idx[0]].prt() == uni10.PRT_ODD:
      template[idx] = -1.0

     else:
      template[idx] = +1.0

    if Sys[9]=="U1":
     if bdi.Qlist()[idx[0]].U1() % 2 ==1 and  bdi1.Qlist()[idx[1]].U1() % 2 ==1:
       template[idx] = -1.0
     else:
      template[idx] = +1.0

  T.setRawElem(template.reshape(-1))
 elif Sys[0]=="Bos":
  T.identity()


 T.permute([1,2,3,4],2)

 return T









def   Init_Q_list(N_x, d_in, d_out,Sys):


 bdi = uni10.Bond(uni10.BD_IN, d_in)
 bdo = uni10.Bond(uni10.BD_OUT, d_out)

 T0=uni10.UniTensor([bdi, bdi, bdi, bdi, bdo])
 T0.setLabel([1,2,3,4,5])
 if Sys[5]=="rand": T0.randomize()

 Q_list=[None]*N_x
 for i in xrange(N_x):
  Q_list[i]=[None]*N_x

 for i in xrange(N_x):
  for j in xrange(N_x):
   blk_qnums = T0.blockQnum()
   for qnum in blk_qnums:

    T_mat=T0.getBlock(qnum)
    U=T_mat.svd()
    if T_mat.row()>=T_mat.col():
     U_mat=U[0]
     T0.putBlock(qnum,U_mat)
    else:
     U_mat=U[2]
     T0.putBlock(qnum,U_mat)
   if Sys[5]=="iden":T0.identity() 
   Q_list[i][j]=T0*1.0 


 if Sys[5]=="part" and Sys[0]=="Fer": 

  T = uni10.Matrix(bdi.dim()*bdi.dim()*bdi.dim()*bdi.dim(),bdo.dim())
  DO=bdo.dim()
  DI=bdi.dim()
  for i in xrange(bdi.dim()):
   for j in xrange(bdi.dim()):
    for m in xrange(bdi.dim()):
     for n in xrange(bdi.dim()):
      for p in xrange(bdo.dim()):
        if V_2(p,DO)==V_1(i,DI)+V_1(j,DI)+V_1(m,DI)+V_1(n,DI):
         norm=0
         for i1 in xrange(bdi.dim()):
          for j1 in xrange(bdi.dim()):
           for m1 in xrange(bdi.dim()):
            for n1 in xrange(bdi.dim()):
             if V_2(p,DO)==V_1(i1,DI)+V_1(j1,DI)+V_1(m1,DI)+V_1(n1,DI):norm=norm+1.0
         T[i*bdo.dim()*bdi.dim()*bdi.dim()*bdi.dim()+j*bdo.dim()*bdi.dim()*bdi.dim()+m*bdo.dim()*bdi.dim()+n*bdo.dim()+p]=1.0*(1./(norm**0.5))
        else:
         T[i*bdo.dim()*bdi.dim()*bdi.dim()*bdi.dim()+j*bdo.dim()*bdi.dim()*bdi.dim()+m*bdo.dim()*bdi.dim()+n*bdo.dim()+p]=0.0
  #print T
  #T1=T*1.0
  #T1.transpose()
  #print T1*T



  T0=uni10.UniTensor([bdi, bdi, bdi, bdi, bdo]) 
  T0.setRawElem(T)

  #print T0
  T0.setLabel([1,2,3,4,5])
  T1=T0*1
  T1.setLabel([1,2,3,4,-5])
  Q_h=T0*T1
  Q_h.permute([-5,5],1)
  #print Q_h

  for i in xrange(N_x):
   for j in xrange(N_x):
     #T0.identity()
     Q_list[i][j]=T0*1.0 







 if Sys[5]=="part" and Sys[0]=="Bos": 

  T = uni10.Matrix(bdi.dim()*bdi.dim()*bdi.dim()*bdi.dim(),bdo.dim())
  DO=bdo.dim()
  DI=bdi.dim()
  for i in xrange(bdi.dim()):
   for j in xrange(bdi.dim()):
    for m in xrange(bdi.dim()):
     for n in xrange(bdi.dim()):
      for p in xrange(bdo.dim()):
        if V_22(p,DO)==V_11(i,DI)+V_11(j,DI)+V_11(m,DI)+V_11(n,DI):
         norm=0
         for i1 in xrange(bdi.dim()):
          for j1 in xrange(bdi.dim()):
           for m1 in xrange(bdi.dim()):
            for n1 in xrange(bdi.dim()):
             if V_22(p,DO)==V_11(i1,DI)+V_11(j1,DI)+V_11(m1,DI)+V_11(n1,DI):norm=norm+1.0
         T[i*bdo.dim()*bdi.dim()*bdi.dim()*bdi.dim()+j*bdo.dim()*bdi.dim()*bdi.dim()+m*bdo.dim()*bdi.dim()+n*bdo.dim()+p]=1.0*(1./(norm**0.5))
        else:
         T[i*bdo.dim()*bdi.dim()*bdi.dim()*bdi.dim()+j*bdo.dim()*bdi.dim()*bdi.dim()+m*bdo.dim()*bdi.dim()+n*bdo.dim()+p]=0.0
  #print T
  #T1=T*1.0
  #T1.transpose()
  #print T1*T



  T0=uni10.UniTensor([bdi, bdi, bdi, bdi, bdo]) 
  T0.setRawElem(T)

  #print T0
  T0.setLabel([1,2,3,4,5])
  T1=T0*1
  T1.setLabel([1,2,3,4,-5])
  Q_h=T0*T1
  Q_h.permute([-5,5],1)
  #print Q_h

  for i in xrange(N_x):
   for j in xrange(N_x):
     #T0.identity()
     Q_list[i][j]=T0*1.0 



 return  Q_list


def V_2(i,D):

  if i<=((D/2)-1):
    return 2*i
  else:
   ip=i-(D/2)
   return 2*ip+1


def V_1(i,D):
  if i<=((D/2)-1):
    return 2*i
  else:
   ip=i-(D/2)
   return 2*ip+1



def V_22(i,D):
 return i

def V_11(i,D):
 return i





def full_make_bond( Model, D, chi_boundry, chi_single, chi_try, d_in, d_out ):

######################### No-symmetry #############################################
 if  Model[0] is "Heis"  or  Model[0] is "ITF" or Model[0] is "Fer" or Model[0] is "FFI" or Model[0] is "Fer_BOS":
  q0_even = uni10.Qnum(0,uni10.PRT_EVEN)
  q_list=[q0_even]
  qchi_list=[q0_even]

  qchi_boundry_list=[q0_even]
  qchi_single_list=[q0_even]
  qchi_try_list=[q0_even]

  q_d_in=[q0_even]*d_in[0]
  q_d_out=[q0_even]*d_out[0]
  q_list_out=[q0_even]


  q_D, q_chi_boundry, q_chi_single, q_chi_try,q_d_in, q_d_out=make_bond(D, q_list, chi_boundry, chi_single, chi_try, qchi_boundry_list, qchi_single_list, qchi_try_list, d_in, d_out, q_list_out )

 ###################### Z(2) ######################################
 if Model[0] is "Heis_Z2" or Model[0] is "ITF_Z2" or Model[0] is "Fer_Z2" or Model[0] is "FFI_Z2" or Model[0] is "Fer_BOS_Z2":
  q0_even = uni10.Qnum(0,uni10.PRT_EVEN)
  q0_odd = uni10.Qnum(0,uni10.PRT_ODD)
  q_list=[q0_even,q0_odd]
  qchi_boundry_list=[q0_even]
  qchi_single_list=[q0_even]
  qchi_try_list=[q0_even]
  q_d_in=[q0_even,q0_odd]
  q_d_out=[q0_even,q0_odd]
  q_list_out=[q0_even,q0_odd]

  q_D, q_chi_boundry, q_chi_single, q_chi_try,q_d_in, q_d_out=make_bond(D, q_list, chi_boundry, chi_single, chi_try, qchi_boundry_list, qchi_single_list, qchi_try_list, d_in, d_out, q_list_out)


######################################
#  c_val=0
#  for i in xrange(len(d_out)):
#   for q in xrange(d_out[i]):
#    c_val=c_val+1

#  q_d_out=[]
#  for i in xrange(c_val):
#   if i%2==0:
#    q_d_out.append(q0_even)
#   else:
#    q_d_out.append(q0_odd)

#######################################
 if Model[0] is "Fer_U1": 

  q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
  q0_odd = uni10.Qnum(0,uni10.PRT_ODD);
  q1_even = uni10.Qnum(1,uni10.PRT_EVEN)
  q2_even = uni10.Qnum(2,uni10.PRT_EVEN)

  q2_odd = uni10.Qnum(2,uni10.PRT_ODD)
  q1_odd = uni10.Qnum(1,uni10.PRT_ODD)

  q_1_odd = uni10.Qnum(-1,uni10.PRT_ODD)
  q_2_odd = uni10.Qnum(-2,uni10.PRT_ODD)

  q2_even = uni10.Qnum(2,uni10.PRT_EVEN);
  q3_even = uni10.Qnum(3,uni10.PRT_EVEN);
  q4_even = uni10.Qnum(4,uni10.PRT_EVEN);
  q5_even = uni10.Qnum(5,uni10.PRT_EVEN);
  q6_even = uni10.Qnum(6,uni10.PRT_EVEN);
  q7_even = uni10.Qnum(7,uni10.PRT_EVEN);

  q_1_even = uni10.Qnum(-1,uni10.PRT_EVEN);
  q_2_even = uni10.Qnum(-2,uni10.PRT_EVEN);
  q_3_even = uni10.Qnum(-3,uni10.PRT_EVEN);
  q_4_even = uni10.Qnum(-4,uni10.PRT_EVEN);
  q_5_even = uni10.Qnum(-5,uni10.PRT_EVEN);
  q_6_even = uni10.Qnum(-6,uni10.PRT_EVEN);
  q_7_even = uni10.Qnum(-7,uni10.PRT_EVEN);

  #qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
  #qchi_list=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
  #qchi_list=[q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even]
  qchi_list=[q0_even]
  qchi_boundry_list=[q0_even]
  qchi_single_list=[q0_even]
  qchi_try_list=[q0_even]

  #qchi_list=[q_4_even,q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even,q4_even]
  #qchi_list=[q_5_even,q_4_even,q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even,q4_even,q5_even]
  #qchi_list=[q_1_even,q1_even]

  #q_list=[q_2_even,q_1_even,q0_odd,q1_even,q2_even]
  q_list=[q_1_even,q0_even,q1_even]
  #q_list=[q_2_odd,q_1_odd,q_1_even,q0_odd,q1_even,q1_odd,q2_odd]
  #q_list=[q_1_even,q0_even,q1_even]
  #q_list=[q0_even,q1_even,q2_even]
  #q_list=[q1_even,q3_even,q2_even]
  #q_list=[q_2_even,q_1_even,q0_even,q1_even, q2_even]
  #q_list=[q_3_even,q_2_even,q_1_even,q0_even,q1_even, q2_even,q3_even]
  #q_list=[q_1_even,q1_even]

  #q_phys=[q_1_even,q1_even]
  q_phys=[q0_even,q0_even,q1_even,q_1_even]
  q_d_in=q_phys
  q_list_out=[q_3_even,q_2_even,q_1_even,q0_even,q1_even,q2_even,q3_even]
#  q_list_out=[q_2_even,q_1_even,q0_even,q1_even,q2_even]
#  q_list_out=[q_1_even,q0_even,q1_even]

  q_D, q_chi_boundry, q_chi_single, q_chi_try,q_d_in_try, q_d_out=make_bond(D, q_list, chi_boundry, chi_single, chi_try, qchi_boundry_list, qchi_single_list, qchi_try_list, d_in, d_out,q_list_out  )





  #q_d_out=q_phys


 return q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d_in, q_d_out


def make_bond(D, q_list, chi_boundry, chi_single, chi_try, qchi_boundry_list, qchi_single_list, qchi_try_list, d_in, d_out, q_list_out):
 q_D=[]
 q_chi_boundry=[]
 q_chi_single=[]
 q_chi_try=[]
 q_d_in=[]
 q_d_out=[]

 for i in xrange(len(D)):
  for q in xrange(D[i]):
   q_D.append(q_list[i])


 for i in xrange(len(d_in)):
  for q in xrange(d_in[i]):
   q_d_in.append(q_list[i])

 for i in xrange(len(d_out)):
  for q in xrange(d_out[i]):
   q_d_out.append(q_list_out[i])



 for i in xrange(len(chi_boundry)):
  for q in xrange(chi_boundry[i]):
   q_chi_boundry.append(qchi_boundry_list[i])


 for i in xrange(len(chi_single)):
  for q in xrange(chi_single[i]):
   q_chi_single.append(qchi_single_list[i])

 for i in xrange(len(chi_try)):
  for q in xrange(chi_try[i]):
   q_chi_try.append(qchi_try_list[i])


 return q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d_in, q_d_out







##@profile
def  Ascend_f_col(Q_list, H_list, N_x,Sys,H_long):

 h_long=H_long[0][0]
 HA_list=[None]*N_x
 for i in xrange(N_x):
  HA_list[i]=[None]*(N_x-1)

 Swap=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1))

 for i in xrange(N_x):
  for j in xrange(N_x):

   if j==N_x-1:
    Q_list[i][j].setLabel([1,2,3,4,5])
    Swap.setLabel([-2,-3,2,3])

    Q_list_ferm=Swap*Q_list[i][j]
    Q_list_ferm.permute([1,-2,-3,4,5],5)

    Q_list_ferm.setLabel([1,2,3,4,5])
    Q_trans=Q_list_ferm*1.0
    Q_trans.transpose()
    Q_trans.setLabel([-1,2,-3,4,-5])


    Q_tem=Q_list[i][j]*1.0
    Q_tem.setLabel([6,7,8,9,10])
    Q_temN=Q_tem*1.0
    Q_temN.transpose()
    Q_temN.setLabel([-10,6,7,8,9])

    H_list[2*i][2*j].setLabel([-1,-3,1,3])
    result=((Q_list_ferm*H_list[2*i][2*j])*Q_trans)*(Q_tem*Q_temN)
    result.permute([-10,-5,10,5],2)
    HA_list[i][j-1]=result+HA_list[i][j-1]

################################################

    Q_list[i][j].setLabel([1,2,3,4,5])
    Swap.setLabel([-3,-4,3,4])
    Q_list_ferm=Swap*Q_list[i][j]
    Q_list_ferm.permute([1,2,-3,-4,5],5)

    Q_list_ferm.setLabel([1,2,3,4,5])
    Q_trans=Q_list_ferm*1.0
    Q_trans.transpose()
    Q_trans.setLabel([1,-2,3,-4,-5])


    Q_tem=Q_list[i][j]*1.0
    Q_tem.setLabel([6,7,8,9,10])
    Q_temN=Q_tem*1.0
    Q_temN.transpose()
    Q_temN.setLabel([-10,6,7,8,9])

    H_list[2*i+1][2*j].setLabel([-2,-4,2,4])

    result=((Q_list_ferm*H_list[2*i+1][2*j])*Q_trans)*(Q_tem*Q_temN)
    result.permute([-10,-5,10,5],2)
    HA_list[i][j-1]=result+HA_list[i][j-1]
   else:

    Q_list[i][j].setLabel([1,2,3,4,5])
    Swap.setLabel([-2,-3,2,3])
    Q_list_ferm=Swap*Q_list[i][j]
    Q_list_ferm.permute([1,-2,-3,4,5],5)

    Q_list_ferm.setLabel([1,2,3,4,5])
    Q_trans=Q_list_ferm*1.0
    Q_trans.transpose()
    Q_trans.setLabel([-1,2,-3,4,-5])

    Q_list[i][j+1].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i][j+1]*1.0
    Q_transN.transpose()
    Q_transN.setLabel([-10,6,7,8,9])

    H_list[2*i][2*j].setLabel([-1,-3,1,3])

    result=((Q_list_ferm*H_list[2*i][2*j])*Q_trans)*(Q_list[i][j+1]*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result*1.0

################################################

    Q_list[i][j].setLabel([1,2,3,4,5])
    Swap.setLabel([-3,-4,3,4])
    Q_list_ferm=Swap*Q_list[i][j]
    Q_list_ferm.permute([1,2,-3,-4,5],5)

    Q_list_ferm.setLabel([1,2,3,4,5])
    Q_trans=Q_list_ferm*1.0
    Q_trans.transpose()
    Q_trans.setLabel([1,-2,3,-4,-5])

    Q_list[i][j+1].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i][j+1]*1.0
    Q_transN.transpose()
    Q_transN.setLabel([-10,6,7,8,9])

    H_list[2*i+1][2*j].setLabel([-2,-4,2,4])

    result=((Q_list_ferm*H_list[2*i+1][2*j])*Q_trans)*(Q_list[i][j+1]*Q_transN)
    result.permute([-5,-10,5,10],2)

    HA_list[i][j]=result+HA_list[i][j]

#####################################################################################

    Swap.setLabel([12,4,6,-4])

    Q_list[i][j].setLabel([1,2,3,-4,5])
    Q_list_fermin=Q_list[i][j]*Swap
    Q_list_fermin.permute([1,2,3,4,5,12,6],7)
    Q_trans=Q_list_fermin*1.0
    Q_trans.transpose()
    Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

    Q_list[i][j+1].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i][j+1]*1.0
    Q_transN.transpose()
    Q_transN.setLabel([-10,-6,7,8,9])

    H_list[2*i][2*j+1].setLabel([-3,-12,3,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*H_list[2*i][2*j+1])*Q_trans)*(Q_list[i][j+1]*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]

    Swap.setLabel([7,6,-7,-6])

    Q_list[i][j].setLabel([1,2,3,4,5])
    Q_trans=Q_list[i][j]*1.0
    Q_trans.transpose()
    Q_trans.setLabel([-5,1,2,3,-4])
    
    
    Q_list[i][j+1].setLabel([-6,-7,8,9,10])
    Q_list_ferminP=Q_list[i][j+1]*Swap
    Q_list_ferminP.permute([6,7,8,9,10],5)
    Q_transN=Q_list_ferminP*1.0
    Q_transN.transpose()
    Q_transN.setLabel([6,-7,8,9,-10])


    H_list[2*i+1][2*j+1].setLabel([-4,-7,4,7])
    result=((Q_list[i][j]*H_list[2*i+1][2*j+1])*Q_trans)*(Q_list_ferminP*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]
    #print 2*i+1, 2*j+1, H_list[2*i+1][2*j+1]


 return HA_list



def  Ascend_f_row(Q_list, H_list, N_x,Sys,H_long):

 h_long=H_long[0][0]

 Swap=fermionicOPT( Sys, Q_list[0][0].bond(1), Q_list[0][0].bond(2) )
 Swap1=fermionicOPT( Sys, Q_list[0][0].bond(1), Q_list[0][0].bond(2) )

 HA_list=[None]*(N_x-1)
 for i in xrange(N_x-1):
  HA_list[i]=[None]*(N_x)

 for i in xrange(N_x):
  for j in xrange(N_x):

   if i==N_x-1:
    Q_list[i][j].setLabel([1,2,3,4,5])
    Q_trans=Q_list[i][j]*1.0
    Q_trans.transpose()
    Q_trans.setLabel([-5,1,2,-3,-4])

    Q_tem=Q_list[i][j]*1.0
    Q_tem.setLabel([6,7,8,9,10])
    Q_temN=Q_tem*1.0
    Q_temN.transpose()
    Q_temN.setLabel([-10,6,7,8,9])

    H_list[2*i][2*j+1].setLabel([-3,-4,3,4])
    result=((Q_list[i][j]*H_list[2*i][2*j+1])*Q_trans)*(Q_tem*Q_temN)
    result.permute([-10,-5,10,5],2)
    HA_list[i-1][j]=result+HA_list[i-1][j]
 ################################################

    Q_list[i][j].setLabel([1,2,3,4,5])
    Q_trans=Q_list[i][j]*1.0
    Q_trans.transpose()
    Q_trans.setLabel([-5,-1,-2,3,4])

    Q_tem=Q_list[i][j]*1.0
    Q_tem.setLabel([6,7,8,9,10])
    Q_temN=Q_tem*1.0
    Q_temN.transpose()
    Q_temN.setLabel([-10,6,7,8,9])

    H_list[2*i][2*j].setLabel([-1,-2,1,2])
    result=((Q_list[i][j]*H_list[2*i][2*j])*Q_trans)*(Q_tem*Q_temN)
    result.permute([-10,-5,10,5],2)
    HA_list[i-1][j]=result+HA_list[i-1][j]
   else:
    #print "i, j", 2*i, 2*j+1, H_list[2*i][2*j+1]
    Q_list[i][j].setLabel([1,2,3,4,5])
    Q_trans=Q_list[i][j]*1.0
    Q_trans.transpose()
    Q_trans.setLabel([-5,1,2,-3,-4])

    Q_list[i+1][j].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i+1][j]*1.0
    Q_transN.transpose()
    Q_transN.setLabel([-10,6,7,8,9])

    H_list[2*i][2*j+1].setLabel([-3,-4,3,4])
    result=((Q_list[i][j]*H_list[2*i][2*j+1])*Q_trans)*(Q_list[i+1][j]*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result*1.0
################################################

    Q_list[i][j].setLabel([1,2,3,4,5])
    Q_trans=Q_list[i][j]*1.0
    Q_trans.transpose()
    Q_trans.setLabel([-5,-1,-2,3,4])

    Q_list[i+1][j].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i+1][j]*1.0
    Q_transN.transpose()
    Q_transN.setLabel([-10,6,7,8,9])

    H_list[2*i][2*j].setLabel([-1,-2,1,2])

    result=((Q_list[i][j]*H_list[2*i][2*j])*Q_trans)*(Q_list[i+1][j]*Q_transN)
    result.permute([-5,-10,5,10],2)

    HA_list[i][j]=result+HA_list[i][j]

##################################################
    Swap.setLabel([15,6,16,-6])
    Swap1.setLabel([16,7,8,-7])

    Q_list[i][j].setLabel([1,2,3,4,5])
    Q_trans=Q_list[i][j]*1.0
    Q_trans.transpose()
    Q_trans.setLabel([-5,1,2,3,-4])

    Q_list[i+1][j].setLabel([-6,-7,8,9,10])
    Q_list_ferminP=Q_list[i+1][j]*(Swap1*Swap)
    Q_list_ferminP.permute([6,7,9,10,15],5)
    Q_transN=Q_list_ferminP*1.0
    Q_transN.transpose()
    Q_transN.setLabel([6,7,9,-10,-15])


    H_list[2*i+1][2*j+1].setLabel([-4,-15,4,15])
    result=((Q_list[i][j]*H_list[2*i+1][2*j+1])*Q_trans)*(Q_list_ferminP*Q_transN)
    result.permute([-5,-10,5,10],2)

    HA_list[i][j]=result+HA_list[i][j]



    Swap.setLabel([12,3,13,-3])
    Swap1.setLabel([13,4,6,-4])

    Q_list[i][j].setLabel([1,2,-3,-4,5])
    Q_list_fermin=Q_list[i][j]*(Swap1*Swap)
    Q_list_fermin.permute([1,2,3,4,5,12,6],7)
    Q_trans=Q_list_fermin*1.0
    Q_trans.transpose()
    Q_trans.setLabel([1,-2,3,4,-5,-12,-6])


    Q_list[i+1][j].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i+1][j]*1.0
    Q_transN.transpose()
    Q_transN.setLabel([-10,-6,7,8,9])

    H_list[2*i+1][2*j].setLabel([-2,-12,2,12])
    result=((Q_list_fermin*H_list[2*i+1][2*j])*Q_trans)*(Q_list[i+1][j]*Q_transN)
    result.permute([-5,-10,5,10],2)

    HA_list[i][j]=result+HA_list[i][j]

 return HA_list




def  Q_cost_val_middle( H_col, H_row, Q_list, i, j , rho_row, rho_col, N_x, Sys,HH_long):


 h_long=HH_long[0][0]
 result=0
 Swap=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )
 Swap1=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )

############################   col   #####################################
 if j==N_x-1:
  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-2,-3,2,3])

  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,-2,-3,4,5],5)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-1,2,-3,4,-5])


  Q_list[i][j-1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j-1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])


  H_col[2*i][2*j].setLabel([-1,-3,1,3])
  rho_col[i][j-1].setLabel([-10,-5,10,5])

  result1=(((Q_list_ferm*H_col[2*i][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)
  result=result1[0]+result
  #print "12",result1
################################################

  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-3,-4,3,4])
  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,2,-3,-4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,-2,3,-4])


  Q_list[i][j-1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j-1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])

  H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
  rho_col[i][j-1].setLabel([-10,-5,10,5])
  result1=(((Q_list_ferm*H_col[2*i+1][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)

  result=result1[0]+result
  #print "11",result1

#######################################################


  Swap.setLabel([12,4,6,-4])

  Q_list[i][j-1].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j-1]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])


  rho_col[i][j-1].setLabel([-5,-10,5,10])
  H_col[2*i][2*j-1].setLabel([-3,-12,3,12])

  result1=((Q_list[i][j]*Q_transN)*rho_col[i][j-1])*(((Q_list_fermin*H_col[2*i][2*j-1])*Q_trans))
  result=result1[0]+result
  #print "10",result1

  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j-1]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],4)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,-7,8,9])

  H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=((Q_list_ferminP*H_col[2*i+1][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_trans)
  result=result1[0]+result



 if j<N_x-1:

  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-2,-3,2,3])

  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,-2,-3,4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,-1,2,-3,4])

  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])
 
  H_col[2*i][2*j].setLabel([-1,-3,1,3])
  rho_col[i][j].setLabel([-5,-10,5,10])
  result1=((Q_list_ferm*H_col[2*i][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]
  result=result1[0]+result
  #print "8", result1[0]

################################################


  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-3,-4,3,4])
  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,2,-3,-4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,-2,3,-4])


  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])

  H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
  rho_col[i][j].setLabel([-5,-10,5,10])
  result1=((Q_list_ferm*H_col[2*i+1][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]
  result=result1[0]+result
  #print "7", result1[0]
#######################################################

  Swap.setLabel([12,4,6,-4])

  Q_list[i][j].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])

  rho_col[i][j].setLabel([-5,-10,5,10])
  #print rho_col[i][j].printDiagram()

  H_col[2*i][2*j+1].setLabel([-3,-12,3,12])
  result1=(Q_list_fermin*(Q_trans*((Q_list[i][j+1]*Q_transN)*rho_col[i][j])))*H_col[2*i][2*j+1]
  result=result1[0]+result
  #print "6", result1


  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j+1].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j+1]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,-7,8,9,-10])


  H_col[2*i+1][2*j+1].setLabel([-4,-7,4,7])
  rho_col[i][j].setLabel([-5,-10,5,10])
  result1=(((Q_list[i][j]*Q_trans)*H_col[2*i+1][2*j+1])*rho_col[i][j])*(Q_list_ferminP*Q_transN)
  result=result1[0]+result
  #print "5",result1

###############################################################################
 if j-1>=0 and j!=N_x-1:

  Swap.setLabel([12,4,6,-4])

  Q_list[i][j-1].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j-1]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])


  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])


  rho_col[i][j-1].setLabel([-5,-10,5,10])

  H_col[2*i][2*j-1].setLabel([-3,-12,3,12])
  result1=((Q_list[i][j]*H_col[2*i][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list_fermin*Q_trans)
  result=result1[0]+result
  #print "outb0", result1[0]
  #print "4",result1

  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j-1]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,-7,8,9,-10])

  H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=((Q_list_ferminP*H_col[2*i+1][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_trans)
  result=result1[0]+result
  #print "3", result1





##################### Row ####################################
# 
# 
 if i==N_x-1:

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,4])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,-8,-9])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  H_row[2*i][2*j+1].setLabel([-8,-9,8,9])
  result1=((Q_list[i][j]*H_row[2*i][2*j+1]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
  result=result1[0]+result
  #print "1", result1


  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,4])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,-7,8,9])

  H_row[2*i][2*j].setLabel([-6,-7,6,7])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=((Q_list[i][j]*H_row[2*i][2*j]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
  result=result1[0]+result
  #print "2",result1





  Swap.setLabel([15,6,16,-6])
  Swap1.setLabel([16,7,8,-7])

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
  Q_list_ferminP.permute([6,7,9,10,15],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,7,9,-10,-15])

  rho_row[i-1][j].setLabel([-5,-10,5,10])

  H_row[2*i-1][2*j+1].setLabel([-4,-15,4,15])
  result1=((Q_list_ferminP*H_row[2*i-1][2*j+1]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
  result=result1[0]+result
  #print result1
  #print "16",result1

  Swap.setLabel([12,3,13,-3])
  Swap1.setLabel([13,4,6,-4])

  Q_list[i-1][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=(Q_list[i-1][j]*Swap)*Swap1
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,-2,3,4,-5,-12,-6])


  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])


  H_row[2*i-1][2*j].setLabel([-2,-12,2,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=((Q_list[i][j]*Q_transN)*rho_row[i-1][j])*((Q_list_fermin*Q_trans)*H_row[2*i-1][2*j])

  result=result1[0]+result
  #print "17",result1





 if i<N_x-1:
 #######################################################
  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,-3,-4])

  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])

  rho_row[i][j].setLabel([-5,-10,5,10])
  H_row[2*i][2*j+1].setLabel([-3,-4,3,4])
  result1=(Q_list[i][j]*H_row[2*i][2*j+1]*Q_trans)*((Q_list[i+1][j]*Q_transN)*rho_row[i][j])
  result=result1[0]+result
  #print "19", result1

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,-1,-2,3,4])

  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])

  H_row[2*i][2*j].setLabel([-1,-2,1,2])
  rho_row[i][j].setLabel([-5,-10,5,10])

  result1=((Q_list[i][j]*H_row[2*i][2*j]*Q_trans)*rho_row[i][j])*(Q_list[i+1][j]*Q_transN)

  result=result1[0]+result
  #print "18", result1

 #######################################################

  Swap.setLabel([15,6,16,-6])
  Swap1.setLabel([16,7,8,-7])

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i+1][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=(Q_list[i+1][j]*Swap)*Swap1
  Q_list_ferminP.permute([6,7,9,10,15],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,7,9,-10,-15])

  rho_row[i][j].setLabel([-5,-10,5,10])
  H_row[2*i+1][2*j+1].setLabel([-4,-15,4,15])
  result1=(Q_list[i][j]*H_row[2*i+1][2*j+1]*Q_trans)*((Q_list_ferminP*Q_transN)*rho_row[i][j])
  result=result1[0]+result
  #print "20",result1


  Swap.setLabel([12,3,13,-3])
  Swap1.setLabel([13,4,6,-4])
  Swapp=Swap*1.0
  Swapp1=Swap1*1.0
  Swapp.setLabel([-12,3,130,-33])
  Swapp1.setLabel([130,4,-6,-44])



  Q_list[i][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=(Q_list[i][j])*(Swap1*Swap)
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,-2,3,4,-5,-12,-6])

  Q_transS=Q_list[i][j]*1.0
  Q_transS.setLabel([1,-2,-33,-44,-5])



  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.transpose()

  Q_transN.setLabel([-10,-6,7,8,9])

  H_row[2*i+1][2*j].setLabel([-2,-12,2,12])
  rho_row[i][j].setLabel([-5,-10,5,10])

  result1=H_row[2*i+1][2*j]*((Q_trans*((Q_list[i+1][j]*Q_transN)*rho_row[i][j]))*Q_list_fermin)


  result=result1[0]+result
  #print "21",result1


 ######################################################################################

 if i-1>=0 and i!=N_x-1:

  Swap.setLabel([15,6,16,-6])
  Swap1.setLabel([16,7,8,-7])

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
  Q_list_ferminP.permute([6,7,9,10,15],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,7,9,-10,-15])


  rho_row[i-1][j].setLabel([-5,-10,5,10])

  H_row[2*i-1][2*j+1].setLabel([-4,-15,4,15])
  result1=((Q_list_ferminP*H_row[2*i-1][2*j+1]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)

  result=result1[0]+result
  #print "22",result1


  Swap.setLabel([12,3,13,-3])
  Swap1.setLabel([13,4,6,-4])

  Q_list[i-1][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=(Q_list[i-1][j]*Swap)*Swap1
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,-2,3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])


  H_row[2*i-1][2*j].setLabel([-2,-12,2,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=((Q_list[i][j]*Q_transN)*rho_row[i-1][j])*((Q_list_fermin*Q_trans)*H_row[2*i-1][2*j])
  result=result1[0]+result
  #print "23",result1



 return result


def  Grad_Q_middle( H_col, H_row, Q_list, i, j , rho_row, rho_col, N_x, Sys,HH_long):


 h_long=HH_long[0][0]

 Q_list[i][j].setLabel([0,1,2,3,4])
 Q_list[i][j].permute([0,1,2,3,4],4)
 Q_t=Q_list[i][j]*1.0
 Q_t.transpose()
 Q_t.setLabel([4,0,1,2,3])
 Q_t.permute([0,1,2,3,4],4)
 result=Q_t*0.0
 
 Swap=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )
 Swap1=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )

 ###########################   col   #####################################
 if j==N_x-1:

  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-2,-3,2,3])

  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,-2,-3,4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,-1,2,-3,4])

  Q_list[i][j-1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j-1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])


  H_col[2*i][2*j].setLabel([-1,-3,1,3])
  rho_col[i][j-1].setLabel([-10,-5,10,5])

  result1=(((H_col[2*i][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)
  result1.permute([1,2,3,4,5],4)

  Swap_t=fermionicOPT( Sys, result1.bond(1), result1.bond(2) )

  Swap_t.setLabel([-2,-3,2,3])
  result1=Swap_t*result1
  result1.permute([1,-2,-3,4,5],4)
  result1.setLabel([1,2,3,4,5])

  result=result1+result

 ################################################


  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-3,-4,3,4])
  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,2,-3,-4,5],5)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,-2,3,-4,-5])


  Q_list[i][j-1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j-1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])

  H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
  rho_col[i][j-1].setLabel([-10,-5,10,5])
  result1=(((H_col[2*i+1][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)
  result1.permute([1,2,3,4,5],4)

  Swap_t=fermionicOPT( Sys, result1.bond(2), result1.bond(3) )

  Swap_t.setLabel([-3,-4,3,4])
  result1=Swap_t*result1
  result1.permute([1,2,-3,-4,5],4)
  result1.setLabel([1,2,3,4,5])

  result=result1+result


#################################################################################
  Swap.setLabel([12,4,6,-4])

  SwapP=Swap*1.0
  SwapP.setLabel([-12,4,-6,-40])


  Q_list[i][j-1].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j-1]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

  Q_transS=Q_list[i][j-1]*1.0
  Q_transS.transpose()
  Q_transS.setLabel([-5,1,2,-3,-40])
  Q_transS.permute([1,2,-3,-40,-5],4)

  SwapP=fermionicOPT( Sys, Q_transS.bond(2), Q_transS.bond(3) )
  SwapP.setLabel([-12,4,-6,-40])


  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])

  rho_col[i][j-1].setLabel([-5,-10,5,10])
  H_col[2*i][2*j-1].setLabel([-3,-12,3,12])

  result1=((((Q_list[i][j-1]*Swap)*H_col[2*i][2*j-1])*(Q_transS*SwapP))*rho_col[i][j-1])*Q_transN

  result1.permute([6,7,8,9,10],4)

  result=result1+result

##################

  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j-1]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,-7,8,9,-10])


  H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])

  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=((H_col[2*i+1][2*j-1])*((Q_list[i][j-1]*Q_trans)*rho_col[i][j-1]))*Q_transN
  result1.permute([6,7,8,9,10],4)

  Swap_t=fermionicOPT( Sys, result1.bond(1), result1.bond(0) )
  Swap_t.setLabel([-12,4,-6,-40])


  Swap_t.setLabel([-7,-6,7,6])
  result1=result1*Swap_t
  result1.permute([-6,-7,8,9,10],4)

  result=result1+result



 if j<N_x-1:


  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-2,-3,2,3])

  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,-2,-3,4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,-1,2,-3,4])


  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_list[i][j+1].permute([6,7,8,9,10],4)
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])

  H_col[2*i][2*j].setLabel([-1,-3,1,3])
  rho_col[i][j].setLabel([-5,-10,5,10])

  result1=((H_col[2*i][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]
  result1.permute([1,2,3,4,5],4)

  Swap_t=fermionicOPT(Sys,result1.bond(1), result1.bond(2))

  Swap_t.setLabel([-2,-3,2,3])
  result1=Swap_t*result1
  result1.permute([1,-2,-3,4,5],4)
  result1.setLabel([1,2,3,4,5])

  result=result1+result

######################################################################

  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-3,-4,3,4])
  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,2,-3,-4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,-2,3,-4])

  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])

  H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
  rho_col[i][j].setLabel([-5,-10,5,10])
  result1=((H_col[2*i+1][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]
  result1.permute([1,2,3,4,5],4)
  Swap_t=fermionicOPT(Sys,result1.bond(2), result1.bond(3))
  Swap_t.setLabel([-3,-4,3,4])
  result1=Swap_t*result1
  result1.permute([1,2,-3,-4,5],4)
  result1.setLabel([1,2,3,4,5])
  result1.permute([1,2,3,4,5],4)
  result=result1+result

#####################################################################


  Swap.setLabel([12,4,6,-4])

  Q_list[i][j].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])

  rho_col[i][j].setLabel([-5,-10,5,10])

  H_col[2*i][2*j+1].setLabel([-3,-12,3,12])

  result1=(H_col[2*i][2*j+1])*(((Q_list[i][j+1]*Q_transN)*rho_col[i][j])*Q_trans)
  result1.permute([1,2,3,4,5,12,6],7)

  Swap_t=fermionicOPT(Sys,result1.bond(5), result1.bond(3))
  Swap_t.setLabel([6,-4,12,4])

  result1=result1*Swap_t
  result1.permute([1,2,3,-4,5],4)
  result=result1+result


  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j+1].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j+1]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,-7,8,9,-10])

  H_col[2*i+1][2*j+1].setLabel([-4,-7,4,7])

  rho_col[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_ferminP*Q_transN)*rho_col[i][j])*H_col[2*i+1][2*j+1])*Q_trans
  result1.permute([1,2,3,4,5],4)

  result=result1+result


##############################################################################
 if j-1>=0 and j!=N_x-1:

  Swap.setLabel([12,4,6,-4])

  Q_list[i][j-1].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j-1]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])


  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])


  rho_col[i][j-1].setLabel([-5,-10,5,10])


  H_col[2*i][2*j-1].setLabel([-3,-12,3,12])
  result1=(((Q_list_fermin*Q_trans)*rho_col[i][j-1])*H_col[2*i][2*j-1])*Q_transN

  result1.permute([6,7,8,9,10],4)
  result=result1+result


  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j-1]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,-7,8,9,-10])

  H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=(((Q_list[i][j-1]*Q_trans)*rho_col[i][j-1])*H_col[2*i+1][2*j-1])*Q_transN
  result1.permute([6,7,8,9,10],4)

  Swap_t=fermionicOPT(Sys,result1.bond(1), result1.bond(0))
  Swap_t.setLabel([-7,-6,7,6])

  result1=result1*Swap_t
  result1.permute([-6,-7,8,9,10],4)

  result=result1+result
####################################################################################

################ Row ####################################

 if i==N_x-1:

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,4])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,-8,-9])

  rho_row[i-1][j].setLabel([-5,-10,5,10])

  H_row[2*i][2*j+1].setLabel([-8,-9,8,9])
  result1=((H_row[2*i][2*j+1]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
  result1.permute([6,7,8,9,10],4)  
  result=result1+result

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,4])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,-7,8,9])

  H_row[2*i][2*j].setLabel([-6,-7,6,7])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=((H_row[2*i][2*j]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
  result1.permute([6,7,8,9,10],4)  
  result=result1+result


  Swap.setLabel([15,6,16,-6])
  Swap1.setLabel([16,7,8,-7])

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
  Q_list_ferminP.permute([6,7,9,10,15],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,7,9,-10,-15])


  rho_row[i-1][j].setLabel([-5,-10,5,10])

  H_row[2*i-1][2*j+1].setLabel([-4,-15,4,15])

  result1=(((Q_list[i-1][j]*Q_trans)*rho_row[i-1][j])*H_row[2*i-1][2*j+1])*Q_transN
  result1.permute([6,7,9,10,15],4)    

  Swap_t=fermionicOPT(Sys,result1.bond(4), result1.bond(0))
  Swap_t1=fermionicOPT(Sys,result1.bond(4), result1.bond(1))

  Swap_t.setLabel([15,-6,16,6])
  Swap_t1.setLabel([16,-7,8,7])


  result1=(result1*Swap_t)*Swap_t1
  result1.permute([-6,-7,8,9,10],4)

  result=result1+result
#############################

  Swap.setLabel([12,3,13,-3])
  Swap1.setLabel([13,4,6,-4])
  Swapp=Swap*1.0
  Swapp1=Swap1*1.0
  Swapp.setLabel([-12,-33,-13,3])
  Swapp1.setLabel([-13,4,-6,-44])


  Q_list[i-1][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=(Q_list[i-1][j]*Swap)*Swap1
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,-2,3,4,-5,-12,-6])

  Q_list[i-1][j].setLabel([1,2,-3,-4,5])
  Q_transS=Q_list[i-1][j]*1.0
  Q_transS.transpose()
  Q_transS.setLabel([-5,1,-2,-33,-44])
  Q_transS.permute([1,-2,-33,-44,-5],4)

  Swapp=fermionicOPT(Sys,Q_trans.bond(5), Q_transS.bond(2))
  Swapp1=fermionicOPT(Sys,Q_transS.bond(1), Q_transS.bond(3))

  Swapp.setLabel([-13,3,-12,-33])
  Swapp1.setLabel([-13,4,-6,-44])


  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])

  H_row[2*i-1][2*j].setLabel([-2,-12,2,12])

  rho_row[i-1][j].setLabel([-5,-10,5,10])

  #result1=(((Q_list_fermin*Q_trans)*rho_row[i-1][j])*H_row[2*i-1][2*j])*Q_transN




  result1=(((((Q_list[i-1][j]*(Swap*Swap1))*H_row[2*i-1][2*j])*((Swapp*Swapp1)))*Q_transS)*rho_row[i-1][j])*Q_transN

  result1.permute([6,7,8,9,10],4)

  result=result1+result



 if i<N_x-1:
 ##################################################
  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,-3,-4])

  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])

  rho_row[i][j].setLabel([-5,-10,5,10])
  H_row[2*i][2*j+1].setLabel([-3,-4,3,4])
  result1=(H_row[2*i][2*j+1]*Q_trans)*((Q_list[i+1][j]*Q_transN)*rho_row[i][j])
  result1.permute([1,2,3,4,5],4)
  result=result1+result

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,-1,-2,3,4])

  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,8,9])

  H_row[2*i][2*j].setLabel([-1,-2,1,2])
  rho_row[i][j].setLabel([-5,-10,5,10])

  result1=(H_row[2*i][2*j])*(((Q_list[i+1][j]*Q_transN)*rho_row[i][j])*Q_trans)
  result1.permute([1,2,3,4,5],4)    
  result=result1+result

 ##################################################
  Swap.setLabel([15,6,16,-6])
  Swap1.setLabel([16,7,8,-7])

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i+1][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i+1][j]*(Swap1*Swap)
  Q_list_ferminP.permute([6,7,9,10,15],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,7,9,-10,-15])



  rho_row[i][j].setLabel([-5,-10,5,10])
  H_row[2*i+1][2*j+1].setLabel([-4,-15,4,15])

  result1=(H_row[2*i+1][2*j+1]*Q_trans)*((Q_list_ferminP*Q_transN)*rho_row[i][j])
  result1.permute([1,2,3,4,5],4)
  result=result1+result


  Swap.setLabel([12,3,13,-3])
  Swap1.setLabel([13,4,6,-4])

  Q_list[i][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=Q_list[i][j]*(Swap1*Swap)
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,-2,3,4,-5,-12,-6])


  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])

  H_row[2*i+1][2*j].setLabel([-2,-12,2,12])

  rho_row[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list[i+1][j]*Q_transN)*rho_row[i][j])*Q_trans)*H_row[2*i+1][2*j]
  result1.permute([1,2,3,4,5,12,6],7)

  Swap_t=fermionicOPT(Sys,result1.bond(5), result1.bond(2))
  Swap_t1=fermionicOPT(Sys,result1.bond(6), result1.bond(3))

  Swap_t.setLabel([13,-3,12,3])
  Swap_t1.setLabel([13,-4,6,4])

  result1=result1*(Swap_t1*Swap_t)
  result1.permute([1,2,-3,-4,5],4)

  result=result1+result
 ################################################################################




 if i-1>=0 and i!=N_x-1:

  Swap.setLabel([15,6,16,-6])
  Swap1.setLabel([16,7,8,-7])

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,-4])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
  Q_list_ferminP.permute([6,7,9,10,15],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.transpose()
  Q_transN.setLabel([6,7,9,-10,-15])

  rho_row[i-1][j].setLabel([-5,-10,5,10])
  H_row[2*i-1][2*j+1].setLabel([-4,-15,4,15])

  result1=(((Q_list[i-1][j]*Q_trans)*rho_row[i-1][j])*H_row[2*i-1][2*j+1])*Q_transN
  result1.permute([6,7,9,10,15],4)

  result1=(result1*Swap)*Swap1
  result1.permute([-6,-7,8,9,10],4)
  result=result1+result

  Swap.setLabel([12,3,13,-3])
  Swap1.setLabel([13,4,6,-4])

  Q_list[i-1][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=(Q_list[i-1][j]*Swap)*Swap1
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,-2,3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,7,8,9])


  H_row[2*i-1][2*j].setLabel([-2,-12,2,12])

  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*Q_trans)*rho_row[i-1][j])*H_row[2*i-1][2*j])*Q_transN
  result1.permute([6,7,8,9,10],4)
  result=result1+result

 return result



# 
# def Q_cost_val_middle( H_col, H_row, Q_list, i, j , rho_row, rho_col, N_x, Sys,HH_long):
#  h_long=HH_long[0][0]
#  result=0
#  #print h_long
#  Swap=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )
#  Swap1=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )
# 
# ############################   col   #####################################
#  if j==N_x-1:
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Swap.setLabel([-2,-3,2,3])
# 
#   Q_list_ferm=Swap*Q_list[i][j]
#   Q_list_ferm.permute([1,-2,-3,4,5],5)
# 
#   Q_list_ferm.setLabel([1,2,3,4,5])
#   Q_trans=Q_list_ferm*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-1,2,-3,4,-5])
# 
# 
#   Q_list[i][j-1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j-1]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,6,7,8,9])
#   H_col[2*i][2*j].setLabel([-1,-3,1,3])
#   rho_col[i][j-1].setLabel([-10,-5,10,5])
# 
#   result1=(((Q_list_ferm*H_col[2*i][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)
#   result=result1[0]+result
#   #print "12",result1
# ################################################
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Swap.setLabel([-3,-4,3,4])
#   Q_list_ferm=Swap*Q_list[i][j]
#   Q_list_ferm.permute([1,2,-3,-4,5],5)
# 
#   Q_list_ferm.setLabel([1,2,3,4,5])
#   Q_trans=Q_list_ferm*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([1,-2,3,-4,-5])
# 
# 
#   Q_list[i][j-1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j-1]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,6,7,8,9])
# 
#   H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
#   rho_col[i][j-1].setLabel([-10,-5,10,5])
#   result1=(((Q_list_ferm*H_col[2*i+1][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)
# 
#   result=result1[0]+result
#   #print "11",result1
# 
# #######################################################
# 
# 
#   Swap.setLabel([12,4,6,-4])
# 
#   Q_list[i][j-1].setLabel([1,2,3,-4,5])
#   Q_list_fermin=Q_list[i][j-1]*Swap
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([1,2,-3,4,-5,-12,-6])
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,-6,7,8,9])
# 
# 
#   rho_col[i][j-1].setLabel([-5,-10,5,10])
#   H_col[2*i][2*j-1].setLabel([-3,-12,3,12])
# 
#   result1=((Q_list[i][j]*Q_transN)*rho_col[i][j-1])*(((Q_list_fermin*H_col[2*i][2*j-1])*Q_trans))
#   result=result1[0]+result
#   #print "10",result1
# 
#   Swap.setLabel([7,6,-7,-6])
# 
#   Q_list[i][j-1].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j-1]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,2,3,-4])
# 
#   Q_list[i][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=Q_list[i][j]*Swap
#   Q_list_ferminP.permute([6,7,8,9,10],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([6,-7,8,9,-10])
# 
#   H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
#   rho_col[i][j-1].setLabel([-5,-10,5,10])
# 
#   result1=((Q_list_ferminP*H_col[2*i+1][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_trans)
#   result=result1[0]+result
# 
#   #print "9",result1
# 
# 
# 
# 
#  if j<N_x-1:
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Swap.setLabel([-2,-3,2,3])
# 
#   Q_list_ferm=Swap*Q_list[i][j]
#   Q_list_ferm.permute([1,-2,-3,4,5],5)
# 
#   Q_list_ferm.setLabel([1,2,3,4,5])
#   Q_trans=Q_list_ferm*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-1,2,-3,4,-5])
# 
#   Q_list[i][j+1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j+1]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,6,7,8,9])
#  
#   H_col[2*i][2*j].setLabel([-1,-3,1,3])
#   rho_col[i][j].setLabel([-5,-10,5,10])
#   result1=((Q_list_ferm*H_col[2*i][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]
#   result=result1[0]+result
#   #print "8", result1[0]
# 
# ################################################
# 
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Swap.setLabel([-3,-4,3,4])
#   Q_list_ferm=Swap*Q_list[i][j]
#   Q_list_ferm.permute([1,2,-3,-4,5],4)
# 
#   Q_list_ferm.setLabel([1,2,3,4,5])
#   Q_trans=Q_list_ferm*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,-2,3,-4])
# 
# 
#   Q_list[i][j+1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j+1]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,6,7,8,9])
# 
#   H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
#   rho_col[i][j].setLabel([-5,-10,5,10])
#   result1=((Q_list_ferm*H_col[2*i+1][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]
#   result=result1[0]+result
#   #print "7", result1[0]
# #######################################################
# 
#   Swap.setLabel([12,4,6,-4])
# 
#   Q_list[i][j].setLabel([1,2,3,-4,5])
#   Q_list_fermin=Q_list[i][j]*Swap
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([1,2,-3,4,-5,-12,-6])
# 
#   Q_list[i][j+1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j+1]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,-6,7,8,9])
# 
#   rho_col[i][j].setLabel([-5,-10,5,10])
# 
#   H_col[2*i][2*j+1].setLabel([-3,-12,3,12])
#   result1=(Q_list_fermin*(Q_trans*((Q_list[i][j+1]*Q_transN)*rho_col[i][j])))*H_col[2*i][2*j+1]
#   result=result1[0]+result
#   #print "6", result1
# 
# 
#   Swap.setLabel([7,6,-7,-6])
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,2,3,-4])
# 
#   Q_list[i][j+1].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=Q_list[i][j+1]*Swap
#   Q_list_ferminP.permute([6,7,8,9,10],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([6,-7,8,9,-10])
# 
# 
#   H_col[2*i+1][2*j+1].setLabel([-4,-7,4,7])
#   rho_col[i][j].setLabel([-5,-10,5,10])
# 
#   result1=(((Q_list[i][j]*Q_trans)*H_col[2*i+1][2*j+1])*rho_col[i][j])*(Q_list_ferminP*Q_transN)
#   result=result1[0]+result
#   #print "5",result1
# 
# 
# 
# ###############################################################################
#  if j-1>=0 and j!=N_x-1:
# 
#   Swap.setLabel([12,4,6,-4])
# 
#   Q_list[i][j-1].setLabel([1,2,3,-4,5])
#   Q_list_fermin=Q_list[i][j-1]*Swap
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([1,2,-3,4,-5,-12,-6])
# 
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,-6,7,8,9])
# 
# 
#   rho_col[i][j-1].setLabel([-5,-10,5,10])
# 
#   H_col[2*i][2*j-1].setLabel([-3,-12,3,12])
#   result1=((Q_list[i][j]*H_col[2*i][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list_fermin*Q_trans)
#   result=result1[0]+result
#   #print "outb0", result1[0]
#   #print "4",result1
# 
#   Swap.setLabel([7,6,-7,-6])
# 
#   Q_list[i][j-1].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j-1]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,2,3,-4])
# 
#   Q_list[i][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=Q_list[i][j]*Swap
#   Q_list_ferminP.permute([6,7,8,9,10],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([6,-7,8,9,-10])
# 
#   H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
#   rho_col[i][j-1].setLabel([-5,-10,5,10])
# 
#   result1=((Q_list_ferminP*H_col[2*i+1][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_trans)
#   result=result1[0]+result
#   #print "3", result1
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##################### Row ####################################
# 
# 
#  if i==N_x-1:
# 
#   Q_list[i-1][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i-1][j]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,2,3,4])
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,6,7,-8,-9])
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   H_row[2*i][2*j+1].setLabel([-8,-9,8,9])
#   result1=((Q_list[i][j]*H_row[2*i][2*j+1]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
#   result=result1[0]+result
#   #print "1", result1
# 
# 
#   Q_list[i-1][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i-1][j]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,2,3,4])
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,-6,-7,8,9])
# 
#   H_row[2*i][2*j].setLabel([-6,-7,6,7])
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   result1=((Q_list[i][j]*H_row[2*i][2*j]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
#   result=result1[0]+result
#   #print "2",result1
# 
# 
# 
# 
# 
#   Swap.setLabel([15,6,16,-6])
#   Swap1.setLabel([16,7,8,-7])
# 
#   Q_list[i-1][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i-1][j]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,2,3,-4])
# 
#   Q_list[i][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
#   Q_list_ferminP.permute([6,7,9,10,15],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([6,7,9,-10,-15])
# 
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   H_row[2*i-1][2*j+1].setLabel([-4,-15,4,15])
#   result1=((Q_list_ferminP*H_row[2*i-1][2*j+1]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
#   result=result1[0]+result
#   #print result1
#   #print "16",result1
# 
#   Swap.setLabel([12,3,13,-3])
#   Swap1.setLabel([13,4,6,-4])
# 
#   Q_list[i-1][j].setLabel([1,2,-3,-4,5])
#   Q_list_fermin=(Q_list[i-1][j]*Swap)*Swap1
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([1,-2,3,4,-5,-12,-6])
# 
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,-6,7,8,9])
# 
# 
#   H_row[2*i-1][2*j].setLabel([-2,-12,2,12])
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   result1=((Q_list[i][j]*Q_transN)*rho_row[i-1][j])*((Q_list_fermin*Q_trans)*H_row[2*i-1][2*j])
# 
#   result=result1[0]+result
#   #print "17",result1
# 
# 
# 
# 
# 
# 
# 
# 
#  if i<N_x-1:
#  #######################################################
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,2,-3,-4])
# 
#   Q_list[i+1][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i+1][j]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,6,7,8,9])
# 
#   rho_row[i][j].setLabel([-5,-10,5,10])
#   H_row[2*i][2*j+1].setLabel([-3,-4,3,4])
#   result1=(Q_list[i][j]*H_row[2*i][2*j+1]*Q_trans)*((Q_list[i+1][j]*Q_transN)*rho_row[i][j])
#   result=result1[0]+result
#   #print "19", result1
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,-1,-2,3,4])
# 
#   Q_list[i+1][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i+1][j]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,6,7,8,9])
# 
#   H_row[2*i][2*j].setLabel([-1,-2,1,2])
#   rho_row[i][j].setLabel([-5,-10,5,10])
# 
#   result1=((Q_list[i][j]*H_row[2*i][2*j]*Q_trans)*rho_row[i][j])*(Q_list[i+1][j]*Q_transN)
# 
#   result=result1[0]+result
#   #print "18", result1
# 
#  #######################################################
# 
#   Swap.setLabel([15,6,16,-6])
#   Swap1.setLabel([16,7,8,-7])
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,2,3,-4])
# 
#   Q_list[i+1][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=(Q_list[i+1][j]*Swap)*Swap1
#   Q_list_ferminP.permute([6,7,9,10,15],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.transpose()
# 
#   Q_transN.setLabel([6,7,9,-10,-15])
# 
#   rho_row[i][j].setLabel([-5,-10,5,10])
#   H_row[2*i+1][2*j+1].setLabel([-4,-15,4,15])
#   result1=(Q_list[i][j]*H_row[2*i+1][2*j+1]*Q_trans)*((Q_list_ferminP*Q_transN)*rho_row[i][j])
#   result=result1[0]+result
#   #print "20",result1
# 
# 
#   Swap.setLabel([12,3,13,-3])
#   Swap1.setLabel([13,4,6,-4])
#   Swapp=Swap*1.0
#   Swapp1=Swap1*1.0
#   Swapp.setLabel([-12,3,130,-33])
#   Swapp1.setLabel([130,4,-6,-44])
# 
# 
# 
#   Q_list[i][j].setLabel([1,2,-3,-4,5])
#   Q_list_fermin=(Q_list[i][j])*(Swap1*Swap)
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([1,-2,3,4,-5,-12,-6])
# 
#   Q_transS=Q_list[i][j]*1.0
#   Q_transS.transpose()
#   Q_transS.setLabel([-5,1,-2,-33,-44])
# 
# 
# 
#   Q_list[i+1][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i+1][j]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,-6,7,8,9])
# 
#   H_row[2*i+1][2*j].setLabel([-2,-12,2,12])
#   rho_row[i][j].setLabel([-5,-10,5,10])
# 
#   result1=H_row[2*i+1][2*j]*((Q_trans*((Q_list[i+1][j]*Q_transN)*rho_row[i][j]))*Q_list_fermin)
# 
# 
#   result=result1[0]+result
#   #print "21",result1
# 
# 
#  ########################################################################################
# 
# 
#  if i-1>=0 and i!=N_x-1:
# 
#   Swap.setLabel([15,6,16,-6])
#   Swap1.setLabel([16,7,8,-7])
# 
#   Q_list[i-1][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i-1][j]*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([-5,1,2,3,-4])
# 
#   Q_list[i][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
#   Q_list_ferminP.permute([6,7,9,10,15],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([6,7,9,-10,-15])
# 
# 
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   H_row[2*i-1][2*j+1].setLabel([-4,-15,4,15])
#   result1=((Q_list_ferminP*H_row[2*i-1][2*j+1]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
# 
#   result=result1[0]+result
#   #print "22",result1
# 
# 
#   Swap.setLabel([12,3,13,-3])
#   Swap1.setLabel([13,4,6,-4])
# 
#   Q_list[i-1][j].setLabel([1,2,-3,-4,5])
#   Q_list_fermin=(Q_list[i-1][j]*Swap)*Swap1
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   Q_trans.transpose()
#   Q_trans.setLabel([1,-2,3,4,-5,-12,-6])
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.transpose()
#   Q_transN.setLabel([-10,-6,7,8,9])
# 
# 
#   H_row[2*i-1][2*j].setLabel([-2,-12,2,12])
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   result1=((Q_list[i][j]*Q_transN)*rho_row[i-1][j])*((Q_list_fermin*Q_trans)*H_row[2*i-1][2*j])
#   result=result1[0]+result
#   #print "23",result1
# 
# 
# 
#  return result
# 
# 
# def  Grad_Q_middle( H_col, H_row, Q_list, i, j , rho_row, rho_col, N_x, Sys,HH_long):
# 
# 
#  h_long=HH_long[0][0]
# 
#  Q_list[i][j].setLabel([0,1,2,3,4])
#  Q_list[i][j].permute([0,1,2,3,4],4)
# 
#  Q_t=Q_list[i][j]*1.0
#  Q_t.transpose()
#  Q_t.setLabel([4,0,1,2,3])
#  Q_t.permute([0,1,2,3,4],4)
#  result=Q_t*0.0
# 
#  Swap=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )
#  Swap1=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )
# 
#  ###########################   col   #####################################
#  if j==N_x-1:
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Swap.setLabel([-2,-3,2,3])
# 
#   Q_list_ferm=Swap*Q_list[i][j]
#   Q_list_ferm.permute([1,-2,-3,4,5],5)
#   Q_list_ferm.setLabel([1,2,3,4,5])
#   Q_trans=Q_list_ferm*1.0
#   #Q_trans.transpose()
#   Q_trans.setLabel([-1,2,-3,4,-5])
# 
#   Q_list[i][j-1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j-1]*1.0
#   #Q_transN.transpose()
#   Q_transN.setLabel([6,7,8,9,-10])
# 
#   H_col[2*i][2*j].setLabel([-1,-3,1,3])
#   rho_col[i][j-1].setLabel([-10,-5,10,5])
# 
#   result1=(((H_col[2*i][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)
#   result1.permute([1,2,3,4,5],4)
# 
#   Swap.setLabel([-2,-3,2,3])
#   result1=Swap*result1
#   result1.permute([1,-2,-3,4,5],4)
#   result1.setLabel([1,2,3,4,5])
# 
#   result=result1+result
# 
#  ################################################
# 
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Swap.setLabel([-3,-4,3,4])
#   Q_list_ferm=Swap*Q_list[i][j]
#   Q_list_ferm.permute([1,2,-3,-4,5],4)
# 
#   Q_list_ferm.setLabel([1,2,3,4,5])
#   Q_trans=Q_list_ferm*1.0
#   #Q_trans.transpose()
#   Q_trans.setLabel([1,-2,3,-4,-5])
# 
# 
#   Q_list[i][j-1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j-1]*1.0
#   #Q_transN.transpose()
#   Q_transN.setLabel([6,7,8,9,-10])
# 
#   H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
#   rho_col[i][j-1].setLabel([-10,-5,10,5])
#   result1=(((H_col[2*i+1][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)
#   result1.permute([1,2,3,4,5],4)
# 
#   Swap.setLabel([-3,-4,3,4])
#   result1=Swap*result1
#   result1.permute([1,2,-3,-4,5],4)
#   result1.setLabel([1,2,3,4,5])
# 
#   result=result1+result
# 
# 
# #################################################################################
#   Swap.setLabel([12,4,6,-4])
# 
#   SwapP=Swap*1.0
#   SwapP.setLabel([-12,4,-6,-40])
# 
# 
#   Q_list[i][j-1].setLabel([1,2,3,-4,5])
#   Q_list_fermin=Q_list[i][j-1]*Swap
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   #Q_trans.transpose()
#   Q_trans.setLabel([1,2,-3,4,-5,-12,-6])
# 
#   Q_transS=Q_list[i][j-1]*1.0
#   #Q_transS.transpose()
#   Q_transS.setLabel([1,2,-3,-40,-5])
# 
# 
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   #Q_transN.transpose()
#   Q_transN.setLabel([-6,7,8,9,-10])
# 
# 
#   rho_col[i][j-1].setLabel([-5,-10,5,10])
#   H_col[2*i][2*j-1].setLabel([-3,-12,3,12])
# 
# 
#   result1=((((Q_list[i][j-1]*Swap)*H_col[2*i][2*j-1])*(Q_transS*SwapP))*rho_col[i][j-1])*Q_transN
# 
#   result1.permute([6,7,8,9,10],4)
# 
#   result=result1+result
# 
# ##################
# 
#   Swap.setLabel([7,6,-7,-6])
# 
#   Q_list[i][j-1].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j-1]*1.0
#   #Q_trans.transpose()
#   Q_trans.setLabel([1,2,3,-4,-5])
# 
#   Q_list[i][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=Q_list[i][j]*Swap
#   Q_list_ferminP.permute([6,7,8,9,10],5)
#   Q_transN=Q_list_ferminP*1.0
#   #Q_transN.transpose()
#   Q_transN.setLabel([6,-7,8,9,-10])
# 
# 
#   H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
# 
#   rho_col[i][j-1].setLabel([-5,-10,5,10])
# 
#   result1=((H_col[2*i+1][2*j-1])*((Q_list[i][j-1]*Q_trans)*rho_col[i][j-1]))*Q_transN
#   result1.permute([6,7,8,9,10],4)
# 
#   Swap.setLabel([7,6,-7,-6])
#   result1=result1*Swap
#   result1.permute([-6,-7,8,9,10],4)
# 
#   result=result1+result
# 
# 
# 
#  if j<N_x-1:
# 
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Swap.setLabel([-2,-3,2,3])
# 
#   Q_list_ferm=Swap*Q_list[i][j]
#   Q_list_ferm.permute([1,-2,-3,4,5],5)
# 
#   Q_list_ferm.setLabel([1,2,3,4,5])
#   Q_trans=Q_list_ferm*1.0
#   #Q_trans.transpose()
#   Q_trans.setLabel([-1,2,-3,4,-5])
# 
# 
#   Q_list[i][j+1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j+1]*1.0
#   #Q_transN.transpose()
#   Q_transN.setLabel([6,7,8,9,-10])
# 
#   H_col[2*i][2*j].setLabel([-1,-3,1,3])
#   rho_col[i][j].setLabel([-5,-10,5,10])
# 
#   result1=((H_col[2*i][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]
#   result1.permute([1,2,3,4,5],4)
# 
#   Swap.setLabel([-2,-3,2,3])
#   result1=Swap*result1
#   result1.permute([1,-2,-3,4,5],4)
#   result1.setLabel([1,2,3,4,5])
# 
# 
#   result=result1+result
# 
# ######################################################################
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Swap.setLabel([-3,-4,3,4])
#   Q_list_ferm=Swap*Q_list[i][j]
#   Q_list_ferm.permute([1,2,-3,-4,5],5)
# 
#   Q_list_ferm.setLabel([1,2,3,4,5])
#   Q_trans=Q_list_ferm*1.0
#   #Q_trans.transpose()
#   Q_trans.setLabel([1,-2,3,-4,-5])
# 
# 
#   Q_list[i][j+1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j+1]*1.0
#   #Q_transN.transpose()
#   Q_transN.setLabel([6,7,8,9,-10])
# 
#   H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
#   rho_col[i][j].setLabel([-5,-10,5,10])
#   result1=((H_col[2*i+1][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]
# 
#   Swap.setLabel([-3,-4,3,4])
#   result1=Swap*result1
#   result1.permute([1,2,-3,-4,5],4)
#   result1.setLabel([1,2,3,4,5])
# 
#   result1.permute([1,2,3,4,5],4)
#   result=result1+result
# 
# #####################################################################
# 
# 
#   Swap.setLabel([12,4,6,-4])
# 
#   Q_list[i][j].setLabel([1,2,3,-4,5])
#   Q_list_fermin=Q_list[i][j]*Swap
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   #Q_trans.transpose()
#   Q_trans.setLabel([1,2,-3,4,-5,-12,-6])
# 
#   Q_list[i][j+1].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j+1]*1.0
#   #Q_transN.transpose()
#   Q_transN.setLabel([-6,7,8,9,-10])
# 
#   rho_col[i][j].setLabel([-5,-10,5,10])
# 
#   H_col[2*i][2*j+1].setLabel([-3,-12,3,12])
# 
#   result1=(H_col[2*i][2*j+1])*(((Q_list[i][j+1]*Q_transN)*rho_col[i][j])*Q_trans)
#   result1.permute([1,2,3,4,5,12,6],4)
# 
# 
#   result1=result1*Swap
#   result1.permute([1,2,3,-4,5],4)
#   result=result1+result
# 
# 
#   Swap.setLabel([7,6,-7,-6])
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j]*1.0
#   #Q_trans.transpose()
#   Q_trans.setLabel([1,2,3,-4,-5])
# 
#   Q_list[i][j+1].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=Q_list[i][j+1]*Swap
#   Q_list_ferminP.permute([6,7,8,9,10],5)
#   Q_transN=Q_list_ferminP*1.0
#   #Q_transN.transpose()
#   Q_transN.setLabel([6,-7,8,9,-10])
# 
#   H_col[2*i+1][2*j+1].setLabel([-4,-7,4,7])
# 
#   rho_col[i][j].setLabel([-5,-10,5,10])
# 
#   result1=(((Q_list_ferminP*Q_transN)*rho_col[i][j])*H_col[2*i+1][2*j+1])*Q_trans
#   result1.permute([1,2,3,4,5],4)
# 
#   result=result1+result
# 
# 
# 
# ##############################################################################
#  if j-1>=0 and j!=N_x-1:
# 
#   Swap.setLabel([12,4,6,-4])
# 
#   Q_list[i][j-1].setLabel([1,2,3,-4,5])
#   Q_list_fermin=Q_list[i][j-1]*Swap
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   #Q_trans.transpose()
#   Q_trans.setLabel([1,2,-3,4,-5,-12,-6])
# 
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   #Q_transN.transpose()
#   Q_transN.setLabel([-6,7,8,9,-10])
# 
# 
#   rho_col[i][j-1].setLabel([-5,-10,5,10])
# 
# 
#   H_col[2*i][2*j-1].setLabel([-3,-12,3,12])
#   result1=(((Q_list_fermin*Q_trans)*rho_col[i][j-1])*H_col[2*i][2*j-1])*Q_transN
# 
#   result1.permute([6,7,8,9,10],4)
#   result=result1+result
# 
# 
#   Swap.setLabel([7,6,-7,-6])
# 
#   Q_list[i][j-1].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j-1]*1.0
#   Q_trans.setLabel([1,2,3,-4,-5])
# 
#   Q_list[i][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=Q_list[i][j]*Swap
#   Q_list_ferminP.permute([6,7,8,9,10],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.setLabel([6,-7,8,9,-10])
# 
#   H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
#   rho_col[i][j-1].setLabel([-5,-10,5,10])
# 
#   result1=(((Q_list[i][j-1]*Q_trans)*rho_col[i][j-1])*H_col[2*i+1][2*j-1])*Q_transN
#   result1.permute([6,7,8,9,10],4)
# 
#   result1=result1*Swap
#   result1.permute([-6,-7,8,9,10],4)
# 
#   result=result1+result
# 
# ##################### Row ####################################
# 
#  if i==N_x-1:
# 
#   Q_list[i-1][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i-1][j]*1.0
#   Q_trans.setLabel([1,2,3,4,-5])
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.setLabel([6,7,-8,-9,-10])
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   H_row[2*i][2*j+1].setLabel([-8,-9,8,9])
#   result1=((H_row[2*i][2*j+1]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
#   result1.permute([6,7,8,9,10],4)  
#   result=result1+result
# 
#   Q_list[i-1][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i-1][j]*1.0
#   Q_trans.setLabel([1,2,3,4,-5])
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.setLabel([-6,-7,8,9,-10])
# 
#   H_row[2*i][2*j].setLabel([-6,-7,6,7])
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   result1=((H_row[2*i][2*j]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
#   result1.permute([6,7,8,9,10],4)  
#   result=result1+result
# 
# 
#   Swap.setLabel([15,6,16,-6])
#   Swap1.setLabel([16,7,8,-7])
# 
#   Q_list[i-1][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i-1][j]*1.0
#   Q_trans.setLabel([1,2,3,-4,-5])
# 
#   Q_list[i][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
#   Q_list_ferminP.permute([6,7,9,10,15],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.setLabel([6,7,9,-10,-15])
# 
# 
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   H_row[2*i-1][2*j+1].setLabel([-4,-15,4,15])
# 
#   result1=(((Q_list[i-1][j]*Q_trans)*rho_row[i-1][j])*H_row[2*i-1][2*j+1])*Q_transN
#   result1.permute([6,7,9,10,15],4)    
# 
#   result1=(result1*Swap)*Swap1
#   result1.permute([-6,-7,8,9,10],4)
# 
#   result=result1+result
# ##################################
# 
#   Swap.setLabel([12,3,13,-3])
#   Swap1.setLabel([13,4,6,-4])
#   Swapp=Swap*1.0
#   Swapp1=Swap1*1.0
#   Swapp.setLabel([-12,-33,-13,3])
#   Swapp1.setLabel([-13,4,-6,-44])
# 
# 
#   Q_list[i-1][j].setLabel([1,2,-3,-4,5])
#   Q_list_fermin=(Q_list[i-1][j]*Swap)*Swap1
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   Q_trans.setLabel([1,-2,3,4,-5,-12,-6])
# 
#   Q_list[i-1][j].setLabel([1,2,-3,-4,5])
#   Q_transS=Q_list[i-1][j]*1.0
#   Q_transS.setLabel([1,-2,-33,-44,-5])
# 
# 
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.setLabel([-6,7,8,9,-10])
# 
#   H_row[2*i-1][2*j].setLabel([-2,-12,2,12])
# 
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
# 
#   result1=(((((Q_list[i-1][j]*(Swap*Swap1))*H_row[2*i-1][2*j])*((Swapp*Swapp1)))*Q_transS)*rho_row[i-1][j])*Q_transN
# 
#   result1.permute([6,7,8,9,10],4)
# 
#   result=result1+result
# 
# 
#  if i<N_x-1:
#  #######################################################
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j]*1.0
#   Q_trans.setLabel([1,2,-3,-4,-5])
# 
#   Q_list[i+1][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i+1][j]*1.0
#   Q_transN.setLabel([6,7,8,9,-10])
# 
#   rho_row[i][j].setLabel([-5,-10,5,10])
#   H_row[2*i][2*j+1].setLabel([-3,-4,3,4])
#   result1=(H_row[2*i][2*j+1]*Q_trans)*((Q_list[i+1][j]*Q_transN)*rho_row[i][j])
#   result1.permute([1,2,3,4,5],4)    
#   result=result1+result
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j]*1.0
#   Q_trans.setLabel([-1,-2,3,4,-5])
# 
#   Q_list[i+1][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i+1][j]*1.0
#   Q_transN.setLabel([6,7,8,9,-10])
# 
#   H_row[2*i][2*j].setLabel([-1,-2,1,2])
#   rho_row[i][j].setLabel([-5,-10,5,10])
# 
#   result1=(H_row[2*i][2*j])*(((Q_list[i+1][j]*Q_transN)*rho_row[i][j])*Q_trans)
#   result1.permute([1,2,3,4,5],4)    
#   result=result1+result
# 
#  #######################################################
#   Swap.setLabel([15,6,16,-6])
#   Swap1.setLabel([16,7,8,-7])
# 
#   Q_list[i][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i][j]*1.0
#   Q_trans.setLabel([1,2,3,-4,-5])
# 
#   Q_list[i+1][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=Q_list[i+1][j]*(Swap1*Swap)
#   Q_list_ferminP.permute([6,7,9,10,15],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.setLabel([6,7,9,-10,-15])
# 
# 
# 
#   rho_row[i][j].setLabel([-5,-10,5,10])
#   H_row[2*i+1][2*j+1].setLabel([-4,-15,4,15])
# 
#   result1=(H_row[2*i+1][2*j+1]*Q_trans)*((Q_list_ferminP*Q_transN)*rho_row[i][j])
#   result1.permute([1,2,3,4,5],4)      
#   result=result1+result
# 
# 
#   Swap.setLabel([12,3,13,-3])
#   Swap1.setLabel([13,4,6,-4])
# 
#   Q_list[i][j].setLabel([1,2,-3,-4,5])
#   Q_list_fermin=Q_list[i][j]*(Swap1*Swap)
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   Q_trans.setLabel([1,-2,3,4,-5,-12,-6])
# 
# 
#   Q_list[i+1][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i+1][j]*1.0
#   Q_transN.setLabel([-6,7,8,9,-10])
# 
#   H_row[2*i+1][2*j].setLabel([-2,-12,2,12])
# 
#   rho_row[i][j].setLabel([-5,-10,5,10])
# 
#   result1=(((Q_list[i+1][j]*Q_transN)*rho_row[i][j])*Q_trans)*H_row[2*i+1][2*j]
#   result1.permute([1,2,3,4,5,12,6],4)
# 
#   result1=result1*(Swap1*Swap)
#   result1.permute([1,2,-3,-4,5],4)
# 
#   result=result1+result
#  ########################################################################################
# 
# 
#  if i-1>=0 and i!=N_x-1:
# 
#   Swap.setLabel([15,6,16,-6])
#   Swap1.setLabel([16,7,8,-7])
# 
#   Q_list[i-1][j].setLabel([1,2,3,4,5])
#   Q_trans=Q_list[i-1][j]*1.0
#   Q_trans.setLabel([1,2,3,-4,-5])
# 
#   Q_list[i][j].setLabel([-6,-7,8,9,10])
#   Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
#   Q_list_ferminP.permute([6,7,9,10,15],5)
#   Q_transN=Q_list_ferminP*1.0
#   Q_transN.setLabel([6,7,9,-10,-15])
# 
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
#   H_row[2*i-1][2*j+1].setLabel([-4,-15,4,15])
# 
#   result1=(((Q_list[i-1][j]*Q_trans)*rho_row[i-1][j])*H_row[2*i-1][2*j+1])*Q_transN
#   result1.permute([6,7,9,10,15],4)
# 
#   result1=(result1*Swap)*Swap1
#   result1.permute([-6,-7,8,9,10],4)
#   result=result1+result
# 
#   Swap.setLabel([12,3,13,-3])
#   Swap1.setLabel([13,4,6,-4])
# 
#   Q_list[i-1][j].setLabel([1,2,-3,-4,5])
#   Q_list_fermin=(Q_list[i-1][j]*Swap)*Swap1
#   Q_list_fermin.permute([1,2,3,4,5,12,6],7)
#   Q_trans=Q_list_fermin*1.0
#   Q_trans.setLabel([1,-2,3,4,-5,-12,-6])
# 
#   Q_list[i][j].setLabel([6,7,8,9,10])
#   Q_transN=Q_list[i][j]*1.0
#   Q_transN.setLabel([-6,7,8,9,-10])
# 
# 
#   H_row[2*i-1][2*j].setLabel([-2,-12,2,12])
# 
#   rho_row[i-1][j].setLabel([-5,-10,5,10])
# 
#   result1=(((Q_list_fermin*Q_trans)*rho_row[i-1][j])*H_row[2*i-1][2*j])*Q_transN
#   result1.permute([6,7,8,9,10],4)
#   result=result1+result
# 
# 
# 
# 
# 
# 
#  return result
# 




















def make_Env_singleLayer( PEPS_listten, Location, mps_boundry, d, chi_boundry, N_x,Sys):

 Peps_ket=[]
 for i in xrange( len(PEPS_listten) ):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)
  Peps_ket.append(A*1.0)


 Peps_bra=[]
 for i in xrange(len(PEPS_listten)):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)

  A.permute([1,2,3,4,5],3)
  A.setLabel([1,2,3,4,5])

  A_conj=A*1.0
  A_conj.transpose()
  A_conj.setLabel([-4,-5,-1,-2,3])
  A_conj.permute([-1,-2,3,-4,-5],5)
  Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
  Swap2.setLabel([-6,7,-4,-5])
  Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
  Swap3.setLabel([8,9,-1,-2])
  A_conj=(A_conj*Swap2)*Swap3
  A_conj.permute([8,9,3,-6,7],5)
  A_conj.setLabel([-1,-2,-3,-4,-5])

  Peps_bra.append(A_conj*1.0)


 chi_boundry_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi_boundry)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_boundry_list.append(dim)

 DX=0
 D=0
 if Location!=(N_x-1):
  DX=Peps_ket[0].bond(3).Qlist()
  D=Peps_bra[0].bond(3).Qlist()
 else:
  DX=Peps_ket[0].bond(0).Qlist()
  D=Peps_bra[0].bond(0).Qlist()

 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 #############Zero#######################
 bdi_1 = uni10.Bond( uni10.BD_IN, 1)
 bdo_1 = uni10.Bond( uni10.BD_IN, 1)
 Tem0=uni10.UniTensor([bdi_1])
 Tem0.identity()
 Tem0.setLabel([0])
 #print Tem0, 

 Tem1=uni10.UniTensor([bdi_1])
 Tem1.identity()
 Tem1.setLabel([1])
 #print Tem1

 Tem11=uni10.UniTensor([bdi_1])
 Tem11.identity()
 Tem11.setLabel([11])
 #print Tem11


 Tem10=uni10.UniTensor([bdi_1])
 Tem10.identity()
 Tem10.setLabel([10])
 #print Tem10


 if Location==0:
  mps_list=[None]*(N_x*2)

  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([1,10,-1,-10])

  Peps_ket[0].setLabel([0,-1,2,3,4])
  Peps_bra[0].setLabel([-10,11,2,-3,-4])

  Tem10.identity()
  Tem1.identity()
  Tem11.identity()

  mps_list[0]=(((Peps_ket[0]*Tem0)*Tem1)*((Peps_bra[0])*Tem11))*SwapX
  mps_list[0].permute([10, -3, 4, 3, -4], 2)
  ##########################################


  bdi = uni10.Bond(uni10.BD_IN, Peps_ket[0].bond(4).Qlist())
  bdo = uni10.Bond(uni10.BD_OUT, Peps_ket[0].bond(4).Qlist())
  IdenX=uni10.UniTensor([bdi, bdo])
  IdenX.setLabel([1,4])
  IdenX.identity()
  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapX.setLabel([6,11,3,-4])

  mps_list[1]=SwapX*IdenX
  mps_list[1].permute([4,3,-4,6,1,11],4)
  #mps_list.append(results)

  #print mps_list[0].printDiagram(), mps_list[1].printDiagram() 

  ####################Middle##############
  for  q  in  xrange(1,len(Peps_ket)):
   Peps_ket[q].setLabel([0,1,2,3,4])
   Peps_bra[q].setLabel([10,11,2,-3,-4])
   Tem0.identity()
   Tem0.setLabel([0])
   Tem10.identity()
   Tem10.setLabel([-10])
   #print "Hi4"

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(1), Peps_bra[q].bond(0))
   SwapX.setLabel([-1,-10,1,10])

   mps_list[2*q]=(((Peps_ket[q]*Tem0))*((Peps_bra[q])*Tem10))*SwapX
   mps_list[2*q].permute([-1,11,-3,4,3,-4],3)

   ######################################################


   bdi = uni10.Bond( uni10.BD_IN, Peps_ket[q].bond(4).Qlist())
   bdo = uni10.Bond( uni10.BD_OUT, Peps_ket[q].bond(4).Qlist())
   IdenX=uni10.UniTensor([bdi, bdo])
   IdenX.setLabel([1,4])
   IdenX.identity()

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapX.setLabel([6,11,3,-4])

   mps_list[2*q+1]=IdenX*SwapX
   mps_list[2*q+1].permute([4,3,-4,6,1,11],4)

  mps_boundry=MPSclass.MPS(2,2,len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "0, Single", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()
  return mps_boundry

  ######################## Add next Layer ###########

 if Location!=0:
  mps_list=[None]*(N_x*2)

  ###############################First-Layer###############################

  mps_boundry[0].setLabel([1,2,3])
  Tem1.setLabel([4])
  SwapX=fermionicOPT(Sys,mps_boundry[0].bond(1), Peps_ket[0].bond(1))
  SwapX.setLabel([5,6,2,4])
  mps_list[0]=(mps_boundry[0]*SwapX)*Tem1
  mps_list[0].permute([1,5,3,6],2)

  mps_boundry[1].setLabel([5,0,6])
  Peps_ket[0].setLabel([0,1,2,3,4])
  Tem1.setLabel([1])
  mps_list[1]=(mps_boundry[1]*(Peps_ket[0]))
  mps_list[1].permute([5,1,2,3,6,4],4)

  ###########################################
 
  #print mps_list[0].printDiagram(), mps_list[1].printDiagram()
   ###############################################################################

  for q in xrange(1,len(Peps_ket)):
   ###########################################

    SwapX=fermionicOPT(Sys,mps_boundry[2*q].bond(1), Peps_ket[q-1].bond(4))
    SwapX.setLabel([10,-4,0,4])
    mps_boundry[2*q].setLabel([6,0,-6])

    mps_list[2*q]=mps_boundry[2*q]*SwapX
    mps_list[2*q].permute([6,4,10,-6,-4],3)



    mps_boundry[2*q+1].setLabel([5,0,6])
    Peps_ket[q].setLabel([0,1,2,3,4])

    mps_list[2*q+1]=mps_boundry[2*q+1]*Peps_ket[q]
    mps_list[2*q+1].permute([5,1,2,3,6,4],4)


  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
   mps_boundry[i]=mps_list[i]*1.0

  #print mps_list[2].printDiagram(), mps_list[3].printDiagram()

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD( sum(chi_boundry_list) )
  #print "Middle", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t


######################### Next Absorption ####################################
  for  q  in  xrange(0, mps_boundry.N, 2):
    mps_boundry[q].setLabel([10,3,2])
    Peps_bra[q/2].setLabel([3,11,1,4,-2])
    mps_list[q]=((mps_boundry[q]*Peps_bra[q/2]))
    mps_list[q].permute([10,11,4,2,1,-2],3)




    mps_boundry[q+1].setLabel([2,1,3,5])
    SwapX=fermionicOPT(Sys, Peps_bra[q/2].bond(4), mps_boundry[q+1].bond(2) )
    SwapX.setLabel([ -5, -3, -2, 3 ])
    mps_list[q+1]=mps_boundry[q+1]*SwapX
    mps_list[q+1].permute([2,1,-2,-3,5,-5],4)




  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Last", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
 
  return mps_boundry



######@profile
def make_Env_singleLayer_right( PEPS_listten, Location, mps_boundry, d, chi_boundry, N_x,Sys):

 Peps_ket=[]
 for i in xrange( len(PEPS_listten) ):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)
  Peps_ket.append(A*1.0)


# Peps_bra=[]
# for i in xrange(len(PEPS_listten)):
#  A=PEPS_listten[i]*1.0
#  A.setLabel([1,2,3,4,5])
#  A.permute([1,2,3,4,5],5)

#  Swap2=fermionicOPT(Sys,A.bond(3), A.bond(4))
#  Swap2.setLabel([-6,7,-4,-5])
#  Swap3=fermionicOPT(Sys,A.bond(0), A.bond(1))
#  Swap3.setLabel([-1,-2,8,9])
#  A.permute([1,2,3,4,5],3)
#  A.setLabel([1,2,3,4,5])

#  A_conj=A*1.0
#  A_conj.transpose()
#  A_conj.setLabel([-4,-5,-1,-2,3])
#  A_conj=(A_conj*Swap2)*Swap3
#  A_conj.permute([8,9,3,-6,7],5)
#  A_conj.setLabel([-1,-2,-3,-4,-5])
#  Peps_bra.append(A_conj*1.0)






 Peps_bra=[]
 for i in xrange(len(PEPS_listten)):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)

  A.permute([1,2,3,4,5],3)
  A.setLabel([1,2,3,4,5])

  A_conj=A*1.0
  A_conj.transpose()
  A_conj.setLabel([-4,-5,-1,-2,3])
  A_conj.permute([-1,-2,3,-4,-5],5)
  Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
  Swap2.setLabel([-6,7,-4,-5])
  Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
  Swap3.setLabel([8,9,-1,-2])
  A_conj=(A_conj*Swap2)*Swap3
  A_conj.permute([8,9,3,-6,7],5)
  A_conj.setLabel([-1,-2,-3,-4,-5])

  Peps_bra.append(A_conj*1.0)





 chi_boundry_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi_boundry)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_boundry_list.append(dim)

 DX=0
 D=0
 if Location!=(N_x-1):
  DX=Peps_ket[0].bond(3).Qlist()
  D=Peps_bra[0].bond(3).Qlist()
 else:
  DX=Peps_ket[0].bond(0).Qlist()
  D=Peps_bra[0].bond(0).Qlist()

 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 #############   Zero   ##############################################
 bdi_1 = uni10.Bond( uni10.BD_IN, 1)
 bdo_1 = uni10.Bond( uni10.BD_IN, 1)
 Tem0=uni10.UniTensor([bdi_1])
 Tem0.identity()
 Tem0.setLabel([0])
 #print Tem0, 

 Tem1=uni10.UniTensor([bdi_1])
 Tem1.identity()
 Tem1.setLabel([1])
 #print Tem1

 Tem11=uni10.UniTensor([bdi_1])
 Tem11.identity()
 Tem11.setLabel([11])
 #print Tem11


 Tem10=uni10.UniTensor([bdi_1])
 Tem10.identity()
 Tem10.setLabel([10])
 #print Tem10

 if Location==N_x-1:
  mps_list=[None]*(2*N_x)

  SwapX=fermionicOPT(Sys, Peps_ket[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([1,10,-1,-10])

  SwapXX=fermionicOPT(Sys, Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapXX.setLabel([5,6,3,-4])

  Peps_ket[0].setLabel([0,-1,2,3,4])
  Peps_bra[0].setLabel([-10,11,2,-3,-4])

  Tem10.identity()
  Tem1.identity()
  Tem11.identity()
  Tem10.setLabel([10])
  Tem1.setLabel([1])
  Tem0.setLabel([0])
  Temr=Tem0*1.0
  Temr.setLabel([-3])
  
  mps_list[0]=(Peps_bra[0]*Temr)*(SwapX*Tem1)
  mps_list[0].permute([11, 10, -4, 2, -1], 2)


  Temr.setLabel([5])
  mps_list[1]=(((Peps_ket[0]*Temr)*SwapXX))
  mps_list[1].permute([ -4, 2, -1,0,6,4], 4)


  ####################   Middle   ##########################################
  for  q  in  xrange(1,len(Peps_ket)):
   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(1), Peps_bra[q].bond(0))
   SwapX.setLabel([1,10,-1,-10])

   SwapXX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapXX.setLabel([5,6,3,-4])

   Peps_ket[q].setLabel([0,-1,2,3,4])
   Peps_bra[q].setLabel([-10,11,2,-3,-4])

   Tem10.identity()
   Tem1.identity()
   Tem11.identity()
   Tem10.setLabel([10])
   Tem1.setLabel([1])
   Tem0.setLabel([0])
   Temr=Tem0*1.0
   Temr.setLabel([-3])
   
   mps_list[2*q]=(Peps_bra[q]*Temr)*(SwapX)
   mps_list[2*q].permute([11,1, 10, -4, 2, -1], 3)

   Temr.setLabel([5])
   mps_list[2*q+1]=(((Peps_ket[q]*Temr)*SwapXX))
   mps_list[2*q+1].permute([ -4, 2, -1,0,6,4], 4)

  mps_boundry=MPSclass.MPS(2,2,len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=copy.copy(mps_list[i])

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  return mps_boundry

  ########################   Add next Layer   ################################

 if Location!=(N_x-1):
  mps_list=[None]*(N_x*2)

  ###############################   First-Layer   ###############################

  mps_boundry[0].setLabel([5,-3,6])

  Peps_bra[0].setLabel([-10,11,2,-3,-4])
  Tem1.setLabel([1])
  Tem0.setLabel([0])
  Temr=Tem0*1.0
  Temr.setLabel([5])

  mps_list[0]=(Temr*mps_boundry[0])*(Peps_bra[0])
  mps_list[0].permute([11,-10,2,6,-4],3)


  mps_boundry[1].setLabel([7,5,8])

  SwapXX=fermionicOPT(Sys,Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapXX.setLabel([5,6,3,-4])

  mps_list[1]=SwapXX*mps_boundry[1]
  mps_list[1].permute([7,-4,3,8,6],3)

###########################################

  for q in xrange(1, len(Peps_ket)):
###########################################
   mps_boundry[2*q].setLabel([5,-3,6])

   Peps_bra[q].setLabel([-10,11,2,-3,-4])


   mps_list[2*q]=(mps_boundry[2*q])*((Peps_bra[q]))

   mps_list[2*q].permute([5,11,-10,2,6,-4],4)


   mps_boundry[2*q+1].setLabel([7,5,8])
   SwapXX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapXX.setLabel([5,6,3,-4])


   mps_list[2*q+1]=SwapXX*mps_boundry[2*q+1]
   mps_list[2*q+1].permute([7,-4,3,8,6],3)

  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Middle", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t



  mps_boundry[0].setLabel([11,-10,2,6])
  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([1,10,-1,-10])
  Tem1.setLabel([1])

  mps_list[0]=(SwapX*Tem1)*mps_boundry[0]
  mps_list[0].permute([11,10,6,2,-1],2)

  mps_boundry[1].setLabel([7,3,8])
  Peps_ket[0].setLabel([0,-1,2,3,4])
  mps_list[1]=mps_boundry[1]*Peps_ket[0]
  mps_list[1].permute([7,2,-1,0,8,4],4)

######################### Next Absorption ####################################
  for  q  in  xrange(1, len(Peps_ket)):
   mps_boundry[2*q].setLabel([11,-10,2,6])
   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(1), Peps_bra[q].bond(0))
   SwapX.setLabel([1,10,-1,-10])
   Tem1.setLabel([1])

   mps_list[2*q]=SwapX*mps_boundry[2*q]
   mps_list[2*q].permute([11,1,10,6,2,-1],3)

   mps_boundry[2*q+1].setLabel([7,3,8])
   Peps_ket[q].setLabel([0,-1,2,3,4])
   mps_list[2*q+1]=mps_boundry[2*q+1]*Peps_ket[q]
   mps_list[2*q+1].permute([7,2,-1,0,8,4],4)


  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Last", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
 
  return mps_boundry


######@profile
def make_Env_singleLayer_down( PEPS_listten, Location, mps_boundry, d, chi_boundry, N_x,Sys):

 Peps_ket=[]
 for i in xrange( len(PEPS_listten) ):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)
  Peps_ket.append(A*1.0)


 Peps_bra=[]
 for i in xrange(len(PEPS_listten)):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)

  A.permute([1,2,3,4,5],3)
  A.setLabel([1,2,3,4,5])

  A_conj=A*1.0
  A_conj.transpose()
  A_conj.setLabel([-4,-5,-1,-2,3])
  A_conj.permute([-1,-2,3,-4,-5],5)
  Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
  Swap2.setLabel([-6,7,-4,-5])
  Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
  Swap3.setLabel([8,9,-1,-2])
  A_conj=(A_conj*Swap2)*Swap3
  A_conj.permute([8,9,3,-6,7],5)
  A_conj.setLabel([-1,-2,-3,-4,-5])

  Peps_bra.append(A_conj*1.0)




 chi_boundry_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi_boundry)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_boundry_list.append(dim)

 DX=0
 D=0
 if Location!=(N_x-1):
  DX=Peps_ket[0].bond(3).Qlist()
  D=Peps_bra[0].bond(3).Qlist()
 else:
  DX=Peps_ket[0].bond(0).Qlist()
  D=Peps_bra[0].bond(0).Qlist()

 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 #############Zero#######################
 bdi_1 = uni10.Bond( uni10.BD_IN, 1)
 bdo_1 = uni10.Bond( uni10.BD_IN, 1)
 Tem0=uni10.UniTensor([bdi_1])
 Tem0.identity()
 Tem0.setLabel([0])
 #print Tem0, 

 Tem1=uni10.UniTensor([bdi_1])
 Tem1.identity()
 Tem1.setLabel([1])
 #print Tem1

 Tem11=uni10.UniTensor([bdi_1])
 Tem11.identity()
 Tem11.setLabel([11])
 #print Tem11


 Tem10=uni10.UniTensor([bdi_1])
 Tem10.identity()
 Tem10.setLabel([10])
 #print Tem10

 Temr=uni10.UniTensor([bdi_1])
 Temr.identity()
 Temr.setLabel([10])

 Temrr=uni10.UniTensor([bdi_1])
 Temrr.identity()
 Temrr.setLabel([10])

 if Location==0:

  mps_list=[None]*(N_x*2)

  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([1,10,-1,-10])
  Peps_ket[0].setLabel([0,-1,2,3,4])
  Tem1.setLabel([1])
  Tem0.setLabel([0])
  
  mps_list[0]=(Peps_ket[0]*Tem0)*(SwapX*Tem1)
  mps_list[0].permute([10, 4, -10,2,3], 2)
  ##########################################

  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapX.setLabel([6,11,3,-4])
  Peps_bra[0].setLabel([-10,-11,2,-3,-4])
  Tem11.setLabel([-11])

  mps_list[1]=SwapX*(Peps_bra[0]*Tem11)
  mps_list[1].permute([-10,2,3,11,-3,6],4)
  #mps_list.append(results)

  ####################Middle##############
  for  q  in  xrange(1,len(Peps_ket)):

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(1), Peps_bra[q].bond(0))
   SwapX.setLabel([1,10,-1,-10])
   Peps_ket[q].setLabel([6,-1,2,3,4])
   Tem1.setLabel([1])
   
   mps_list[2*q]=((Peps_ket[q]))*(SwapX*Tem1)
   mps_list[2*q].permute([10,6, 4, -10,2,3], 3)
   ##########################################

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapX.setLabel([6,11,3,-4])
   Peps_bra[q].setLabel([-10,-11,2,-3,-4])
   Tem11.setLabel([-11])

   mps_list[2*q+1]=SwapX*(Peps_bra[q]*Tem11)
   mps_list[2*q+1].permute([-10,2,3,11,-3,6],4)


  mps_boundry=MPSclass.MPS(2,2,len(mps_list))
  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "0, Single", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()
  return mps_boundry

  ########################  Add next Layer  ###########

 if Location!=0:
  mps_list=[None]*(N_x*2)
  ###############################First-Layer###############################
  #print  mps_boundry[0].printDiagram()
  mps_boundry[0].setLabel([1,2,3])
  Tem1.setLabel([4])
  SwapX=fermionicOPT(Sys,mps_boundry[0].bond(1), Peps_bra[0].bond(0))
  SwapX.setLabel([5,6,2,4])
  Tem1.setLabel([6])


  mps_list[0]=(mps_boundry[0]*SwapX)*Tem1
  mps_list[0].permute([1,5,3,4],2)

  mps_boundry[1].setLabel([3,0,6])
  Peps_bra[0].setLabel([4,0,2,7,8])
  mps_list[1]=(mps_boundry[1]*Peps_bra[0])
  mps_list[1].permute([3,4,2,8,6,7],4)

  ###########################################
 

   ###############################################################################

  for q in xrange(1,len(Peps_ket)):
   ###########################################

    SwapX=fermionicOPT(Sys,mps_boundry[2*q].bond(1), Peps_bra[q].bond(0))
    SwapX.setLabel([5,6,2,4])
    mps_boundry[2*q].setLabel([1,2,3])

    mps_list[2*q]=mps_boundry[2*q]*SwapX
    mps_list[2*q].permute([1,6,5,3,4],3)



    mps_boundry[2*q+1].setLabel([3,0,6])
    Peps_bra[q].setLabel([4,0,2,7,8])

    mps_list[2*q+1]=mps_boundry[2*q+1]*Peps_bra[q]
    mps_list[2*q+1].permute([3,4,2,8,6,7],4)


  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Middle", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()


######################### Next Absorption ####################################

  mps_boundry[0].setLabel([6,5,4])
  #print "Peps_ket", Peps_ket[0].printDiagram()
  Peps_ket[0].setLabel([ 3, 5, 2, -4, -2 ])
  Tem1.setLabel([3])
  mps_list[0]=mps_boundry[0]*(Peps_ket[0]*Tem1)
  mps_list[0].permute([6,-2,4,2,-4],2)

  
  mps_boundry[1].setLabel([4,2,8,6])
  SwapX=fermionicOPT(Sys, Peps_ket[0].bond(3), mps_boundry[1].bond(2) )
  SwapX.setLabel([ -4, -2, 3, 8 ])
  mps_list[1]=mps_boundry[1]*SwapX
  mps_list[1].permute([4,2,3,-2,6,-4],4)

  #print "Hi", mps_list[0].printDiagram(), mps_list[1].printDiagram()

  for  q  in  xrange(1, len(Peps_ket)):

   mps_boundry[2*q].setLabel([6,5,4])
   Peps_ket[q].setLabel([ 3, 5, 2, -4, -2 ])
   mps_list[2*q]=mps_boundry[2*q]*Peps_ket[q]
   mps_list[2*q].permute([6,3,-2,4,2,-4],3)

   mps_boundry[2*q+1].setLabel([4,2,8,6])
   SwapX=fermionicOPT(Sys, Peps_ket[q].bond(3), mps_boundry[2*q+1].bond(2) )
   SwapX.setLabel( [ -4, -2, 3, 8 ] )
   mps_list[2*q+1]=mps_boundry[2*q+1]*SwapX
   mps_list[2*q+1].permute([4,2,3,-2,6,-4],4)




  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))
  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Last", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()

  return mps_boundry




######@profile
def make_Env_singleLayer_up( PEPS_listten, Location, mps_boundry, d, chi_boundry, N_x,Sys):

 Peps_ket=[]
 for i in xrange( len(PEPS_listten) ):
  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)
  Peps_ket.append(A*1.0)


 Peps_bra=[]
 for i in xrange(len(PEPS_listten)):

  A=PEPS_listten[i]*1.0
  A.setLabel([1,2,3,4,5])
  A.permute([1,2,3,4,5],5)

  A.permute([1,2,3,4,5],3)
  A.setLabel([1,2,3,4,5])

  A_conj=A*1.0
  A_conj.transpose()
  A_conj.setLabel([-4,-5,-1,-2,3])
  A_conj.permute([-1,-2,3,-4,-5],5)
  Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
  Swap2.setLabel([-6,7,-4,-5])
  Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
  Swap3.setLabel([8,9,-1,-2])
  A_conj=(A_conj*Swap2)*Swap3
  A_conj.permute([8,9,3,-6,7],5)
  A_conj.setLabel([-1,-2,-3,-4,-5])
  Peps_bra.append(A_conj*1.0)

 chi_boundry_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi_boundry)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_boundry_list.append(dim)

 DX=0
 D=0
 if Location!=(N_x-1):
  DX=Peps_ket[0].bond(3).Qlist()
  D=Peps_bra[0].bond(3).Qlist()
 else:
  DX=Peps_ket[0].bond(0).Qlist()
  D=Peps_bra[0].bond(0).Qlist()

 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdiX=uni10.Bond(uni10.BD_IN, DX)
 bdoX=uni10.Bond(uni10.BD_OUT, DX)

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 #############Zero#######################
 bdi_1 = uni10.Bond( uni10.BD_IN, 1)
 bdo_1 = uni10.Bond( uni10.BD_IN, 1)
 Tem0=uni10.UniTensor([bdi_1])
 Tem0.identity()
 Tem0.setLabel([0])
 #print Tem0, 

 Tem1=uni10.UniTensor([bdi_1])
 Tem1.identity()
 Tem1.setLabel([1])
 #print Tem1

 Tem11=uni10.UniTensor([bdi_1])
 Tem11.identity()
 Tem11.setLabel([11])
 #print Tem11

 Tem10=uni10.UniTensor([bdi_1])
 Tem10.identity()
 Tem10.setLabel([10])
 #print Tem10

 Temr=uni10.UniTensor([bdi_1])
 Temr.identity()
 Temr.setLabel([10])

 Temrr=uni10.UniTensor([bdi_1])
 Temrr.identity()
 Temrr.setLabel([10])

 if Location==N_x-1:

  mps_list=[None]*(N_x*2)

  SwapX=fermionicOPT(Sys, Peps_ket[0].bond(1), Peps_bra[0].bond(0) )
  SwapX.setLabel( [1,10,-1,-10] )
  Peps_ket[0].setLabel( [0,-1,2,3,4] )
  Tem1.setLabel( [10] )
  Temr.setLabel( [4] )

  mps_list[0]=(Peps_ket[0]*Temr)*(SwapX*Tem1)
  mps_list[0].permute( [0, 1, -10, 2 , 3] , 2 )

  ##########################################

  SwapX=fermionicOPT(Sys,Peps_ket[0].bond(3), Peps_bra[0].bond(4))
  SwapX.setLabel([6,11,3,-4])
  Peps_bra[0].setLabel([-10,-11,2,-3,-4])
  Tem11.setLabel([11])

  mps_list[1]=(SwapX*Tem11)*Peps_bra[0]
  mps_list[1].permute([-10,2,3,-11,-3,6],4)
  #mps_list.append(results)

  ####################Middle##############
  for  q  in  xrange(1,len(Peps_ket)):

   SwapX=fermionicOPT(Sys, Peps_ket[q].bond(1), Peps_bra[q].bond(0) )
   SwapX.setLabel([1,10,-1,-10])
   Peps_ket[q].setLabel([6,-1,2,3,4])
   Tem1.setLabel([4])
   
   mps_list[2*q]=((Peps_ket[q]))*(SwapX*Tem1)
   mps_list[2*q].permute([ 10, 6, 1, -10,2, 3], 3)
   ##########################################

   SwapX=fermionicOPT(Sys,Peps_ket[q].bond(3), Peps_bra[q].bond(4))
   SwapX.setLabel([6,11,3,-4])
   Peps_bra[q].setLabel([-10,-11,2,-3,-4])
   Tem11.setLabel([11])

   mps_list[2*q+1]=(SwapX*Tem11)*(Peps_bra[q])
   mps_list[2*q+1].permute([-10,2,3,-11,-3,6],4)


  mps_boundry=MPSclass.MPS(2,2,len(mps_list))
  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "0, Single", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  #print  mps_boundry[0].printDiagram(), mps_boundry[1].printDiagram(), mps_boundry[2].printDiagram(),mps_boundry[4].printDiagram()
  return mps_boundry
 ######################### Next Absorption ####################################
 if Location!=(N_x-1):

  mps_list=[None]*(N_x*2)
  mps_boundry[0].setLabel([6,5,4])
  Peps_ket[0].setLabel([ 3, -2, 2, -4, 5 ])
  Tem1.setLabel([3])
  mps_list[0]=mps_boundry[0]*(Peps_ket[0]*Tem1)
  mps_list[0].permute([6,-2,2,-4,4],3)

  mps_boundry[1].setLabel([4,-2,6])
  SwapX=fermionicOPT(Sys, Peps_ket[0].bond(3), mps_boundry[1].bond(1) )
  SwapX.setLabel([ -4, 8, 3, -2 ])

  mps_list[1]=mps_boundry[1]*SwapX
  mps_list[1].permute([3,4,8,-4,6],3)

  for  q  in  xrange(1, len(Peps_ket)):

   mps_boundry[2*q].setLabel([6,5,4])
   Peps_ket[q].setLabel( [ 3, -2, 2, -4, 5 ] )
   mps_list[2*q]=mps_boundry[2*q]*Peps_ket[q]
   mps_list[2*q].permute( [3,6,-2,2,-4,4], 4 )

   mps_boundry[2*q+1].setLabel( [4,-2,6] )
   SwapX=fermionicOPT(Sys, Peps_ket[q].bond(3), mps_boundry[2*q+1].bond(1) )
   SwapX.setLabel( [ -4, 8, 3, -2 ] )
   mps_list[2*q+1]=mps_boundry[2*q+1]*SwapX
   mps_list[2*q+1].permute([3,4,8,-4,6],3)

  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))
  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

#   A_t=mps_boundry.norm()
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))

########################  Add next Layer  ###########

###############################First-Layer###############################
  mps_boundry[0].setLabel([3,-2,2,-4])
  SwapX=fermionicOPT(Sys, mps_boundry[0].bond(1), Peps_bra[0].bond(0) )
  SwapX.setLabel([20,6,-2,4])
  Tem1.setLabel([6])
  mps_list[0]=(mps_boundry[0]*SwapX)*Tem1
  mps_list[0].permute([3,20,4,2,-4],2)

  mps_boundry[1].setLabel([3,8,-4])
  Peps_bra[0].setLabel([4,0,2,7,8])
  mps_list[1]=(mps_boundry[1]*Peps_bra[0])
  mps_list[1].permute([4,2,3,0,7,-4],4)
###########################################
###########################################

  for q in xrange(1,len(Peps_ket)):
   ###########################################
   mps_boundry[2*q].setLabel([3,-2,2,-4])
   SwapX=fermionicOPT(Sys, mps_boundry[2*q].bond(1), Peps_bra[q].bond(0) )
   SwapX.setLabel([20,6,-2,4])
   mps_list[2*q]=(mps_boundry[2*q]*SwapX)
   mps_list[2*q].permute([6,3,20,4,2,-4],3)

   mps_boundry[2*q+1].setLabel([3,8,-4])
   Peps_bra[q].setLabel([4,0,2,7,8])
   mps_list[2*q+1]=(mps_boundry[2*q+1]*Peps_bra[q])
   mps_list[2*q+1].permute([4,2,3,0,7,-4],4)

  mps_boundry=MPSclass.MPS( 2, 2, len(mps_list))

  for i in xrange(len(mps_list)):
    mps_boundry[i]=mps_list[i]*1.0

  #A_t=mps_boundry.norm()
  #mps_boundry_t=mps_boundry*1.0
  mps_boundry=mps_boundry.appSVD(sum(chi_boundry_list))
  #print "Middle", mps_boundry.norm(), A_t , (A_t-mps_boundry.norm()) / A_t
  return mps_boundry




def   setTruncation1(theta, chi):
    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge1(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0),theta.bond(1),theta.bond(2),bdo_mid])
    GB.assign([bdi_mid,theta.bond(3),theta.bond(4),theta.bond(5)])
    LA.assign([bdi_mid, b_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim)  )


def   sv_merge1(svs, bidxs, bidx, sv_mat, chi, len_qn):
    if(len(svs)):
        length = len(svs) + sv_mat.elemNum()
        length = length if length < chi else chi
        ori_svs = svs
        ori_bidxs = bidxs
        svs = [0] * length
        bidxs = [0] * length
        svs = []
        bidxs = []
        cnt  = 0
        cur1 = 0
        cur2 = 0
        while cnt < length:
            if(cur1 < len(ori_svs)) and cur2 < sv_mat.elemNum():
                if ori_svs[cur1] >= sv_mat[cur2]:
                    if (ori_svs[cur1] > -0.01):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                    cur1 += 1
                else:
                    if (sv_mat[cur2] > -0.01):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                    cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > -0.01):
                     svs.append(sv_mat[cur2]) 
                     bidxs.append(bidx) 
                    cur2 += 1
                break
            else:
                for i in xrange(cur1, len(ori_svs)):
                 svs.append(ori_svs[i])
                 bidxs.append(ori_bidxs[i]) 
                break
            cnt += 1
    else:
       if (len_qn is 1):
        bidxs = [bidx] * chi  
        svs = [sv_mat[i] for i in xrange(chi)]
       elif (sv_mat[0] > -0.01):
        bidxs = [bidx] * sv_mat.elemNum()
        svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
       else: bidxs = [bidx];  svs = [sv_mat[0]];  
    return svs, bidxs


def setTruncation3(theta, chi):
    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge1(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0),theta.bond(1), bdo_mid])
    GB.assign([bdi_mid, theta.bond(2), theta.bond(3)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim)  )
    return GA, GB, LA

#####@profile
def C_i(dim):
  Mat=uni10.Matrix(dim, dim)
  Mat.set_zero()
  for i in xrange(dim):
   for j in xrange(dim):
     if j==i+1:
      Mat[i*dim+j]=j**(0.5) 
  iden_m=Mat*1.0
  iden_m.identity()
  Mat2=Mat*1.0
  Mat2.transpose()
  commutator=Mat2*Mat+(-1.0)*Mat*Mat2
  #print Mat, Mat2 ,commutator, Mat2*Mat, Mat*Mat2
  return Mat , Mat2, iden_m






def V_ext( i, j, L, N_x):
 perid=1.0
 return 0.0
 return ((((i+1)*(L/(N_x+1)))-(L/2.0))**2)   +   ((((j+1)*(L/(N_x+1)))-(L/2.0))**2)
 #return +((np.cos((perid*np.pi*(i+1))/(N_x+1))**2))  +  ((np.cos((perid*np.pi*(j+1))/(N_x+1))**2))


def Mat_np_to_Uni(Mat_np):
 d0=np.size(Mat_np,0)
 d1=np.size(Mat_np,1)
 Mat_uni=uni10.Matrix(d0,d1)
 for i in xrange(d0):
  for j in xrange(d1):
   Mat_uni[i*d1+j]=Mat_np[i,j]
 return  Mat_uni


 
def Mat_uni_to_np(Mat_uni):
 dim0=int(Mat_uni.row())
 dim1=int(Mat_uni.col())
 Mat_np=np.zeros((dim0,dim1))
 for i in xrange(dim0):
  for j in xrange(dim1):
   Mat_np[i,j]=Mat_uni[i*dim1+j]
 return  Mat_np





#def Ham_mat( J, g1, g2, J1, J2, h1, h2):
# dim=4
# template = np.zeros([dim, dim, dim, dim])
# for idx in np.ndindex(template.shape):
#  if idx[0]==0 and idx[1]==1 and idx[2]==2 and idx[3]==3: template[idx] = +J
#  if idx[0]==0 and idx[1]==1 and idx[2]==3 and idx[3]==2: template[idx] = -J
#  if idx[0]==0 and idx[1]==1 and idx[2]==0 and idx[3]==1: template[idx] = g2+2*J2
#  if idx[0]==0 and idx[1]==2 and idx[2]==2 and idx[3]==0: template[idx] = +J
#  if idx[0]==0 and idx[1]==2 and idx[2]==0 and idx[3]==2: template[idx] = +J2+h2
#  if idx[0]==0 and idx[1]==3 and idx[2]==3 and idx[3]==0: template[idx] = +J
#  if idx[0]==0 and idx[1]==3 and idx[2]==0 and idx[3]==3: template[idx] = +J2-h2
#  if idx[0]==1 and idx[1]==0 and idx[2]==3 and idx[3]==2: template[idx] = -J
#  if idx[0]==1 and idx[1]==0 and idx[2]==2 and idx[3]==3: template[idx] = +J
#  if idx[0]==1 and idx[1]==0 and idx[2]==1 and idx[3]==0: template[idx] = g1+2*J1
#  if idx[0]==1 and idx[1]==1 and idx[2]==1 and idx[3]==1: template[idx] =g1+g2+2.0*J1+2.0*J2
#  if idx[0]==1 and idx[1]==2 and idx[2]==2 and idx[3]==1: template[idx] =-J
#  if idx[0]==1 and idx[1]==2 and idx[2]==1 and idx[3]==2: template[idx] =g1+2*J1+J2+h2
#  if idx[0]==1 and idx[1]==3 and idx[2]==3 and idx[3]==1: template[idx] =-J
#  if idx[0]==1 and idx[1]==3 and idx[2]==1 and idx[3]==3: template[idx] =g1+2*J1+J2-h2
#  if idx[0]==2 and idx[1]==0 and idx[2]==0 and idx[3]==2: template[idx] =J
#  if idx[0]==2 and idx[1]==0 and idx[2]==2 and idx[3]==0: template[idx] =+J1+h1
#  if idx[0]==2 and idx[1]==1 and idx[2]==1 and idx[3]==2: template[idx] =-J
#  if idx[0]==2 and idx[1]==1 and idx[2]==2 and idx[3]==1: template[idx] =g2+J1+2*J2+h1
#  if idx[0]==2 and idx[1]==2 and idx[2]==2 and idx[3]==2: template[idx] =J1+J2+h1+h2
#  if idx[0]==2 and idx[1]==3 and idx[2]==0 and idx[3]==1: template[idx] =J
#  if idx[0]==2 and idx[1]==3 and idx[2]==1 and idx[3]==0: template[idx] =J
#  if idx[0]==2 and idx[1]==3 and idx[2]==2 and idx[3]==3: template[idx] =J1+J2+h1-h2
#  if idx[0]==3 and idx[1]==0 and idx[2]==0 and idx[3]==3: template[idx] =J
#  if idx[0]==3 and idx[1]==0 and idx[2]==3 and idx[3]==0: template[idx] =+J1-h1
#  if idx[0]==3 and idx[1]==1 and idx[2]==1 and idx[3]==3: template[idx] =-J
#  if idx[0]==3 and idx[1]==1 and idx[2]==3 and idx[3]==1: template[idx] =g2+J1+2*J2-h1
#  if idx[0]==3 and idx[1]==2 and idx[2]==1 and idx[3]==0: template[idx] =-J
#  if idx[0]==3 and idx[1]==2 and idx[2]==0 and idx[3]==1: template[idx] =-J
#  if idx[0]==3 and idx[1]==2 and idx[2]==3 and idx[3]==2: template[idx] =J1+J2-h1+h2
#  if idx[0]==3 and idx[1]==3 and idx[2]==3 and idx[3]==3: template[idx] =J1+J2-h1-h2


# #print template

# return template.reshape(-1)



def Ham_mat( J, g1, g2, J1, J2, h1, h2):
  dim=4
  Mat=uni10.Matrix(dim*dim, dim*dim)
  Mat.set_zero()
  for i in xrange(dim):
   for j in xrange(dim):
    for m in xrange(dim):
     for n in xrange(dim):
 
      if i==0 and j==1  and m==2 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J
      if i==2 and j==3  and m==0 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J
      if i==0 and j==1  and m==3 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==3 and j==2 and m==0 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J


      if i==1 and j==0  and m==3 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==3 and j==2  and m==1 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==1 and j==0  and m==2 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J
      if i==2 and j==3  and m==1 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J

      if i==0 and j==2 and m==2 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] =+J
      if i==2 and j==0  and m==0 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J

      if i==0 and j==3  and m==3 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J
      if i==3 and j==0  and m==0 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J

      if i==1 and j==2  and m==2 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==2 and j==1  and m==1 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J

      if i==1 and j==3  and m==3 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J
      if i==3 and j==1  and m==1 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = -J


      if i==0 and j==1  and m==0 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g2+(2*J2)
      if i==0 and j==2 and m==0 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] =+J2+h2
      if i==0 and j==3  and m==0 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J2-h2
      if i==1 and j==0  and m==1 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g1+(2*J1)
      if i==1 and j==1  and m==1 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g1+g2+(2.0*J1)+(2.0*J2)
      if i==1 and j==2  and m==1 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g1+(2*J1)+J2+h2
      if i==1 and j==3  and m==1 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g1+(2*J1)+J2-h2
      if i==2 and j==0  and m==2 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J1+h1
      if i==2 and j==1  and m==2 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g2+J1+(2*J2)+h1
      if i==2 and j==2  and m==2 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J1+J2+h1+h2
      if i==2 and j==3  and m==2 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J1+J2+h1-h2
      if i==3 and j==0  and m==3 and n==0:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = +J1-h1
      if i==3 and j==1  and m==3 and n==1:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = g2+J1+(2*J2)-h1
      if i==3 and j==2 and m==3 and n==2:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J1+J2-h1+h2
      if i==3 and j==3 and m==3 and n==3:
        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J1+J2-h1-h2



  #print Mat
  #Mat_np=Mat_uni_to_np(Mat)
  #print Mat_np
  return Mat




# 
# def Ham_mat( J, g1, g2, J1, J2, h1, h2):
#  dim=2
#  Mat=uni10.Matrix(dim*dim, dim*dim)
#  Mat.set_zero()
#  for i in xrange(dim):
#   for j in xrange(dim):
#    for m in xrange(dim):
#     for n in xrange(dim):
# 
#      if i==0 and j==1  and m==1 and n==0:
#        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J
#      if i==1 and j==0  and m==0 and n==1:
#        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J
# 
#      if i==0 and j==1 and m==0 and n==1:
#        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] =J2
#      if i==1 and j==0 and m==1 and n==0:
#        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] =J1
# 
#      if i==1 and j==1  and m==1 and n==1:
#        Mat[m*dim*dim*dim+n*dim*dim+i*dim+j] = J1+J2
# 
#  Mat_np=Mat_uni_to_np(Mat)
#  return Mat



# def Ham_mat( J, g1, g2, J1, J2, h1, h2):
#  dim=2
#  c_iu, c_iu_dag, iden=C_i(dim)
#  ham=J*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
#  #print c_iu, c_iu_dag
#  ham=g1*uni10.otimes(c_iu_dag*c_iu,iden)+ham
#  ham=g2*uni10.otimes(iden,c_iu_dag*c_iu)+ham
#  ham=-1.0*h2*uni10.otimes(iden,c_iu_dag*c_iu)+ham
#  ham=-1.0*h1*uni10.otimes(c_iu_dag*c_iu,iden)+ham
#  ham=1.0*J2*uni10.otimes(iden,c_iu_dag*c_iu)+ham
#  ham=1.0*J1*uni10.otimes(c_iu_dag*c_iu,iden)+ham
# 
#  #print ham
#  return ham



def make_H_col( N_x, d_phys, h, Model, L):

 

 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*(N_x-1) 



 if Model[0]=="ITF_Z2"  or Model[0]=="ITF":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
  sx = matSx()
  sy = matSy()
  sz = matSz()
  iden = matIden()
  ham=uni10.otimes(sx,sx)
  for i in xrange(N_x):
   for j in xrange(N_x-1):
    ham1=0
    if j==0:
     ham1=(uni10.otimes(iden,sz)+2.0*uni10.otimes(sz,iden))
    elif j==N_x-2:
     ham1=(2.0*uni10.otimes(iden,sz)+uni10.otimes(sz,iden))
    else:
     ham1=(uni10.otimes(iden,sz)+uni10.otimes(sz,iden))

    ham_m=h[0]*(-1)*ham+h[1]*ham1*(-0.25)
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0



 if Model[0]=="Heis_Z2"  or Model[0]=="Heis":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
  sx = matSx()
  sy = matSy()
  sz = matSz()
  iden = matIden()
  ham=(uni10.otimes(sx,sx)+(-1.0)*uni10.otimes(sy,sy)+uni10.otimes(sz,sz))*(0.25)
  for i in xrange(N_x):
   for j in xrange(N_x-1):
    ham_m=ham*1.0
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0



 if Model[0]=="FFI_Z2" or Model[0]=="FFI":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "FFI_Z2")
  c_iu, c_iu_dag, iden=C_i(len(d_phys))
  c_id, c_id_dag, iden=C_i(len(d_phys))
  c_id=c_id*0
  c_id_dag=c_id_dag*0
  for i in xrange(N_x):
   for j in xrange(N_x-1):
    if j==0:
     ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
     ham=ham+h[1]*uni10.otimes(c_iu_dag*c_iu,c_iu_dag*c_iu)
    elif j==N_x-2:
     ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
     ham=ham+h[1]*uni10.otimes(c_iu_dag*c_iu,c_iu_dag*c_iu)
    else:
     ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
     ham=ham+h[1]*uni10.otimes(c_iu_dag*c_iu,c_iu_dag*c_iu)

    ham_m=ham*1.0
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0




 if Model[0]=="Fer_Z2" or Model[0]=="Fer_U1" or Model[0]=="Fer" or Model[0]=="Fer_BOS" or Model[0]=="Fer_BOS_Z2":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")

  ham=Ham_mat( 1.0, 0., 0., 0., 0., 0., 0.)

  for i in xrange(N_x):
   for j in xrange(N_x-1):
    ham1=0
    if j==0:
     ham1=Ham_mat( 0.0, 0.0, 0.0, 2.0, 1.0, 0., 0.)
    elif j==N_x-2 :
     ham1=Ham_mat( 0.0, 0.0, 0.0, 1.0, 2.0, 0., 0.)
    else:
     ham1=Ham_mat( 0.0, 0.0, 0.0, 1.0, 1.0, 0., 0.)
    ham_m=h[0]*ham+h[2]*ham1*0.25+h[1]*ham1*0.25
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0



#############################################################################
  for i in xrange(N_x):
   for j in xrange(N_x-1):
    ham1=0
    if j==0:
     ham1=Ham_mat( 0.0, 2.0, 1.0, 0., 0.0, 0., 0.)
    elif j==N_x-2 :
     ham1=Ham_mat( 0.0, 1.0, 2.0, 0.0, 0., 0., 0.)
    else:
     ham1=Ham_mat( 0.0, 1.0, 1.0, 0, 0, 0., 0.)

    ham_m=h[3]*ham1*0.25
    H.setRawElem(ham_m)
    H_list[i][j]=H+H_list[i][j]

########################################################################
  for i in xrange(N_x):
   for j in xrange(N_x-1):
    ham1=0
    if j==0:
     ham1=Ham_mat( 0.0, 0., 0.0, 0., 0.0, 2.0, 1.)
    elif j==N_x-2 :
     ham1=Ham_mat( 0.0, 0.0, 0., 0.0, 0., 1., 2.0)
    else:
     ham1=Ham_mat( 0.0, 0.0, 0.0, 0, 0, 1., 1.)
    ham_m=h[4]*ham1*0.25
    H.setRawElem(ham_m)
    H_list[i][j]=H+H_list[i][j]




 return H_list


def make_H_row(N_x, d_phys,h, Model,L):

 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heis")
 #H_list=[None]*(N_x-1)

 H_list=[None]*(N_x-1)
 for i in xrange(N_x-1):
  H_list[i]=[None]*N_x 


 if Model[0]=="Heis_Z2"  or Model[0]=="Heis":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
  sx = matSx()
  sy = matSy()
  sz = matSz()
  iden = matIden()
  ham=(uni10.otimes(sx,sx)+(-1.0)*uni10.otimes(sy,sy)+uni10.otimes(sz,sz))*(0.25)
  for i in xrange(N_x-1):
   for j in xrange(N_x):
    ham_m=ham*1.0
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0



 if Model[0]=="ITF_Z2"  or Model[0]=="ITF" :
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
  sx = matSx()
  sy = matSy()
  sz = matSz()
  iden = matIden()
  ham=uni10.otimes(sx,sx)
  for i in xrange(N_x-1):
   for j in xrange(N_x):
    ham1=0
    if i==0:
     ham1=(uni10.otimes(iden,sz)+2.0*uni10.otimes(sz,iden))
    elif i==N_x-2:
     ham1=(2.0*uni10.otimes(iden,sz)+uni10.otimes(sz,iden))
    else:
     ham1=(uni10.otimes(iden,sz)+uni10.otimes(sz,iden))

    ham_m=h[0]*(-1)*ham+h[1]*ham1*(-0.25)
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0

##########################################

 if Model[0]=="FFI_Z2" or Model[0]=="FFI":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "FFI_Z2")
  c_iu, c_iu_dag, iden=C_i(len(d_phys))
  c_id, c_id_dag, iden=C_i(len(d_phys))
  c_id=c_id*0
  c_id_dag=c_id_dag*0
  for i in xrange(N_x-1):
   for j in xrange(N_x):
    if i==0:
     ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
     ham=ham+h[1]*uni10.otimes(c_iu_dag*c_iu,c_iu_dag*c_iu)
    elif i==N_x-2:
     ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
     ham=ham+h[1]*uni10.otimes(c_iu_dag*c_iu,c_iu_dag*c_iu)
    else:
     ham=h[0]*(uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu))
     ham=ham+h[1]*uni10.otimes(c_iu_dag*c_iu,c_iu_dag*c_iu)

    ham_m=ham*1.0
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0



 if Model[0]=="Fer_Z2"  or Model[0]=="Fer_U1" or Model[0]=="Fer" or Model[0]=="Fer_BOS" or Model[0]=="Fer_BOS_Z2":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
  ham=Ham_mat( 1.0, 0, 0, 0, 0, 0., 0.)

  for i in xrange(N_x-1):
   for j in xrange(N_x):
    ham1=0
    if i==0:
     ham1=Ham_mat( 0.0, 0.0, 0.0, 2.0, 1.0, 0., 0.)
    elif i==N_x-2 :
     ham1=Ham_mat( 0.0, 0.0, 0.0, 1.0, 2.0, 0., 0.)
    else:
     ham1=Ham_mat( 0.0, 0.0, 0.0, 1.0, 1.0, 0., 0.)

    ham_m=h[0]*ham+h[1]*ham1*0.25+h[2]*ham1*0.25
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0


###########################################################################


  for i in xrange(N_x-1):
   for j in xrange(N_x):
    ham1=0
    if i==0:
     ham1=Ham_mat( 0.0, 2.0, 1.0, 0.0, 0.0, 0., 0.)
    elif i==N_x-2 :
     ham1=Ham_mat( 0.0, 1.0, 2.0, 0.0, 0.0, 0., 0.)
    else:
     ham1=Ham_mat( 0.0, 1.0, 1.0, 0.0, 0.0, 0., 0.)

    ham_m=h[3]*ham1*0.25
    H.setRawElem(ham_m)
    H_list[i][j]=H+H_list[i][j]


##############################################################################
  for i in xrange(N_x-1):
   for j in xrange(N_x):
    ham1=0
    if i==0:
     ham1=Ham_mat( 0.0, 0.0, 0.0, 0.0, 0.0, 2., 1.)
    elif i==N_x-2 :
     ham1=Ham_mat( 0.0, 0.0, 0.0, 0.0, 0.0, 1., 2.)
    else:
     ham1=Ham_mat( 0.0, 0.0, 0.0, 0.0, 0.0, 1., 1.)

    ham_m=h[4]*ham1*0.25
    H.setRawElem(ham_m)
    H_list[i][j]=H+H_list[i][j]


 return H_list


def  make_H_long(N_x, d_phys,h, Model):

 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heis")
 #H_list=[None]*(N_x-1)

 
 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*N_x 


###########################################################
 if Model[0]=="ITF_Z2" or Model[0]=="ITF" or Model[0]=="FFI_Z2" or Model[0]=="FFI" or Model[0]=="Heis" or Model[0]=="Heis_Z2":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [ bdi, bdi, bdo, bdo], "Heis")
  sx = matSx()
  sy = matSy()
  sz = matSz()
  iden = matIden()
  ham=uni10.otimes(sx,sx)*(-1)
  for i in xrange(N_x):
   for j in xrange(N_x):
    ham1=0
    if i==0:
     ham1=h[1]*(uni10.otimes(iden,sz)+2.0*uni10.otimes(sz,iden))
    elif i==N_x-2:
     ham1=h[1]*(2.0*uni10.otimes(iden,sz)+uni10.otimes(sz,iden))
    else:
     ham1=h[1]*(uni10.otimes(iden,sz)+uni10.otimes(sz,iden))

    ham_m=0*ham+h[1]*ham1*(-0.0)
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0

########################################################



 if Model[0]=="Fer_Z2"  or Model[0]=="Fer_U1" or Model[0]=="Fer" or Model[0]=="FFI_Z2" or Model[0]=="Fer_BOS" or Model[0]=="Fer_BOS_Z2":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
#  c_iu, c_iu_dag, iden=C_i(len(d_phys))
#  c_id, c_id_dag, iden=C_i(len(d_phys))
#  c_id=0*c_id
#  c_id_dag=0*c_id_dag


#  c_iu, c_iu_dag, iden=C_i_spinUP()
#  c_id, c_id_dag, iden=C_i_spinDOWN()




#  ham=uni10.otimes(c_iu,c_iu_dag)+uni10.otimes(c_iu_dag,c_iu)
#  ham=ham+uni10.otimes(c_id,c_id_dag)+uni10.otimes(c_id_dag,c_id)
  ham=Ham_mat( 1.0, 0, 0, 0, 0, 0., 0.)

  
  for i in xrange(N_x):
   for j in xrange(N_x):
    ham_m=ham*1.0
    H.setRawElem(ham_m)
    H_list[i][j]=H*1.0*h[2]
    H_list[i][j]=H*0.0*h[2]

 #print H_list
 return H_list






def make_Sz_col( N_x, d_phys, h, Model):

 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*(N_x-1) 

 bdi = uni10.Bond( uni10.BD_IN, d_phys)
 bdo = uni10.Bond( uni10.BD_OUT, d_phys)
 H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")




 for i in xrange(N_x):
  for j in xrange(N_x-1):
   if j==0:
    ham=Ham_mat( 0.0,  0.0, 0.0, 0., 0., 1.0*2., 1.0)

   elif j==N_x-2 :
    ham=Ham_mat( 0.0,  0.0, 0.0, 0., 0.,1.0, 1.0*2.)

   else:
    ham=Ham_mat( 0.0,  0.0, 0.0, 0., 0., 1.0, 1.0)

   H.setRawElem(ham*0.25)
   H_list[i][j]=H*1.0

 return H_list








def  make_Sz_row(N_x, d_phys,h, Model):

 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heis")

 H_list=[None]*(N_x-1)
 for i in xrange(N_x-1):
  H_list[i]=[None]*N_x 

 bdi = uni10.Bond( uni10.BD_IN, d_phys)
 bdo = uni10.Bond( uni10.BD_OUT, d_phys)
 H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")

 for i in xrange(N_x-1):
  for j in xrange(N_x):

   if i==0 :
    ham=Ham_mat( 0.0,  0.0, 0.0,  0., 0., 1.0*2., 1.0)

   elif i==N_x-2:
    ham=Ham_mat( 0.0,  0.0, 0.0, 0., 0.,1.0, 1.0*2.)

   else:
    ham=Ham_mat( 0.0,  0.0, 0.0, 0., 0.,1.0, 1.0)

   H.setRawElem(ham*0.25)
   H_list[i][j]=H*1.0


 return H_list





def make_Magz_col( N_x, d_phys, h, Model):

 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*(N_x) 

 bdi = uni10.Bond( uni10.BD_IN, d_phys)
 bdo = uni10.Bond( uni10.BD_OUT, d_phys)
 H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")

 for i in xrange(N_x/2):
  for j in xrange(N_x/2):

   ham=Ham_mat( 0.0,  0.0, 0.0, 0.0, 0.0 , 1., 0.)
   H.setRawElem(ham)
   H_list[2*i][2*j]=H*1.0
   ham=Ham_mat( 0.0,  0.0, 0.0, 0.0, 0.0 , 1., 0.)
   H.setRawElem(ham)
   H_list[2*i+1][2*j]=H*1.0

   ham=Ham_mat( 0.0,  0.0, 0.0 ,0.0, 0.0, 0., 1.)
   H.setRawElem(ham)
   H_list[2*i][2*j+1]=H*1.0

   ham=Ham_mat( 0.0, 0.0, 0.0,0.0, 0.0 , 0., 1.)
   H.setRawElem(ham)
   H_list[2*i+1][2*j+1]=H*1.0


 return H_list








def make_Magz_col_direct( N_x, d_phys, h, Model):

 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*(N_x) 

 bdi = uni10.Bond( uni10.BD_IN, d_phys)
 bdo = uni10.Bond( uni10.BD_OUT, d_phys)
 H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")

 for i in xrange(N_x):
  for j in xrange(N_x):
   if j==N_x-1:
    ham=Ham_mat( 0.0,  0.0, 0.0, 0.0, 0.0 , 0., 1.)
    H.setRawElem(ham)
    H_list[i][j]=H*1.0
   
   else:
    ham=Ham_mat( 0.0,  0.0, 0.0, 0.0, 0.0 , 1., 0.)
    H.setRawElem(ham)
    H_list[i][j]=H*1.0

 return H_list







def make_N_col( N_x, d_phys, h, Model):

 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*(N_x-1) 

 bdi = uni10.Bond( uni10.BD_IN, d_phys)
 bdo = uni10.Bond( uni10.BD_OUT, d_phys)
 H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
 c_iu, c_iu_dag, iden=C_i(len(d_phys))
 c_id, c_id_dag, iden=C_i(len(d_phys))
 c_id=0*c_id
 c_id_dag=0*c_id_dag


# c_iu, c_iu_dag, iden=C_i_spinUP()
# c_id, c_id_dag, iden=C_i_spinDOWN()

 for i in xrange(N_x):
  for j in xrange(N_x-1):
   if j==0:
    ham=Ham_mat( 0.0,  0.0, 0.0,1.0*2., 1.0 , 0., 0.)

   elif j==N_x-2 :
    ham=Ham_mat( 0.0,  0.0, 0.0,1.0, 1.0*2. , 0., 0.)

   else:
    ham=Ham_mat( 0.0,  0.0, 0.0,1.0, 1.0 , 0., 0.)

   H.setRawElem(ham*0.25)
   H_list[i][j]=H*1.0

 return H_list

def make_particle_col( N_x, d_phys, h, Model):

 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*(N_x) 

 bdi = uni10.Bond( uni10.BD_IN, d_phys)
 bdo = uni10.Bond( uni10.BD_OUT, d_phys)
 H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
 c_iu, c_iu_dag, iden=C_i(len(d_phys))
 c_id, c_id_dag, iden=C_i(len(d_phys))
 c_id=0*c_id
 c_id_dag=0*c_id_dag

 for i in xrange(N_x/2):
  for j in xrange(N_x/2):

   ham=Ham_mat( 0.0,  0.0, 0.0,1.0, 0.0 , 0., 0.)

   H.setRawElem(ham)
   H_list[2*i][2*j]=H*1.0

   ham=Ham_mat( 0.0,  0.0, 0.0,1.0, 0.0 , 0., 0.)

   H.setRawElem(ham)
   H_list[2*i+1][2*j]=H*1.0

   ham=Ham_mat( 0.0,  0.0, 0.0 ,0.0, 1.0, 0., 0.)
   H.setRawElem(ham)
   H_list[2*i][2*j+1]=H*1.0

   ham=Ham_mat( 0.0, 0.0, 0.0,0.0, 1.0 , 0., 0.)
   H.setRawElem(ham)
   H_list[2*i+1][2*j+1]=H*1.0


 return H_list








def make_particle_col_direct( N_x, d_phys, h, Model):

 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*(N_x) 

 bdi = uni10.Bond( uni10.BD_IN, d_phys)
 bdo = uni10.Bond( uni10.BD_OUT, d_phys)
 H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")

 for i in xrange(N_x):
  for j in xrange(N_x):
   if j==N_x-1:
    ham=Ham_mat( 0.0,  0.0, 0.0, 0.0, 1.0 , 0., 0.)
    H.setRawElem(ham)
    H_list[i][j]=H*1.0
   else:
    ham=Ham_mat( 0.0,  0.0, 0.0, 1.0, 0.0 , 0., 0.)
    H.setRawElem(ham)
    H_list[i][j]=H*1.0
 return H_list



def  make_N_row(N_x, d_phys,h, Model):

 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heis")

 H_list=[None]*(N_x-1)
 for i in xrange(N_x-1):
  H_list[i]=[None]*N_x 

 bdi = uni10.Bond( uni10.BD_IN, d_phys)
 bdo = uni10.Bond( uni10.BD_OUT, d_phys)
 H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
 #c_i, c_i_dag, iden=C_i(len(d_phys))

 c_iu, c_iu_dag, iden=C_i(len(d_phys))
 c_id, c_id_dag, iden=C_i(len(d_phys))
 c_id=0*c_id
 c_id_dag=0*c_id_dag


# c_iu, c_iu_dag, iden=C_i_spinUP()
# c_id, c_id_dag, iden=C_i_spinDOWN()

 for i in xrange(N_x-1):
  for j in xrange(N_x):

   if i==0 :
#    ham = (0.25)*(2.0*uni10.otimes(c_iu_dag*c_iu,iden)+uni10.otimes(iden,c_iu_dag*c_iu))
#    ham =ham+ (0.25)*(2.0*uni10.otimes(c_id_dag*c_id,iden)+uni10.otimes(iden,c_id_dag*c_id))
    ham=Ham_mat( 0.0,  0.0, 0.0,1.0*2., 1.0 , 0., 0.)

   elif i==N_x-2:
#    ham = (0.25)*(uni10.otimes(c_iu_dag*c_iu,iden)+2.0*uni10.otimes(iden,c_iu_dag*c_iu))
#    ham =ham+ (0.25)*(uni10.otimes(c_id_dag*c_id,iden)+2.0*uni10.otimes(iden,c_id_dag*c_id))
    ham=Ham_mat( 0.0,  0.0, 0.0,1.0, 1.0*2. , 0., 0.)

   else:
#    ham = (0.25)*(uni10.otimes(c_iu_dag*c_iu,iden)+uni10.otimes(iden,c_iu_dag*c_iu))
#    ham =ham+ (0.25)*(uni10.otimes(c_id_dag*c_id,iden)+uni10.otimes(iden,c_id_dag*c_id))
    ham=Ham_mat( 0.0,  0.0, 0.0,1.0, 1.0 , 0., 0.)

   H.setRawElem(ham*0.25)
   H_list[i][j]=H*1.0


 return H_list



def  make_N_long(N_x, d_phys,h, Model):

 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heis")

 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*N_x 

 bdi = uni10.Bond( uni10.BD_IN, d_phys)
 bdo = uni10.Bond( uni10.BD_OUT, d_phys)
 H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
# c_i, c_i_dag, iden=C_i(len(d_phys))

 c_iu, c_iu_dag, iden=C_i(len(d_phys))
 c_id, c_id_dag, iden=C_i(len(d_phys))
 c_id=0*c_id
 c_id_dag=0*c_id_dag

# c_iu, c_iu_dag, iden=C_i_spinUP()
# c_id, c_id_dag, iden=C_i_spinDOWN()

 for i in xrange(N_x):
  for j in xrange(N_x):

   if i==0 :
    ham = (0.25)*(2.0*uni10.otimes(c_iu_dag*c_iu,iden)+uni10.otimes(iden,c_iu_dag*c_iu))
   elif i==N_x-2:
    ham = (0.25)*(uni10.otimes(c_iu_dag*c_iu,iden)+2.0*uni10.otimes(iden,c_iu_dag*c_iu))
   else:
    ham = (0.25)*(uni10.otimes(c_iu_dag*c_iu,iden)+uni10.otimes(iden,c_iu_dag*c_iu))

   H.setRawElem(ham)
   H_list[i][j]=H*0.0


 return H_list



def  make_U_col(H_col, N_x, delta):

 U_ham_list=[None]*(N_x)

 for i in xrange(N_x):
  U_ham_list[i]=[None]*(N_x-1) 

 for i in xrange(N_x):
  for j in xrange(N_x-1):
   U_ham = uni10.UniTensor(H_col[i][j].bond(), "U_ham")
   blk_qnums = H_col[i][j].blockQnum()
   for qnum in blk_qnums:
    U_ham.putBlock(qnum, uni10.takeExp(-delta, H_col[i][j].getBlock(qnum)))
   U_ham_list[i][j]=U_ham*1.0

 return  U_ham_list


def  make_U_row(H_row, N_x, delta):
 U_ham_list=[None]*(N_x-1)

 for i in xrange(N_x-1):
  U_ham_list[i]=[None]*(N_x) 

 for i in xrange(N_x-1):
  for j in xrange(N_x):
   U_ham = uni10.UniTensor(H_row[i][j].bond(), "U_ham")
   blk_qnums = H_row[i][j].blockQnum()
   for qnum in blk_qnums:   
    U_ham.putBlock(qnum,uni10.takeExp(-delta, H_row[i][j].getBlock(qnum)))
   U_ham_list[i][j]=U_ham*1.0

 return  U_ham_list










#######@profile
def  TEBD_Full( H_col, H_row, N_x, PEPS_listten, E_iter_list, D, accuracy, N_tebd, i, d, chi_single, chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model, count_list, E_iter_list1, Sys):

 PEPS_listtenU=[None]*N_x
 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_temp=[None]*N_x
 mps_boundry_up=[None]*(N_x+1)
 mps_boundry_down=[None]*N_x

 Error_try=[]
 E_tebd=[]

 for i in xrange(N_x):
  PEPS_listtenU[i]=Init_PEPS( N_x, D, d, i)

 Fidel_val=1
 E_coulmn=[]
 E_row=[]
 E_mag_coulmn=[]
 E_mag_row=[]
 E_0=1.0
 E_1=1.0


 E_00=1.0
 E_11=1.0
 Energy_val=2.0

 mps_I=MPSclass.MPS(1,1,2*N_x)
 for i_ind in xrange(mps_I.N):
  mps_I[i_ind].identity()

 List_delN=Short_TrotterSteps(N_tebd)

 print List_delN
 E_00_f=10.0
 E_min=1
 E_0=1
 E_1=100
 count_iter=0
 for delta, N_iter in List_delN:
  #delta=0
  print delta, N_iter

  U_ham_col=make_U_col(H_col, N_x, delta)
  U_ham_row=make_U_row(H_row, N_x, delta)

  PEPS_listten, norm_val, count_val=Normalize_PEPS(PEPS_listten, N_x, D, chi_try, d, threshold, interval,Sys)
  Energy_val=Energy_cal( PEPS_listten, d, chi_single, N_x, D, H_col, H_row, Sys)
  E_0=Energy_val
  #E_00=E_00_f*1.0
  E_11=Energy_val*1.0
  print "E_0", Energy_val
  E_iter_list1.append("E_0, PEPS")
  E_iter_list1.append(Energy_val)

  PEPS_listtenU=copy_f(PEPS_listten,N_x, PEPS_listtenU)
  PEPS_listten=copy_f( PEPS_listtenU, N_x, PEPS_listten)

  for q_iter in xrange(N_iter[0]):
   PEPS_listten, norm_val, count_val=Normalize_PEPS(PEPS_listten, N_x, D, chi_try, d, threshold, interval,Sys)
   #Energy_val=Energy_cal(PEPS_listten, d, chi_single, N_x, D, Model, h_coupling,Sys)
   #print "E_0", Energy_val

##########################   Col   #############################
   for Location in reversed(xrange(1,N_x)):
    mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_single, N_x,Sys)
#############
   E_coulmn_t=[]
   for i_ind in xrange(N_x):
    if i_ind==0:
     update_twotensor(mps_I, mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_x, U_ham_col[i_ind], H_col[i_ind], D, E_coulmn_t,threshold, interval,Sys)
    elif i_ind==N_x-1:
     update_twotensor(mps_boundry_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_x,U_ham_col[i_ind], H_col[i_ind], D, E_coulmn_t,threshold, interval,Sys)
    else:
     update_twotensor(mps_boundry_left[i_ind-1], mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_x,U_ham_col[i_ind], H_col[i_ind], D, E_coulmn_t,threshold, interval,Sys)

    if  i_ind==(N_x-1):break
    mps_boundry_left[i_ind]=make_Env_singleLayer(PEPS_listten[i_ind], i_ind, mps_boundry_left[i_ind-1], d, chi_single, N_x,Sys)


###########################   Row   ###############################

   for Location in reversed(xrange(1,N_x)):

    peps_l=[] 

    for i in xrange(N_x): 
     peps_l.append(PEPS_listten[i][Location]) 
    mps_boundry_up[Location]=make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], d, chi_single, N_x ,Sys)

   E_row_t=[]

   for i_ind in xrange(N_x):
    peps_l=[]
    H_l=[] 
    U_l=[] 
    for i in xrange(N_x):
     peps_l.append(PEPS_listten[i][i_ind])
    for i in xrange(N_x-1): 
     H_l.append(H_row[i][i_ind]) 
    for i in xrange(N_x-1): 
     U_l.append(U_ham_row[i][i_ind]) 


    if i_ind==0:
     update_twotensor_row(mps_I, mps_boundry_up[i_ind+1], peps_l, N_x, U_l, H_l, D, E_row_t,threshold, interval,Sys)
    elif i_ind==N_x-1:
     update_twotensor_row(mps_boundry_down[i_ind-1], mps_I, peps_l, N_x,U_l, H_l, D, E_row_t, threshold, interval,Sys)
    else:
     update_twotensor_row(mps_boundry_down[i_ind-1], mps_boundry_up[i_ind+1], peps_l, N_x, U_l, H_l, D, E_row_t,threshold, interval,Sys)

    for i in xrange(N_x):
     PEPS_listten[i][i_ind]=peps_l[i]

    if  i_ind==(N_x-1):break
    mps_boundry_down[i_ind]=make_Env_singleLayer_down( peps_l, i_ind, mps_boundry_down[i_ind-1], d, chi_single, N_x,Sys)






   E_00=E_11*1.0
   E_11=(sum(E_coulmn_t)+sum(E_row_t))
   #E_11=sum(E_coulmn_t)
   #E_11=sum(E_row_t)
   #print  sum(E_coulmn_t),  sum(E_row_t)
   print "E_t", q_iter, E_11, E_00, sum(E_coulmn_t),  sum(E_row_t), (E_11-E_00) / E_11 
   E_iter_list1.append(E_11)


   Error_try.append((E_11-E_00) / E_11)
   E_tebd.append(E_11)

   file = open("TEBD.txt", "w")
   for index in range(len(Error_try)):
    file.write(str(index) + " " + str(Error_try[index])+" " + str(E_tebd[index])+" "+ "\n")
   file.close()

   if N_iter[1]=="off": 
    Store_f(PEPS_listten, N_x)


   break_criteria="off"
   if N_iter[1]=="on":
    if E_11 < E_00 and q_iter>0:
     Store_f(PEPS_listten, N_x)
     PEPS_listtenU=copy_f(PEPS_listten,N_x, PEPS_listtenU)


    if abs((E_11-E_00)/E_11)<accuracy or E_11>E_00:
     if abs((E_11-E_00)/E_11)<accuracy:
      print "loop_finished:reached_levelofaccuracy"
      E_iter_list1.append("loop_finished:reached_levelofaccuracy")
     else: 
      print "Didnot_get_lowered";
      E_iter_list1.append("Didnot_get_lowered")
     PEPS_listten=copy_f( PEPS_listtenU, N_x, PEPS_listten)
     break_criteria="on"
     break

   if break_criteria=="on": break




# 
# 
# 
#    E_00=E_11*1.0
#    E_11=(sum(E_coulmn_t)+sum(E_row_t))
#    #E_11=sum(E_coulmn_t)
#    #E_11=sum(E_row_t)
#    #print  sum(E_coulmn_t),  sum(E_row_t)
#    print "E_t", q_iter, E_11, E_00, sum(E_coulmn_t),  sum(E_row_t), (E_11-E_00) / E_11 
#    E_iter_list1.append(E_11)
# 
#    if E_11 < E_00 and q_iter>0:
#     Store_f(PEPS_listten, N_x)
#     PEPS_listtenU=copy_f(PEPS_listten,N_x, PEPS_listtenU)
#     #E_00_f=E_11*1.0
# 
# 
#    if abs((E_11-E_00)/E_11)<accuracy or E_11>E_00:
#     if abs((E_11-E_00)/E_11)<accuracy:
#      print "loop_finished:reached_levelofaccuracy"
#      E_iter_list1.append("loop_finished:reached_levelofaccuracy")
#     else: 
#      print "Didnot_get_lowered";
#      E_iter_list1.append("Didnot_get_lowered")
# 
#     PEPS_listten=copy_f( PEPS_listtenU, N_x, PEPS_listten)
#     
#     break







  Energy_val=Energy_cal( PEPS_listten, d, chi_single, N_x, D, H_col, H_row, Sys)
  print "E_0", Energy_val

















#######@profile
def  TEBD_Full_RG( H_col, H_row, N_x, PEPS_listten, E_iter_list, D, accuracy, N_tebd, i, d, chi_single, chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model, count_list, E_iter_list1, Sys, N_part, h_coupling, N_col, N_row, Sz_col, Sz_row):

 PEPS_listtenU=[None]*N_x
 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_temp=[None]*N_x
 mps_boundry_up=[None]*(N_x+1)
 mps_boundry_down=[None]*N_x

 for i in xrange(N_x):
  PEPS_listtenU[i]=Init_PEPS( N_x, D, d, i)

 Fidel_val=1
 E_coulmn=[]
 E_row=[]
 E_mag_coulmn=[]
 E_mag_row=[]
 E_0=1.0
 E_1=1.0

 E_00=1.0
 E_11=1.0
 Energy_val=2.0

 mps_I=MPSclass.MPS(1,1,2*N_x)
 for i_ind in xrange(mps_I.N):
  mps_I[i_ind].identity()

 List_delN=Short_TrotterSteps(N_tebd)
 N_list=[]
 Sz_list=[]
 E_list_try=[]

 Error_try=[]
 E_tebd=[]

 E_E0=[]

 print   List_delN
 E_00_f=10.0
 E_min=1
 E_0=1
 E_1=100
 count_iter=0
 for delta, N_iter in List_delN:
  #delta=0
  print delta, N_iter

  U_ham_col=make_U_col(H_col, N_x, delta)
  U_ham_row=make_U_row(H_row, N_x, delta)

  PEPS_listten, norm_val, count_val=Normalize_PEPS(PEPS_listten, N_x, D, chi_try, d, threshold, interval,Sys)

  rho_row, rho_col=make_density_matrix_sinlgeLayer( PEPS_listten, N_x, chi_single, d, D, Sys)
  N_0=Energy_from_Density(rho_row, rho_col, N_col, N_row, N_x)
  print "N", N_0

  Sz_0=Energy_from_Density(rho_row, rho_col, Sz_col, Sz_row, N_x)
  print "Sz", Sz_0

  Energy_val=Energy_cal( PEPS_listten, d, chi_single, N_x, D, H_col, H_row, Sys)
  E_0=Energy_val
  #E_00=E_00_f*1.0
  E_11=Energy_val*1.0
  print "E_0", Energy_val
  print "E_final",(Energy_val-h_coupling[2]*N_part)*0.5
  print "E_final",(Energy_val-h_coupling[2]*N_0)*0.5

  #E_iter_list1.append("E_0, PEPS")
  #E_iter_list1.append(Energy_val)

  PEPS_listtenU=copy_f(PEPS_listten,N_x, PEPS_listtenU)
  PEPS_listten=copy_f( PEPS_listtenU, N_x, PEPS_listten)

  for q_iter in xrange(N_iter[0]):

   rho_row, rho_col=make_density_matrix_sinlgeLayer( PEPS_listten, N_x, chi_single, d,D,Sys)
   N_0=Energy_from_Density(rho_row, rho_col, N_col, N_row, N_x)
   Sz_0=Energy_from_Density(rho_row, rho_col, Sz_col, Sz_row, N_x)
   print "Sz", Sz_0
   Energy_val=Energy_cal( PEPS_listten, d, chi_single, N_x, D, H_col, H_row, Sys)
   print "N", N_0
   print "E_0", Energy_val
   print "E_final",(Energy_val-h_coupling[2]*N_part)*0.5
   print "E_final",(Energy_val-h_coupling[2]*N_0)*0.5

   E_iter_list.append(Energy_val-h_coupling[2]*N_0)
   E_list_try.append(Energy_val-h_coupling[2]*N_part)
   N_list.append(N_0)
   Sz_list.append(Sz_0)
   E_E0.append(Energy_val)
   file = open("EnergyRG.txt", "w")
   for index in range(len(E_iter_list)):
    file.write(str(index) + " " + str(E_iter_list[index])+" " + str(E_list_try[index])+" "+str(E_E0[index])+" "+str(N_list[index])+" "+str(Sz_list[index])+" "+ "\n")
   file.close()

##########################   Col   #############################
   for Location in reversed(xrange(1,N_x)):
    mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_single, N_x,Sys)
#############
   E_coulmn_t=[]
   for i_ind in xrange(N_x):
    if i_ind==0:
     update_twotensor(mps_I, mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_x, U_ham_col[i_ind], H_col[i_ind], D, E_coulmn_t,threshold, interval,Sys)
    elif i_ind==N_x-1:
     update_twotensor(mps_boundry_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_x,U_ham_col[i_ind], H_col[i_ind], D, E_coulmn_t,threshold, interval,Sys)
    else:
     update_twotensor(mps_boundry_left[i_ind-1], mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_x,U_ham_col[i_ind], H_col[i_ind], D, E_coulmn_t,threshold, interval,Sys)

    if  i_ind==(N_x-1):break
    mps_boundry_left[i_ind]=make_Env_singleLayer(PEPS_listten[i_ind], i_ind, mps_boundry_left[i_ind-1], d, chi_single, N_x,Sys)

###########################   Row   ###############################

   for Location in reversed(xrange(1,N_x)):

    peps_l=[] 

    for i in xrange(N_x): 
     peps_l.append(PEPS_listten[i][Location]) 
    mps_boundry_up[Location]=make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], d, chi_single, N_x ,Sys)

   E_row_t=[]

   for i_ind in xrange(N_x):
    peps_l=[]
    H_l=[] 
    U_l=[] 
    for i in xrange(N_x):
     peps_l.append(PEPS_listten[i][i_ind])
    for i in xrange(N_x-1): 
     H_l.append(H_row[i][i_ind]) 
    for i in xrange(N_x-1): 
     U_l.append(U_ham_row[i][i_ind]) 

    if i_ind==0:
     update_twotensor_row(mps_I, mps_boundry_up[i_ind+1], peps_l, N_x, U_l, H_l, D, E_row_t,threshold, interval,Sys)
    elif i_ind==N_x-1:
     update_twotensor_row(mps_boundry_down[i_ind-1], mps_I, peps_l, N_x,U_l, H_l, D, E_row_t, threshold, interval,Sys)
    else:
     update_twotensor_row(mps_boundry_down[i_ind-1], mps_boundry_up[i_ind+1], peps_l, N_x, U_l, H_l, D, E_row_t,threshold, interval,Sys)

    for i in xrange(N_x):
     PEPS_listten[i][i_ind]=peps_l[i]

    if  i_ind==(N_x-1):break
    mps_boundry_down[i_ind]=make_Env_singleLayer_down( peps_l, i_ind, mps_boundry_down[i_ind-1], d, chi_single, N_x,Sys)

   E_00=E_11*1.0
   E_11=(sum(E_coulmn_t)+sum(E_row_t))
   #E_11=sum(E_coulmn_t)
   #E_11=sum(E_row_t)
   #print  sum(E_coulmn_t),  sum(E_row_t)
   print "E_t", q_iter, E_11, E_00, sum(E_coulmn_t),  sum(E_row_t), (E_11-E_00) / E_11 
   E_iter_list1.append(E_11)

   Error_try.append((E_11-E_00) / E_11)
   E_tebd.append(E_11)

   file = open("TEBD.txt", "w")
   for index in range(len(Error_try)):
    file.write(str(index) + " " + str(Error_try[index])+" " + str(E_tebd[index])+" "+ "\n")
   file.close()

   if N_iter[1]=="off": 
    Store_fRG(PEPS_listten, N_x)

   break_criteria="off"
   if N_iter[1]=="on":
    if E_11 < E_00 and q_iter>0:
     Store_fRG(PEPS_listten, N_x)
     PEPS_listtenU=copy_f(PEPS_listten,N_x, PEPS_listtenU)

    if abs((E_11-E_00)/E_11)<accuracy or E_11>E_00:
     if abs((E_11-E_00)/E_11)<accuracy:
      print "loop_finished:reached_levelofaccuracy"
      E_iter_list1.append("loop_finished:reached_levelofaccuracy")
     else: 
      print "Didnot_get_lowered";
      E_iter_list1.append("Didnot_get_lowered")
     PEPS_listten=copy_f( PEPS_listtenU, N_x, PEPS_listten)
     break_criteria="on"
     break

   if break_criteria=="on": break


  Energy_val=Energy_cal( PEPS_listten, d, chi_single, N_x, D, H_col, H_row, Sys)
  print "E_0", Energy_val






def  renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval):
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 norm_val=Norm_h[0]
 #print "Zero_order_local",  norm_val

 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold[1]:
   r_up=r_up*alpha
   r_dp=r_dp*alpha
   l_up=l_up*alpha
   l_dp=l_dp*alpha
   Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
   norm_val=Norm_h[0]
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<(threshold[0]):
   r_up=r_up*alpha
   r_dp=r_dp*alpha
   l_up=l_up*alpha
   l_dp=l_dp*alpha
   Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
   norm_val=Norm_h[0]

 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 norm_val=Norm_h[0]
 #print "Fixednorm_local", abs(norm_val)
 return r_up, l_up





def increase_physicalbond( Q_list, PEPS_list, N_x, d_out):
 bdi = uni10.Bond( uni10.BD_IN, Q_list[0][0].bond(4).Qlist())
 bdo = uni10.Bond( uni10.BD_OUT, d_out)
 Increase_tensor = uni10.UniTensor( [bdi, bdo], "Increase_tensor")
 Increase_tensor.setLabel([-1, -2])
 Increase_tensor.identity()


 bdo = uni10.Bond( uni10.BD_OUT, PEPS_list[0][0].bond(2).Qlist())
 bdi = uni10.Bond( uni10.BD_IN, d_out)
 Increase_tensor1 = uni10.UniTensor( [bdi, bdo], "Increase_tensor")
 Increase_tensor1.setLabel([-2, -1])
 Increase_tensor1.identity()

 for i in xrange(N_x):
  for j in xrange(N_x):

   Increase_tensor.setLabel([-1, -2])
   Q_list[i][j].setLabel([0,1,2,3,-1])
   Q_list[i][j]=Q_list[i][j]*Increase_tensor
   Q_list[i][j].permute([0,1,2,3,-2],4)

   
   PEPS_list[i][j].setLabel([1,3,-1,2,4])
   PEPS_list[i][j]=PEPS_list[i][j]*Increase_tensor1
   PEPS_list[i][j].permute([1,3,-2,2,4],3)











def Q_feed_p( N_p, PEPS_listten, Q_list):
 
 for i in xrange(1,N_p-1):
  for j in xrange(1,N_p-1):
   #print i,j 
   A=uni10.UniTensor("Store/StoreP/a" + str(i)+"P"+str(j))
   #print A.printDiagram()
   PEPS_listten[i][j]=A*1.0

 for i in xrange(1,N_p-1):
  for j in xrange(1,N_p-1):
   #print i,j 
   A=uni10.UniTensor("StoreQ/StoreQP/Q" + str(i)+"P"+str(j))
   #print A.printDiagram()
   Q_list[i][j]=A*1.0
 

 return Q_list, PEPS_listten






def Store_fRG(PEPS_listten, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   PEPS_listten[i][j].save("StoreRG/a" + str(i)+"P"+str(j))


def Reload_fRG(PEPS_listten, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   PEPS_listten[i][j]=uni10.UniTensor("StoreRG/a" + str(i)+"P"+str(j))





def Store_f(PEPS_listten, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   PEPS_listten[i][j].save("Store/a" + str(i)+"P"+str(j))


def Reload_f(PEPS_listten, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   PEPS_listten[i][j]=uni10.UniTensor("Store/a" + str(i)+"P"+str(j))


def Store_Q_list(Q_list, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   Q_list[i][j].save("StoreQ/Q" + str(i)+"P"+str(j))

def Reload_Q_list(Q_list, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   Q_list[i][j]=uni10.UniTensor("StoreQ/Q" + str(i)+"P"+str(j))


def   copy_f( PEPS_listten, N_x, PEPS_listtenU):
 for i in xrange(N_x):
  for j in xrange(N_x):
   PEPS_listtenU[i][j]=PEPS_listten[i][j]*1.0
 return  PEPS_listtenU







def  cal_rowcol(T):
 blk_qnums = T.blockQnum()
 row_list=[]
 col_list=[]
 for qnum in blk_qnums:
  M_tem=T.getBlock(qnum)
  col_list.append(M_tem.col())
  row_list.append(M_tem.row())
 return  sum(row_list),  sum(col_list)




def  density_matrix_list_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x, D,Location_y,rho_row ,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[2*i+1]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[2*i+1]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=((mps_boundry_up[2*i+1]*Swap1))*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[2*i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[2*i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[2*i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[2*i]*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)


 for i in xrange(N_x-1):
  rho_row[i][Location_y]=Density_local_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, D, Sys)



#####@profile
def  Density_local_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, D, Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)


 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)


 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],2)
 r_u.setLabel([-100,3,9])
 r_u.permute([-100,3,9],2)

 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)



 #####################################################
 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)

 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys, A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys, A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel( [8,9,-1,-2] )
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 Swap1=fermionicOPT(Sys,A_conj.bond(2), A_conj.bond(4))
 Swap1.setLabel([3,10,-3,-10])
 A_conj=Swap1*A_conj
 A_conj.permute([-6,-7,3,-9,10],3)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 A_conj.permute([-6,-7,-10,-3,-9],3)

 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  q_d, V, s=TU.setTruncation(A_conj, row)
 else:
  q_d, V, s=TU.setTruncation(A_conj, colm)


 s.setLabel([1,0])
 V.setLabel([0,-3,-9])
 r_d=V*s
 r_d.permute([1,-3,-9],2)

 q_d.setLabel([-6,-7,-10,-200])
 q_d.permute([-6,-7,-10,-200],4)

 r_d.setLabel([-200,-3,-9])
 r_d.permute([-200,-3,-9],2)



 A=Peps_2*1.0
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],2)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 l_u.setLabel([9,13,-300])
 l_u.permute([9,13,-300],2)


 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)

 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys, A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys, A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel( [8,9,-1,-2] )
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 A_conj.setLabel([-9,-10,-13,-14,-15])
 A_conj.permute([-9,-13,-10,-14,-15],2)


 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  U, qq_d, s=TU.setTruncation(A_conj, row)
 else:
  U, qq_d, s=TU.setTruncation(A_conj, colm)

 s.setLabel([0,1])
 U.setLabel([-9,-13,0])
 l_d=U*s
 l_d.permute([-9,-13,1],2)

 qq_d.setLabel([-400,-10,-14,-15])
 qq_d.permute([-400,-10,-14,-15],4)


 l_d.setLabel([-9,-13,-400])
 l_d.permute([-9,-13,-400],2)


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 Swap4.setLabel([23,400,10,-400])


######################################################################
 mps_boundry_up[Location*2].setLabel([16,10,26])
 mps_boundry_up[Location*2+1].setLabel([26,25,27])
 mps_boundry_up[Location*2+2].setLabel([27,15,28])
 mps_boundry_up[Location*2+3].setLabel([28,24,21])

 mps_boundry_down[Location*2].setLabel([18,22,31])
 mps_boundry_down[Location*2+1].setLabel([31,-7,32])
 mps_boundry_down[Location*2+2].setLabel([32,23,33])
 mps_boundry_down[Location*2+3].setLabel([33,-10,19])
######################################################

 A=E_left*mps_boundry_up[Location*2]
 A=A*(Swap1*q)
 A=A*mps_boundry_down[2*Location]
 A=A*mps_boundry_up[2*Location+1]
 A=A*(Swap3*q_d)
 A=A*mps_boundry_down[2*Location+1]


 B=E_right*mps_boundry_up[2*Location+3]
 B=B*(Swap2*qq_d)
 B=B*mps_boundry_down[2*Location+3]

 B=B*mps_boundry_up[2*Location+2]
 B=B*(Swap4*qq)
 B=B*mps_boundry_down[2*Location+2]
 
 N_ten=A*B
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])
######################################################
 N_ten=(((r_u*r_d)*N_ten)*(l_u*l_d))
 N_ten.permute([-3,-13,3,13],2)

 return  N_ten



def density_matrix_list_col( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x, D, Location_x, rho, Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([13])
   mps_left=mps_boundry_left[2*i+1]*iden
   iden.setLabel([15])
   mps_right=mps_boundry_right[2*i+1]*iden

   iden.setLabel([5])
   iden1=iden*1.0
   iden1.setLabel([11])


   E_list_up[i]=Peps_list*iden
   E_list_up[i]=E_list_up[i]*(Swap1*iden1)
   E_list_up[i]=E_list_up[i]*mps_left
   E_list_up[i]=E_list_up[i]*mps_right
   E_list_up[i]=E_list_up[i]*mps_boundry_left[2*i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=E_list_up[i]*mps_boundry_right[2*i]
   E_list_up[i]=E_list_up[i]*((Peps_list_conj*Swap3)*Swap2)

   E_list_up[i].permute([12,10,9,14],4)
  else:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])
   E_list_up[i+1].setLabel([13,5,11,15])

   E_list_up[i]=((Peps_list*Swap1))*E_list_up[i+1]
   E_list_up[i]=mps_boundry_left[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_right[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_left[2*i]*E_list_up[i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=mps_boundry_right[2*i]*E_list_up[i]
   E_list_up[i]=((Peps_list_conj*Swap3)*Swap2)*E_list_up[i]
   E_list_up[i].permute([12,10,9,14],4)





 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[2*i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[2*i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])


   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])


   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)

 for i in xrange(N_x-1):
   #E_val=local_energy( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H_orig,D,Sys)
   rho[Location_x][i]=Density_local_col( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, D, Sys)




def  Density_local_col( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location ,D, Sys):
 bdi_mid=uni10.Bond( uni10.BD_IN, 1)
 iden=uni10.UniTensor( [bdi_mid] )
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],2)

 r_u.setLabel([-100,3,10])

 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


#####################################################

 A=Peps_1*1.0

 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],3)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 A_conj.permute([-6,-7,-9,-3,-10],3)



 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  q_d, V, s=TU.setTruncation(A_conj, row)
 else:
  q_d, V, s=TU.setTruncation(A_conj, colm)


 s.setLabel([1,0])
 V.setLabel([0,-3,-10])
 r_d=V*s
 r_d.permute([1,-3,-10],3)
 
 q_d.setLabel([-6,-7,-9,-200])
 q_d.permute([-6,-7,-9,-200],4)

 r_d.setLabel([-200,-3,-10])



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)

 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])


 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 Swap1=fermionicOPT(Sys,A_conj.bond(2), A_conj.bond(0))
 Swap1.setLabel([11,10,3,8])
 A_conj=A_conj*Swap1

 A_conj.permute([10,9,11,-6,7],3)
 A_conj.setLabel([-11,-10,-13,-14,-15])
 A_conj.permute([-10,-13,-11,-14,-15],2)


 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  U, qq_d, s=TU.setTruncation(A_conj, row)
 else:
  U, qq_d, s=TU.setTruncation(A_conj, colm)

 s.setLabel([0,1])
 U.setLabel([-10,-13,0])
 l_d=U*s
 l_d.permute([-10,-13,1],3)

 qq_d.setLabel([-400,-11,-14,-15])
 qq_d.permute([-400,-11,-14,-15],4)


 l_d.setLabel([-10,-13,-400])


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([ 90, 200, 9, -200 ])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])


######################################################################

 mps_boundry_left[Location*2].setLabel([16,-60,18])
 mps_boundry_left[Location*2+1].setLabel([18,6,21])
 mps_boundry_left[Location*2+2].setLabel([21,-110,22])
 mps_boundry_left[Location*2+3].setLabel([22,11,25])

 mps_boundry_right[Location*2].setLabel([17,-9,19])
 mps_boundry_right[Location*2+1].setLabel([19,90,20])
 mps_boundry_right[Location*2+2].setLabel([20,-14,23])
 mps_boundry_right[Location*2+3].setLabel([23,140,24])


######################################################

 A=E_left*mps_boundry_left[Location*2]
 A=A*(Swap1*q_d)
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_left[2*Location+1]
 A=A*(Swap3*q)
 A=A*mps_boundry_right[2*Location+1]


 B=E_right*mps_boundry_left[2*Location+3]
 B=B*(Swap2*qq)
 B=B*mps_boundry_right[2*Location+3]

 B=B*mps_boundry_left[2*Location+2]
 B=B*(Swap4*qq_d)
 B=B*mps_boundry_right[2*Location+2]
 

 N_ten=A*B
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 l_d.setLabel([-10,-13,-400])
 l_u.setLabel([10,13,-300])
 r_u.setLabel([-100,3,10])
 r_d.setLabel([-200,-3,-10])

 N_ten=((((N_ten*l_d)*l_u)*r_u)*r_d)
 N_ten.permute([-3,-13,3,13],2)


 return  N_ten




def make_density_matrix_sinlgeLayer( PEPS_listten, N_x, chi_single,d,D,Sys):

 rho_row=[None]*(N_x-1)
 for i in xrange(N_x-1):
  rho_row[i]=[None]*N_x

 rho_col=[None]*(N_x)
 for i in xrange(N_x):
  rho_col[i]=[None]*(N_x-1)


 mps_boundry_left=[None]*(N_x)
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_temp=[None]*(N_x)

 mps_boundry_up=[None]*(N_x+1)
 mps_boundry_down=[None]*N_x


 mps_I=MPSclass.MPS(1,1,2*N_x)
 for i_ind in xrange(mps_I.N):
  mps_I[i_ind].identity()


##################   Col   #####################

 for Location in reversed(xrange(1,N_x)):
  mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_single, N_x,Sys)

 for Location in xrange(N_x-1):
  mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_single, N_x,Sys)
#############################################################


 for i_ind in xrange(N_x):
  if i_ind==0:
   density_matrix_list_col( mps_I,mps_boundry_right[i_ind+1],PEPS_listten[i_ind], N_x, D, i_ind,rho_col,Sys)
  elif i_ind==N_x-1:
   density_matrix_list_col( mps_boundry_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_x, D, i_ind,rho_col,Sys)
  else:
   density_matrix_list_col( mps_boundry_left[i_ind-1], mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_x, D,i_ind,rho_col,Sys)


##################   row   #####################

##################   Env   #####################

 for Location in xrange(N_x):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_down[Location]=make_Env_singleLayer_down( peps_l, Location, mps_boundry_down[Location-1], d, chi_single, N_x,Sys)


 for Location in reversed(xrange(1,N_x)):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_up[Location]=make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], d, chi_single, N_x ,Sys)

#############################################################


 for i_ind in xrange(N_x):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])

  if i_ind==0:
   density_matrix_list_row( mps_I,mps_boundry_up[i_ind+1],peps_l, N_x, D, i_ind,rho_row,Sys)
  elif i_ind==N_x-1:
   density_matrix_list_row( mps_boundry_down[i_ind-1], mps_I, peps_l, N_x, D, i_ind,rho_row,Sys)
  else:
   density_matrix_list_row( mps_boundry_down[i_ind-1], mps_boundry_up[i_ind+1], peps_l, N_x, D,i_ind,rho_row,Sys)


 return rho_row, rho_col


def  Energy_from_Density(rho_row, rho_col, HA_col, HA_row,N_x):
 E_c=0
 E_c_list=[]
 for i in xrange(N_x):
  for j in xrange(N_x-1):
   HA_col[i][j].setLabel([-1,-2,1,2])
   rho_col[i][j].setLabel([-1,-2,1,2])
   val=HA_col[i][j]*rho_col[i][j]
   Iden=HA_col[i][j]*1.0
   Iden.identity()
   Iden.setLabel([-1,-2,1,2])
   norm=rho_col[i][j]*Iden
   E_c=(val[0]/norm[0])+E_c
   E_c_list.append(val[0]/norm[0])
   #print "dens", i,j, val[0]
 E_r=0
 E_r_list=[]
 for i in xrange(N_x-1):
  for j in xrange(N_x):
   HA_row[i][j].setLabel([-1,-2,1,2])
   rho_row[i][j].setLabel([-1,-2,1,2])
   val=HA_row[i][j]*rho_row[i][j]
   Iden=HA_row[i][j]*1.0
   Iden.identity()
   Iden.setLabel([-1,-2,1,2])
   norm=rho_row[i][j]*Iden
   E_r=(val[0]/norm[0])+E_r
   E_r_list.append(val[0]/norm[0])
 #print E_c_list, E_r_list
 return E_c+E_r








def   particle_val( H_col, Q_list, i, j, rho_col, Sys,P_f,N_x):

  Swap=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )
  Swap1=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )


  if j==N_x-1:
   A_iden=rho_col[i][j]*1.0
   A_iden.identity()
   A_iden.setLabel([-10,-5,10,5])
   rho_col[i][j].setLabel([-10,-5,10,5])
   normlization=A_iden*rho_col[i][j]
   
   
   
   Q_list[i][j].setLabel([1,2,3,4,5])
   Swap.setLabel([-2,-3,2,3])

   Q_list_ferm=Swap*Q_list[i][j]
   Q_list_ferm.permute([1,-2,-3,4,5],5)

   Q_list_ferm.setLabel([1,2,3,4,5])
   Q_trans=Q_list_ferm*1.0
   Q_trans.transpose()
   Q_trans.setLabel([-1,2,-3,4,-5])

   Q_list[i][j].setLabel([6,7,8,9,10])
   Q_transN=Q_list[i][j]*1.0
   Q_transN.transpose()
   Q_transN.setLabel([-10,6,7,8,9])
  
   H_col[2*i][2*j].setLabel([-1,-3,1,3])
   rho_col[i][j].setLabel([-10,-5,10,5])
   val=((Q_list_ferm*H_col[2*i][2*j]*Q_trans)*(Q_list[i][j]*Q_transN))*rho_col[i][j]
   P_f[2*i][2*j]=val[0]/normlization[0]
 ################################################


   Q_list[i][j].setLabel([1,2,3,4,5])
   Swap.setLabel([-3,-4,3,4])
   Q_list_ferm=Swap*Q_list[i][j]
   Q_list_ferm.permute([1,2,-3,-4,5],4)

   Q_list_ferm.setLabel([1,2,3,4,5])
   Q_trans=Q_list_ferm*1.0
   Q_trans.transpose()
   Q_trans.setLabel([-5,1,-2,3,-4])


   Q_list[i][j].setLabel([6,7,8,9,10])
   Q_transN=Q_list[i][j]*1.0
   Q_transN.transpose()
   Q_transN.setLabel([-10,6,7,8,9])

   H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
   rho_col[i][j].setLabel([-10,-5,10,5])
   val=((Q_list_ferm*H_col[2*i+1][2*j]*Q_trans)*(Q_list[i][j]*Q_transN))*rho_col[i][j]
   P_f[2*i+1][2*j]=val[0]/normlization[0]
 #######################################################


   Q_list[i][j].setLabel([1,2,3,4,5])
   Swap.setLabel([-2,-3,2,3])

   Q_list_ferm=Swap*Q_list[i][j]
   Q_list_ferm.permute([1,-2,-3,4,5],4)

   Q_list_ferm.setLabel([1,2,3,4,5])
   Q_trans=Q_list_ferm*1.0
   Q_trans.transpose()
   Q_trans.setLabel([-5,-1,2,-3,4])

   Q_list[i][j].setLabel([6,7,8,9,10])
   Q_transN=Q_list[i][j]*1.0
   Q_transN.transpose()
   Q_transN.setLabel([-10,6,7,8,9])
  
   H_col[2*i][2*j+1].setLabel([-1,-3,1,3])
   rho_col[i][j].setLabel([-10,-5,10,5])
   val=((Q_list_ferm*H_col[2*i][2*j+1]*Q_trans)*(Q_list[i][j]*Q_transN))*rho_col[i][j]
   P_f[2*i][2*j+1]=val[0]/normlization[0]
 ################################################


   Q_list[i][j].setLabel([1,2,3,4,5])
   Swap.setLabel([-3,-4,3,4])
   Q_list_ferm=Swap*Q_list[i][j]
   Q_list_ferm.permute([1,2,-3,-4,5],4)

   Q_list_ferm.setLabel([1,2,3,4,5])
   Q_trans=Q_list_ferm*1.0
   Q_trans.transpose()
   Q_trans.setLabel([-5,1,-2,3,-4])


   Q_list[i][j].setLabel([6,7,8,9,10])
   Q_transN=Q_list[i][j]*1.0
   Q_transN.transpose()
   Q_transN.setLabel([-10,6,7,8,9])

   H_col[2*i+1][2*j+1].setLabel([-2,-4,2,4])
   rho_col[i][j].setLabel([-10,-5,10,5])
   val=((Q_list_ferm*H_col[2*i+1][2*j+1]*Q_trans)*(Q_list[i][j]*Q_transN))*rho_col[i][j]
   P_f[2*i+1][2*j+1]=val[0]/normlization[0]



  else:
   A_iden=rho_col[i][j]*1.0
   A_iden.identity()
   A_iden.setLabel([-10,-5,10,5])
   rho_col[i][j].setLabel([-10,-5,10,5])
   normlization=A_iden*rho_col[i][j]
   




   Q_list[i][j].setLabel([1,2,3,4,5])
   Swap.setLabel([-2,-3,2,3])

   Q_list_ferm=Swap*Q_list[i][j]
   Q_list_ferm.permute([1,-2,-3,4,5],5)

   Q_list_ferm.setLabel([1,2,3,4,5])
   Q_trans=Q_list_ferm*1.0
   Q_trans.transpose()
   Q_trans.setLabel([-1,2,-3,4,-5])

   Q_list[i][j].setLabel([6,7,8,9,10])
   Q_transN=Q_list[i][j]*1.0
   Q_transN.transpose()
   Q_transN.setLabel([-10,6,7,8,9])
  
   H_col[2*i][2*j].setLabel([-1,-3,1,3])
   rho_col[i][j].setLabel([-5,-10,5,10])
   val=((Q_list_ferm*H_col[2*i][2*j]*Q_trans)*(Q_list[i][j]*Q_transN))*rho_col[i][j]
   P_f[2*i][2*j]=val[0]/normlization[0]
 ################################################


   Q_list[i][j].setLabel([1,2,3,4,5])
   Swap.setLabel([-3,-4,3,4])
   Q_list_ferm=Swap*Q_list[i][j]
   Q_list_ferm.permute([1,2,-3,-4,5],4)

   Q_list_ferm.setLabel([1,2,3,4,5])
   Q_trans=Q_list_ferm*1.0
   Q_trans.transpose()
   Q_trans.setLabel([-5,1,-2,3,-4])


   Q_list[i][j].setLabel([6,7,8,9,10])
   Q_transN=Q_list[i][j]*1.0
   Q_transN.transpose()
   Q_transN.setLabel([-10,6,7,8,9])

   H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
   rho_col[i][j].setLabel([-5,-10,5,10])
   val=((Q_list_ferm*H_col[2*i+1][2*j]*Q_trans)*(Q_list[i][j]*Q_transN))*rho_col[i][j]
   P_f[2*i+1][2*j]=val[0]/normlization[0]
 #######################################################


   Q_list[i][j].setLabel([1,2,3,4,5])
   Swap.setLabel([-2,-3,2,3])

   Q_list_ferm=Swap*Q_list[i][j]
   Q_list_ferm.permute([1,-2,-3,4,5],5)

   Q_list_ferm.setLabel([1,2,3,4,5])
   Q_trans=Q_list_ferm*1.0
   Q_trans.transpose()
   Q_trans.setLabel([-1,2,-3,4,-5])

   Q_list[i][j].setLabel([6,7,8,9,10])
   Q_transN=Q_list[i][j]*1.0
   Q_transN.transpose()
   Q_transN.setLabel([-10,6,7,8,9])
  
   H_col[2*i][2*j+1].setLabel([-1,-3,1,3])
   rho_col[i][j].setLabel([-5,-10,5,10])
   val=((Q_list_ferm*H_col[2*i][2*j+1]*Q_trans)*(Q_list[i][j]*Q_transN))*rho_col[i][j]
   P_f[2*i][2*j+1]=val[0]/normlization[0]
 ################################################


   Q_list[i][j].setLabel([1,2,3,4,5])
   Swap.setLabel([-3,-4,3,4])
   Q_list_ferm=Swap*Q_list[i][j]
   Q_list_ferm.permute([1,2,-3,-4,5],4)

   Q_list_ferm.setLabel([1,2,3,4,5])
   Q_trans=Q_list_ferm*1.0
   Q_trans.transpose()
   Q_trans.setLabel([-5,1,-2,3,-4])


   Q_list[i][j].setLabel([6,7,8,9,10])
   Q_transN=Q_list[i][j]*1.0
   Q_transN.transpose()
   Q_transN.setLabel([-10,6,7,8,9])

   H_col[2*i+1][2*j+1].setLabel([-2,-4,2,4])
   rho_col[i][j].setLabel([-5,-10,5,10])
   val=((Q_list_ferm*H_col[2*i+1][2*j+1]*Q_trans)*(Q_list[i][j]*Q_transN))*rho_col[i][j]
   P_f[2*i+1][2*j+1]=val[0]/normlization[0]







def  Particle_from_Density(rho_col, particle_col,N_x,Q_list,Sys):

 P_f=[None]*N_x
 for i in xrange(N_x):
  P_f[i]=[None]*(N_x)

 rho_col1=[None]*(N_x/2)
 for i in xrange(N_x/2):
  rho_col1[i]=[None]*(N_x/2)
 for i in xrange(N_x/2):
  for j in xrange((N_x/2)):
   if j==(N_x/2)-1:
    rho_col1[i][j]=rho_col[i][j-1]*1.0
   else:
    rho_col1[i][j]=rho_col[i][j]*1.0

 E_c=0
 E_c_list=[]
 for i in xrange(N_x/2):
  for j in xrange((N_x/2)):
   particle_val( particle_col, Q_list, i, j, rho_col1,  Sys,P_f, N_x/2)

 #print "Hi", P_f
 return P_f










def  Sz_from_Density_direct(rho_col, Magz_col, N_x, Sys):

 P_f=[None]*N_x
 for i in xrange(N_x):
  P_f[i]=[None]*(N_x)

 rho_col1=[None]*(N_x)
 for i in xrange(N_x):
  rho_col1[i]=[None]*(N_x)
 for i in xrange(N_x):
  for j in xrange(N_x):
   if j==(N_x)-1:
    rho_col1[i][j]=rho_col[i][j-1]*1.0
   else:
    rho_col1[i][j]=rho_col[i][j]*1.0

 E_c=0
 E_c_list=[]
 for i in xrange(N_x):
  for j in xrange((N_x)):
   A_iden=rho_col1[i][j]*1.0
   A_iden.identity()
   A_iden.setLabel([-10,-5,10,5])
   rho_col1[i][j].setLabel([-10,-5,10,5])
   normlization=A_iden*rho_col1[i][j]

   Magz_col[i][j].setLabel([-10,-5,10,5])
   rho_col1[i][j].setLabel([-10,-5,10,5])
   val=Magz_col[i][j]*rho_col1[i][j]
   P_f[i][j]=val[0]/normlization[0]


 #print "Hi", P_f
 return P_f



##@profile
def  Q_update( H_col, H_row, N_x, PEPS_listten, E_iter_list, D, accuracy, d, chi_single,chi_try, Q_list, Sys, H_long):

 if Sys[1] is "single" or Sys[1] is "simple":
  rho_row, rho_col=make_density_matrix_sinlgeLayer( PEPS_listten, N_x, chi_single, d, D,Sys)

 if Sys[1] is "double" or Sys[1] is "simple":
  rho_row, rho_col=make_density_matrix_double( PEPS_listten, N_x, chi_single, d, D,Sys)

 HA_col=Ascend_f_col( Q_list, H_col, N_x,Sys,H_long)
 HA_row=Ascend_f_row( Q_list, H_row, N_x,Sys,H_long)
 E_0=Energy_from_Density(rho_row, rho_col, HA_col, HA_row,N_x)
 print "E_form_rho",E_0 
 
 #Energy_val=Energy_cal(PEPS_listten, d, chi_single, N_x, D, HA_col, HA_row)
 #print "E_init", Energy_val
 
 
 HH_col=[None]*(2*N_x)
 for i in xrange(2*N_x):
   HH_col[i]=[None]*(2*N_x-1)

 for i in xrange(2*N_x):
  for j in xrange(2*N_x-1):

   HH_col[i][j]=H_col[i][j]*1.0
   list_energy1=[]
   blk_qnums = HH_col[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_col[i][j].getBlock(qnum)
    e=h_ham.eigh()
    N_col=e[0].col()
    for iter in xrange(N_col):
      list_energy1.append(e[0][iter])

   list_energy1=sorted(list_energy1, key=float)
   blk_qnums = HH_col[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_col[i][j].getBlock(qnum)
    iden_m=h_ham*1.0
    iden_m.identity()
    h_ham=h_ham+(-1.0)*iden_m*abs(max(list_energy1,key=abs))
    HH_col[i][j].putBlock(qnum,h_ham)
   #print "i,j",i,j,HH_col[i][j] 

 HH_row=[None]*(2*N_x-1)
 for i in xrange(2*N_x-1):
   HH_row[i]=[None]*(2*N_x)

 for i in xrange(2*N_x-1):
  for j in xrange(2*N_x):
   list_energy1=[]
   HH_row[i][j]=H_row[i][j]*1.0
   blk_qnums = HH_row[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_row[i][j].getBlock(qnum)
    e=h_ham.eigh()
    N_col=e[0].col()
    for iter in xrange(N_col):
      list_energy1.append(e[0][iter])

   list_energy1=sorted(list_energy1, key=float)
   blk_qnums = HH_row[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_row[i][j].getBlock(qnum)
    iden_m=h_ham*1.0
    iden_m.identity()
    h_ham=h_ham+(-1.0)*iden_m*abs(max(list_energy1,key=abs))
    HH_row[i][j].putBlock(qnum,h_ham)



 HH_long=[None]*(2*N_x)
 for i in xrange(2*N_x):
   HH_long[i]=[None]*(2*N_x)

 for i in xrange(2*N_x):
  for j in xrange(2*N_x):
   list_energy1=[]
   HH_long[i][j]=H_long[i][j]*1.0
   blk_qnums = HH_long[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_long[i][j].getBlock(qnum)
    e=h_ham.eigh()
    N_col=e[0].col()
    for iter in xrange(N_col):
      list_energy1.append(e[0][iter])

   list_energy1=sorted(list_energy1, key=float)
   blk_qnums = HH_long[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_long[i][j].getBlock(qnum)
    iden_m=h_ham*1.0
    iden_m.identity()
    h_ham=h_ham+(-1.0)*iden_m*abs(max(list_energy1,key=abs))
    HH_long[i][j].putBlock(qnum,h_ham)


 #print "hi", H_long[0][0]
 E_iter_list.append("E_Q")

 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  for j_ind in xrange(N_x):

   HA_col=Ascend_f_col( Q_list, H_col, N_x,Sys,H_long)
   HA_row=Ascend_f_row( Q_list, H_row, N_x,Sys,H_long)
   E_0=Energy_from_Density(rho_row, rho_col, HA_col, HA_row,N_x)
   print "E_Q", i_ind, j_ind, E_0
   E_iter_list.append(E_0)

   Q_init=Q_list[i_ind][j_ind]*1.0
   E2=10
   for iter in xrange(Sys[6]):
     E1=Q_cost_val_middle(HH_col, HH_row, Q_list, i_ind, j_ind , rho_row, rho_col, N_x,Sys,HH_long)
     if abs(E1)< 1.0e-14:
       #print "E1< 1.0e-14"
       break
     #print iter, E1, E2, abs((E2-E1)/E1)
     if abs(E1)>abs(E2) or iter is 0:
      Q_init=Q_list[i_ind][j_ind]*1.0
      if abs((E2-E1)/E1) < 1.0e-11:
       #print E2, E1, abs((E2-E1)/E1), i
       break
     else:
      print 'Notoptimized=i, E1, E2=', iter,'  ', E1, '   ', E2
      Q_list[i_ind][j_ind]=Q_init*1.0
      break
     Y=Grad_Q_middle( HH_col, HH_row, Q_list, i_ind, j_ind, rho_row, rho_col, N_x,Sys,HH_long)

     Y.setLabel([1,2,3,4,5])
     Y.permute([1,2,3,4,5],4)
     Y.transpose()
     Y.setLabel([5,1,2,3,4])
     Y.permute([1,2,3,4,5],4)

     blk_qnums = Y.blockQnum()
     for qnum in blk_qnums:
      Y_m=Y.getBlock(qnum)
      svd=Y_m.svd()
      temporary_matrix=svd[0]*svd[2]
      temporary_update=+1.0*temporary_matrix
      Q_list[i_ind][j_ind].putBlock(qnum,temporary_update)
     Q_list[i_ind][j_ind].setLabel([1,2,3,4,5])
     E2=copy.copy(E1)


 Store_Q_list(Q_list, N_x)
#####################################################################################









def Q_cost_val_middle_single( H_col, H_row, Q_list, i, j, N_x, Sys, H_long):
  h_long=H_long[0][0]
  result=0
  #print h_long
  Swap=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )
  Swap1=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )

############################   col   #####################################
  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-2,-3,2,3])

  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,-2,-3,4,5],5)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-1,2,-3,4,-5])


  H_col[2*i][2*j].setLabel([-1,-3,1,3])
  bdi_I=uni10.Bond(uni10.BD_IN, Q_trans.bond(4).Qlist())
  bdo_O=uni10.Bond(uni10.BD_OUT, Q_trans.bond(4).Qlist())

  Iden_t=uni10.UniTensor([bdi_I,bdo_O])
  Iden_t.permute([0,1],1)
  Iden_t.identity()
  Iden_t.setLabel([-5,5])

  result1=((Q_list_ferm*H_col[2*i][2*j])*Q_trans)*Iden_t
  result=result1[0]+result
################################################

  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-3,-4,3,4])
  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,2,-3,-4,5],5)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([1,-2,3,-4,-5])



  H_col[2*i+1][2*j].setLabel([-2,-4,2,4])

  result1=(((Q_list_ferm*H_col[2*i+1][2*j])*Q_trans))*Iden_t

  result=result1[0]+result
  #print "11",result1

#######################################################

#################### Row ####################################




  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,-8,-9])

  H_row[2*i][2*j+1].setLabel([-8,-9,8,9])

  result1=((Q_list[i][j]*H_row[2*i][2*j+1]*Q_transN))*Iden_t
  result=result1[0]+result


  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,2,3,4])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,-7,8,9])

  H_row[2*i][2*j].setLabel([-6,-7,6,7])

  result1=((Q_list[i][j]*H_row[2*i][2*j]*Q_transN))*(Iden_t)
  result=result1[0]+result

  return result





def  Grad_Q_middle_single( H_col, H_row, Q_list, i, j, N_x, Sys, H_long):

  h_long=H_long[0][0]

  Q_list[i][j].setLabel([0,1,2,3,4])
  Q_list[i][j].permute([0,1,2,3,4],4)
  result=Q_list[i][j]*0.0
  Swap=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )
  Swap1=fermionicOPT( Sys, Q_list[0][0].bond(0), Q_list[0][0].bond(1) )

  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-2,-3,2,3])

  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,-2,-3,4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,-1,2,-3,4])
  Q_trans.permute([-1,2,-3,4,-5],4)


  #print Q_trans.printDiagram()

  H_col[2*i][2*j].setLabel([-1,-3,1,3])

  bdi_I=uni10.Bond(uni10.BD_IN, Q_trans.bond(4).Qlist())
  bdo_O=uni10.Bond(uni10.BD_OUT, Q_trans.bond(4).Qlist())

  Iden_t=uni10.UniTensor([bdi_I,bdo_O])
  Iden_t.permute([0,1],1)
  Iden_t.identity()
  Iden_t.setLabel([-5,5])

  result1=(((H_col[2*i][2*j])*Q_trans))*(Iden_t)
  result1.permute([1,2,3,4,5],4)
  #print result1.printDiagram()


  Swap_t=fermionicOPT( Sys, result1.bond(1), result1.bond(2) )
  Swap_t.setLabel([-2,-3,2,3])

  result1=Swap_t*result1
  result1.permute([1,-2,-3,4,5],4)
  result1.setLabel([1,2,3,4,5])
  #print result1.printDiagram(), result.printDiagram()
  result=result1

 ################################################


  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-3,-4,3,4])
  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,2,-3,-4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.transpose()
  Q_trans.setLabel([-5,1,-2,3,-4])



  H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
  result1=((H_col[2*i+1][2*j])*Q_trans)*Iden_t
  result1.permute([1,2,3,4,5],4)

  Swap_t=fermionicOPT( Sys, result1.bond(1), result1.bond(2) )
  Swap_t.setLabel([-3,-4,3,4])

#  Swap.setLabel([-3,-4,3,4])
#  print result1.printDiagram()


  result1=Swap_t*result1
  result1.permute([1,2,-3,-4,5],4)
  result1.setLabel([1,2,3,4,5])
#  print result1.printDiagram()

  result=result1+result

##################### Row ####################################

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,6,7,-8,-9])
  Iden_t.setLabel([-10,10])

  H_row[2*i][2*j+1].setLabel([-8,-9,8,9])
  result1=((H_row[2*i][2*j+1]*Q_transN))*Iden_t
  result1.permute([6,7,8,9,10],4)  
  result=result1+result

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.transpose()
  Q_transN.setLabel([-10,-6,-7,8,9])

  H_row[2*i][2*j].setLabel([-6,-7,6,7])

  result1=(H_row[2*i][2*j]*Q_transN)*Iden_t
  result1.permute([6,7,8,9,10],4)
  result=result1+result

  return result



def  Q_update_init(H_col, H_row, N_x, D, accuracy, d, Q_list, Sys, H_long):

 HH_col=[None]*(2*N_x)
 for i in xrange(2*N_x):
   HH_col[i]=[None]*(2*N_x-1)

 for i in xrange(2*N_x):
  for j in xrange(2*N_x-1):

   HH_col[i][j]=H_col[i][j]*1.0
   list_energy1=[]
   blk_qnums = HH_col[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_col[i][j].getBlock(qnum)
    e=h_ham.eigh()
    N_col=e[0].col()
    for iter in xrange(N_col):
      list_energy1.append(e[0][iter])

   list_energy1=sorted(list_energy1, key=float)
   blk_qnums = HH_col[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_col[i][j].getBlock(qnum)
    iden_m=h_ham*1.0
    iden_m.identity()
    h_ham=h_ham+(-1.0)*iden_m*abs(max(list_energy1,key=abs))
    HH_col[i][j].putBlock(qnum,h_ham)
   #print "i,j",i,j,HH_col[i][j] 

 HH_row=[None]*(2*N_x-1)
 for i in xrange(2*N_x-1):
   HH_row[i]=[None]*(2*N_x)

 for i in xrange(2*N_x-1):
  for j in xrange(2*N_x):
   list_energy1=[]
   HH_row[i][j]=H_row[i][j]*1.0
   blk_qnums = HH_row[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_row[i][j].getBlock(qnum)
    e=h_ham.eigh()
    N_col=e[0].col()
    for iter in xrange(N_col):
      list_energy1.append(e[0][iter])

   list_energy1=sorted(list_energy1, key=float)
   blk_qnums = HH_row[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_row[i][j].getBlock(qnum)
    iden_m=h_ham*1.0
    iden_m.identity()
    h_ham=h_ham+(-1.0)*iden_m*abs(max(list_energy1,key=abs))
    HH_row[i][j].putBlock(qnum,h_ham)


 HH_long=[None]*(2*N_x)
 for i in xrange(2*N_x):
   HH_long[i]=[None]*(2*N_x)

 for i in xrange(2*N_x):
  for j in xrange(2*N_x):
   list_energy1=[]
   #print i,j, H_long[i][j]
   HH_long[i][j]=H_long[i][j]*1.0
   blk_qnums = HH_long[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_long[i][j].getBlock(qnum)
    e=h_ham.eigh()
    N_col=e[0].col()
    for iter in xrange(N_col):
      list_energy1.append(e[0][iter])

   list_energy1=sorted(list_energy1, key=float)
   blk_qnums = HH_long[i][j].blockQnum()
   for qnum in blk_qnums:
    h_ham=HH_long[i][j].getBlock(qnum)
    iden_m=h_ham*1.0
    iden_m.identity()
    h_ham=h_ham+(-1.0)*iden_m*abs(max(list_energy1,key=abs))
    HH_long[i][j].putBlock(qnum,h_ham)



 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  for j_ind in xrange(N_x):

   Q_init=Q_list[i_ind][j_ind]*1.0
   E2=10
   for iter in xrange(20):
     E1=Q_cost_val_middle_single(HH_col, HH_row, Q_list, i_ind, j_ind , N_x,Sys,HH_long)
     if abs(E1)< 1.0e-14:
       #print "E1< 1.0e-14"
       break
     #print iter, E1, E2, abs((E2-E1)/E1)
     if abs(E1)>abs(E2) or iter is 0:
      Q_init=Q_list[i_ind][j_ind]*1.0
      if abs((E2-E1)/E1) < 1.0e-11:
       #print E2, E1, abs((E2-E1)/E1), i
       break
     else:
      print 'Notoptimized=i, E1, E2=', iter,'  ', E1, '   ', E2
      Q_list[i_ind][j_ind]=Q_init*1.0
      break
     Y=Grad_Q_middle_single( HH_col, HH_row, Q_list, i_ind, j_ind, N_x,Sys,HH_long)
     #print Y.printDiagram()
     Y.setLabel([1,2,3,4,5])
     Y.transpose()
     Y.permute([1,2,3,4,5],4)
     blk_qnums = Y.blockQnum()
     for qnum in blk_qnums:
      Y_m=Y.getBlock(qnum)
      svd=Y_m.svd()
      temporary_matrix=svd[0]*svd[2]
      temporary_update=+1.0*temporary_matrix
      Q_list[i_ind][j_ind].putBlock(qnum,temporary_update)
     Q_list[i_ind][j_ind].setLabel([1,2,3,4,5])
     E2=copy.copy(E1)



 return Q_list

















def Energy_cal(PEPS_listten, d, chi_single, N_x, D, H_col, H_row, Sys):

 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_temp=[None]*N_x
 mps_boundry_up=[None]*(N_x+1)
 mps_boundry_down=[None]*N_x


 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_x*2)
 for i_ind in xrange(N_x*2):
  T.identity()
  mps_I[i_ind]=T*1.0

##################   Col   #####################

 for Location in reversed(xrange(1,N_x)):
  mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_single, N_x,Sys)

 for Location in xrange(N_x-1):
   mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_single, N_x,Sys)

 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  if i_ind==0:
   energy_col(mps_I, mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_x,  H_col[i_ind], D, E_coulmn_t,Sys)
  elif i_ind==N_x-1:
   energy_col(mps_boundry_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_x, H_col[i_ind], D, E_coulmn_t,Sys)
  else:
   energy_col(mps_boundry_left[i_ind-1], mps_boundry_right[i_ind+1], PEPS_listten[i_ind], N_x, H_col[i_ind], D, E_coulmn_t,Sys)

  if  i_ind==(N_x-1):break


 E_c=sum(E_coulmn_t)
 E_coulmn=[   E_coulmn_t[i]*1.0   for i in xrange(len(E_coulmn_t))   ]


 #print "\n"
#################   Row   #################
 for Location in xrange(N_x):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_down[Location]=make_Env_singleLayer_down( peps_l, Location, mps_boundry_down[Location-1], d, chi_single, N_x,Sys)


 for Location in reversed(xrange(1,N_x)):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][Location])
  mps_boundry_up[Location]=make_Env_singleLayer_up( peps_l, Location, mps_boundry_up[Location+1], d, chi_single, N_x ,Sys)


 E_row_t=[]
 for i_ind in xrange(N_x):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])
  H_l=[]
  for i in xrange(N_x-1):
   H_l.append(H_row[i][i_ind])

  if i_ind==0:
   energy_row(mps_I, mps_boundry_up[i_ind+1], peps_l,N_x, H_l, D, E_row_t,Sys)
  elif i_ind==N_x-1:
   energy_row(mps_boundry_down[i_ind-1], mps_I, peps_l, N_x, H_l, D, E_row_t,Sys)
  else:
   energy_row(mps_boundry_down[i_ind-1], mps_boundry_up[i_ind+1], peps_l, N_x, H_l, D, E_row_t,Sys)


  if  i_ind==(N_x-1):break

 E_r=sum(E_row_t)
 #print "E=", E_c, E_r
 E_1=E_c+E_r
 
 return E_1


########@profile
def  energy_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x, H , D, E_coulmn,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[2*i+1]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[2*i+1]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=((mps_boundry_up[2*i+1]*Swap1))*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[2*i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[2*i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[2*i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[2*i]*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

# for i in xrange(len(PEPS_listten)-1):
#  E_list_left[i].setLabel([1,2,3,4])
#  E_list_right[i+1].setLabel([1,2,3,4])
#  A=E_list_left[i]*E_list_right[i+1]
#  print "Yoho1", A[0]

 for i in xrange(N_x-1):


  E_val=local_energy_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, H[i], D,Sys)
  E_coulmn.append(E_val)



########@profile
def  local_energy_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, H_orig, D,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)


 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)


 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],2)
 r_u.setLabel([-100,3,9])
 r_u.permute([-100,3,9],2)

 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)



 #####################################################
 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)

 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys, A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys, A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel( [8,9,-1,-2] )
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 Swap1=fermionicOPT(Sys,A_conj.bond(2), A_conj.bond(4))
 Swap1.setLabel([3,10,-3,-10])
 A_conj=Swap1*A_conj
 A_conj.permute([-6,-7,3,-9,10],3)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 A_conj.permute([-6,-7,-10,-3,-9],3)

 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  q_d, V, s=TU.setTruncation(A_conj, row)
 else:
  q_d, V, s=TU.setTruncation(A_conj, colm)


 s.setLabel([1,0])
 V.setLabel([0,-3,-9])
 r_d=V*s
 r_d.permute([1,-3,-9],2)

 q_d.setLabel([-6,-7,-10,-200])
 q_d.permute([-6,-7,-10,-200],4)

 r_d.setLabel([-200,-3,-9])
 r_d.permute([-200,-3,-9],2)



 A=Peps_2*1.0
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],2)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 l_u.setLabel([9,13,-300])
 l_u.permute([9,13,-300],2)


 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)

 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys, A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys, A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel( [8,9,-1,-2] )
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 A_conj.setLabel([-9,-10,-13,-14,-15])
 A_conj.permute([-9,-13,-10,-14,-15],2)


 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  U, qq_d, s=TU.setTruncation(A_conj, row)
 else:
  U, qq_d, s=TU.setTruncation(A_conj, colm)

 s.setLabel([0,1])
 U.setLabel([-9,-13,0])
 l_d=U*s
 l_d.permute([-9,-13,1],2)

 qq_d.setLabel([-400,-10,-14,-15])
 qq_d.permute([-400,-10,-14,-15],4)


 l_d.setLabel([-9,-13,-400])
 l_d.permute([-9,-13,-400],2)


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 #print Swap4.printDiagram()
 Swap4.setLabel([23,400,10,-400])


######################################################################

 mps_boundry_up[Location*2].setLabel([16,10,26])
 mps_boundry_up[Location*2+1].setLabel([26,25,27])
 mps_boundry_up[Location*2+2].setLabel([27,15,28])
 mps_boundry_up[Location*2+3].setLabel([28,24,21])

 mps_boundry_down[Location*2].setLabel([18,22,31])
 mps_boundry_down[Location*2+1].setLabel([31,-7,32])
 mps_boundry_down[Location*2+2].setLabel([32,23,33])
 mps_boundry_down[Location*2+3].setLabel([33,-10,19])


# mps_boundry_left[Location].setLabel([16,6,-60,18])
# mps_boundry_left[Location+1].setLabel([18,11,-110,25])


# mps_boundry_right[Location].setLabel([17,90,-9,19])
# mps_boundry_right[Location+1].setLabel([19,140,-14,24])

######################################################

 A=E_left*mps_boundry_up[Location*2]
 A=A*(Swap1*q)
 A=A*mps_boundry_down[2*Location]
 A=A*mps_boundry_up[2*Location+1]
 A=A*(Swap3*q_d)
 A=A*mps_boundry_down[2*Location+1]


 B=E_right*mps_boundry_up[2*Location+3]
 B=B*(Swap2*qq_d)
 B=B*mps_boundry_down[2*Location+3]

 B=B*mps_boundry_up[2*Location+2]
 B=B*(Swap4*qq)
 B=B*mps_boundry_down[2*Location+2]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 #N_ten=N_Positiv(N_ten)
######################################################



######################################################

# l_up=copy.copy(l_u)
# r_up=copy.copy(r_u)
# l_dp=copy.copy(l_up)
# r_dp=copy.copy(r_up)
# l_dp.transpose()
# r_dp.transpose()
# l_dp.setLabel([-400,-10,-13])
# r_dp.setLabel([-10,-200,-3 ])

 
 iden_h=copy.copy(H_orig)
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])
 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])

 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig

 #print "Fnit_norm", Norm_h

 #print  Norm_h[0], h_h[0]/Norm_h[0]

 return  h_h[0]/Norm_h[0]




########@profile
def  energy_col( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x, H , D, E_coulmn,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([13])
   mps_left=mps_boundry_left[2*i+1]*iden
   iden.setLabel([15])
   mps_right=mps_boundry_right[2*i+1]*iden

   iden.setLabel([5])
   iden1=iden*1.0
   iden1.setLabel([11])


   E_list_up[i]=Peps_list*iden
   E_list_up[i]=E_list_up[i]*(Swap1*iden1)
   E_list_up[i]=E_list_up[i]*mps_left
   E_list_up[i]=E_list_up[i]*mps_right
   E_list_up[i]=E_list_up[i]*mps_boundry_left[2*i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=E_list_up[i]*mps_boundry_right[2*i]
   E_list_up[i]=E_list_up[i]*((Peps_list_conj*Swap3)*Swap2)

   E_list_up[i].permute([12,10,9,14],4)
  else:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])
   E_list_up[i+1].setLabel([13,5,11,15])

   E_list_up[i]=((Peps_list*Swap1))*E_list_up[i+1]
   E_list_up[i]=mps_boundry_left[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_right[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_left[2*i]*E_list_up[i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=mps_boundry_right[2*i]*E_list_up[i]
   E_list_up[i]=((Peps_list_conj*Swap3)*Swap2)*E_list_up[i]
   E_list_up[i].permute([12,10,9,14],4)





 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[2*i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[2*i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])


   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])


   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)



# for i in xrange(len(PEPS_listten)-1):
#  E_list_down[i].setLabel([1,2,3,4])
#  E_list_up[i+1].setLabel([1,2,3,4])
#  A=E_list_down[i]*E_list_up[i+1]
#  print "Yoho", A[0]



 for i in xrange(N_x-1):


  E_val=local_energy( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H[i],D,Sys)
  E_coulmn.append(E_val)

########@profile
def  local_energy( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, H_orig, D,Sys):

 bdi_mid=uni10.Bond( uni10.BD_IN, 1)
 iden=uni10.UniTensor( [bdi_mid] )
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],2)


 r_u.setLabel([-100,3,10])


 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


#####################################################

 A=Peps_1*1.0

 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],3)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 A_conj.permute([-6,-7,-9,-3,-10],3)



 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  q_d, V, s=TU.setTruncation(A_conj, row)
 else:
  q_d, V, s=TU.setTruncation(A_conj, colm)


 s.setLabel([1,0])
 V.setLabel([0,-3,-10])
 r_d=V*s
 r_d.permute([1,-3,-10],3)
 
 q_d.setLabel([-6,-7,-9,-200])
 q_d.permute([-6,-7,-9,-200],4)

 r_d.setLabel([-200,-3,-10])



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)

 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])


 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 Swap1=fermionicOPT(Sys,A_conj.bond(2), A_conj.bond(0))
 Swap1.setLabel([11,10,3,8])
 A_conj=A_conj*Swap1

 A_conj.permute([10,9,11,-6,7],3)
 A_conj.setLabel([-11,-10,-13,-14,-15])
 A_conj.permute([-10,-13,-11,-14,-15],2)


 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  U, qq_d, s=TU.setTruncation(A_conj, row)
 else:
  U, qq_d, s=TU.setTruncation(A_conj, colm)

 s.setLabel([0,1])
 U.setLabel([-10,-13,0])
 l_d=U*s
 l_d.permute([-10,-13,1],3)

 qq_d.setLabel([-400,-11,-14,-15])
 qq_d.permute([-400,-11,-14,-15],4)


 l_d.setLabel([-10,-13,-400])


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([ 90, 200, 9, -200 ])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])


######################################################################

 mps_boundry_left[Location*2].setLabel([16,-60,18])
 mps_boundry_left[Location*2+1].setLabel([18,6,21])
 mps_boundry_left[Location*2+2].setLabel([21,-110,22])
 mps_boundry_left[Location*2+3].setLabel([22,11,25])

 mps_boundry_right[Location*2].setLabel([17,-9,19])
 mps_boundry_right[Location*2+1].setLabel([19,90,20])
 mps_boundry_right[Location*2+2].setLabel([20,-14,23])
 mps_boundry_right[Location*2+3].setLabel([23,140,24])


######################################################

 A=E_left*mps_boundry_left[Location*2]
 A=A*(Swap1*q_d)
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_left[2*Location+1]
 A=A*(Swap3*q)
 A=A*mps_boundry_right[2*Location+1]


 B=E_right*mps_boundry_left[2*Location+3]
 B=B*(Swap2*qq)
 B=B*mps_boundry_right[2*Location+3]

 B=B*mps_boundry_left[2*Location+2]
 B=B*(Swap4*qq_d)
 B=B*mps_boundry_right[2*Location+2]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 #N_ten=N_Positiv(N_ten)
######################################################
# l_up=copy.copy(l_u)
# r_up=copy.copy(r_u)
# l_dp=copy.copy(l_up)
# r_dp=copy.copy(r_up)
# l_dp.transpose()
# r_dp.transpose()
# l_dp.setLabel([-400,-10,-13])
# r_dp.setLabel([-10,-200,-3 ])

 
 iden_h=copy.copy(H_orig)
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])
 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])

 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig

 #print "Fnit_norm", Norm_h

 #print   Norm_h[0], h_h[0]/Norm_h[0]

 return  h_h[0]/Norm_h[0]



def  update_env_LR(PEPS_listten, N_x, d, chi_single,mps_boundry_left, mps_boundry_right, mps_boundry_temp):
 for Location in xrange(N_x-1):
  PEPS_ten=rotate(PEPS_listten[N_x-Location-1])
  mps_boundry_temp[Location]=make_Env_singleLayer(PEPS_ten, Location, mps_boundry_temp[Location-1],d,chi_single,N_x)
  mps_boundry_right[N_x-Location-1]=copy.copy(mps_boundry_temp[Location])
  for Location in xrange(N_x-1):
   mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_single, N_x)


 return mps_boundry_right, mps_boundry_left



def Normalize_PEPS(PEPS_listten, N_x, D, chi_try, d, threshold, interval,Sys):

 norm_val=Cal_norm( PEPS_listten, N_x, D, chi_try, d,Sys)
 #print "Zero_order", norm_val

 count=0

 while count<100:
  #print norm_val, count
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold[1]:
   for i in xrange(N_x):
    for j in xrange(N_x):
     PEPS_listten[i][j]=PEPS_listten[i][j]*alpha
   norm_val=Cal_norm(PEPS_listten,N_x,D,chi_try,d,Sys)

  alpha=1.0+interval
  if abs(norm_val)<(threshold[0]):
   for i in xrange(N_x):
    for j in xrange(N_x):
     PEPS_listten[i][j]=PEPS_listten[i][j]*alpha
   norm_val=Cal_norm(PEPS_listten,N_x,D,chi_try,d,Sys)

 #PEPS_listten=All_dist(PEPS_listten,N_x, D)
 norm_val=Cal_norm(PEPS_listten,N_x,D,chi_try,d,Sys)
 #print "Fixed norm", abs(norm_val)
 return PEPS_listten, abs(norm_val), count


def Cal_norm( PEPS_listten, N_x, D, chi_try, d,Sys):

 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*(N_x+1)
 mps_boundry_temp=[None]*N_x

 for Location in xrange(N_x-1):
  mps_boundry_left[Location]=make_Env_singleLayer(PEPS_listten[Location], Location, mps_boundry_left[Location-1], d, chi_try, N_x,Sys)

 Location=N_x-1
 mps_boundry_right[Location]=make_Env_singleLayer_right( PEPS_listten[Location], Location, mps_boundry_right[Location+1], d, chi_try, N_x,Sys)


 norm_val=mps_boundry_left[N_x-2].product_nonsymm(mps_boundry_right[N_x-1])

 return  norm_val






def MaxAbs(c):
 blk_qnums = c.blockQnum()
 max_list=[]
 for qnum in blk_qnums:
    c_mat=c.getBlock(qnum)
    max_list.append(c_mat.absMax())
 max_list_f=[abs(x) for x in max_list]
 return max(max_list_f)

def max_ten(a):
 Max_val=MaxAbs(a)
 if ( Max_val < 0.5e-1) or (Max_val > 0.5e+1)   :

  if Max_val >= 1:
   print ">1",Max_val
   a=a*(1.00/Max_val)
  if Max_val < 1: 
   print "<1",Max_val
   a=a*(1.00/Max_val)

#  if Max_val >= 1:
#   print ">1",Max_val, Max_val**(1./2.)
#   a=a*(1.00/(Max_val**(1./2.)))
#  if Max_val < 1: 
#   print "<1",Max_val
#   a=a*(1.00/Max_val)

 else: a=a;
 return a

def All_dist(PEPS_listten,N_x, N_y, D):

 for i in xrange(N_x-1):
  for j in xrange(N_y):
   PEPS_listten[i][j], PEPS_listten[i+1][j]=equall_dis_H(PEPS_listten[i][j], PEPS_listten[i+1][j],D)
   PEPS_listten[i][j]=max_ten(PEPS_listten[i][j])
   PEPS_listten[i+1][j]=max_ten(PEPS_listten[i+1][j])
 
 for i in xrange(N_x):
  for j in xrange(N_y-1):
   PEPS_listten[i][j], PEPS_listten[i][j+1]=equall_dis_V(PEPS_listten[i][j], PEPS_listten[i][j+1],D)
   PEPS_listten[i][j]=max_ten(PEPS_listten[i][j])
   PEPS_listten[i][j+1]=max_ten(PEPS_listten[i][j+1])
 
 
 
 return PEPS_listten

def equall_dis_H(PEPS_1, PEPS_2, D):

 A=copy.copy(PEPS_1)
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,5,3,4],3)
 
 #q, r= qr_parity1(A) 
 q,s,V=svd_parity5(A)
 V.setLabel([0,3,4])
 s.setLabel([-100,0])
 V=s*V
 V.permute([-100,3,4],2)
 r=V*1.0
 r.setLabel([-100,3,4])
 q.setLabel([1,2,5,-100])



 A=copy.copy(PEPS_2)
 A.setLabel([4,5,6,7,8])
 A.permute([4, 6, 5,8,7],2)
 U,s,qq=svd_parity6(A) 
 U.setLabel([4,6,0])
 s.setLabel([0,-200])
 U=U*s
 U.permute([4,6,-200],2)
 l=U*1.0
 qq.setLabel([-200,5,8,7])
 l.setLabel([4,6,-200])

 Teta=l*r
 Teta.permute([-100,3,-200,6],2)
 U, V, s= setTruncation3(Teta, D)

 #U,s,V=svd_parity2(Teta)

 U.setLabel([-100,3,17])
 s.setLabel([17,-17])
 V.setLabel([-17,-200,6])
 s=Sqrt(s)
 s.setLabel([17,-17])
 
 U=U*s
 V=s*V

 U.permute([-100,3,-17],1)
 U.setLabel([-100,3,4])
 V.permute([17,-200,6],1)
 V.setLabel([4,-200,6])
 
 PEPS_1=q*U
 PEPS_2=qq*V

 PEPS_1.permute([1,2,3,4,5],3)
 PEPS_2.permute([4,5,6,7,8],3)

 return PEPS_1, PEPS_2


def equall_dis_V(PEPS_1, PEPS_2, D):

 A=copy.copy(PEPS_1)
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,4,3,5],3)
 
 #q, r= qr_parity1(A) 
 q,s,V=svd_parity5(A)
 V.setLabel([0,3,5])
 s.setLabel([-100,0])
 V=s*V
 V.permute([-100,3,5],2)
 r=V*1.0
 r.setLabel([-100,3,5])
 q.setLabel([1,2,4,-100])



 A=copy.copy(PEPS_2)
 A.setLabel([6,5,7,8,9])
 A.permute([5, 7, 6,8,9],2)
 U,s,qq=svd_parity6(A) 
 U.setLabel([5,7,0])
 s.setLabel([0,-200])
 U=U*s
 U.permute([5,7,-200],2)
 l=U*1.0
 qq.setLabel([-200,6,8,9])
 l.setLabel([5,7,-200])

 Teta=l*r
 Teta.permute([-100,3,-200,7],2)
 U, V, s= setTruncation3(Teta, D)

 #U,s,V=svd_parity2(Teta)

 U.setLabel([-100,3,17])
 s.setLabel([17,-17])
 V.setLabel([-17,-200,7])
 s=Sqrt(s)
 s.setLabel([17,-17])
 
 U=U*s
 V=s*V

 U.permute([-100,3,-17],1)
 U.setLabel([-100,3,5])
 V.permute([17,-200,7],1)
 V.setLabel([5,-200,7])
 
 PEPS_1=q*U
 PEPS_2=qq*V

 PEPS_1.permute([1,2,3,4,5],3)
 PEPS_2.permute([6,5,7,8,9],3)

 return PEPS_1, PEPS_2

def  norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H):
 val1=((((r_up*r_dp)*iden_h)*N_ten)*l_up)*l_dp
 #print val1
 #val0=((((r_u*r_d)*H)*N_ten)*l_u)*l_d


 val2=((((r_up*r_d)*Ham)*N_ten)*l_up)*l_d
 val3=((((r_u*r_dp)*Ham)*N_ten)*l_u)*l_dp

 #print "2, 3", val2, val3
 
 return val1[0]-val3[0]-val2[0]#+val0[0]

#################################@profile#####################################
def  optimum_0(N_ten, l_u, r_u, r_d, l_d , l_up, r_up,l_dp, r_dp, Ham,iden_h,H):

 #l_up.setLabel()

 Env_r=((l_dp)*(l_up*iden_h))*N_ten
 Env_r.permute([ -200,-3,-10, -100,3,10],3)

 Env_r1=copy.copy(Env_r)
 Env_r1.transpose()
 Env_r=Env_r+Env_r1



 Env_s=(r_u*N_ten)*((l_u*Ham)*l_dp)
 Env_s.permute([-200,-3,-10],3)

 Env_s1=((((r_d)*Ham)*N_ten)*l_up)*l_d

 Env_s1.permute([-100,3,10],0)
 Env_s1.transpose()
 Env_s=Env_s+Env_s1




 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([5,6,7,8])
 V.setLabel([1,2,3,4])
 S.setLabel([4,5])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,6,7,8],3)
 A2_inv.setLabel([-100,3,10,-200,-3,-10])

 
 r_up=A2_inv*Env_s
 r_up.permute([-100,3,10],3)
 r_up.setLabel([-100,3,10])

 r_dp=r_up*1.0
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])
 return r_up, r_dp


#########@profile
def  optimum_1(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H):
  
 Env_r=((N_ten)*(r_up*r_dp))*iden_h
 Env_r.permute([-10,-13,-400,10,13,-300],3)

 Env_r1=copy.copy(Env_r)
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=(((r_dp*r_u)*N_ten)*(l_u))*Ham
 Env_s.permute([-10,-13,-400],3)

 Env_s1=((((r_up*r_d)*Ham)*N_ten))*l_d
 Env_s1.permute([10,13,-300],0)
 Env_s1.transpose()
 Env_s=Env_s+Env_s1


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)

 U.transpose()
 V.transpose()
 S=inverse(S)


 U.setLabel([5,6,7,8])
 V.setLabel([1,2,3,4])
 S.setLabel([4,5])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,6,7,8],3)
 A2_inv.setLabel([10,13,-300,-10,-13,-400])


 l_up=A2_inv*Env_s
 l_up.permute([10,13,-300],3)

 l_dp=l_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])

 return l_up, l_dp


def  optimum_11(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H):

 Env_r=((N_ten)*(r_up*r_dp))*iden_h
 Env_r.permute([-10,13,-300,10,-13,-400],3)

 Env_s=(((r_dp*r_u)*N_ten)*(l_u))*Ham
 Env_s.permute([-10,13,-300],2)

 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)

# U, S, V=svd_parityrl(Env_r)

 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([5,6,7,8])
 V.setLabel([1,2,3,4])
 S.setLabel([4,5])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,6,7,8],3)
 A2_inv.setLabel([10,-13,-400,-10,13,-300])

 l_up=A2_inv*Env_s
 l_up.permute([10,-13,-400],3)
 l_dp=copy.copy(l_up)
 l_dp.transpose()
 l_dp.setLabel([-10,13,-300])

 return l_up, l_dp











def  update_twotensor( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x, U_ham, H , D, E_coulmn,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([13])
   mps_left=mps_boundry_left[2*i+1]*iden
   iden.setLabel([15])
   mps_right=mps_boundry_right[2*i+1]*iden

   iden.setLabel([5])
   iden1=iden*1.0
   iden1.setLabel([11])

   E_list_up[i]=Peps_list*iden
   E_list_up[i]=E_list_up[i]*(Swap1*iden1)
   E_list_up[i]=E_list_up[i]*mps_left
   E_list_up[i]=E_list_up[i]*mps_right
   E_list_up[i]=E_list_up[i]*mps_boundry_left[2*i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=E_list_up[i]*mps_boundry_right[2*i]
   E_list_up[i]=E_list_up[i]*((Peps_list_conj*Swap3)*Swap2)

   E_list_up[i].permute([12,10,9,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])
   E_list_up[i+1].setLabel([13,5,11,15])

   E_list_up[i]=((Peps_list*Swap1))*E_list_up[i+1]
   E_list_up[i]=mps_boundry_left[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_right[2*i+1]*E_list_up[i]
   E_list_up[i]=mps_boundry_left[2*i]*E_list_up[i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=mps_boundry_right[2*i]*E_list_up[i]
   E_list_up[i]=((Peps_list_conj*Swap3)*Swap2)*E_list_up[i]
   E_list_up[i].permute([12,10,9,14],4)





 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[2*i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[2*i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])


   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])


   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)


# for i in xrange(len(PEPS_listten)-1):
#  E_list_down[i].setLabel([1,2,3,4])
#  E_list_up[i+1].setLabel([1,2,3,4])
#  A=E_list_down[i]*E_list_up[i+1]
#  print "Inside", A[0]


 for i in xrange(N_x-1):
  #print "i", i

#  PEPS_f, PEPS_s, E_val=Update_twotensor_local( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, U_ham[i], H[i],D,threshold, interval,Sys)


  if Sys[2]=="QR":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, U_ham[i], H[i], D, threshold, interval, Sys)
  if Sys[2]=="Inv":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local_inv( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, U_ham[i], H[i], D, threshold, interval, Sys)
  if Sys[2]=="Grad":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local_grad( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, U_ham[i], H[i], D, threshold, interval, Sys)



  E_coulmn.append(E_val)

  PEPS_listten[i]=PEPS_f*1.0
  PEPS_listten[i+1]=PEPS_s*1.0

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[2*i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[2*i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])

   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[2*i].setLabel([12,-8,130])
   mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[2*i].setLabel([14,-6,150])
   mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])

   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)





#######@profile
def  Update_twotensor_local(PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, Ham, H_orig, D,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=copy.copy(PEPS_listten[Location])
 Peps_2=copy.copy(PEPS_listten[Location+1])

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],3)


 r_u.setLabel([-100,3,10])

 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,9,200])
 q_d.permute([6,7,9,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-9,-200,9,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-9,-200],4)

#####################################################

 A=Peps_2*1.0
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])


 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,11,14,15])
 qq_d.permute([400,11,14,15],4)
 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel( [-400, -11, 400, 11] )
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-11,-14,-15],4)



######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([90,200,9,-200])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])


######################################################################

 mps_boundry_left[Location*2].setLabel([16,-60,18])
 mps_boundry_left[Location*2+1].setLabel([18,6,21])
 mps_boundry_left[Location*2+2].setLabel([21,-110,22])
 mps_boundry_left[Location*2+3].setLabel([22,11,25])

 mps_boundry_right[Location*2].setLabel([17,-9,19])
 mps_boundry_right[Location*2+1].setLabel([19,90,20])
 mps_boundry_right[Location*2+2].setLabel([20,-14,23])
 mps_boundry_right[Location*2+3].setLabel([23,140,24])

######################################################

 A=E_left*mps_boundry_left[Location*2]
 A=(A*Swap1)*q_d
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_left[2*Location+1]
 A=A*(Swap3*q)
 A=A*mps_boundry_right[2*Location+1]


 B=E_right*mps_boundry_left[2*Location+3]
 B=B*(Swap2*qq)
 B=B*mps_boundry_right[2*Location+3]

 B=B*mps_boundry_left[2*Location+2]
 B=B*(Swap4*qq_d)
 B=B*mps_boundry_right[2*Location+2]
 
 N_ten=A*B
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])



 N_ten=N_Positiv(N_ten)

######################################################

###############simple_update###########################
 A=r_u*1.0
 A.setLabel([-10,2,3])
 B=l_u*1.0
 B.setLabel([3,6,-30])
 Ham.setLabel([-2,-6,2,6])

 Teta=(A*B)*Ham
 Teta.permute([-10,-2,-6,-30],2)
 
 D_dim=l_u.bond(0).dim()

 D_list=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_list.append(dim)

 D_dim=sum(D_list)


 
 row, colm=cal_rowcol(Teta)
 if (row<=colm and row<=D_dim):
  U,V,S=TU.setTruncation(Teta,row)
 elif (row<=colm and row>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 elif (row>colm and colm<=D_dim):
  U,V,S=TU.setTruncation(Teta,colm)
 elif (row>colm and colm>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)

 U.setLabel([-10,-2,-3])
 V.setLabel([-3,-6,-30])
 S=Sqrt(S)
 S.setLabel([-3,3])
 r_up=U*S
 r_up.permute([-10,-2,3],3)
 r_up.setLabel([-100,3,10])
 S.setLabel([3,-3])
 l_up=V*S
 l_up.permute([3,-6,-30],3)
 l_up.setLabel([10,13,-300])
###########################################



# l_up=l_u*1.0
# r_up=r_u*1.0
 l_dp=l_up*1.0
 r_dp=r_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])


 Ham.setLabel([51,52,53,54])
 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([-3,-13,3,13])
 Ham.setLabel([-3,-13,3,13])
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])

 r_up_init=r_up*1.0
 r_dp_init=r_dp*1.0
 l_up_init=l_up*1.0
 l_dp_init=l_dp*1.0

 valf=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up, l_dp, r_dp, Ham, iden_h, H)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   r_up=r_up_init*1.0
   r_dp=r_dp_init*1.0
   l_up=l_up_init*1.0
   l_dp=l_dp_init*1.0
   break
  else:
   r_up_init=r_up*1.0
   r_dp_init=r_dp*1.0
   l_up_init=l_up*1.0
   l_dp_init=l_dp*1.0

  r_up, r_dp=optimum_0( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
  l_up, l_dp=optimum_1( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)

 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-12 or abs(valf)>1.0e+12:
  print "warning_norm_in_optimization",  abs(valf), "count", count


 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 h_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*H_orig
# print "E_2", h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<(threshold[0]):
  r_up,l_up=renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval)



 PEPS_1=r_up*q
 PEPS_1.permute([6,7,3,9,10],3)

 A=PEPS_1*1.0
 A.permute([6,7,3,9,10],5)

 Swap1=fermionicOPT(Sys,A.bond(2), A.bond(3))
 Swap1.setLabel([-3,-9,3,9])
 A=Swap1*A
 A.permute([6,7,-3,-9,10],3)
 PEPS_1=A*1.0


 PEPS_2=l_up*qq
 PEPS_2.permute([11,10,13,14,15],3)


 

 return  PEPS_1,  PEPS_2,  h_h[0]/Norm_h[0]




########@profile
def update_twotensor_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x, U_ham, H , D, E_coulmn,threshold, interval,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):

  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[2*i+1]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[2*i+1]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=((mps_boundry_up[2*i+1]*Swap1))*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[2*i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[2*i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[2*i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[2*i]*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)



# for i in xrange(len(PEPS_listten)-1):
#  E_list_left[i].setLabel([1,2,3,4])
#  E_list_right[i+1].setLabel([1,2,3,4])
#  A=E_list_left[i]*E_list_right[i+1]
#  print "Y1", A[0]



 for i in xrange(N_x-1):
  #print "i", i

  #PEPS_f, PEPS_s, E_val=Update_twotensor_local_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, U_ham[i], H[i],D,threshold, interval,Sys)


  if Sys[2]=="QR":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, U_ham[i], H[i],D,threshold, interval,Sys)
  elif Sys[2]=="Inv":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local_row_inv( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, U_ham[i], H[i],D,threshold, interval,Sys)
  elif Sys[2]=="Grad":
   PEPS_f, PEPS_s, E_val=Update_twotensor_local_row_grad( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, U_ham[i], H[i],D,threshold, interval,Sys)


  E_coulmn.append(E_val)


  PEPS_listten[i]=PEPS_f*1.0
  PEPS_listten[i+1]=PEPS_s*1.0

  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[2*i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[2*i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[2*i].setLabel([12,5,130])
   mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[2*i].setLabel([14,10,150])
   mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[2*i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[2*i]*E_list_left[i]
   E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(mps_boundry_up[2*i+1]*Swap1)
   E_list_left[i].permute([13,6,-6,15],4)


########@profile
def  Update_twotensor_local_row(PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, Ham, H_orig, D,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)


 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)


 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],3)
 r_u.setLabel([-100,3,9])

 r_u.setLabel([-100,3,10])   #new
 r_u.permute([-100,3,10],3)



 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)



 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,10,200])
 q_d.permute([6,7,10,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-10,-200,10,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-10,-200],4)



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)

 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],3)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,10,14,15])
 qq_d.permute([400,10,14,15],4)

 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel([-400, -10, 400, 10])
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-10,-14,-15],4)


 l_u.setLabel([9,13,-300])

 l_u.setLabel([10,13,-300])  #new
 l_u.permute([10,13,-300],3)

 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])  #new


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 Swap4.setLabel([23,400,10,-400])


######################################################################

 mps_boundry_up[Location*2].setLabel([16,10,26])
 mps_boundry_up[Location*2+1].setLabel([26,25,27])
 mps_boundry_up[Location*2+2].setLabel([27,15,28])
 mps_boundry_up[Location*2+3].setLabel([28,24,21])

 mps_boundry_down[Location*2].setLabel([18,22,31])
 mps_boundry_down[Location*2+1].setLabel([31,-7,32])
 mps_boundry_down[Location*2+2].setLabel([32,23,33])
 mps_boundry_down[Location*2+3].setLabel([33,-10,19])



######################################################

 A=E_left*mps_boundry_up[Location*2]
 A=A*(Swap1*q)
 A=A*mps_boundry_down[2*Location]
 A=A*mps_boundry_up[2*Location+1]
 A=A*(Swap3*q_d)
 A=A*mps_boundry_down[2*Location+1]


 B=E_right*mps_boundry_up[2*Location+3]
 B=B*(Swap2*qq_d)
 B=B*mps_boundry_down[2*Location+3]

 B=B*mps_boundry_up[2*Location+2]
 B=B*(Swap4*qq)
 B=B*mps_boundry_down[2*Location+2]
 
 N_ten=A*B
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 N_ten=N_Positiv(N_ten)
######################################################



######################################################

##############simple_update###########################
 A=r_u*1.0
 A.setLabel([-10,2,3])
 B=l_u*1.0
 B.setLabel([3,6,-30])
 Ham.setLabel([-2,-6,2,6])

 Teta=(A*B)*Ham
 Teta.permute([-10,-2,-6,-30],2)
 
 
 D_dim=l_u.bond(0).dim()
 D_list=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_list.append(dim)

 D_dim=sum(D_list)

 row, colm=cal_rowcol(Teta)
 if (row<=colm and row<=D_dim):
  U,V,S=TU.setTruncation(Teta,row)
 elif (row<=colm and row>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 elif (row>colm and colm<=D_dim):
  U,V,S=TU.setTruncation(Teta,colm)
 elif (row>colm and colm>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)

 U.setLabel([-10,-2,-3])
 V.setLabel([-3,-6,-30])
 S=Sqrt(S)
 S.setLabel([-3,3])
 r_up=U*S
 r_up.permute([-10,-2,3],3)
 r_up.setLabel([-100,3,10])
 S.setLabel([3,-3])
 l_up=V*S
 l_up.permute([3,-6,-30],3)
 l_up.setLabel([10,13,-300])
##########################################


 
# l_up=l_u*1.0
# r_up=r_u*1.0
 l_dp=l_up*1.0
 r_dp=r_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])  #new
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])



 Ham.setLabel([51,52,53,54])

 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([-3,-13,3,13])

 Ham.setLabel([-3,-13,3,13])
 
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])

 r_up_init=r_up*1.0
 r_dp_init=r_dp*1.0
 l_up_init=l_up*1.0
 l_dp_init=l_dp*1.0

 valf=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up, l_dp, r_dp, Ham, iden_h, H)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12 or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   r_up=r_up_init*1.0
   r_dp=r_dp_init*1.0
   l_up=l_up_init*1.0
   l_dp=l_dp_init*1.0
   break
  else:
   r_up_init=r_up*1.0
   r_dp_init=r_dp*1.0
   l_up_init=l_up*1.0
   l_dp_init=l_dp*1.0

  r_up, r_dp=optimum_0( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
  l_up, l_dp=optimum_1( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)

 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-11 or abs(valf)>1.0e+11:
   print "warning_norm_in_optimization",  abs(valf), "count", count



 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 h_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*H_orig
# print "E_2", h_h[0]/Norm_h[0]
# #print Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<(threshold[0]):
   r_up,l_up=renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval)


 r_up.setLabel([-100,3,9])   #new
 l_up.setLabel([9,13,-300])  #new


 PEPS_1=r_up*q
 PEPS_1.permute([6,7,3,9,10],3)


 PEPS_2=l_up*qq
 A=PEPS_2*1.0
 A.permute([9,10,13,14,15],5)
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([9,-10,-13,14,15],3)
 PEPS_2=A*1.0
 PEPS_2.permute([9,-10,-13,14,15],3)


 return PEPS_1, PEPS_2, h_h[0]/Norm_h[0]
















def Mat_uni_to_np(Mat_uni):
 dim0=int(Mat_uni.row())
 dim1=int(Mat_uni.col())
 Mat_np=np.zeros((dim0,dim1))
 for i in xrange(dim0):
  for j in xrange(dim1):
   Mat_np[i,j]=Mat_uni[i*dim1+j]
 return  Mat_np


def condition_number(N):
 N_mat=N.getBlock()
 A_np=Mat_uni_to_np(N_mat)
 norm_val=npLA.norm(A_np)
 #print norm_val
 A_np=A_np*(1.0/norm_val)
 val=npLA.cond(A_np) 
 return val

































def   rotate_all(PEPS_listten, N_x):

 PEPS_listten_t=[None]*N_x
 for i in xrange(N_x):
  PEPS_listten_t[i]=[None]*N_x

 for i in xrange(N_x):
  for j in xrange(N_x):
   ten=copy.copy(PEPS_listten[i][j])
   ten.setLabel([0,1,2,3,4])
   ten.permute([1,0,2,4,3],3)
   PEPS_listten_t[j][i]=ten*1.0


 return PEPS_listten_t









def sqrt_general(N2):
  N_init=copy.copy(N2)
  blk_qnums = N2.blockQnum()
  for qnum in blk_qnums:
   M=N2.getBlock(qnum)
   eig=M.eigh()
   
   e=Sqrt_mat(eig[0])
   U_trans=copy.copy(eig[1])
   U_trans.transpose()
   M=U_trans*e*eig[1]
   N_init.putBlock(qnum,M)
  return N_init
def Sqrt_mat(e):
 d=int(e.row())
 
 for q in xrange(d):
   ##print e[q] 
   if e[q] > 0:  
    e[q]=((e[q])**(1.00/2.00))
   else:  
    e[q]=0.0 
 return e  
 
def N_Positiv(N):
 N.setLabel([-200,-400,-100,-300])
 N1=copy.copy(N)
 N1.transpose()
 N=(N+N1)*(1.00/2.00)
 N1=copy.copy(N)
 N1.setLabel([2,3,-100,-300])
 N.setLabel([-200,-400,2,3])

 N=N*N1
 N.permute([-200,-400,-100,-300],2)
 N_final=sqrt_general(N)
 N_final.setLabel([-200,-400,-100,-300])
 N_final.permute([-200,-400,-100,-300], 2)
 return N_final              










#####@profile
def update_energy_eff(PEPS_listmps, U_eval, U_eval0, U_evalN, d, D, Location, N_x, chi_express):

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)
 N_y=PEPS_listmps.N
 B_list=[None]*N_y
 if Location == 0:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi1,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
 elif Location==N_x-1:
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo1,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
 else:
  A_list=[None]*N_y
  for i in xrange(N_y):
   if i == 0:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_first")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   elif i ==(N_y-1):
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_Last")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A
   else:
     A=uni10.UniTensor([bdi,PEPS_listmps[i].bond(0),bdiphy,bdo,PEPS_listmps[i].bond(2)], "A_middle")
     A.setLabel([0,1,2,3,4])
     A.permute([1,3,2,0,4],4)
     A.putBlock(PEPS_listmps[i].getBlock())
     B_list[i]=A


 B_list_u=[None]*N_y
 for i in xrange(N_y):
  B_list[i].setLabel([1,3,2,0,4])
  B_list_u[i]=copy.copy(B_list[i])
  B_list_u[i].permute([1,3,2,0,4],6)
  B_list_u[i].combineBond([3,2])
  B_list_u[i].combineBond([3,0])
  B_list_u[i].permute([1,3,4],2)

 list_bond=[]
 for q in xrange(len(B_list_u)):
   list_bond.append(B_list_u[q].bond(2).dim())

 mps_b=MPSclass.MPS(B_list_u[1].bond(1).dim(),max(list_bond)*4,len(B_list_u))

 for i in xrange(len(B_list_u)):
   mps_b[i]=copy.copy(B_list_u[i])


 Norm_init=mps_b.norm()
 #print Norm_init

 for i in xrange(N_y-1):
  row=B_list_u[i].bond(0).dim()*B_list_u[i].bond(1).dim()
  colm=B_list_u[i].bond(2).dim()
  if (row<=colm):
   U,V,s=MPSclass.setTruncation1(B_list_u[i],row)
  else:
   U,V,s=MPSclass.setTruncation1(B_list_u[i],colm)
  B_list_u[i]=copy.copy(U)
  s.setLabel([1,2])
  V.setLabel([2,3])
  V=s*V
  V.permute([1,3],1)
  B_list_u[i+1].setLabel([3,2,4])
  B_list_u[i+1]=V*B_list_u[i+1]
  B_list_u[i+1].permute([1,2,4],2)
  #A_l_can[i+1]=copy.copy(A_l[i+1])

 test_norm=B_list_u[N_y-1]*B_list_u[N_y-1]
 test_norm1=B_list_u[N_y-2]*B_list_u[N_y-2]

 #print test_norm[0]#,test_norm1[0]

#U_eval, U_eval0, U_evalN

 for i in xrange(N_y-1):
  A=uni10.UniTensor([ B_list_u[N_y-1-i].bond(0), B_list[N_y-1-i].bond(1), B_list[N_y-1-i].bond(2), B_list[N_y-1-i].bond(3), B_list_u[N_y-1-i].bond(2)])
  A.putBlock(B_list_u[N_y-1-i].getBlock())
  A.setLabel([1,3,2,0,4])

  B=uni10.UniTensor([ B_list_u[N_y-2-i].bond(0), B_list[N_y-2-i].bond(1), B_list[N_y-2-i].bond(2), B_list[N_y-2-i].bond(3), B_list_u[N_y-2-i].bond(2)])
  B.putBlock(B_list_u[N_y-2-i].getBlock())
  B.setLabel([-11,-3,-2,-10,1])

  if i==0:
    H_eff=U_evalN
    H_eff.setLabel([20,30,-2,2])
  elif i==N_y-2:
    H_eff=U_eval0
    H_eff.setLabel([20,30,-2,2])
  else:
    H_eff=U_eval
    H_eff.setLabel([20,30,-2,2])

  B.setLabel([-11,-3,-2,-10,1])
  B.permute([-11,-3,-10,-2,1],3)
  q,s,V=svd_parity5(B)
  s.setLabel([0,10])
  V.setLabel([0,2,3])
  r_u=V*s
  r_u.permute([10,2,3],2)
  r_u.setLabel([10,-2,1])
  
  q.setLabel([-11,-3,-10,10])
  r_u.setLabel([10,-2,1])

  A.setLabel([1,3,2,0,4])
  A.permute([1,2,3,0,4],2)
  U,s,qq=svd_parity6(A)
  s.setLabel([0,20])
  U.setLabel([1,2,0])
  l_u=U*s
  l_u.permute([1,2,20],2)

  qq.setLabel([20,3,0,4])
  l_u.setLabel([1,2,20])

 ##############QR_update###########################
  A=copy.copy(r_u)
  A.setLabel([-10,2,3])
  B=copy.copy(l_u)
  B.setLabel([3,6,-30])
  H_eff.setLabel([-2,-6,2,6])

  Teta=(A*B)*H_eff
  Teta.permute([-10,-2,-6,-30],2)
  Teta_mat=Teta.getBlock()
  col1=Teta_mat.col()
  row1=Teta_mat.row()
  dim_mat=1
  if col1<row1:
   dim_mat=col1
  else:
   dim_mat=row1
  #print "Info", dim_mat, chi_express
  if chi_express<dim_mat:
   U, V, S= setTruncation3(Teta, chi_express)
  else:
   U, V, S= setTruncation3(Teta, dim_mat)
  U.setLabel([-10,-2,-3])
  V.setLabel([-3,-6,-30])
  S=Sqrt(S)
  S.setLabel([3,-3])
  r_up=U*S
  r_up.permute([-10,-2,3],2)
  r_up.setLabel([10,-2,1])
  l_up=V*S
  l_up.permute([3,-6,-30],2)
  l_up.setLabel([1,2,20])

  U=q*r_up
  U.permute([-11,-3,-2,-10,1],4)
  V=qq*l_up
  V.permute([1,3,2,0,4],4)

  U.setLabel([-11, -3, 20, -10, 1])
  U.permute([-11, -3, 20, -10, 1],4)
  U.combineBond([-3,20])
  U.combineBond([-3,-10])
  U.permute([-11, -3, 1],2)

  V.setLabel([1, 3, 30, 0, 4])
  V.combineBond([3,30])
  V.combineBond([3,0])
  V.permute([1, 3, 4],2)



#  Theta=(A*H_eff)*B
#  #print Theta.printDiagram()
#  Theta.permute([-11, -3, 20, -10, 3, 30, 0, 4],4)
#  U, V, S=setTruncation2(Theta, D)
#  U.setLabel([-11, -3, 20, -10, -1])
#  S.setLabel([-1,1])
#  U=U*S
#  V.setLabel([1, 3, 30, 0, 4])

##################
#  Theta1=A*B
#  Theta2=U*V
#  Theta2.permute([-11, -3, 20, -10, 3, 30, 0, 4],4)
#  Theta2.setLabel([-11, -3, -2, -10, 3, 2, 0, 4])
##################

#  U.permute([-11, -3, 20, -10, 1],4)
#  U.combineBond([-3,20])
#  U.combineBond([-3,-10])
#  U.permute([-11, -3, 1],2)

#  V.combineBond([3,30])
#  V.combineBond([3,0])
#  V.permute([1, 3, 4],2)

#################################
#  norm_test1=Theta1*Theta1
#  norm_test2=Theta2*Theta2
#  norm_testf=Theta1*Theta2
#  print "i", i, norm_test1[0], norm_test2[0], norm_testf[0]/((norm_test2[0]*norm_test1[0])**(0.5))
##################################################
  B_list_u[N_y-1-i]=copy.copy(V)
  B_list_u[N_y-2-i]=copy.copy(U)

  B_list_u[N_y-1-i].setLabel([1,2,3])
  B_list_u[N_y-1-i].permute([2,3,1],2)


  row=B_list_u[N_y-1-i].bond(0).dim()*B_list_u[N_y-1-i].bond(1).dim()
  colm=B_list_u[N_y-1-i].bond(2).dim()
 
  if (row<=colm):
   U,V,s=MPSclass.setTruncation1(B_list_u[N_y-1-i],row)
  else:
   U,V,s=MPSclass.setTruncation1(B_list_u[N_y-1-i],colm)


  U.setLabel([2,3,1])
  U.permute([1,2,3],2)
  B_list_u[N_y-1-i]=copy.copy(U)
  s.setLabel([3,4])
  V.setLabel([4,5])
  V_f=(s*V)
  V_f.permute([3,5],1)
  B_list_u[N_y-2-i].setLabel([-1,2,5])
  B_list_u[N_y-2-i]=V_f*B_list_u[N_y-2-i]
  B_list_u[N_y-2-i].permute([-1,2,3],2)


 list_bond=[]
 for  q  in  xrange(len(B_list_u)):
   list_bond.append(B_list_u[q].bond(2).dim())

 mps_b=MPSclass.MPS( B_list_u[1].bond(1).dim(), max(list_bond), len(B_list_u))

 for  i  in  xrange(len(B_list_u)):
   mps_b[i]=copy.copy(B_list_u[i])


 #mps_b=mps_b.normalize()

 return mps_b










def matSx():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [0.0, 1.0, 1.00, 0.0])
  return Mat 

def matSz():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [1.0, 0, 0, -1.0]);
  return Mat 

def matSy():
  spin = 0.5
  dim = int(spin * 2 + 1)
  Mat=(1.0)*uni10.Matrix(dim, dim, [0.0, -1.00, 1.00, 0.00]);
  return Mat 


def matIden():
    spin_t=0.5
    dimT = int(2*spin_t + 1)
    Mat=uni10.Matrix(dimT, dimT,[1,0,0,1])
    return Mat



def  make_density_matrix_double(PEPS_listten, N_x, chi_boundry, d, D,Sys):

 rho_row=[None]*(N_x-1)
 for i in xrange(N_x-1):
  rho_row[i]=[None]*N_x

 rho_col=[None]*(N_x)
 for i in xrange(N_x):
  rho_col[i]=[None]*(N_x-1)



 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_x)
 for i_ind in xrange(mps_I.N):
  T.identity()
  mps_I[i_ind]=T*1.0


 MPO_Ten=make_boundry_MPO(PEPS_listten, N_x,Sys)
 Env_left, Env_right=make_ENVMPS(MPO_Ten, N_x, chi_boundry)


 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  if i_ind==0:
   density_double_col(mps_I, Env_right[i_ind+1], PEPS_listten[i_ind], N_x, rho_col, D, i_ind,Sys)
  elif i_ind==N_x-1:
   density_double_col(Env_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_x, rho_col, D, i_ind,Sys)
  else:
   density_double_col(Env_left[i_ind-1], Env_right[i_ind+1], PEPS_listten[i_ind], N_x, rho_col, D, i_ind,Sys)

 E_c=sum(E_coulmn_t)
 E_coulmn=[ E_coulmn_t[i]*1.0   for i in xrange(len(E_coulmn_t))  ]

#################   Row   #################


 MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,Sys)
 Env_up, Env_down=make_ENV_updown( MPO_Ten, N_x, chi_boundry)

 E_row_t=[]
 for i_ind in xrange(N_x):
  peps_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])

  if i_ind==0:
   density_double_row( mps_I, Env_up[i_ind+1], peps_l, N_x, rho_row, D, i_ind,Sys)
  elif i_ind==N_x-1:
   density_double_row( Env_down[i_ind-1], mps_I, peps_l, N_x, rho_row, D, i_ind,Sys)
  else:
   density_double_row( Env_down[i_ind-1], Env_up[i_ind+1], peps_l, N_x, rho_row, D, i_ind,Sys)





 return rho_row, rho_col



def  density_double_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x,  rho_row , D, Location_x,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):

  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[i]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[i]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   #E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=Swap1*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[i]*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)



# for i in xrange(len(PEPS_listten)-1):
#  E_list_left[i].setLabel([1,2,3,4])
#  E_list_right[i+1].setLabel([1,2,3,4])
#  A=E_list_left[i]*E_list_right[i+1]
#  print "Y1", A[0]



 for i in xrange(N_x-1):



  rho_row[i][Location_x]=density_row_double_local( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i,D,Sys)



########@profile
def  density_row_double_local( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, D,Sys ):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)


 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)


 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],3)
 r_u.setLabel([-100,3,9])

 r_u.setLabel([-100,3,10])   #new
 r_u.permute([-100,3,10],3)



 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)



 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,10,200])
 q_d.permute([6,7,10,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-10,-200,10,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-10,-200],4)



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)

 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],3)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,10,14,15])
 qq_d.permute([400,10,14,15],4)

 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel([-400, -10, 400, 10])
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-10,-14,-15],4)


 l_u.setLabel([9,13,-300])

 l_u.setLabel([10,13,-300])  #new
 l_u.permute([10,13,-300],3)

 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])  #new


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 Swap4.setLabel([23,400,10,-400])


######################################################################

 mps_boundry_up[Location].setLabel([16,10,25,27])
 #mps_boundry_up[Location*2+1].setLabel([26,25,27])
 mps_boundry_up[Location+1].setLabel([27,15,24,21])
 #mps_boundry_up[Location*2+3].setLabel([28,24,21])

 mps_boundry_down[Location].setLabel([18,22,-7,32])
 #mps_boundry_down[Location*2+1].setLabel([31,-7,32])
 mps_boundry_down[Location+1].setLabel([32,23,-10,19])
 #mps_boundry_down[Location*2+3].setLabel([33,-10,19])



######################################################

 A=E_left*mps_boundry_up[Location]
 A=A*(Swap1*q)
 A=A*mps_boundry_down[Location]
 A=A*mps_boundry_up[Location+1]
 A=A*(Swap3*q_d)
 A=A*mps_boundry_down[Location+1]


 B=E_right#*mps_boundry_up[2*Location+3]
 B=(Swap2*qq_d)*B
 #B=B*mps_boundry_down[2*Location+3]

 #B=B*mps_boundry_up[2*Location+2]
 B=B*(Swap4*qq)
 #B=B*mps_boundry_down[2*Location+2]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])


 N_ten=(((r_u*r_d)*N_ten)*(l_u*l_d))
 N_ten.permute([-3,-13,3,13],2)


 return   N_ten








########@profile
def density_double_col( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x, rho_col , D, loacation_y,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)

 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])

   iden.setLabel([5])
   Peps_list=Peps_list*iden
   iden.setLabel([11])
   Peps_list_conj=Peps_list_conj*iden

   PEP_com=((Peps_list*Swap1)*Swap4)*((Peps_list_conj*Swap3)*Swap2)
   
   iden.setLabel([13])
   mps_boundry_leftA=mps_boundry_left[i]*iden
   iden.setLabel([15])
   mps_boundry_rightA=mps_boundry_right[i]*iden

   result=(PEP_com*mps_boundry_rightA)*mps_boundry_leftA

   result.permute([12,10,9,14],4)
   E_list_up[i]=result
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   PEP_com=((Peps_list*Swap1)*Swap4)*((Peps_list_conj*Swap3)*Swap2)


   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   E_list_up[i+1].setLabel([13,5,11,15])
   result=((((PEP_com*E_list_up[i+1])*mps_boundry_left[i])*mps_boundry_right[i]))
   result.permute([12,10,9,14],4)
   E_list_up[i]=result



 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   PEP_com=((Peps_list*Swap1)*Swap4)*((Peps_list_conj*Swap3)*Swap2)

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   iden.setLabel([10])
   PEP_com=PEP_com*iden
   iden.setLabel([9])
   PEP_com=PEP_com*iden

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[i]*iden

   result=((PEP_com*mps_leftA)*mps_rightA)

   result.permute([13,5,11,15],4)
   E_list_down[i]=result

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   PEP_com=((Peps_list*Swap1)*((Peps_list_conj*Swap3)*Swap2))*Swap4

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   E_list_down[i-1].setLabel([12,10,9,14])

   result=((((PEP_com*E_list_down[i-1])*mps_boundry_left[i])*mps_boundry_right[i]))
   result.permute([13,5,11,15],4)
   E_list_down[i]=result

#  for i in xrange(len(PEPS_listten)-1):
#   E_list_down[i].setLabel([1,2,3,4])
#   E_list_up[i+1].setLabel([1,2,3,4])
#   A=E_list_down[i]*E_list_up[i+1]
#   print "double" , A[0]


 for i in xrange(N_x-1):


  rho_col[loacation_y][i]=density_local_col_double( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i,D,Sys)



########@profile
def  density_local_col_double(PEPS_listten,E_list_down,E_list_up,mps_boundry_left,mps_boundry_right, Location,D,Sys):

 bdi_mid=uni10.Bond( uni10.BD_IN, 1)
 iden=uni10.UniTensor( [bdi_mid] )
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],2)


 r_u.setLabel([-100,3,10])


 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


#####################################################

 A=Peps_1*1.0

 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],3)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 A_conj.permute([-6,-7,-9,-3,-10],3)



 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  q_d, V, s=TU.setTruncation(A_conj, row)
 else:
  q_d, V, s=TU.setTruncation(A_conj, colm)


 s.setLabel([1,0])
 V.setLabel([0,-3,-10])
 r_d=V*s
 r_d.permute([1,-3,-10],3)
 
 q_d.setLabel([-6,-7,-9,-200])
 q_d.permute([-6,-7,-9,-200],4)

 r_d.setLabel([-200,-3,-10])



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)

 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])


 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 Swap1=fermionicOPT(Sys,A_conj.bond(2), A_conj.bond(0))
 Swap1.setLabel([11,10,3,8])
 A_conj=A_conj*Swap1

 A_conj.permute([10,9,11,-6,7],3)
 A_conj.setLabel([-11,-10,-13,-14,-15])
 A_conj.permute([-10,-13,-11,-14,-15],2)


 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  U, qq_d, s=TU.setTruncation(A_conj, row)
 else:
  U, qq_d, s=TU.setTruncation(A_conj, colm)

 s.setLabel([0,1])
 U.setLabel([-10,-13,0])
 l_d=U*s
 l_d.permute([-10,-13,1],3)

 qq_d.setLabel([-400,-11,-14,-15])
 qq_d.permute([-400,-11,-14,-15],4)


 l_d.setLabel([-10,-13,-400])


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([ 90, 200, 9, -200 ])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])

######################################################################
 mps_boundry_left[Location].setLabel([16,6,-60,18])
 mps_boundry_left[Location+1].setLabel([18,11,-110,25])


 mps_boundry_right[Location].setLabel([17,90,-9,19])
 mps_boundry_right[Location+1].setLabel([19,140,-14,24])

######################################################

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*q)
 A=A*(Swap3*q_d)
 A=A*mps_boundry_right[Location]

 B=E_right*mps_boundry_left[Location+1]
 B=B*(Swap4*qq)
 B=B*(Swap2*qq_d)
 B=B*mps_boundry_right[Location+1]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 #N_ten=N_Positiv(N_ten)
######################################################


 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))
 Norm_h.permute([-3,-13,3,13],2)

 return  Norm_h









def  Norm_based_on_LR(MPS_R_right, MPS_R_left, D):

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)
 E_a=0
 A_l=[]
 for i in xrange(MPS_R_left.N):

  Ten_L=uni10.UniTensor([MPS_R_left[i].bond(0), bdi, bdi, MPS_R_left[i].bond(2)], "Ten_R")
  Ten_L.putBlock(MPS_R_left[i].getBlock())
  Ten_L.setLabel([-1,-2,0,-3])

  Ten_R=uni10.UniTensor([MPS_R_right[i].bond(0), bdi, bdi, MPS_R_right[i].bond(2)], "Ten_R")
  Ten_R.putBlock(MPS_R_right[i].getBlock())
  Ten_R.setLabel([1,2,0,3])

  A_ten=Ten_R*Ten_L
  A_ten.permute([1,-1,2,-2,3,-3],4)
  A_ten.combineBond([1,-1])
  A_ten.combineBond([2,-2])
  A_ten.combineBond([3,-3])
  A_ten.permute([1,2,3],2)
  
  if i == 0:
    #print A_l[i]
    A_ten.setLabel([-1,-2,1])
    A_l_dag=copy.copy(A_ten)
    A_l_dag.setLabel([-1,-2,2])
    E_a=A_l_dag*A_ten
    E_a.permute([1,2],1)
    E_a.setLabel([-3,-4])
  elif i == (MPS_R_left.N-1):
    A_ten.setLabel([-3,-2,1])
    A_l_dag=copy.copy(A_ten)
    A_l_dag.setLabel([-4,-2,1])
    #print A_l[i].printDiagram()
    E_a=A_l_dag*(E_a*A_ten)
  else:
    A_ten.setLabel([-3,-2,1])
    A_l_dag=copy.copy(A_ten)
    A_l_dag.setLabel([-4,-2,2])
    E_a=A_l_dag*(E_a*A_ten)
    E_a.permute([1,2],1)
    E_a.setLabel([-3,-4])

# list_bond=[]
# for q in xrange(MPS_R_left.N):
#   list_bond.append(A_list[q].bond(2).dim())
# 
# print "hi", max(list_bond)
# 
# mps_R=MPSclass.MPS(A_list[1].bond(1).dim(),max(list_bond),MPS_R_left.N)

# for i in xrange(MPS_R_left.N):
#   mps_R[i]=copy.copy(A_list[i])
    
 return   E_a[0]

def make_ENV_updown(MPO_Ten, N_x, chi):

 chi_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_list.append(dim)

 chi_dim=sum(chi_list)
 #print "chi_dim", chi_dim
 bdi1 = uni10.Bond(uni10.BD_IN, 1)
 Ten1=uni10.UniTensor([bdi1])
 Ten2=uni10.UniTensor([bdi1])
 Ten1.identity()
 Ten2.identity()

###################################################################
 Env_up=[None]*N_x
 cont_list=[None]*N_x
 for j in xrange(N_x):
  MPO_Ten[j][N_x-1].setLabel([1,-1,2,-2,4,-4,5,-5])
  Ten1.setLabel([5])
  Ten2.setLabel([-5])

  A=(MPO_Ten[j][N_x-1]*Ten1)*Ten2
  A.permute([1,-1,2,-2,4,-4],4)
  cont_list[j]=copy.copy(A)

 mps_A=MPSclass.MPS( 2, 2, N_x )

 for i in xrange(N_x):
  mps_A[i]=cont_list[i]*1.0

 Env_up[N_x-1]=mps_A.appSVD(chi_dim)


 for q in xrange(N_x-1):
  for j in xrange(N_x):
   A=Env_up[N_x-1-q][j]*1.0
   A.setLabel([10,5,-5,20])
   B=MPO_Ten[j][N_x-2-q]*1.0
   B.setLabel([1,-1,2,-2,4,-4,5,-5])
   Result=A*B
   Result.permute([10,1,-1,2,-2,20,4,-4],5)
   cont_list[j]=Result*1.0


  mps_A=MPSclass.MPS( 2, 2, N_x)
  for i in xrange(N_x):
   mps_A[i]=cont_list[i]*1.0


  Env_up[N_x-2-q]=mps_A.appSVD(sum(chi_list))
#############################################################


###################################################################
 Env_down=[None]*N_x
 cont_list=[None]*N_x
 for j in xrange(N_x):
  MPO_Ten[j][0].setLabel([1,-1,2,-2,4,-4,5,-5])
  Ten1.setLabel([2])
  Ten2.setLabel([-2])

  A=(MPO_Ten[j][0]*Ten1)*Ten2
  A.permute([1,-1,5,-5,4,-4],4)
  cont_list[j]=A*1.0

 mps_A=MPSclass.MPS( 2, 2, N_x )

 for i in xrange(N_x):
  mps_A[i]=cont_list[i]*1.0

 Env_down[0]=mps_A.appSVD(chi_dim)


 for q in xrange(N_x-1):
  for j in xrange(N_x):
   A=Env_down[q][j]*1.0
   A.setLabel([10,2,-2,20])
   B=MPO_Ten[j][q+1]*1.0
   B.setLabel([1,-1,2,-2,4,-4,5,-5])
   Result=A*B
   Result.permute([10,1,-1,5,-5,20,4,-4],5)
   cont_list[j]=Result*1.0


  mps_A=MPSclass.MPS( 2, 2, N_x)
  for i in xrange(N_x):
   mps_A[i]=cont_list[i]*1.0


  Env_down[q+1]=mps_A.appSVD(sum(chi_list))


 return Env_up, Env_down



def make_ENVMPS(MPO_Ten, N_x, chi):

 chi_list=[]
 bdi = uni10.Bond(uni10.BD_IN, chi)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  chi_list.append(dim)

 chi_dim=sum(chi_list)
 #print "chi_dim", chi_dim

 bdi1 = uni10.Bond(uni10.BD_IN, 1)
 Ten1=uni10.UniTensor([bdi1])
 Ten2=uni10.UniTensor([bdi1])
 Ten1.identity()
 Ten2.identity()

 Env_left=[None]*N_x
 cont_list=[None]*N_x
 for j in xrange(N_x):
  MPO_Ten[0][j].setLabel([1,-1,2,-2,4,-4,5,-5])
  Ten1.setLabel([1])
  Ten2.setLabel([-1])

  A=(MPO_Ten[0][j]*Ten1)*Ten2
  A.permute([2,-2,4,-4,5,-5],4)
  cont_list[j]=copy.copy(A)

 mps_A=MPSclass.MPS( 2, 2, N_x)

 for i in xrange(N_x):
  mps_A[i]=copy.copy(cont_list[i])

 Env_left[0]=mps_A.appSVD(chi_dim)

 for q in xrange(N_x-1):
  for j in xrange(N_x):
   A=copy.copy(Env_left[q][j])
   A.setLabel([10,1,-1,20])
   B=MPO_Ten[q+1][j]*1.0
   B.setLabel([1,-1,2,-2,3,-3,4,-4])
   Result=A*B
   Result.permute([10,2,-2,3,-3,20,4,-4],5)
   cont_list[j]=Result*1.0

  mps_A=MPSclass.MPS( 2, 2, N_x)
  for i in xrange(N_x):
   mps_A[i]=copy.copy(cont_list[i])

  Env_left[q+1]=mps_A.appSVD(sum(chi_list))


##################################################

###############################

 Env_right=[None]*N_x
 cont_list=[None]*N_x
 for j in xrange(N_x):
   MPO_Ten[N_x-1][j].setLabel([1,-1,2,-2,3,-3,4,-4])
   Ten1.setLabel([3])
   Ten2.setLabel([-3])

   A=(MPO_Ten[N_x-1][j]*Ten1)*Ten2
   A.permute([2,-2,1,-1,4,-4],4)
   cont_list[j]=copy.copy(A)


 mps_A=MPSclass.MPS( 2, 2, N_x)

 for i in xrange(N_x):
  mps_A[i]=cont_list[i]*1.0

 #print mps_A.norm()
 Env_right[N_x-1]=mps_A.appSVD(sum(chi_list))
#  print N_x-1,  Env_right[N_x-1].norm()#, ( Env_right[N_x-1].norm()-mps_A.norm()) / mps_A.norm()


 for q in xrange(N_x-1):
  for j in xrange(N_x):
   A=copy.copy(Env_right[N_x-1-q][j])
   A.setLabel([10,3,-3,20])
   B=copy.copy(MPO_Ten[N_x-2-q][j])
   B.setLabel([1,-1,2,-2,3,-3,4,-4])
   Result=A*B
   Result.permute([10,2,-2,1,-1,20,4,-4],5)
   cont_list[j]=copy.copy(Result)


  mps_A=MPSclass.MPS( 2, 2, N_x)
  for i in xrange(N_x):
   mps_A[i]=copy.copy(cont_list[i])

  Env_right[N_x-2-q]=mps_A.appSVD(sum(chi_list))
#   print N_x-2-q, Env_right[N_x-2-q].norm()#, mps_A.norm(), ( Env_right[N_x-2-q].norm()-mps_A.norm() ) / mps_A.norm()


 return Env_left, Env_right




def make_boundry_MPO( PEPS_listten, N_x,Sys):

 MPO=[None]*N_x
 for i in xrange(N_x):
   MPO[i]=[None]*N_x

 for i in xrange(N_x):
  for j in xrange(N_x):
   A=copy.copy(PEPS_listten[i][j])
   A.setLabel([1,2,3,4,5])
   A.permute([1,2,3,4,5],3)
   A_t=A*1.0
   A_t.transpose()
   A_t.setLabel([-4,-5,-1,-2,3])
   A_t.permute([-1,-2,3,-4,-5],5)
   A.permute([1,2,3,4,5],5)

   #print A.printDiagram(), Swap1.printDiagram(), Swap2.printDiagram() , A_t.printDiagram(), Swap4.printDiagram(), Swap3.printDiagram()

   Swap2=fermionicOPT(Sys,A_t.bond(3), A_t.bond(4))
   Swap2.setLabel([-6,7,-4,-5])

   Swap3=fermionicOPT(Sys,A_t.bond(0), A_t.bond(1))
   Swap3.setLabel([8,9,-1,-2])

   Swap4=fermionicOPT(Sys,A_t.bond(0), A.bond(1))
   Swap4.setLabel([-8,10,8,2])

   Swap1=fermionicOPT(Sys,A.bond(3), A_t.bond(4))
   Swap1.setLabel([6,11,4,7])

   
   A=(A*Swap1)*(((A_t*Swap2)*Swap3)*Swap4)
   A.permute([1,-8,10,9,6,-6,5,11],4)
   #A.combineBond([1,-1])
   #A.combineBond([2,-2])
   #A.combineBond([4,-4])
   #A.combineBond([5,-5])
   #A.permute([1,2,4,5],2)
   MPO[i][j]=A*1.0

 return   MPO



















def Init_PEPS( N_x, D, d, q):
 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdiphy=uni10.Bond(uni10.BD_IN, d)
 bdophy=uni10.Bond(uni10.BD_OUT, d)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)
 A_list=[None]*N_x
 B_list=[None]*N_x

 if q == 0:
  for i in xrange(N_x):
   if i == 0:
     A=uni10.UniTensor([bdi1,bdi1,bdiphy,bdo,bdo], "A_first")
     A.randomize()
     A.setLabel([0,1,2,3,4])
     A.permute([0,1,2,3,4],3)
     A_list[i]=A
   elif i ==(N_x-1):
     A=uni10.UniTensor([bdi1,bdi,bdiphy,bdo,bdo1], "A_Last")
     A.randomize()
     A.setLabel([0,1,2,3,4])
     A.permute([0,1,2,3,4],3)
     A_list[i]=A
   else:
     A=uni10.UniTensor([bdi1,bdi,bdiphy,bdo,bdo], "A_middle")
     A.randomize()
     A.setLabel([0,1,2,3,4])
     A.permute([0,1,2,3,4],3)
     A_list[i]=A
 elif q==N_x-1:
  A_list=[None]*N_x
  for i in xrange(N_x):
   if i == 0:
     A=uni10.UniTensor([bdi,bdi1,bdiphy,bdo1,bdo], "A_first")
     A.randomize()
     A.setLabel([0,1,2,3,4])
     A.permute([0,1,2,3,4],3)
     A_list[i]=A
   elif i ==(N_x-1):
     A=uni10.UniTensor([bdi,bdi,bdiphy,bdo1,bdo1], "A_Last")
     A.randomize()
     A.setLabel([0,1,2,3,4])
     A.permute([0,1,2,3,4],3)
     A_list[i]=A
   else:
     A=uni10.UniTensor([bdi,bdi,bdiphy,bdo1,bdo], "A_middle")
     A.randomize()
     A.setLabel([0,1,2,3,4])
     A.permute([0,1,2,3,4],3)
     A_list[i]=A
 else:
  A_list=[None]*N_x
  for i in xrange(N_x):
   if i == 0:
     A=uni10.UniTensor([bdi,bdi1,bdiphy,bdo,bdo], "A_first")
     A.randomize()
     A.setLabel([0,1,2,3,4])
     A.permute([0,1,2,3,4],3)
     A_list[i]=A
   elif i ==(N_x-1):
     A=uni10.UniTensor([bdi,bdi,bdiphy,bdo,bdo1], "A_Last")
     A.randomize()
     A.setLabel([0,1,2,3,4])
     A.permute([0,1,2,3,4],3)
     A_list[i]=A
   else:
     A=uni10.UniTensor([bdi,bdi,bdiphy,bdo,bdo], "A_middle")
     A.randomize()
     A.setLabel([0,1,2,3,4])
     A.permute([0,1,2,3,4],3)
     A_list[i]=A

 return A_list





def inverse(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 blk_qnums=Landa2.blockQnum()
 for qnum in blk_qnums:
  D=int(Landa2.getBlock(qnum).row())
  D1=int(Landa2.getBlock(qnum).col())
  invL2 = uni10.Matrix(D, D1,True)
  invLt = uni10.Matrix(D, D1,True)
  invLt=Landa2.getBlock(qnum,True)
  #print invLt[0], invLt[1], invLt[2], invLt[3]
  for i in xrange(D):
      invL2[i] = 0 if ((invLt[i].real) < 1.0e-10) else (1.00 / (invLt[i].real))

  invLanda2.putBlock(qnum,invL2)
 return invLanda2


#########  prerequisite functions  #############
def   setTruncation1(theta, chi):
    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge1(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0),theta.bond(1),theta.bond(2),bdo_mid])
    GB.assign([bdi_mid,theta.bond(3),theta.bond(4),theta.bond(5)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim)  )
def   setTruncationMPS(theta, chi):
    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_mergemps(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0), theta.bond(1),bdo_mid])
    GB.assign([bdi_mid,  theta.bond(2)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim)  )
    return GA, GB, LA

def   sv_mergemps(svs, bidxs, bidx, sv_mat, chi, len_qn):
    if(len(svs)):
        length = len(svs) + sv_mat.elemNum()
        length = length if length < chi else chi
        ori_svs = svs
        ori_bidxs = bidxs
        svs = [0] * length
        bidxs = [0] * length
        svs = []
        bidxs = []
        cnt  = 0
        cur1 = 0
        cur2 = 0
        while cnt < length:
            if(cur1 < len(ori_svs)) and cur2 < sv_mat.elemNum():
                if ori_svs[cur1] >= sv_mat[cur2]:
                    if (ori_svs[cur1] > -0.01):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                    cur1 += 1
                else:
                    if (sv_mat[cur2] > -0.01):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                    cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > -0.01):
                     svs.append(sv_mat[cur2]) 
                     bidxs.append(bidx) 
                    cur2 += 1
                break
            else:
                for i in xrange(cur1, len(ori_svs)):
                 svs.append(ori_svs[i])
                 bidxs.append(ori_bidxs[i]) 
                break
            cnt += 1
    else:
       if (len_qn is 1):
        bidxs = [bidx] * chi  
        svs = [sv_mat[i] for i in xrange(chi)]
       elif (sv_mat[0] > -0.01):
        bidxs = [bidx] * sv_mat.elemNum()
        svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
       else: bidxs = [bidx];  svs = [sv_mat[0]];  
    return svs, bidxs











def   sv_merge1(svs, bidxs, bidx, sv_mat, chi, len_qn):
    if(len(svs)):
        length = len(svs) + sv_mat.elemNum()
        length = length if length < chi else chi
        ori_svs = svs
        ori_bidxs = bidxs
        svs = [0] * length
        bidxs = [0] * length
        svs = []
        bidxs = []
        cnt  = 0
        cur1 = 0
        cur2 = 0
        while cnt < length:
            if(cur1 < len(ori_svs)) and cur2 < sv_mat.elemNum():
                if ori_svs[cur1] >= sv_mat[cur2]:
                    if (ori_svs[cur1] > -0.01):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                    cur1 += 1
                else:
                    if (sv_mat[cur2] > -0.01):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                    cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > -0.01):
                     svs.append(sv_mat[cur2]) 
                     bidxs.append(bidx) 
                    cur2 += 1
                break
            else:
                for i in xrange(cur1, len(ori_svs)):
                 svs.append(ori_svs[i])
                 bidxs.append(ori_bidxs[i]) 
                break
            cnt += 1
    else:
       if (len_qn is 1):
        bidxs = [bidx] * chi  
        svs = [sv_mat[i] for i in xrange(chi)]
       elif (sv_mat[0] > -0.01):
        bidxs = [bidx] * sv_mat.elemNum()
        svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
       else: bidxs = [bidx];  svs = [sv_mat[0]];  
    return svs, bidxs



    
    
    
  
  
 
 

def svd_parityrl(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
    bd2=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
    bd3=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
    LA=uni10.UniTensor([bd1,bd2,bd3,theta.bond(3),theta.bond(4),theta.bond(5)])
    GB=uni10.UniTensor([bd1,bd2,bd3,theta.bond(3),theta.bond(4),theta.bond(5)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB

def svd_parity(theta):

    bd1=uni10.Bond(uni10.BD_IN,theta.bond(1).Qlist())

    GA=uni10.UniTensor([theta.bond(0),theta.bond(1)])
    LA=uni10.UniTensor([bd1,theta.bond(1)])
    GB=uni10.UniTensor([bd1,theta.bond(1)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB


def svd_parity1(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim())
    #print bdi1
    GA=uni10.UniTensor([theta.bond(0),theta.bond(1), bdo1])
    LA=uni10.UniTensor([bdi1,bdo1])
    GB=uni10.UniTensor([bdi1,theta.bond(2),theta.bond(3)])

    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        GA.putBlock(qnum, svds[qnum][0])
        LA.putBlock(qnum, svds[qnum][1])
        GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB


def svd_parity2(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(3).dim()*theta.bond(4).dim()*theta.bond(5).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(3).dim()*theta.bond(4).dim()*theta.bond(5).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())

    if bdi.dim()<=bdi1.dim():
     #print theta.printDiagram()
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(3),theta.bond(4), theta.bond(5)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
    else:
     #print theta.printDiagram()
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(3),theta.bond(4), theta.bond(5)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB




def svd_parity5(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(3).dim()*theta.bond(4).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(3).dim()*theta.bond(4).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim())

    if bdi.dim()<=bdi1.dim():
     #print theta.printDiagram()
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(3),theta.bond(4)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
    else:
     #print theta.printDiagram()
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(3),theta.bond(4)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB


def svd_parity6(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim())

    if bdi.dim()<=bdi1.dim():
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(2),theta.bond(3),theta.bond(4)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
    else:
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(2),theta.bond(3),theta.bond(4)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB







def setTruncation2(theta, chi):
    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge1(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3), bdo_mid])
    GB.assign([bdi_mid, theta.bond(4), theta.bond(5), theta.bond(6), theta.bond(7)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim)  )
    return GA, GB, LA


def setTruncation3(theta, chi):
    LA=uni10.UniTensor(theta.bond())
    GA=uni10.UniTensor(theta.bond())
    GB=uni10.UniTensor(theta.bond())
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
        dim_svd.append(int(svds[qnum][1].col()))
    svs = []
    bidxs = []
    for bidx in xrange(len(blk_qnums)):
        svs, bidxs = sv_merge1(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi,len(blk_qnums))
    dims = [0] * len(blk_qnums)
    for bidx in bidxs:
        dims[bidx] += 1  
    qnums = []
    for bidx in xrange(len(blk_qnums)):
        qnums += [blk_qnums[bidx]] * dims[bidx]
    bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
    bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
    GA.assign([theta.bond(0),theta.bond(1), bdo_mid])
    GB.assign([bdi_mid, theta.bond(2), theta.bond(3)])
    LA.assign([bdi_mid, bdo_mid])
    degs = bdi_mid.degeneracy()
    for qnum, dim in degs.iteritems():
        if qnum not in svds:
            raise Exception("In setTruncaton(): Fatal error.")
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
        GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
        LA.putBlock(qnum, svd[1].resize(dim, dim)  )
    return GA, GB, LA

def sv_merge1(svs, bidxs, bidx, sv_mat, chi, len_qn):
    if(len(svs)):
        length = len(svs) + sv_mat.elemNum()
        length = length if length < chi else chi
        ori_svs = svs
        ori_bidxs = bidxs
        svs = [0] * length
        bidxs = [0] * length
        svs = []
        bidxs = []
        cnt  = 0
        cur1 = 0
        cur2 = 0
        while cnt < length:
            if(cur1 < len(ori_svs)) and cur2 < sv_mat.elemNum():
                if ori_svs[cur1] >= sv_mat[cur2]:
                    if (ori_svs[cur1] > -0.01):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                    cur1 += 1
                else:
                    if (sv_mat[cur2] > -0.01):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                    cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > -0.01):
                     svs.append(sv_mat[cur2]) 
                     bidxs.append(bidx) 
                    cur2 += 1
                break
            else:
                for i in xrange(cur1, len(ori_svs)):
                 svs.append(ori_svs[i])
                 bidxs.append(ori_bidxs[i]) 
                break
            cnt += 1
    else:
       if (len_qn is 1):
        bidxs = [bidx] * chi  
        svs = [sv_mat[i] for i in xrange(chi)]
       elif (sv_mat[0] > -0.01):
        bidxs = [bidx] * sv_mat.elemNum()
        svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
       else: bidxs = [bidx];  svs = [sv_mat[0]];  
    return svs, bidxs


def svd_parity3(theta):

    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(4).dim()*theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(4).dim()*theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim())

    if bdi.dim()<=bdi1.dim():
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(4),theta.bond(5), theta.bond(6),theta.bond(7)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
    else:
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(4),theta.bond(5), theta.bond(6),theta.bond(7)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB

def svd_parity4(theta):
    bdi1=uni10.Bond(uni10.BD_IN,theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim()*theta.bond(8).dim()*theta.bond(9).dim())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(5).dim()*theta.bond(6).dim()*theta.bond(7).dim()*theta.bond(8).dim()*theta.bond(9).dim())

    bdi=uni10.Bond(uni10.BD_IN,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(0).dim()*theta.bond(1).dim()*theta.bond(2).dim()*theta.bond(3).dim()*theta.bond(4).dim())

    if bdi.dim()<=bdi1.dim():
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3),theta.bond(4), bdo])
     LA=uni10.UniTensor([bdi,bdo])
     GB=uni10.UniTensor([bdi,theta.bond(5),theta.bond(6), theta.bond(7),theta.bond(8),theta.bond(9)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]

     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])
    else:
     #print theta.printDiagram()
     GA=uni10.UniTensor([theta.bond(0),theta.bond(1), theta.bond(2),theta.bond(3),theta.bond(4), bdo1])
     LA=uni10.UniTensor([bdi1,bdo1])
     GB=uni10.UniTensor([bdi1,theta.bond(5),theta.bond(6), theta.bond(7),theta.bond(8),theta.bond(9)])

     svds = {}
     blk_qnums = theta.blockQnum()
     dim_svd=[]
     for qnum in blk_qnums:
         svds[qnum] = theta.getBlock(qnum).svd()
         #print   svds[qnum][0].row(), svds[qnum][0].col()
         GA.putBlock(qnum, svds[qnum][0])
         LA.putBlock(qnum, svds[qnum][1])
         GB.putBlock(qnum, svds[qnum][2])

    return GA, LA, GB


# def svd_parity3(theta):
# 
#  bd1=uni10.Bond(uni10.BD_IN,theta.bond(3).Qlist())
#  bd2=uni10.Bond(uni10.BD_IN,theta.bond(4).Qlist())
#  bd3=uni10.Bond(uni10.BD_IN,theta.bond(5).Qlist())
#  #bd4=uni10.Bond(uni10.BD_IN,theta.bond(7).Qlist())
# 
#  GA=uni10.UniTensor([theta.bond(0),theta.bond(1),theta.bond(2),theta.bond(3),theta.bond(4),theta.bond(5)])
#  LA=uni10.UniTensor([bd1,bd2,bd3,theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7)])
#  GB=uni10.UniTensor([bd1,bd2,bd3,bd4,theta.bond(4),theta.bond(5),theta.bond(6),theta.bond(7)])
# 
#  svds = {}
#  blk_qnums = theta.blockQnum()
#  dim_svd=[]
#  for qnum in blk_qnums:
#      svds[qnum] = theta.getBlock(qnum).svd()
#      GA.putBlock(qnum, svds[qnum][0])
#      LA.putBlock(qnum, svds[qnum][1])
#      GB.putBlock(qnum, svds[qnum][2])
# 
# #    print LA
#  return GA, LA, GB

def Sqrt(Landa):
  Landa_cp=copy.copy(Landa)
  blk_qnums=Landa.blockQnum()
  for qnum in blk_qnums:
   D=int(Landa_cp.getBlock(qnum).col())
   Landa_cpm=Landa_cp.getBlock(qnum)
   Landam=Landa_cp.getBlock(qnum)
   for i in xrange(D):
    for j in xrange(D):
     if Landam[i*D+j] > 1.0e-12:
      Landa_cpm[i*D+j]=Landam[i*D+j]**(1.00/2.00)
     else:
      Landa_cpm[i*D+j]=0
   Landa_cp.putBlock(qnum,Landa_cpm)
  return Landa_cp 


def Sqrt_minor(Landa):
  Landa_cp=copy.copy(Landa)
  blk_qnums=Landa.blockQnum()
  for qnum in blk_qnums:
   D=int(Landa_cp.getBlock(qnum).col())
   Landa_cpm=Landa_cp.getBlock(qnum)
   Landam=Landa_cp.getBlock(qnum)
   for i in xrange(D):
    for j in xrange(D):
     if Landam[i*D+j] > 1.0e-10:
      Landa_cpm[i*D+j]=Landam[i*D+j]**(1.00/2.00)
     else:
      Landa_cpm[i*D+j]=0
   Landa_cp.putBlock(qnum,Landa_cpm)
  return Landa_cp 
























def inv_landa_col_row( Landa_col, Landa_row, N_x): 

 Landa_col_inv=[None]*N_x
 for i in xrange(N_x):
  Landa_col_inv[i]=[None]*(N_x+1)

 Landa_row_inv=[None]*(N_x+1)
 for i in xrange(N_x+1):
  Landa_row_inv[i]=[None]*(N_x)

 for i in xrange(N_x):
  for j in xrange(N_x+1):
    Landa_col_inv[i][j]=inverse(Landa_col[i][j])

 for i in xrange(N_x+1):
  for j in xrange(N_x):
    Landa_row_inv[i][j]=inverse(Landa_row[i][j])

 return   Landa_col_inv,   Landa_row_inv

def Landa_f_col_iden(Landa_col,N_x):

 Landa_col1=[None]*N_x
 for i in xrange(N_x):
  Landa_col1[i]=[None]*(N_x+1)

 for i in xrange(N_x):
  for j in xrange(N_x+1):
   Landa_col[i][j].identity()
   Landa_col1[i][j]=Landa_col[i][j]

 return   Landa_col1


def Landa_f_row_iden(Landa_row,N_x):



 Landa_row1=[None]*(N_x+1)
 for i in xrange(N_x+1):
  Landa_row1[i]=[None]*(N_x)
 
 for i in xrange(N_x+1):
  for j in xrange(N_x):
    Landa_row[i][j].identity()
    Landa_row1[i][j]=Landa_row[i][j]

 return Landa_row1


def Landa_f_row_rebonding(PEPS_listtenRG, Landa_row, N_x):

 for i in xrange(N_x+1):
  for j in xrange(N_x):
   if i<N_x:
    bdi=uni10.Bond(uni10.BD_IN, PEPS_listtenRG[i][j].bond(0).Qlist())
    bdo=uni10.Bond(uni10.BD_OUT, PEPS_listtenRG[i][j].bond(0).Qlist())
    Landa_row[i][j]=uni10.UniTensor([bdi, bdo])
    Landa_row[i][j].identity()
   if i==N_x:
    bdi=uni10.Bond(uni10.BD_IN, PEPS_listtenRG[i-1][j].bond(3).Qlist())
    bdo=uni10.Bond(uni10.BD_OUT, PEPS_listtenRG[i-1][j].bond(3).Qlist())
    Landa_row[i][j]=uni10.UniTensor([bdi, bdo])
    Landa_row[i][j].identity()


 #print Landa_row[0][0]


def Landa_f_col_rebonding(PEPS_listtenRG, Landa_col, N_x):

 for i in xrange(N_x):
  for j in xrange(N_x+1):
   if j<N_x:
    bdi=uni10.Bond(uni10.BD_IN, PEPS_listtenRG[i][j].bond(1).Qlist())
    bdo=uni10.Bond(uni10.BD_OUT, PEPS_listtenRG[i][j].bond(1).Qlist())
    Landa_col[i][j]=uni10.UniTensor([bdi, bdo])
    Landa_col[i][j].identity()
   if j==N_x:
    bdi=uni10.Bond(uni10.BD_IN, PEPS_listtenRG[i][j-1].bond(4).Qlist())
    bdo=uni10.Bond(uni10.BD_OUT, PEPS_listtenRG[i][j-1].bond(4).Qlist())
    Landa_col[i][j]=uni10.UniTensor([bdi, bdo])
    Landa_col[i][j].identity()











def Landa_f_col(D, N_x):

 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q0_odd = uni10.Qnum(0,uni10.PRT_ODD);

 bdi = uni10.Bond(uni10.BD_IN, D)
 bdo = uni10.Bond(uni10.BD_OUT, D)
 Landa=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa.identity()
 Landa.randomize()
 #Landa=Landa*.0001
 #Landa.orthoRand()
 Landa=sparce_init(Landa)

 bdi = uni10.Bond(uni10.BD_IN, 1)
 bdo = uni10.Bond(uni10.BD_OUT, 1)
 Landa_iden=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa_iden.identity()


 Landa_col=[None]*N_x
 for i in xrange(N_x):
  Landa_col[i]=[None]*(N_x+1)

 for i in xrange(N_x):
  for j in xrange(N_x+1):
   if j==0 or j==N_x:
    Landa_col[i][j]=Landa_iden*1.0
   else:
    Landa_col[i][j]=Landa*1.0

 return   Landa_col


def Landa_f_row(D, N_x):

 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q0_odd = uni10.Qnum(0,uni10.PRT_ODD);

 bdi = uni10.Bond(uni10.BD_IN, D)
 bdo = uni10.Bond(uni10.BD_OUT, D)
 Landa=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa.randomize()
 #Landa.identity()
 #Landa.orthoRand()
 Landa=sparce_init(Landa)

 bdi = uni10.Bond(uni10.BD_IN, 1)
 bdo = uni10.Bond(uni10.BD_OUT, 1)
 Landa_iden=uni10.UniTensor([bdi,bdo],"Landa_1")
 Landa_iden.identity()


 Landa_row=[None]*(N_x+1)
 for i in xrange(N_x+1):
  Landa_row[i]=[None]*(N_x)
 
 for i in xrange(N_x+1):
  for j in xrange(N_x):
   if i==0 or i==N_x:
    Landa_row[i][j]=Landa_iden*1.0
   else:
    Landa_row[i][j]=Landa*1.0

 return Landa_row

def sparce_init(Landa):
 q0_even = uni10.Qnum(0,uni10.PRT_EVEN);
 q0_odd = uni10.Qnum(0,uni10.PRT_ODD);

 blk_qnums = Landa.blockQnum()
 for qnum in blk_qnums:
  M=Landa.getBlock(qnum)
  if qnum == q0_even:
   M.randomize()
   M=M*0.1
   M[0]=1.00
  else: 
   M.randomize()
   M=M*0.1
   #M[0]=.200
  Landa.putBlock(qnum,M)

 return Landa

def norm_Symmetry(LA):
 norm=0
 blk_qnums = LA.blockQnum()
 for qnum in blk_qnums:
  M=LA.getBlock(qnum)  
  norm=norm+(M.norm()*M.norm())
 norm=norm**(1.00/2.00)
 return norm

def inverse(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 blk_qnums=Landa2.blockQnum()
 for qnum in blk_qnums:
  D=int(Landa2.getBlock(qnum).row())
  D1=int(Landa2.getBlock(qnum).col())
  invL2 = uni10.Matrix(D, D1,True)
  invLt = uni10.Matrix(D, D1,True)
  invLt=Landa2.getBlock(qnum,True)

  for i in xrange(D):
   invL2[i] = 0 if ((invLt[i].real) < 1.0e-12) else (1.00 / (invLt[i].real))
  invLanda2.putBlock(qnum,invL2)
 return invLanda2

#############@profile
def update_row( PEPS_A, PEPS_B, U_ham, i_x, j_y, Landa_col, Landa_row, D,Sys, N_iterF):

 bdi = uni10.Bond(uni10.BD_IN, D)
 D_dim=[]
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_dim.append(dim)

 PEPS_A.setLabel([1,2,3,4,5])
 PEPS_B.setLabel([6,7,8,9,10])

 B=PEPS_B*1.0
 B.setLabel([6,-7,8,-9,-10])
 B.permute([6,-7,8,-9,-10],5)
 Swap=fermionicOPT(Sys,B.bond(1), B.bond(2))
 Swap.setLabel([7,-8,-7,8])
 B=Swap*B
 B.permute([6,7,-8,-9,-10],3)
 B.setLabel([6,7,8,9,10])

 U_ham.setLabel([11,12,3,8])

 Landa_row[i_x][j_y].setLabel([-1,1])
 Landa_row[i_x+1][j_y].setLabel([4,6])
 Landa_row[i_x+2][j_y].setLabel([9,-9])
 Landa_col[i_x][j_y].setLabel([-2,2])
 Landa_col[i_x][j_y+1].setLabel([5,-5])
 Landa_col[i_x+1][j_y].setLabel([-7,7])
 Landa_col[i_x+1][j_y+1].setLabel([10,-10])

 PEPS_AA=((((PEPS_A*Landa_row[i_x][j_y])*Landa_row[i_x+1][j_y])*Landa_col[i_x][j_y])*Landa_col[i_x][j_y+1])
 PEPS_BB=(((B*Landa_row[i_x+2][j_y])*Landa_col[i_x+1][j_y])*Landa_col[i_x+1][j_y+1])

 PEPS_AA.permute([-1,-2,3,6,-5],3)
 PEPS_BB.permute([6,-7,8,-9,-10],3)


 if N_iterF[1]=="full":
  Theta=(PEPS_AA*U_ham)*PEPS_BB
  Theta.permute([-1,-2,11,-5,-7,12,-9,-10],4)

  U, V, LA=TU.setTruncation(Theta, sum(D_dim))
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
 ###########################################################################

  blk_qnums = LA.blockQnum()
  Landa_row[i_x+1][j_y].assign(LA.bond()) 
  for qnum in blk_qnums:
   Landa_row[i_x+1][j_y].putBlock(qnum,LA.getBlock(qnum))

  U.setLabel([-1,-2,11,-5,6])
  U.permute([-1,-2,11,6,-5],3)

  V.setLabel([6,-7,12,-9,-10])
  V.permute([6,12,-7,-9,-10],3)

  invLanda1=inverse_ten(Landa_row[i_x][j_y])
  invLanda2=inverse_ten(Landa_col[i_x][j_y])
  invLanda3=inverse_ten(Landa_col[i_x][j_y+1])

  invLanda1.setLabel([1,-1])
  invLanda2.setLabel([2,-2])
  invLanda3.setLabel([-5,5])

  U=((U*invLanda1)*invLanda2)*invLanda3
  U.permute([1,2,11,6,5],3)

  invLanda8=inverse_ten(Landa_row[i_x+2][j_y])
  invLanda9=inverse_ten(Landa_col[i_x+1][j_y])
  invLanda10=inverse_ten(Landa_col[i_x+1][j_y+1])

  invLanda8.setLabel([-9,9])
  invLanda9.setLabel([7,-7])
  invLanda10.setLabel([-10,10])

  V=((V*invLanda8)*invLanda9)*invLanda10
  V.permute([6,7,12,9,10],3)

  B=V*1.0
  B.setLabel([6,-7,8,-9,-10])
  B.permute([6,-7,8,-9,-10],5)
  Swap=fermionicOPT( Sys, B.bond(1), B.bond(2))
  Swap.setLabel([7,-8,-7,8])
  B=Swap*B
  B.permute([6,7,-8,-9,-10],3)
  B.setLabel([6,7,8,9,10])
  
 elif N_iterF[1]=="QR":

  A=copy.copy(PEPS_AA)
  A.setLabel([-1,-2,3,6,-5])
  A.permute([-1,-2,-5,3,6],3)


  row, colm=cal_rowcol(A)
  if (row<=colm):
   q, V, s=TU.setTruncation(A, row)
  else:
   q, V, s=TU.setTruncation(A, colm)


  s.setLabel([1,0])
  V.setLabel([0,3,6])
  r_u=V*s
  r_u.permute([1,3,6],1)
  q.setLabel([-1,-2,-5,-100])
  r_u.setLabel([-100,3,6])


  A=copy.copy(PEPS_BB)
  A.permute([6,8,-7,-9,-10],2)

  row, colm=cal_rowcol(A)
  if (row<=colm):
   U, qq, s=TU.setTruncation(A, row)
  else:
   U, qq, s=TU.setTruncation(A, colm)

  s.setLabel([0,1])
  U.setLabel([6,8,0])
  l_u=U*s
  l_u.permute([6,8,1],2)

  qq.setLabel([-400,-7,-9,-10])
  l_u.setLabel([6,8,-400])

  Theta=(l_u*U_ham)*r_u
  Theta.permute([-100,11,-400,12],2)

  U, V, LA=TU.setTruncation(Theta, sum(D_dim))
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
 ##########################################################################################

  blk_qnums = LA.blockQnum()
  Landa_row[i_x+1][j_y].assign(LA.bond()) 
  for qnum in blk_qnums:
   Landa_row[i_x+1][j_y].putBlock(qnum,LA.getBlock(qnum))

  U.setLabel([-100,11,6])
  U=U*q
  U.permute([-1,-2,11,6,-5],3)

  V.setLabel([6,-400,12])
  V=V*qq
  #V.permute([6,-7,12,-9,-10],3)

  V.permute([6,12,-7,-9,-10],3)

 # Swap=fermionicOPT(Sys,V.bond(1), V.bond(2))
 # Swap.setLabel([12,-7,-12,7])
 # #print Swap.printDiagram(), A.printDiagram(), Swap
 # V=Swap*V
 # V.permute([6,7,-12,-9,-10],3)
 # V.setLabel([6,-7,12,-9,-10])


  invLanda1=inverse_ten(Landa_row[i_x][j_y])
  invLanda2=inverse_ten(Landa_col[i_x][j_y])
  invLanda3=inverse_ten(Landa_col[i_x][j_y+1])

  invLanda1.setLabel([1,-1])
  invLanda2.setLabel([2,-2])
  invLanda3.setLabel([-5,5])

  U=((U*invLanda1)*invLanda2)*invLanda3
  U.permute([1,2,11,6,5],3)

  invLanda8=inverse_ten(Landa_row[i_x+2][j_y])
  invLanda9=inverse_ten(Landa_col[i_x+1][j_y])
  invLanda10=inverse_ten(Landa_col[i_x+1][j_y+1])

  invLanda8.setLabel([-9,9])
  invLanda9.setLabel([7,-7])
  invLanda10.setLabel([-10,10])

  V=((V*invLanda8)*invLanda9)*invLanda10
  V.permute([6,7,12,9,10],3)

  B=V*1.0
  B.setLabel([6,-7,8,-9,-10])
  B.permute([6,-7,8,-9,-10],5)
  Swap=fermionicOPT(Sys,B.bond(1), B.bond(2))
  Swap.setLabel([7,-8,-7,8])
  B=Swap*B
  B.permute([6,7,-8,-9,-10],3)
  B.setLabel([6,7,8,9,10])


 return U, B



#############@profile
def update_col(PEPS_A, PEPS_B, U_ham, i_x, j_y, Landa_col, Landa_row, D, Sys, N_iterF):
 bdi = uni10.Bond(uni10.BD_IN, D)
 D_dim=[]
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_dim.append(dim)

 PEPS_A.setLabel([1,2,3,4,5])
 A=PEPS_A*1.0
 A.setLabel([ 1, 2, 3, 4, 5])
 A.permute([ 1, 2, 3, 4, 5], 6)
 Swap=fermionicOPT(Sys,A.bond(2), A.bond(3))
 Swap.setLabel([-3,-4,3,4])
 #print A.printDiagram(), Swap.printDiagram()
 A=Swap*A
 A.permute([ 1, 2, -3, -4, 5], 3)
 A.setLabel([ 1, 2, 3, 4, 5])


 PEPS_B.setLabel([6,7,8,9,10])
 U_ham.setLabel([11,12,3,8])

 Landa_row[i_x][j_y].setLabel([-1,1])
 Landa_row[i_x+1][j_y].setLabel([4,-4])
 Landa_row[i_x][j_y+1].setLabel([-6,6])
 Landa_row[i_x+1][j_y+1].setLabel([9,-9])

 Landa_col[i_x][j_y].setLabel([-2,2])
 Landa_col[i_x][j_y+1].setLabel([5,7])
 Landa_col[i_x][j_y+2].setLabel([10,-10])

 PEPS_AA=((((A*Landa_row[i_x][j_y])*Landa_row[i_x+1][j_y])*Landa_col[i_x][j_y])*Landa_col[i_x][j_y+1])
 PEPS_BB=(((PEPS_B*Landa_row[i_x][j_y+1])*Landa_row[i_x+1][j_y+1])*Landa_col[i_x][j_y+2])

 PEPS_AA.permute([-1,-2,3,-4,7],3)
 PEPS_BB.permute([-6,7,8,-9,-10],3)

 if N_iterF[1]=="full":
  Theta=(PEPS_AA*U_ham)*PEPS_BB
  Theta.permute([-1,-2,11,-4,-6,12,-9,-10],4)

  U, V, LA=TU.setTruncation(Theta, sum(D_dim))
  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
 ##############################################################################

  blk_qnums = LA.blockQnum()
  Landa_col[i_x][j_y+1].assign(LA.bond()) 
  for qnum in blk_qnums:
   Landa_col[i_x][j_y+1].putBlock(qnum,LA.getBlock(qnum))

  U.setLabel([-1,-2,11,-4,7])
  U.permute([-1,-2,11,-4,7],3)

  V.setLabel([7,-6,12,-9,-10])
  V.permute([-6,7,12,-9,-10],3)

  invLanda1=inverse_ten(Landa_row[i_x][j_y])
  invLanda2=inverse_ten(Landa_col[i_x][j_y])
  invLanda3=inverse_ten(Landa_row[i_x+1][j_y])

  invLanda1.setLabel([1,-1])
  invLanda2.setLabel([2,-2])
  invLanda3.setLabel([-4,4])

  U=((U*invLanda1)*invLanda2)*invLanda3

  U.permute([1,2,11,4,7],3)
  U.setLabel([1,2,3,4,5])
  U=U*1.0
  U.setLabel([ 1, 2, 3, 4, 5])
  U.permute([ 1, 2, 3, 4, 5], 6)
  Swap=fermionicOPT(Sys,U.bond(2), U.bond(3))
  Swap.setLabel([-3,-4,3,4])
  U=Swap*U
  U.permute([ 1, 2, -3, -4, 5], 3)
  U.setLabel([ 1, 2, 3, 4, 5])


  invLanda8=inverse_ten(Landa_row[i_x+1][j_y+1])
  invLanda9=inverse_ten(Landa_row[i_x][j_y+1])
  invLanda10=inverse_ten(Landa_col[i_x][j_y+2])

  invLanda8.setLabel([-9,9])
  invLanda9.setLabel([6,-6])
  invLanda10.setLabel([-10,10])
  V=((V*invLanda8)*invLanda9)*invLanda10
  V.permute([6,7,12,9,10],3)

 elif N_iterF[1]=="QR":

  A=PEPS_AA*1.0
  A.permute([ -1, -2, -4, 3, 7], 3)

  row, colm=cal_rowcol(A)
  if (row<=colm):
   q, V, s=TU.setTruncation(A, row)
  else:
   q, V, s=TU.setTruncation(A, colm)

  s.setLabel([1,0])
  V.setLabel([ 0, 3, 7])
  r_u=V*s
  r_u.permute([ 1, 3, 7],1)
  q.setLabel([ -1, -2, -4, -100])
  r_u.setLabel([ -100, 3, 7])

  A=copy.copy(PEPS_BB)
  A.setLabel([-6,7,8,-9,-10])
  A.permute([7,8,-6,-9,-10],2)

  row, colm=cal_rowcol(A)
  if (row<=colm):
   U, qq, s=TU.setTruncation(A, row)
  else:
   U, qq, s=TU.setTruncation(A, colm)

  s.setLabel([0,1])
  U.setLabel([7,8,0])
  l_u=U*s
  l_u.permute([7,8,1],2)

  qq.setLabel([-400,-6,-9,-10])
  l_u.setLabel([7,8,-400])

  Theta=(l_u*U_ham)*r_u
  Theta.permute([-100,11,-400,12],2)

  #print Theta.printDiagram(), sum(D_dim)
  U, V, LA=TU.setTruncation(Theta, sum(D_dim))

  norm=norm_Symmetry(LA)
  LA=LA*(1.00/norm)
  #print LA

 ##########################################################################################

  blk_qnums = LA.blockQnum()
  Landa_col[i_x][j_y+1].assign(LA.bond()) 
  for qnum in blk_qnums:
   Landa_col[i_x][j_y+1].putBlock(qnum,LA.getBlock(qnum))

  U.setLabel([-100,11,7])
  U=U*q
  U.permute([-1,-2,11,-4,7],3)

  V.setLabel([7,-400,12])
  V=V*qq
  V.permute([-6,7,12,-9,-10],3)

  invLanda1=inverse_ten(Landa_row[i_x][j_y])
  invLanda2=inverse_ten(Landa_col[i_x][j_y])
  invLanda3=inverse_ten(Landa_row[i_x+1][j_y])

  invLanda1.setLabel([1,-1])
  invLanda2.setLabel([2,-2])
  invLanda3.setLabel([-4,4])

  U=((U*invLanda1)*invLanda2)*invLanda3

  U.permute([1,2,11,4,7],3)
  U.setLabel([1,2,3,4,5])
  U=U*1.0
  U.setLabel([ 1, 2, 3, 4, 5])
  U.permute([ 1, 2, 3, 4, 5], 6)
  Swap=fermionicOPT(Sys,U.bond(2), U.bond(3))
  Swap.setLabel([-3,-4,3,4])
  U=Swap*U
  U.permute([ 1, 2, -3, -4, 5], 3)
  U.setLabel([ 1, 2, 3, 4, 5])

  invLanda8=inverse_ten(Landa_row[i_x+1][j_y+1])
  invLanda9=inverse_ten(Landa_row[i_x][j_y+1])
  invLanda10=inverse_ten(Landa_col[i_x][j_y+2])

  invLanda8.setLabel([-9,9])
  invLanda9.setLabel([6,-6])
  invLanda10.setLabel([-10,10])
  V=((V*invLanda8)*invLanda9)*invLanda10
  V.permute([6,7,12,9,10],3)

 return   U,  V



#############@profile########################
def simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iterF, H_col,H_row, d, Model, N_x, D,chi_try, threshold, interval, Sys):


 #print H, H0, HN

 for i in xrange(1,600):

  delta=start_itebd/pow(division_itebd,i) 

  if delta>1.0e-1:
   N_iter=N_iterF
  if delta<1.0e-1 and delta>1.0e-3:
   N_iter=N_iterF
  if delta<1.0e-3  and delta>1.0e-7:
   N_iter=N_iterF
  if delta<1.0e-8:
   break



  U_evolv=1.0
  print 'delta =', delta
  print "N_iterF=", N_iterF



#  for i_x in xrange(N_x):
#   for j_y in xrange(N_x):
#    PEPS_listten[i_x][j_y]=max_ten(PEPS_listten[i_x][j_y])
    #norm=norm_Symmetry(PEPS_listten[i_x][j_y])
    #PEPS_listten[i_x][j_y]=PEPS_listten[i_x][j_y]*(1.00/norm)

#  PEPS_listten=make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x)
#  PEPS_listten, norm_val, count_val=Normalize_PEPS( PEPS_listten, N_x, D, chi_try, d, threshold, interval)
#  Landa_col_inv,Landa_row_inv=inv_landa_col_row( Landa_col, Landa_row, N_x)
#  PEPS_listten=make_PEPS_tensors(PEPS_listten, Landa_row_inv, Landa_col_inv,N_x)


  for q in xrange(N_iterF[0]):

   for i_x in xrange(N_x-1):
    for j_y in xrange(N_x):
     #print i_x, j_y, H_row[i_x][j_y]
     U_ham = uni10.UniTensor( H_row[i_x][j_y].bond(), "U")
     A = H_row[i_x][j_y].blockQnum()
     blk_qnums = H_row[i_x][j_y].blockQnum()
     for qnum in blk_qnums:
       U_ham.putBlock(qnum, uni10.takeExp(-delta, H_row[i_x][j_y].getBlock(qnum)))

     Peps_A, Peps_B=update_row(PEPS_listten[i_x][j_y], PEPS_listten[i_x+1][j_y], U_ham, i_x, j_y, Landa_col, Landa_row, D,Sys,N_iterF)
     PEPS_listten[i_x][j_y]=Peps_A*1.0
     PEPS_listten[i_x+1][j_y]=Peps_B*1.0

   for j_y in xrange(N_x-1):
    for i_x in xrange(N_x):

     U_ham = uni10.UniTensor( H_col[i_x][j_y].bond(), "U")
     blk_qnums = H_col[i_x][j_y].blockQnum()
     for qnum in blk_qnums:
       U_ham.putBlock(qnum, uni10.takeExp(-delta, H_col[i_x][j_y].getBlock(qnum)))

     Peps_A, Peps_B=update_col(PEPS_listten[i_x][j_y], PEPS_listten[i_x][j_y+1], U_ham, i_x, j_y, Landa_col, Landa_row, D,Sys,N_iterF)
     PEPS_listten[i_x][j_y]=Peps_A*1.0
     PEPS_listten[i_x][j_y+1]=Peps_B*1.0

# for i_x in xrange(N_x):
#  for j_y in xrange(N_x):
#   #PEPS_listten[i_x][j_y]=max_ten(PEPS_listten[i_x][j_y])
#   norm=norm_Symmetry(PEPS_listten[i_x][j_y])
#   PEPS_listten[i_x][j_y]=PEPS_listten[i_x][j_y]*(1.00/norm)


 return PEPS_listten, Landa_col, Landa_row


def make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x): 

 for i_x in xrange(N_x):
  for j_y in xrange(N_x):
   PEPS_listten[i_x][j_y].setLabel([1,2,3,4,5]) 
   Landa_col[i_x][j_y].setLabel([-2,2])
   Landa_row[i_x][j_y].setLabel([-1,1])
   result=(PEPS_listten[i_x][j_y]*Landa_row[i_x][j_y])*Landa_col[i_x][j_y]
   result.permute([-1,-2,3,4,5],3)
   PEPS_listten[i_x][j_y]=result*1.0

 return PEPS_listten



def  Store_Landa_row(Landa_row, N_x):
 for j in xrange(N_x):
  for i in xrange(N_x+1):
   Landa_row[i][j].save("StoreGamma/row" + str(i)+"P"+str(j))


def  Store_Landa_col( Landa_col, N_x):
 for j in xrange(N_x+1):
  for i in xrange(N_x):
   Landa_col[i][j].save("StoreGamma/col" + str(i)+"P"+str(j))








def  Store_Landa_rowRG(Landa_row, N_x):
 for j in xrange(N_x):
  for i in xrange(N_x+1):
   Landa_row[i][j].save("StoreGammaRG/row" + str(i)+"P"+str(j))


def  Store_Landa_colRG( Landa_col, N_x):
 for j in xrange(N_x+1):
  for i in xrange(N_x):
   Landa_col[i][j].save("StoreGammaRG/col" + str(i)+"P"+str(j))







def  Reload_Landa_row( Landa_row, N_x):
 for j in xrange(N_x):
  for i in xrange(N_x+1):
   Landa_row[i][j]=uni10.UniTensor("StoreGamma/row" + str(i)+"P"+str(j))

def  Reload_Landa_col( Landa_col, N_x):
 for j in xrange(N_x+1):
  for i in xrange(N_x):
   Landa_col[i][j]=uni10.UniTensor("StoreGamma/col" + str(i)+"P"+str(j))



def  Store_Gamma( PEPS_listten, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   PEPS_listten[i][j].save("StoreGamma/Gamma" + str(i)+"P"+str(j))


def  Reload_Gamma( PEPS_listten, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   PEPS_listten[i][j]=uni10.UniTensor("StoreGamma/Gamma" + str(i)+"P"+str(j))





def  Store_GammaRG( PEPS_listten, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   PEPS_listten[i][j].save("StoreGammaRG/Gamma" + str(i)+"P"+str(j))


def  Reload_GammaRG( PEPS_listten, N_x):
 for i in xrange(N_x):
  for j in xrange(N_x):
   PEPS_listten[i][j]=uni10.UniTensor("StoreGammaRG/Gamma" + str(i)+"P"+str(j))









def inverse_ten(Landa2):
 invLanda2=uni10.UniTensor(Landa2.bond())
 blk_qnums=Landa2.blockQnum()
 for qnum in blk_qnums:
  D=int(Landa2.getBlock(qnum).row())
  D1=int(Landa2.getBlock(qnum).col())
  invL2 = uni10.Matrix(D, D1,True)
  invLt = uni10.Matrix(D, D1,True)
  invLt=Landa2.getBlock(qnum,True)
  #print invLt[0], invLt[1], invLt[2], invLt[3]
  for i in xrange(D):
      invL2[i] = 0 if ((invLt[i].real) < 1.0e-10) else (1.00 / (invLt[i].real))

  invLanda2.putBlock(qnum,invL2)
 return invLanda2

#####@profile
def  cal_energy_double(PEPS_listten, N_x, H_col, H_row, d, Model, chi_boundry, D,Sys):

 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_x)
 for i_ind in xrange(mps_I.N):
  T.identity()
  mps_I[i_ind]=T*1.0


 MPO_Ten=make_boundry_MPO(PEPS_listten, N_x,Sys)
 Env_left, Env_right=make_ENVMPS(MPO_Ten, N_x, chi_boundry)


 E_coulmn_t=[]
 for i_ind in xrange(N_x):
  if i_ind==0:
   energy_double_col(mps_I, Env_right[i_ind+1], PEPS_listten[i_ind], N_x, H_col[i_ind], D, E_coulmn_t,Sys)
  elif i_ind==N_x-1:
   energy_double_col(Env_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_x, H_col[i_ind], D, E_coulmn_t,Sys)
  else:
   energy_double_col(Env_left[i_ind-1], Env_right[i_ind+1], PEPS_listten[i_ind], N_x, H_col[i_ind], D, E_coulmn_t,Sys)

 E_c=sum(E_coulmn_t)
 E_coulmn=[ E_coulmn_t[i]*1.0   for i in xrange(len(E_coulmn_t))  ]

#################   Row   #################


 MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,Sys)
 Env_up, Env_down=make_ENV_updown( MPO_Ten, N_x, chi_boundry)

 E_row_t=[]
 for i_ind in xrange(N_x):
  peps_l=[]
  H_l=[]
  for i in xrange(N_x):
   peps_l.append(PEPS_listten[i][i_ind])
  for i in xrange(N_x-1):
   H_l.append(H_row[i][i_ind])

  if i_ind==0:
   energy_double_row( mps_I, Env_up[i_ind+1], peps_l, N_x, H_l, D, E_row_t,Sys)
  elif i_ind==N_x-1:
   energy_double_row( Env_down[i_ind-1], mps_I, peps_l, N_x, H_l, D, E_row_t,Sys)
  else:
   energy_double_row( Env_down[i_ind-1], Env_up[i_ind+1], peps_l, N_x, H_l, D, E_row_t,Sys)



 E_r=sum(E_row_t)
 #print "E=", E_r, E_c, (E_c+E_r)/(N_x*N_x)
 E_row=[ E_row_t[i]*1.0   for i in xrange(len(E_row_t))  ]
 #E_0=E_1*1.0
 E_1=(E_c+E_r)
 #print E_row_t, E_coulmn_t

 return E_1



def  energy_double_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x,  H , D, E_coulmn,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):

  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[i]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[i]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   #E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=Swap1*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[i]*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)



# for i in xrange(len(PEPS_listten)-1):
#  E_list_left[i].setLabel([1,2,3,4])
#  E_list_right[i+1].setLabel([1,2,3,4])
#  A=E_list_left[i]*E_list_right[i+1]
#  print "Y1", A[0]



 for i in xrange(N_x-1):



  E_val=energy_row_double_local( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, H[i],D,Sys)
  E_coulmn.append(E_val)



########@profile
def  energy_row_double_local( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, H_orig, D,Sys ):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)


 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)


 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],3)
 r_u.setLabel([-100,3,9])

 r_u.setLabel([-100,3,10])   #new
 r_u.permute([-100,3,10],3)



 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)



 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,10,200])
 q_d.permute([6,7,10,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-10,-200,10,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-10,-200],4)



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)

 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],3)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,10,14,15])
 qq_d.permute([400,10,14,15],4)

 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel([-400, -10, 400, 10])
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-10,-14,-15],4)


 l_u.setLabel([9,13,-300])

 l_u.setLabel([10,13,-300])  #new
 l_u.permute([10,13,-300],3)

 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])  #new


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 Swap4.setLabel([23,400,10,-400])


######################################################################

 mps_boundry_up[Location].setLabel([16,10,25,27])
 #mps_boundry_up[Location*2+1].setLabel([26,25,27])
 mps_boundry_up[Location+1].setLabel([27,15,24,21])
 #mps_boundry_up[Location*2+3].setLabel([28,24,21])

 mps_boundry_down[Location].setLabel([18,22,-7,32])
 #mps_boundry_down[Location*2+1].setLabel([31,-7,32])
 mps_boundry_down[Location+1].setLabel([32,23,-10,19])
 #mps_boundry_down[Location*2+3].setLabel([33,-10,19])



######################################################

 A=E_left*mps_boundry_up[Location]
 A=A*(Swap1*q)
 A=A*mps_boundry_down[Location]
 A=A*mps_boundry_up[Location+1]
 A=A*(Swap3*q_d)
 A=A*mps_boundry_down[Location+1]


 B=E_right#*mps_boundry_up[2*Location+3]
 B=(Swap2*qq_d)*B
 #B=B*mps_boundry_down[2*Location+3]

 #B=B*mps_boundry_up[2*Location+2]
 B=B*(Swap4*qq)
 #B=B*mps_boundry_down[2*Location+2]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])



 H_orig.setLabel([-3,-13,3,13])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig
 #print "E_1", Norm_h[0], h_h[0]/Norm_h[0]



 return h_h[0]/Norm_h[0]








########@profile
def energy_double_col( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x,  H , D, E_coulmn,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 bdi=uni10.Bond(uni10.BD_IN, D)
 bdo=uni10.Bond(uni10.BD_OUT, D)

 bdi1=uni10.Bond(uni10.BD_IN, 1)
 bdo1=uni10.Bond(uni10.BD_OUT, 1)

 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])

   iden.setLabel([5])
   Peps_list=Peps_list*iden
   iden.setLabel([11])
   Peps_list_conj=Peps_list_conj*iden

   PEP_com=((Peps_list*Swap1)*Swap4)*((Peps_list_conj*Swap3)*Swap2)
   
   iden.setLabel([13])
   mps_boundry_leftA=mps_boundry_left[i]*iden
   iden.setLabel([15])
   mps_boundry_rightA=mps_boundry_right[i]*iden

   result=(PEP_com*mps_boundry_rightA)*mps_boundry_leftA

   result.permute([12,10,9,14],4)
   E_list_up[i]=result
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   PEP_com=((Peps_list*Swap1)*Swap4)*((Peps_list_conj*Swap3)*Swap2)


   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   E_list_up[i+1].setLabel([13,5,11,15])
   result=((((PEP_com*E_list_up[i+1])*mps_boundry_left[i])*mps_boundry_right[i]))
   result.permute([12,10,9,14],4)
   E_list_up[i]=result



 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:

   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   PEP_com=((Peps_list*Swap1)*Swap4)*((Peps_list_conj*Swap3)*Swap2)

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   iden.setLabel([10])
   PEP_com=PEP_com*iden
   iden.setLabel([9])
   PEP_com=PEP_com*iden

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[i]*iden

   result=((PEP_com*mps_leftA)*mps_rightA)

   result.permute([13,5,11,15],4)
   E_list_down[i]=result

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)

   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   PEP_com=((Peps_list*Swap1)*((Peps_list_conj*Swap3)*Swap2))*Swap4

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])


   E_list_down[i-1].setLabel([12,10,9,14])

   result=((((PEP_com*E_list_down[i-1])*mps_boundry_left[i])*mps_boundry_right[i]))
   result.permute([13,5,11,15],4)
   E_list_down[i]=result

#  for i in xrange(len(PEPS_listten)-1):
#   E_list_down[i].setLabel([1,2,3,4])
#   E_list_up[i+1].setLabel([1,2,3,4])
#   A=E_list_down[i]*E_list_up[i+1]
#   print "double" , A[0]


 for i in xrange(N_x-1):


  E_val=double_local_energy_col( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H[i],D,Sys)
  E_coulmn.append(E_val)



########@profile
def  double_local_energy_col(PEPS_listten,E_list_down,E_list_up,mps_boundry_left,mps_boundry_right, Location, H_orig,D,Sys):

 bdi_mid=uni10.Bond( uni10.BD_IN, 1)
 iden=uni10.UniTensor( [bdi_mid] )
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],2)


 r_u.setLabel([-100,3,10])


 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


#####################################################

 A=Peps_1*1.0

 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)


 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])
 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],3)
 A_conj.setLabel([-6,-7,-3,-9,-10])

 A_conj.permute([-6,-7,-9,-3,-10],3)



 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  q_d, V, s=TU.setTruncation(A_conj, row)
 else:
  q_d, V, s=TU.setTruncation(A_conj, colm)


 s.setLabel([1,0])
 V.setLabel([0,-3,-10])
 r_d=V*s
 r_d.permute([1,-3,-10],3)
 
 q_d.setLabel([-6,-7,-9,-200])
 q_d.permute([-6,-7,-9,-200],4)

 r_d.setLabel([-200,-3,-10])



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 A=Peps_2*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],3)
 A_conj=A*1.0
 A_conj.transpose()
 A_conj.setLabel([-4,-5,-1,-2,3])
 A_conj.permute([-1,-2,3,-4,-5],5)

 Swap2=fermionicOPT(Sys,A_conj.bond(3), A_conj.bond(4))
 Swap2.setLabel([-6,7,-4,-5])
 Swap3=fermionicOPT(Sys,A_conj.bond(0), A_conj.bond(1))
 Swap3.setLabel([8,9,-1,-2])


 A_conj=(A_conj*Swap2)*Swap3
 A_conj.permute([8,9,3,-6,7],5)
 Swap1=fermionicOPT(Sys,A_conj.bond(2), A_conj.bond(0))
 Swap1.setLabel([11,10,3,8])
 A_conj=A_conj*Swap1

 A_conj.permute([10,9,11,-6,7],3)
 A_conj.setLabel([-11,-10,-13,-14,-15])
 A_conj.permute([-10,-13,-11,-14,-15],2)


 row, colm=cal_rowcol(A_conj)
 if (row<=colm):
  U, qq_d, s=TU.setTruncation(A_conj, row)
 else:
  U, qq_d, s=TU.setTruncation(A_conj, colm)

 s.setLabel([0,1])
 U.setLabel([-10,-13,0])
 l_d=U*s
 l_d.permute([-10,-13,1],3)

 qq_d.setLabel([-400,-11,-14,-15])
 qq_d.permute([-400,-11,-14,-15],4)


 l_d.setLabel([-10,-13,-400])


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([ 90, 200, 9, -200 ])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])

######################################################################
 mps_boundry_left[Location].setLabel([16,6,-60,18])
 mps_boundry_left[Location+1].setLabel([18,11,-110,25])


 mps_boundry_right[Location].setLabel([17,90,-9,19])
 mps_boundry_right[Location+1].setLabel([19,140,-14,24])

######################################################

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*q)
 A=A*(Swap3*q_d)
 A=A*mps_boundry_right[Location]

 B=E_right*mps_boundry_left[Location+1]
 B=B*(Swap4*qq)
 B=B*(Swap2*qq_d)
 B=B*mps_boundry_right[Location+1]
 
 N_ten=A*B
 #print N_ten.printDiagram()
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 #N_ten=N_Positiv(N_ten)
######################################################


 iden_h=copy.copy(H_orig)
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])
 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])

 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig
 #print "Fnit_norm", Norm_h[0], h_h[0]/Norm_h[0]

 return  h_h[0]/Norm_h[0]


def  TEBD_Full_double(H_col, H_row, N_x, PEPS_listten, E_iter_list, D, accuracy, N_tebd, i, d, chi_boundry, chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model,count_list,E_iter_list1,Sys):

 PEPS_mps_leftU=[None]*N_x
 PEPS_mps_rightU=[None]*N_x
 PEPS_listtenU=[None]*N_x
 mps_boundry_left=[None]*N_x
 mps_boundry_right=[None]*N_x
 mps_boundry_temp=[None]*N_x

 for i in xrange(N_x):
  PEPS_listtenU[i]=Init_PEPS( N_x, D, d, i)

 Fidel_val=1
 E_coulmn=[]
 E_row=[]
 E_mag_coulmn=[]
 E_mag_row=[]
 E_0=1.0
 E_1=1.0

 E_00=1.0
 E_11=1.0


 bdi=uni10.Bond(uni10.BD_IN,1)
 bdo=uni10.Bond(uni10.BD_OUT,1) 
 T=uni10.UniTensor([bdi,bdi,bdi,bdo])
 mps_I=MPSclass.MPS(1,1,N_x)
 for i_ind in xrange(mps_I.N):
  T.identity()
  mps_I[i_ind]=T*1.0

 List_delN=Short_TrotterSteps(N_tebd)

 print List_delN
 E_00_f=10.0
 E_min=1
 E_0=1
 E_1=100
 count_iter=0
 for delta, N_iter in List_delN:
  #delta=0
  print delta, N_iter

  U_ham_col=make_U_col(H_col, N_x, delta)
  U_ham_row=make_U_row(H_row, N_x, delta)


  PEPS_listten, norm_val, count_val=Normalize_PEPS(PEPS_listten, N_x, D, chi_try, d, threshold, interval,Sys)
  Energy_val=cal_energy_double(PEPS_listten, N_x, H_col,H_row, d, Model, chi_boundry, D,Sys)
  E_0=Energy_val
  E_11=Energy_val*1.0
  print "E_0", Energy_val

  PEPS_listtenU=copy_f(PEPS_listten,N_x, PEPS_listtenU)
  PEPS_listten=copy_f( PEPS_listtenU, N_x, PEPS_listten)

  for q_iter in xrange(N_iter):
   #Energy_val=cal_energy_double(PEPS_listten, N_x, H_col,H_row, d, Model, chi_boundry, D,Sys)
   #print q_iter, "E_0", Energy_val

###########################   Col   #############################
###   make_right_Env   ######
   MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,Sys)
   Env_left, Env_right=make_ENVMPS( MPO_Ten, N_x, chi_boundry)
#############
   E_coulmn_t=[]
   for i_ind in xrange(N_x):
    if i_ind==0:
     update_local_double(mps_I, Env_right[i_ind+1], PEPS_listten[i_ind], N_x, U_ham_col[i_ind], H_col[i_ind], D, E_coulmn_t,threshold, interval,Sys)
    elif i_ind==N_x-1:
     update_local_double(Env_left[i_ind-1], mps_I, PEPS_listten[i_ind], N_x, U_ham_col[i_ind], H_col[i_ind], D, E_coulmn_t,threshold, interval,Sys)
    else:
     update_local_double(Env_left[i_ind-1], Env_right[i_ind+1], PEPS_listten[i_ind], N_x, U_ham_col[i_ind], H_col[i_ind], D, E_coulmn_t,threshold, interval,Sys)

    if  i_ind==(N_x-1):break
    MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,Sys)
    Env_left, Env_right=make_ENVMPS( MPO_Ten, N_x, chi_boundry)

##########################   Row   ###############################
   #PEPS_listten=rotate_all(PEPS_listten, N_x)


   MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,Sys)
   Env_up, Env_down=make_ENV_updown( MPO_Ten, N_x, chi_boundry)

   #Env_left, Env_right=make_ENVMPS( MPO_Ten, N_x, chi_boundry)

   E_row_t=[]
   for i_ind in xrange(N_x):

    peps_l=[]
    H_l=[]
    U_l=[]

    for i in xrange(N_x):
     peps_l.append(PEPS_listten[i][i_ind])

    for i in xrange(N_x-1):
     H_l.append(H_row[i][i_ind])

    for i in xrange(N_x-1):
     U_l.append(U_ham_row[i][i_ind])


    if i_ind==0:
     update_local_double_row( mps_I, Env_up[i_ind+1], peps_l, N_x, U_l, H_l, D, E_row_t, threshold, interval,Sys)
    elif i_ind==N_x-1:
     update_local_double_row( Env_down[i_ind-1], mps_I, peps_l, N_x, U_l, H_l, D, E_row_t, threshold, interval,Sys)
    else:
     update_local_double_row( Env_down[i_ind-1], Env_up[i_ind+1], peps_l, N_x, U_l, H_l, D, E_row_t, threshold, interval,Sys)

    for i in xrange(N_x):
     PEPS_listten[i][i_ind]=peps_l[i]

    if i_ind==(N_x-1):break
    MPO_Ten=make_boundry_MPO( PEPS_listten, N_x,Sys)
    Env_up, Env_down=make_ENV_updown( MPO_Ten, N_x, chi_boundry)


   E_00=E_11*1.0
   E_11=(sum(E_coulmn_t)+sum(E_row_t))
   #print E_coulmn_t, E_row_t

   print "E_t", q_iter, E_11, E_00, sum(E_coulmn_t), sum(E_row_t),(E_11-E_00) / E_11 
   E_iter_list1.append(E_11)

   if E_11 < E_00 and q_iter>0:
    Store_f(PEPS_listten, N_x)
    PEPS_listtenU=copy_f(PEPS_listten,N_x, PEPS_listtenU)
    #E_00_f=E_11*1.0


   if abs((E_11-E_00)/E_11)<accuracy or E_11>E_00:
    if abs((E_11-E_00)/E_11)<accuracy:print "loop_finished:reached_levelofaccuracy"
    else: 
     print "Didnot_get_lowered"
    PEPS_listten=copy_f( PEPS_listtenU, N_x, PEPS_listten)
    Store_f(PEPS_listten, N_x)
    break

  Energy_val=cal_energy_double(PEPS_listten, N_x, H_col,H_row, d, Model, chi_boundry, D,Sys)
  print "E_f", Energy_val

#  return PEPS_listten

#####@profile
def update_local_double( mps_boundry_left, mps_boundry_right, PEPS_listten, N_x, U_ham, H, D, E_coulmn,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_up=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):
  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])

   iden.setLabel([13])
   mps_left=mps_boundry_left[i]*iden
   iden.setLabel([15])
   mps_right=mps_boundry_right[i]*iden

   iden.setLabel([5])
   iden1=iden*1.0
   iden1.setLabel([11])

   E_list_up[i]=Peps_list*iden
   E_list_up[i]=E_list_up[i]*(Swap1*iden1)
   E_list_up[i]=E_list_up[i]*mps_left
   E_list_up[i]=E_list_up[i]*mps_right
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=E_list_up[i]*((Peps_list_conj*Swap3)*Swap2)

   E_list_up[i].permute([12,10,9,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   mps_boundry_right[i].setLabel([14,6,-6,15])

   E_list_up[i+1].setLabel([13,5,11,15])

   E_list_up[i]=((Peps_list*Swap1))*E_list_up[i+1]
   E_list_up[i]=mps_boundry_left[i]*E_list_up[i]
   E_list_up[i]=mps_boundry_right[i]*E_list_up[i]
   E_list_up[i]=E_list_up[i]*Swap4
   E_list_up[i]=((Peps_list_conj*Swap3)*Swap2)*E_list_up[i]
   E_list_up[i].permute([12,10,9,14],4)

 E_list_down=[None]*N_x
 for i in xrange(len(PEPS_listten)):
  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   #mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[i].setLabel([14,6,-6,15])
   #mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])

   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   #E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   #E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   #mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[i].setLabel([14,6,-6,15])
   #mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])

   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[i]
   #E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   #E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)


# for i in xrange(len(PEPS_listten)-1):
#  E_list_down[i].setLabel([1,2,3,4])
#  E_list_up[i+1].setLabel([1,2,3,4])
#  A=E_list_down[i]*E_list_up[i+1]
#  print "Inside", A[0]


 for i in xrange(N_x-1):

#  PEPS_f, PEPS_s, E_val=double_local_update( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H[i], U_ham[i], D, threshold, interval,Sys)

  if Sys[2]=="QR":
   PEPS_f, PEPS_s, E_val=double_local_update( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H[i], U_ham[i], D, threshold, interval,Sys)
  elif Sys[2]=="Inv":
   PEPS_f, PEPS_s, E_val=double_local_update_inv( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, i, H[i], U_ham[i], D, threshold, interval,Sys)




  E_coulmn.append(E_val)
  PEPS_listten[i]=PEPS_f*1.0
  PEPS_listten[i+1]=PEPS_s*1.0

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   #mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[i].setLabel([14,6,-6,15])
   #mps_boundry_right[2*i+1].setLabel([150,6,15])

   iden.setLabel([12])
   mps_leftA=mps_boundry_left[i]*iden
   iden.setLabel([14])
   mps_rightA=mps_boundry_right[i]*iden

   iden.setLabel([10])
   iden1=iden*1.0
   iden1.setLabel([9])


   E_list_down[i]=(((Peps_list_conj*(iden1*Swap3))*Swap2))*mps_leftA
   E_list_down[i]=(E_list_down[i]*(Swap4*iden))*mps_rightA
   #E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   #E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*(Peps_list*Swap1)
   E_list_down[i].permute([13,5,11,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_left[i].setLabel([12,1,-8,13])
   #mps_boundry_left[2*i+1].setLabel([130,1,13])

   mps_boundry_right[i].setLabel([14,6,-6,15])
   #mps_boundry_right[2*i+1].setLabel([150,6,15])

   E_list_down[i-1].setLabel([12,10,9,14])


   E_list_down[i]=(E_list_down[i-1]*Swap4)*Swap3
   E_list_down[i]=E_list_down[i]*(Peps_list_conj*Swap2)
   E_list_down[i]=E_list_down[i]*mps_boundry_left[i]
   E_list_down[i]=E_list_down[i]*mps_boundry_right[i]
   #E_list_down[i]=E_list_down[i]*mps_boundry_right[2*i+1]
   E_list_down[i]=E_list_down[i]*Swap1
   E_list_down[i]=E_list_down[i]*Peps_list
   #E_list_down[i]=E_list_down[i]*mps_boundry_left[2*i+1]
   E_list_down[i].permute([13,5,11,15],4)



#####@profile
def double_local_update( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, H_orig, Ham ,D, threshold, interval,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=copy.copy(PEPS_listten[Location])
 Peps_2=copy.copy(PEPS_listten[Location+1])

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.permute([1,2,3,4,5],5)

 Swap1=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap1.setLabel([6,7,4,3])
 A=Swap1*A
 A.permute([1,2,7,6,5],3)

 A.setLabel([6,7,3,9,10])
 A.permute([6,7,9,3,10],3)



 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)



 s.setLabel([1,0])
 V.setLabel([0,3,10])
 r_u=V*s
 r_u.permute([1,3,10],3)


 r_u.setLabel([-100,3,10])

 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])


 q.setLabel([6,7,9,-100])
 q.permute([6,7,9,-100],4)


 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,9,200])
 q_d.permute([6,7,9,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-9,-200,9,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-9,-200],4)

#####################################################

 A=Peps_2*1.0
 A.setLabel([11,10,13,14,15])
 A.permute([10,13,11,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)


 s.setLabel([0,1])
 U.setLabel([10,13,0])
 l_u=U*s
 l_u.permute([10,13,1],3)

 qq.setLabel([-300,11,14,15])
 qq.permute([-300,11,14,15],4)

 l_u.setLabel([10,13,-300])


 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])


 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,11,14,15])
 qq_d.permute([400,11,14,15],4)
 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel( [-400, -11, 400, 11] )
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-11,-14,-15],4)



######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([70,-60,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([140,-150,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(2), q_d.bond(3) )
 Swap3.setLabel([90,200,9,-200])

 Swap4=fermionicOPT(Sys, qq.bond(0), qq_d.bond(1) )
 Swap4.setLabel([300,-110,-300,-11])


######################################################################

 mps_boundry_left[Location].setLabel([16,6,-60,21])
 mps_boundry_left[Location+1].setLabel([21,11,-110,25])

 mps_boundry_right[Location].setLabel([17,90,-9,20])
 mps_boundry_right[Location+1].setLabel([20,140,-14,24])


######################################################

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*q_d)
 A=A*mps_boundry_right[Location]
 A=A*(Swap3*q)


 B=mps_boundry_left[Location+1]*E_right
 B=B*(Swap2*qq)
 B=B*mps_boundry_right[Location+1]
 B=B*(Swap4*qq_d)


 N_ten=A*B
 N_ten.permute([200,-400,-100,300],2)
 N_ten.setLabel([-200,-400,-100,-300])

 H_orig.setLabel([-3,-13,3,13])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig
 #print "E_1", Norm_h[0], h_h[0]/Norm_h[0]


 N_ten=N_Positiv(N_ten)

######################################################

###############simple_update###########################
 A=r_u*1.0
 A.setLabel([-10,2,3])
 B=l_u*1.0
 B.setLabel([3,6,-30])
 Ham.setLabel([-2,-6,2,6])

 Teta=(A*B)*Ham
 Teta.permute([-10,-2,-6,-30],2)
 
 D_dim=l_u.bond(0).dim()
 #chi_dim=sum(D)
 D_list=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_list.append(dim)

 D_dim=sum(D_list)
 
   
 row, colm=cal_rowcol(Teta)
 if (row<=colm and row<=D_dim):
  U,V,S=TU.setTruncation(Teta,row)
 elif (row<=colm and row>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 elif (row>colm and colm<=D_dim):
  U,V,S=TU.setTruncation(Teta,colm)
 elif (row>colm and colm>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)

 U.setLabel([-10,-2,-3])
 V.setLabel([-3,-6,-30])
 S=Sqrt(S)
 S.setLabel([-3,3])
 r_up=U*S
 r_up.permute([-10,-2,3],3)
 r_up.setLabel([-100,3,10])
 S.setLabel([3,-3])
 l_up=V*S
 l_up.permute([3,-6,-30],3)
 l_up.setLabel([10,13,-300])
###########################################



# l_up=l_u*1.0
# r_up=r_u*1.0
 l_dp=l_up*1.0
 r_dp=r_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])


 Ham.setLabel([51,52,53,54])
 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([-3,-13,3,13])
 Ham.setLabel([-3,-13,3,13])
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])

 r_up_init=r_up*1.0
 r_dp_init=r_dp*1.0
 l_up_init=l_up*1.0
 l_dp_init=l_dp*1.0

 valf=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up, l_dp, r_dp, Ham, iden_h, H)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   r_up=r_up_init*1.0
   r_dp=r_dp_init*1.0
   l_up=l_up_init*1.0
   l_dp=l_dp_init*1.0
   break
  else:
   r_up_init=r_up*1.0
   r_dp_init=r_dp*1.0
   l_up_init=l_up*1.0
   l_dp_init=l_dp*1.0

  r_up, r_dp=optimum_0( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
  l_up, l_dp=optimum_1( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)



 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-12 or abs(valf)>1.0e+12:
  print "warning_norm_in_optimization",  abs(valf), "count", count


 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 h_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*H_orig
 #print "E_2c", Norm_h[0], h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<(threshold[0]):
  r_up,l_up=renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval)


 PEPS_1=r_up*q
 PEPS_1.permute([6,7,3,9,10],3)

 A=PEPS_1*1.0
 A.permute([6,7,3,9,10],5)

 Swap1=fermionicOPT(Sys,A.bond(2), A.bond(3))
 Swap1.setLabel([-3,-9,3,9])
 A=Swap1*A
 A.permute([6,7,-3,-9,10],3)
 PEPS_1=A*1.0


 PEPS_2=l_up*qq
 PEPS_2.permute([11,10,13,14,15],3)


 

 return PEPS_1, PEPS_2, h_h[0]/Norm_h[0]







#####@profile
def update_local_double_row( mps_boundry_down, mps_boundry_up, PEPS_listten, N_x, U_ham, H, D, E_coulmn,threshold, interval,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()


 E_list_right=[None]*N_x
 for i in reversed(xrange(len(PEPS_listten))):

  if i==(len(PEPS_listten)-1):
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([13])
   mps_up=mps_boundry_up[i]*iden
   iden.setLabel([15])
   mps_down=mps_boundry_down[i]*iden

   iden.setLabel([6])
   iden1=iden*1.0
   iden1.setLabel([-6])

   E_list_right[i]=(Swap1*iden)*mps_up
   E_list_right[i]=E_list_right[i]*(Swap2*iden1)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   E_list_right[i]=E_list_right[i]*mps_down
   #E_list_right[i]=E_list_right[i]*mps_boundry_up[2*i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i]
   E_list_right[i].permute([12,1,-8,14],4)
  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_right[i+1].setLabel([13,6,-6,15])

   E_list_right[i]=Swap1*E_list_right[i+1]
   E_list_right[i]=E_list_right[i]*(Swap2)
   E_list_right[i]=E_list_right[i]*((Peps_list_conj*Swap3))
   #E_list_right[i]=E_list_right[i]*mps_boundry_down[2*i+1]
   E_list_right[i]=E_list_right[i]*mps_boundry_up[i]
   E_list_right[i]=E_list_right[i]*Peps_list
   E_list_right[i]=E_list_right[i]*Swap4
   E_list_right[i]=E_list_right[i]*mps_boundry_down[i]
   E_list_right[i].permute([12,1,-8,14],4)





 E_list_left=[None]*N_x
 for i in xrange(len(PEPS_listten)):

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[i]*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)



# for i in xrange(len(PEPS_listten)-1):
#  E_list_left[i].setLabel([1,2,3,4])
#  E_list_right[i+1].setLabel([1,2,3,4])
#  A=E_list_left[i]*E_list_right[i+1]
#  print "Y1", A[0]



 for i in xrange(N_x-1):
#  PEPS_f, PEPS_s, E_val=Update_twotensor_double_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, U_ham[i], H[i], D, threshold, interval,Sys)




  if Sys[2]=="QR":
   PEPS_f, PEPS_s, E_val=Update_twotensor_double_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, U_ham[i], H[i],D,threshold, interval,Sys)
  elif Sys[2]=="Inv":
   PEPS_f, PEPS_s, E_val=Update_twotensor_double_row_inv( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, i, U_ham[i], H[i],D,threshold, interval,Sys)



  E_coulmn.append(E_val)

  PEPS_listten[i]=PEPS_f*1.0
  PEPS_listten[i+1]=PEPS_s*1.0

  if i==0:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])

   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   iden.setLabel([12])
   mps_upA=mps_boundry_up[i]*iden
   iden.setLabel([14])
   mps_downA=mps_boundry_down[i]*iden

   iden.setLabel([1])
   iden1=iden*1.0
   iden1.setLabel([-8])

   E_list_left[i]=mps_upA
   E_list_left[i]=E_list_left[i]*(iden*Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4*iden1)
   E_list_left[i]=mps_downA*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)

  else:
   Peps_list=PEPS_listten[i]*1.0
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list.permute([1,2,3,4,5],3)
   Peps_list.setLabel([1,2,3,4,5])
   Peps_list_conj=Peps_list*1.0
   Peps_list_conj.transpose()
   Peps_list_conj.setLabel([-4,-5,-1,-2,3])
   Peps_list_conj.permute([-1,-2,3,-4,-5],5)
   Peps_list.permute([1,2,3,4,5],5)
   Swap1=fermionicOPT(Sys,Peps_list.bond(3), Peps_list_conj.bond(4))
   Swap1.setLabel([6,11,4,7])
   Swap2=fermionicOPT(Sys,Peps_list_conj.bond(3), Peps_list_conj.bond(4))
   Swap2.setLabel([-6,7,-4,-5])
   Swap3=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list_conj.bond(1))
   Swap3.setLabel([8,9,-1,-2])
   Swap4=fermionicOPT(Sys,Peps_list_conj.bond(0), Peps_list.bond(1))
   Swap4.setLabel([-8,10,8,2])


   mps_boundry_up[i].setLabel([12,5,11,13])
   #mps_boundry_up[2*i+1].setLabel([130,11,13])

   mps_boundry_down[i].setLabel([14,10,9,15])
   #mps_boundry_down[2*i+1].setLabel([150,9,15])

   E_list_left[i-1].setLabel( [12,1,-8,14] )

   E_list_left[i]=mps_boundry_up[i]*E_list_left[i-1]
   E_list_left[i]=E_list_left[i]*(Peps_list)
   E_list_left[i]=E_list_left[i]*(Swap4)
   E_list_left[i]=mps_boundry_down[i]*E_list_left[i]
   #E_list_left[i]=E_list_left[i]*mps_boundry_down[2*i+1]
   E_list_left[i]=E_list_left[i]*((Peps_list_conj*Swap3)*Swap2)
   E_list_left[i]=E_list_left[i]*(Swap1)
   E_list_left[i].permute([13,6,-6,15],4)


#####@profile
def Update_twotensor_double_row( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, Ham, H_orig, D, threshold, interval, Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])

 A=Peps_1*1.0
 A.setLabel([1,2,3,4,5])
 A.setLabel([6,7,3,9,10])
 A.permute([6,7,10,3,9],3)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  q, V, s=TU.setTruncation(A, row)
 else:
  q, V, s=TU.setTruncation(A, colm)

 s.setLabel([1,0])
 V.setLabel([0,3,9])
 r_u=V*s

 r_u.permute([1,3,9],3)
 r_u.setLabel([-100,3,9])

 r_u.setLabel([-100,3,10])
 r_u.permute([-100,3,10],3)

 r_d=r_u*1.0
 r_d.transpose()
 r_d.setLabel([-200,-3,-10])

 q.setLabel([6,7,10,-100])
 q.permute([6,7,10,-100],4)

 q_d=q*1.0
 q_d.transpose()
 q_d.setLabel([6,7,10,200])
 q_d.permute([6,7,10,200],4)
 Swap=fermionicOPT(Sys,q_d.bond(0), q_d.bond(1))
 Swap.setLabel([-6,-7,6,7])
 Swap1=fermionicOPT(Sys,q_d.bond(2), q_d.bond(3))
 Swap1.setLabel([-10,-200,10,200])
 q_d=(q_d*Swap1)*Swap
 q_d.permute([-6,-7,-10,-200],4)



 A=copy.copy(Peps_2)
 A.setLabel([11,10,13,14,15])
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([11,-10,-13,14,15],3)
 A.setLabel([9,10,13,14,15])

 A.permute([9,13,10,14,15],2)

 row, colm=cal_rowcol(A)
 if (row<=colm):
  U, qq, s=TU.setTruncation(A, row)
 else:
  U, qq, s=TU.setTruncation(A, colm)

 s.setLabel([0,1])
 U.setLabel([9,13,0])
 l_u=U*s
 l_u.permute([9,13,1],3)

 qq.setLabel([-300,10,14,15])
 qq.permute([-300,10,14,15],4)

 qq_d=qq*1.0
 qq_d.transpose()
 qq_d.setLabel([400,10,14,15])
 qq_d.permute([400,10,14,15],4)

 Swap=fermionicOPT(Sys, qq_d.bond(0), qq_d.bond(1) )
 Swap.setLabel([-400, -10, 400, 10])
 Swap1=fermionicOPT(Sys, qq_d.bond(2), qq_d.bond(3) )
 Swap1.setLabel([-14,-15,14,15])
 qq_d=(qq_d*Swap1)*Swap
 qq_d.permute([-400,-10,-14,-15],4)


 l_u.setLabel([9,13,-300])

 l_u.setLabel([10,13,-300])
 l_u.permute([10,13,-300],3)

 l_d=l_u*1.0
 l_d.transpose()
 l_d.setLabel([-10,-13,-400])


######################################################################
 Swap1=fermionicOPT(Sys, q.bond(1), q_d.bond(0) )
 Swap1.setLabel([22,17,7,-6])

 Swap2=fermionicOPT(Sys, qq.bond(2), qq_d.bond(3) )
 Swap2.setLabel([20,24,14,-15])


 Swap3=fermionicOPT(Sys, q.bond(3), q_d.bond(2) )
 Swap3.setLabel([100,25,-100,-10])

 Swap4=fermionicOPT(Sys, qq.bond(1), qq_d.bond(0) )
 Swap4.setLabel([23,400,10,-400])

######################################################################

 mps_boundry_up[Location].setLabel([16,10,25,27])
 mps_boundry_up[Location+1].setLabel([27,15,24,21])

 mps_boundry_down[Location].setLabel([18,22,-7,32])
 mps_boundry_down[Location+1].setLabel([32,23,-10,19])

######################################################

 A=E_left*mps_boundry_up[Location]
 A=A*(q*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*(Swap3*q_d)

 B=E_right*mps_boundry_up[Location+1]
 B=(Swap2*qq_d)*B

 B=B*mps_boundry_down[Location+1]
 B=B*(Swap4*qq)
 N_ten=A*B
 
 N_ten.permute([-200,400,100,-300],2)
 N_ten.setLabel([-200,-400,-100,-300])



 H_orig.setLabel([-3,-13,3,13])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*iden_h
 h_h=(((r_u*r_d)*N_ten)*(l_u*l_d))*H_orig
 #print "E_1", Norm_h[0], h_h[0]/Norm_h[0]


 N_ten=N_Positiv(N_ten)
######################################################



######################################################

#############simple_update###########################
 A=r_u*1.0
 A.setLabel([-10,2,3])
 B=l_u*1.0
 B.setLabel([3,6,-30])
 Ham.setLabel([-2,-6,2,6])
 
 
 Teta=(A*B)*Ham
 Teta.permute([-10,-2,-6,-30],2)
 
 
 D_dim=l_u.bond(0).dim()

 D_list=[]
 bdi = uni10.Bond(uni10.BD_IN, D)
 degs = bdi.degeneracy()
 for qnum, dim in degs.iteritems():
  D_list.append(dim)

 D_dim=sum(D_list)
 
 
 row, colm=cal_rowcol(Teta)
 if (row<=colm and row<=D_dim):
  U,V,S=TU.setTruncation(Teta,row)
 elif (row<=colm and row>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 elif (row>colm and colm<=D_dim):
  U,V,S=TU.setTruncation(Teta,colm)
 elif (row>colm and colm>D_dim):
  U,V,S=TU.setTruncation(Teta,D_dim)
 
 U.setLabel([-10,-2,-3])
 V.setLabel([-3,-6,-30])
 S=Sqrt(S)
 S.setLabel([-3,3])
 r_up=U*S
 r_up.permute([-10,-2,3],3)
 r_up.setLabel([-100,3,10])
 S.setLabel([3,-3])
 l_up=V*S
 l_up.permute([3,-6,-30],3)
 l_up.setLabel([10,13,-300])
#########################################


 
# l_up=l_u*1.0
# r_up=r_u*1.0
 l_dp=l_up*1.0
 r_dp=r_up*1.0
 l_dp.transpose()
 l_dp.setLabel([-10,-13,-400])  #new
 r_dp.transpose()
 r_dp.setLabel([-200,-3,-10])



 Ham.setLabel([51,52,53,54])

 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([-3,-13,3,13])

 Ham.setLabel([-3,-13,3,13])
 
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()

 iden_h.setLabel([-3,-13,3,13])

 r_up_init=r_up*1.0
 r_dp_init=r_dp*1.0
 l_up_init=l_up*1.0
 l_dp_init=l_dp*1.0

 valf=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val( N_ten, l_u, r_u, r_d, l_d, l_up, r_up, l_dp, r_dp, Ham, iden_h, H)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12 or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   r_up=r_up_init*1.0
   r_dp=r_dp_init*1.0
   l_up=l_up_init*1.0
   l_dp=l_dp_init*1.0
   break
  else:
   r_up_init=r_up*1.0
   r_dp_init=r_dp*1.0
   l_up_init=l_up*1.0
   l_dp_init=l_dp*1.0

  r_up, r_dp=optimum_0( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)
  l_up, l_dp=optimum_1( N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham, iden_h, H)

 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-11 or abs(valf)>1.0e+11:
   print "warning_norm_in_optimization",  abs(valf), "count", count



 
 H_orig.setLabel([-3,-13,3,13])
 iden_h.setLabel([-3,-13,3,13])
 Norm_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*iden_h
 h_h=(((r_up*r_dp)*N_ten)*(l_up*l_dp))*H_orig
 #print "E_2r", Norm_h[0], h_h[0]/Norm_h[0]
 #print H_orig


 if Norm_h[0]> threshold[1] or Norm_h[0]<(threshold[0]):
   r_up,l_up=renormalize_subtensor(r_up,l_up,N_ten,iden_h,r_dp,l_dp,threshold,interval)


 r_up.setLabel([-100,3,9])
 l_up.setLabel([9,13,-300])


 PEPS_1=r_up*q
 PEPS_1.permute([6,7,3,9,10],3)


 PEPS_2=l_up*qq
 A=PEPS_2*1.0
 A.permute([9,10,13,14,15],5)
 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(2))
 Swap1.setLabel([-10,-13,10,13])
 A=Swap1*A
 A.permute([9,-10,-13,14,15],3)
 PEPS_2=A*1.0
 PEPS_2.permute([9,-10,-13,14,15],3)


 return PEPS_1, PEPS_2, h_h[0]/Norm_h[0]










######@profile
def norm_f_val_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left):

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 A=A*iden_h

 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*((Peps_2p*Swap11)*mps_boundry_left[Location+1])
 B=((B*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8))

 val1=A*B

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1*Swap6)
 A=A*Ham

 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*((Peps_2*Swap11)*mps_boundry_left[Location+1])
 B=((B*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8))

 val2=A*B

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 A=A*Ham

 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*((Peps_2p*Swap11)*mps_boundry_left[Location+1])
 B=((B*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))

 val3=A*B
 #print val1, val2, val3
 return val1[0]-val2[0]-val3[0]






######@profile
def norm_f_val_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left):

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1p*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 A=A*iden_h

 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]

 B=B*((Peps_2pd*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)

 val1=A*B

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1p*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 A=A*Ham

 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)

 val2=A*B

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 A=A*Ham

 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]

 B=B*((Peps_2pd*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2*Swap11)

 val3=A*B
 #print val1, val2, val3
 return val1[0]-val2[0]-val3[0]




######@profile
def optimum_0_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left):

 Env_r=mps_boundry_right[Location+1]*E_right
 Env_r=(Env_r*Swap12)*((Peps_2p*Swap11)*mps_boundry_left[Location+1])
 Env_r=((Env_r*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8))
 Env_r=Env_r*iden_h

 Env_r=(Env_r*mps_boundry_right[Location])
 Env_r=(Env_r*Swap4)*(Swap3)
 Env_r=Env_r*Swap5
 Env_r=Env_r*Swap6

 A=(mps_boundry_left[Location]*Swap1)
 A=A*(E_left*Swap2)
 Env_r=Env_r*A
 Env_r.permute([ 53,56,39, 54,57,6,35,37,55,34],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=mps_boundry_right[Location+1]*E_right
 Env_s=(Env_s*Swap12)*((Peps_2p*Swap11)*mps_boundry_left[Location+1])
 Env_s=((Env_s*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))
 Env_s=Env_s*Ham
 Env_s=((Env_s))*(Swap6)

 A=E_left*(mps_boundry_left[Location]*Swap1)
 A=A*((Peps_1d*Swap2)*Swap3)
 A=A*mps_boundry_right[Location]
 A=(A*Swap4)*Swap5
 Env_s=Env_s*A
 Env_s.permute([6,35,37,55,34],5)



 Env_s1=mps_boundry_right[Location+1]*E_right
 Env_s1=(Env_s1*Swap12)*((Peps_2*Swap11)*mps_boundry_left[Location+1])
 Env_s1=((Env_s1*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8))
 Env_s1=Env_s1*Ham
 Env_s1=((Env_s1)*(Swap4*Swap5))

 A=E_left*(mps_boundry_left[Location]*Swap1)
 A=A*(Peps_1*(Swap6))
 Env_s1=Env_s1*A
 Env_s1=Env_s1*mps_boundry_right[Location]

 Env_s1=Env_s1*((Swap2)*Swap3)

 Env_s1.permute([53,56,39,54,57],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([6,35,37,55,34,53,56,39,54,57])

 
 Peps_1p=A2_inv*Env_s
 Peps_1p.permute([6,35,37,55,34],5)
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.setLabel([53,56,39,54,57])

 return Peps_1p, Peps_1pd





######@profile
def optimum_0_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left):


 Env_r=E_right*mps_boundry_up[Location+1]
 Env_r=Swap7*Env_r
 Env_r=Env_r*mps_boundry_down[Location+1]

 Env_r=Env_r*((Peps_2pd*Swap6)*Swap5)
 Env_r=(Env_r*Swap4)*Swap12
 Env_r=Env_r*(Peps_2p*Swap11)
 Env_r=Env_r*iden_h

 Env_r=Env_r*(Swap10*mps_boundry_up[Location])

 Env_r=(Env_r)*(Swap9)
 Env_r=(Env_r)*(Swap8)
 Env_r=Env_r*Swap3


 A=E_left*mps_boundry_down[Location]
 A=A*Swap1
 A=A*Swap2


 Env_r=(Env_r)*(A)


 Env_r.permute([ 32,33,37, 35,34,6,31,58,59,10],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=E_right*mps_boundry_up[Location+1]
 Env_s=Swap7*Env_s
 Env_s=Env_s*mps_boundry_down[Location+1]

 Env_s=Env_s*((Peps_2d*Swap6)*Swap5)
 Env_s=(Env_s*Swap4)*Swap12
 Env_s=Env_s*(Peps_2p*Swap11)
 Env_s=Env_s*Ham


 Env_s=(Env_s*(Swap8*Swap9))*Swap10


 Env_s=(Env_s*((Peps_1d*Swap3)*Swap2))*mps_boundry_down[Location]
 Env_s=(Env_s)*(E_left*Swap1)

 Env_s=Env_s*mps_boundry_up[Location]

 Env_s.permute([6,31,58,59,10],5)

 Env_s1=E_right*mps_boundry_up[Location+1]
 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*mps_boundry_down[Location+1]

 Env_s1=Env_s1*((Peps_2pd*Swap6)*Swap5)
 Env_s1=(Env_s1*Swap4)*Swap12
 Env_s1=Env_s1*(Peps_2*Swap11)
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*(mps_boundry_up[Location]*Swap10)
 Env_s1=(Env_s1*(Swap8*Swap9))
 Env_s1=(Env_s1*(Swap3))
 Env_s1=Env_s1*Peps_1
 Env_s1=Env_s1*E_left
 Env_s1=Env_s1*((mps_boundry_down[Location]*Swap1)*Swap2)
 Env_s1.permute([32,33,37,35,34],0)
 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([6,31,58,59,10,32,33,37,35,34])

 Peps_1p=A2_inv*Env_s
 Peps_1p.permute([6,31,58,59,10],5)
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.setLabel([32,33,37,35,34])

 return Peps_1p, Peps_1pd


######@profile
def optimum_1_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left):

 Env_r=E_left*mps_boundry_up[Location]
 Env_r=Env_r*(Peps_1p*Swap1)
 Env_r=Env_r*mps_boundry_down[Location]
 Env_r=Env_r*((Swap3*Peps_1pd)*Swap2)
 Env_r=(Env_r*(Swap8*Swap9))*Swap10
 Env_r=Env_r*iden_h

 Env_r=Env_r*Swap4
 Env_r=Env_r*mps_boundry_down[Location+1]
 Env_r=Env_r*Swap5
 Env_r=Env_r*Swap12
 Env_r=Env_r*Swap11

 A=mps_boundry_up[Location+1]*E_right
 A=Swap7*A
 A=Swap6*A

 Env_r=(Env_r*A)

 Env_r.permute([43,42,47,45,44,54,52,51,53,15],5)
 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1

 Env_s=E_left*mps_boundry_up[Location]
 Env_s=Env_s*(Peps_1p*Swap1)
 Env_s=Env_s*mps_boundry_down[Location]
 Env_s=Env_s*((Swap3*Peps_1d)*Swap2)
 Env_s=(Env_s*(Swap8*Swap9))*Swap10
 Env_s=Env_s*Ham
 Env_s=((Swap11*Env_s)*Swap12)*Swap4
 Env_s=Env_s*mps_boundry_down[Location+1]
 Env_s=Env_s*((Swap5*Peps_2d)*Swap6)
 Env_s=Env_s*E_right
 Env_s=Swap7*Env_s

 Env_s=Env_s*mps_boundry_up[Location+1]

 Env_s.permute([ 54,52,51,53,15],5)

 Env_s1=E_left*mps_boundry_up[Location]
 Env_s1=Env_s1*(Peps_1*Swap1)
 Env_s1=Env_s1*mps_boundry_down[Location]
 Env_s1=Env_s1*((Swap3*Peps_1pd)*Swap2)
 Env_s1=(Env_s1*(Swap8*Swap9))*Swap10
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*((Peps_2*Swap11))
 
 Env_s1=(Env_s1*Swap12)*Swap4

 Env_s1=Env_s1*mps_boundry_up[Location+1]

 Env_s1=Env_s1*mps_boundry_down[Location+1]

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*E_right
 Env_s1=Env_s1*(Swap6*Swap5)
 Env_s1.permute([ 43,42,47,45,44],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([54,52,51,53,15,43,42,47,45,44])
 
 
 Peps_2p=A2_inv*Env_s
 Peps_2p.permute([54,52,51,53,15],5)
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.setLabel([43,42,47,45,44])

 return Peps_2p, Peps_2pd



######@profile
def optimum_1_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left):

 Env_r=E_left*mps_boundry_left[Location]
 Env_r=Env_r*(Swap1*((Peps_1pd*Swap3)*Swap2))
 Env_r=Env_r*mps_boundry_right[Location]
 Env_r=((Env_r*Swap4)*Swap5)*(Peps_1p*Swap6)
 Env_r=Env_r*iden_h
 Env_r=((Env_r*Swap10)*Swap11)
 Env_r=((Env_r*Swap9)*Swap8)
 Env_r=Env_r*mps_boundry_left[Location+1]
 A=mps_boundry_right[Location+1]*Swap12
 A=A*Swap7
 A=A*E_right
 Env_r=Env_r*A

 Env_r.permute([ 48,58,46, 49,57,11,31,32,33,15],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1

 Env_s=E_left*mps_boundry_left[Location]
 Env_s=Env_s*(Swap1*((Peps_1pd*Swap3)*Swap2))
 Env_s=Env_s*mps_boundry_right[Location]
 Env_s=((Env_s*Swap4)*Swap5)*(Peps_1*Swap6)
 Env_s=Env_s*Ham
 Env_s=((Env_s*Swap10)*Swap11)
 Env_s=(Env_s*mps_boundry_left[Location+1])*Peps_2
 Env_s=Env_s*E_right
 Env_s=Env_s*(mps_boundry_right[Location+1]*Swap12)
 Env_s=Env_s*Swap7
 Env_s=(Env_s)*(Swap9*Swap8)
 Env_s.permute([48,58,46,49,57],5)



 Env_s1=E_left*mps_boundry_left[Location]
 Env_s1=Env_s1*(Swap1*((Peps_1d*Swap3)*Swap2))
 Env_s1=Env_s1*mps_boundry_right[Location]
 Env_s1=((Env_s1*Swap4)*Swap5)*(Peps_1p*Swap6)
 Env_s1=Env_s1*Ham
 A=((mps_boundry_right[Location+1]*Swap12)*(Swap7*(Swap8*Peps_2d)))*E_right
 Env_s1=A*(Env_s1*Swap9)
 Env_s1=((Env_s1*Swap10)*Swap11)
 Env_s1=(Env_s1*mps_boundry_left[Location+1])
 Env_s1.permute([11,31,32,33,15],0)

 Env_s1.transpose()
 Env_s=Env_s+Env_s1


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([11,31,32,33,15,48,58,46, 49,57])
 
 
 Peps_2p=A2_inv*Env_s
 Peps_2p.permute([11,31,32,33,15],5)
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.setLabel([48,58,46, 49,57])

 return Peps_2p, Peps_2pd






######@profile
def double_local_update_inv( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, H_orig, Ham ,D, threshold, interval,Sys):


 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=copy.copy(PEPS_listten[Location])
 Peps_2=copy.copy(PEPS_listten[Location+1])

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 A=Peps_1*1.0
 Peps_1.setLabel([6,35,37,55,34])
 Peps_1.permute([6,35,37,55,34],5)
 A.setLabel([6,35,37,55,34])
 A.permute([6,35,37,55,34],5)

 Swap1=fermionicOPT(Sys,A.bond(1), A.bond(0))
 Swap1.setLabel([70,-60,35,36])
 Swap2=Swap1*1.0
 Swap2.setLabel([-7,36,56,53])
 Swap3=fermionicOPT(Sys,A.bond(4), A.bond(3))
 Swap3.setLabel([57,54,40,-9])
 Swap4=Swap3*1.0
 Swap4.setLabel([40,58,41,90])
 Swap6=fermionicOPT(Sys,A.bond(2), A.bond(3))
 Swap6.setLabel([37,55,42,38])
 Swap5=fermionicOPT(Sys,A.bond(3), A.bond(2))
 Swap5.setLabel([38,39,58,44])


#####################################################

 A1=Peps_2*1.0
 A1.setLabel([11,31,32,33,15])
 A1.permute([11,31,32,33,15],5)
 Peps_2.setLabel([11,31,32,33,15])
 Peps_2.permute([11,31,32,33,15],5)

 Swap7=fermionicOPT(Sys,A1.bond(3), A1.bond(4))
 Swap7.setLabel([49,57,-14,50])
 Swap8=fermionicOPT(Sys,A1.bond(0), A1.bond(1))
 Swap8.setLabel([48,58,47,41])
 Swap9=fermionicOPT(Sys,A1.bond(0), A1.bond(2))
 Swap9.setLabel([51,45,47,46])
 Swap10=fermionicOPT(Sys,A1.bond(0), A1.bond(2))
 Swap10.setLabel([52,32,51,43])
 Swap11=fermionicOPT(Sys,A1.bond(0), A1.bond(1))
 Swap11.setLabel([-110,34,52,31])
 Swap12=fermionicOPT(Sys,A1.bond(3), A1.bond(4))
 Swap12.setLabel([33,50,140,-150])



######################################################################
 Peps_1d=Peps_1*1.0
 Peps_1d.setLabel([53,56,39,54,57])
 Peps_2d=Peps_2*1.0
 Peps_2d.setLabel([48,58,46,49,57])

 mps_boundry_left[Location].setLabel([16,6,-60,21])
 mps_boundry_left[Location+1].setLabel([21,11,-110,25])

 mps_boundry_right[Location].setLabel([17,90,-9,20])
 mps_boundry_right[Location+1].setLabel([20,140,-14,24])


######################################################

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1*Swap6)


 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*(Peps_2*mps_boundry_left[Location+1])
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2d*Swap7)*Swap8))


 N_ten=A*B
 #print N_ten
 #N_ten.permute([200,-400,-100,300],2)
 #N_ten.setLabel([-200,-400,-100,-300])

 H_orig.setLabel([44,45,42,43])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_1c", Norm_h[0], h_h[0]/Norm_h[0]


 #N_ten=N_Positiv(N_ten)

######################################################



# l_up=l_u*1.0
# r_up=r_u*1.0
 Peps_1p=Peps_1*1.0
 Peps_2p=Peps_2*1.0
 Peps_1pd=Peps_1d*1.0
 Peps_2pd=Peps_2d*1.0

 Peps_2d.setLabel([48,58,46,49,57])
 Peps_1d.setLabel([53,56,39,54,57])


 Ham.setLabel([51,52,53,54])
 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([44,45,42,43])
 Ham.setLabel([44,45,42,43])
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])

 Peps_1p_init=Peps_1p*1.0
 Peps_1pd_init=Peps_1pd*1.0
 Peps_2p_init=Peps_2p*1.0
 Peps_2pd_init=Peps_2pd*1.0

 valf=norm_f_val_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   Peps_1p=Peps_1p_init*1.0
   Peps_1pd=Peps_1pd_init*1.0
   Peps_2p=Peps_2p_init*1.0
   Peps_2pd=Peps_2pd_init*1.0
   break
  else:
   Peps_1p_init=Peps_1p*1.0
   Peps_1pd_init=Peps_1pd*1.0
   Peps_2p_init=Peps_2p*1.0
   Peps_2pd_init=Peps_2pd*1.0


  Peps_1p, Peps_1pd=optimum_0_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)

  Peps_2p, Peps_2pd=optimum_1_inv( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)



 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count

 if abs(valf)<1.0e-12 or abs(valf)>1.0e+12:
  print "warning_norm_in_optimization",  abs(valf), "count", count







 Peps_1d=Peps_1*1.0
 Peps_1d.setLabel([53,56,39,54,57])
 Peps_2d=Peps_2*1.0
 Peps_2d.setLabel([48,58,46,49,57])

 mps_boundry_left[Location].setLabel([16,6,-60,21])
 mps_boundry_left[Location+1].setLabel([21,11,-110,25])

 mps_boundry_right[Location].setLabel([17,90,-9,20])
 mps_boundry_right[Location+1].setLabel([20,140,-14,24])


######################################################

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)


 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))

 N_ten=A*B

 H_orig.setLabel([44,45,42,43])
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_2c", Norm_h[0], h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<(threshold[0]):
  Peps_1p, Peps_2p=renormalize_tensor_c(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,threshold,interval)


 return Peps_1p, Peps_2p, h_h[0]/Norm_h[0]



######@profile
def renormalize_tensor_c(Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,threshold,interval):

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "norm_val", norm_val
 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold[1]:
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*mps_boundry_left[Location]
   A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
   A=A*mps_boundry_right[Location]
   A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
   B=mps_boundry_right[Location+1]*E_right
   B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
   B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
   N_ten=A*B
   iden_h.setLabel([44,45,42,43])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<(threshold[0]):
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*mps_boundry_left[Location]
   A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
   A=A*mps_boundry_right[Location]
   A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
   B=mps_boundry_right[Location+1]*E_right
   B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
   B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
   N_ten=A*B
   iden_h.setLabel([44,45,42,43])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])

 A=E_left*mps_boundry_left[Location]
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 B=mps_boundry_right[Location+1]*E_right
 B=(B*Swap12)*(Peps_2p*mps_boundry_left[Location+1])
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "Fixednorm_local", abs(norm_val)
 return Peps_1p, Peps_2p


######@profile
def  Update_twotensor_double_row_inv( PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, Ham, H_orig, D, threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])



 Peps_1.setLabel([6,31,58,59,10])
 Peps_1d=Peps_1*1.0
 Peps_1d.setLabel([32,33,37,35,34])
 Peps_2.setLabel([54,52,51,53,15])
 Peps_2d=Peps_2*1.0
 Peps_2d.setLabel([43,42,47,45,44])


 Swap1=fermionicOPT(Sys,Peps_1.bond(1), Peps_1.bond(0))
 Swap1.setLabel([31,17,22,30])
 Swap2=Swap1*1.0
 Swap2.setLabel([33,32,-7,30])
 Swap3=fermionicOPT(Sys,Peps_1.bond(4), Peps_1.bond(3))
 Swap3.setLabel([34,35,38,36])
 Swap4=fermionicOPT(Sys,Peps_2.bond(1), Peps_2.bond(0))
 Swap4.setLabel([48,36,23,41])
 Swap5=Swap4*1.0
 Swap5.setLabel([42,43,-10,41])
 Swap6=fermionicOPT(Sys,Peps_2.bond(4), Peps_2.bond(3))
 Swap6.setLabel([44,45,46,-14])

 Swap7=fermionicOPT(Sys,Peps_2.bond(4), Peps_2.bond(3))
 Swap7.setLabel([46,53,24,20])

 Swap8=fermionicOPT(Sys,Peps_1.bond(2), Peps_1.bond(4))
 Swap8.setLabel([37,38,40,39])

 Swap9=Swap8*1.0
 Swap9.setLabel([58,39,56,57])

 Swap10=Swap3*1.0
 Swap10.setLabel([57,59,25,54])

 Swap11=fermionicOPT(Sys,Peps_2.bond(2), Peps_2.bond(1))
 Swap11.setLabel([55,50,51,52])

 Swap12=Swap11*1.0
 Swap12.setLabel([47,48,49,50])



######################################################################

 mps_boundry_up[Location].setLabel([16,10,25,27])
 mps_boundry_up[Location+1].setLabel([27,15,24,21])

 mps_boundry_down[Location].setLabel([18,22,-7,322])
 mps_boundry_down[Location+1].setLabel([322,23,-10,19])



######################################################

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10

 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2*Swap11)

 N_ten=A*B
 
 N_ten.permute([40,49,56,55],2)



 H_orig.setLabel([40,49,56,55])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_1r", Norm_h[0], h_h[0]/Norm_h[0]



 Peps_1p=Peps_1*1.0
 Peps_2p=Peps_2*1.0
 
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.setLabel([32,33,37,35,34])
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.setLabel([43,42,47,45,44])



 Ham.setLabel([51,52,53,54])

 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([40,49,56,55])

 Ham.setLabel([40,49,56,55])
 
 iden_h=Ham*1.0
 iden_h.setLabel([40,49,56,55])
 iden_h.identity()

 iden_h.setLabel([40,49,56,55])

 Peps_1p_init=Peps_1p*1.0
 Peps_1pd_init=Peps_1pd*1.0
 Peps_2p_init=Peps_2p*1.0
 Peps_2pd_init=Peps_2pd*1.0

 valf=norm_f_val_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(14):
  count=count+1
  E_2=E_1*1.0
  val=norm_f_val_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   Peps_1p=Peps_1p_init*1.0
   Peps_1pd=Peps_1pd_init*1.0
   Peps_2p=Peps_2p_init*1.0
   Peps_2pd=Peps_2pd_init*1.0
   break
  else:
   Peps_1p_init=Peps_1p*1.0
   Peps_1pd_init=Peps_1pd*1.0
   Peps_2p_init=Peps_2p*1.0
   Peps_2pd_init=Peps_2pd*1.0


  Peps_1p, Peps_1pd=optimum_0_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left)

  Peps_2p, Peps_2pd=optimum_1_inv_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left)








 if abs(valf)<1.0e-11 or abs(valf)>1.0e+11:
   print "warning_norm_in_optimization",  abs(valf), "count", count



 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1p*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]
 B=B*((Peps_2pd*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)
 N_ten=A*B
 #N_ten.permute([40,49,56,55],2)


 H_orig.setLabel([40,49,56,55])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_2r", Norm_h[0], h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<(threshold[0]):
  Peps_1p, Peps_2p=renormalize_tensor_r(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,threshold,interval)


 return Peps_1p, Peps_2p, h_h[0]/Norm_h[0]




######@profile
def renormalize_tensor_r(Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,threshold,interval):

 A=E_left*mps_boundry_up[Location]
 A=A*(Peps_1p*Swap1)
 A=A*mps_boundry_down[Location]
 A=A*((Swap3*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 B=E_right*mps_boundry_up[Location+1]
 B=Swap7*B
 B=B*mps_boundry_down[Location+1]
 B=B*((Peps_2pd*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)
 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "norm_val", norm_val
 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold[1]:
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*mps_boundry_up[Location]
   A=A*(Peps_1p*Swap1)
   A=A*mps_boundry_down[Location]
   A=A*((Swap3*Peps_1pd)*Swap2)
   A=(A*(Swap8*Swap9))*Swap10
   B=E_right*mps_boundry_up[Location+1]
   B=Swap7*B
   B=B*mps_boundry_down[Location+1]
   B=B*((Peps_2pd*Swap6)*Swap5)
   B=(B*Swap4)*Swap12
   B=B*(Peps_2p*Swap11)
   N_ten=A*B
   iden_h.setLabel([40,49,56,55])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<(threshold[0]):
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*mps_boundry_up[Location]
   A=A*(Peps_1p*Swap1)
   A=A*mps_boundry_down[Location]
   A=A*((Swap3*Peps_1pd)*Swap2)
   A=(A*(Swap8*Swap9))*Swap10
   B=E_right*mps_boundry_up[Location+1]
   B=Swap7*B
   B=B*mps_boundry_down[Location+1]
   B=B*((Peps_2pd*Swap6)*Swap5)
   B=(B*Swap4)*Swap12
   B=B*(Peps_2p*Swap11)
   N_ten=A*B
   iden_h.setLabel([40,49,56,55])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])

  A=E_left*mps_boundry_up[Location]
  A=A*(Peps_1p*Swap1)
  A=A*mps_boundry_down[Location]
  A=A*((Swap3*Peps_1pd)*Swap2)
  A=(A*(Swap8*Swap9))*Swap10
  B=E_right*mps_boundry_up[Location+1]
  B=Swap7*B
  B=B*mps_boundry_down[Location+1]
  B=B*((Peps_2pd*Swap6)*Swap5)
  B=(B*Swap4)*Swap12
  B=B*(Peps_2p*Swap11)
  N_ten=A*B
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "Fixednorm_local", abs(norm_val)
 return Peps_1p, Peps_2p





def optimum_0_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p):

 Env_r=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_r=(Env_r*Swap12)*((Peps_2p*Swap11p)*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))


 Env_r=((Env_r*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))
 Env_r=Env_r*iden_h

 Env_r=(Env_r*(mps_boundry_right[2*Location+1]*mps_boundry_right[2*Location]))

 Env_r=(Env_r*Swap4p)*(Swap3p)
 Env_r=Env_r*Swap5
 Env_r=Env_r*Swap6

 A=(mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1

 A=A*(E_left*Swap2)
 Env_r=Env_r*A
 Env_r.permute([ 53,56,39, 54,57,6,35,37,55,34],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_s=(Env_s*Swap12)*(Peps_2p*(Swap11p*(mps_boundry_left[2*Location+3]*mps_boundry_left[2*Location+2])))

 Env_s=((Env_s*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))
 Env_s=Env_s*Ham
 Env_s=((Env_s))*(Swap6)

 A=E_left*((mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1)

 A=A*((Peps_1d*Swap2)*Swap3)
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_right[2*Location+1]
 A=(A*Swap4)*Swap5
 Env_s=Env_s*A
 Env_s.permute([6,35,37,55,34],5)



 Env_s1=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_s1=(Env_s1*Swap12)*((Peps_2)*((mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])*Swap11))

 Env_s1=((Env_s1*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))
 Env_s1=Env_s1*Ham
 Env_s1=((Env_s1)*(Swap4p*Swap5))

 A=E_left*((mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1)

 A=A*(Peps_1*(Swap6))
 Env_s1=Env_s1*A
 Env_s1=Env_s1*mps_boundry_right[2*Location]
 Env_s1=Env_s1*mps_boundry_right[2*Location+1]

 Env_s1=Env_s1*((Swap2)*Swap3p)

 Env_s1.permute([53,56,39,54,57],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([6,35,37,55,34,53,56,39,54,57])

 
 Peps_1p=A2_inv*Env_s
 Peps_1p.permute([6,35,37,55,34],3)
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.transpose()
 Peps_1pd.setLabel([54,57,53,56,39])

 return Peps_1p, Peps_1pd



def optimum_1_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p):

 Env_r=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])

 Env_r=Env_r*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 Env_r=Env_r*(mps_boundry_right[2*Location]*mps_boundry_right[2*Location+1])


 Env_r=((Env_r*Swap4p)*Swap5)*(Peps_1p*Swap6)
 Env_r=Env_r*iden_h
 Env_r=((Env_r*Swap10)*Swap11p)
 Env_r=((Env_r*Swap9)*Swap8p)
 Env_r=Env_r*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])

 A=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12
 A=A*Swap7
 A=A*E_right

 Env_r=Env_r*A
 Env_r.permute([ 48,58,46, 49,57,11,31,32,33,15],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1

 Env_s=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])

 Env_s=Env_s*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 Env_s=Env_s*(mps_boundry_right[2*Location]*mps_boundry_right[2*Location+1])

 Env_s=((Env_s*Swap4p)*Swap5)*(Peps_1*Swap6)
 Env_s=Env_s*Ham
 Env_s=((Env_s*Swap10)*Swap11)
 Env_s=(Env_s*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))*Peps_2

 Env_s=Env_s*E_right
 Env_s=Env_s*((mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12)

 Env_s=Env_s*Swap7
 Env_s=(Env_s)*(Swap9*Swap8p)
 Env_s.permute([48,58,46,49,57],5)



 Env_s1=E_left*mps_boundry_left[2*Location]
 Env_s1=Env_s1*mps_boundry_left[2*Location+1]
 Env_s1=Env_s1*(Swap1*((Peps_1d*Swap3)*Swap2))
 Env_s1=Env_s1*mps_boundry_right[2*Location]
 Env_s1=Env_s1*mps_boundry_right[2*Location+1]

 Env_s1=((Env_s1*Swap4)*Swap5)*(Peps_1p*Swap6)
 Env_s1=Env_s1*Ham
 A=(((mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12)*(Swap7*(Swap8*Peps_2d)))*E_right

 Env_s1=A*(Env_s1*Swap9)
 Env_s1=((Env_s1*Swap10)*Swap11p)
 Env_s1=(Env_s1*mps_boundry_left[2*Location+2])
 Env_s1=(Env_s1*mps_boundry_left[2*Location+3])

 Env_s1.permute([11,31,32,33,15],0)

 Env_s1.transpose()
 Env_s=Env_s+Env_s1


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([11,31,32,33,15,48,58,46, 49,57])
 
 
 Peps_2p=A2_inv*Env_s
 Peps_2p.permute([11,31,32,33,15],3)
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.transpose()
 Peps_2pd.setLabel([ 49,57,48,58,46])

 return Peps_2p, Peps_2pd





#@profile
def norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p):

#  A=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])
#  A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
#  A=A*mps_boundry_right[2*Location]
#  A=A*mps_boundry_right[2*Location+1]
# 
#  A=((A*Swap4)*Swap5)*(Peps_1*Swap6)
#  A=A*H
# 
#  B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right
# 
#  B=(B*Swap12)*((Peps_2*Swap11)*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))
#  B=((B*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))
# 
#  val=A*B

###############
 A=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])
 A=A*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_right[2*Location+1]

 A=((A*Swap4p)*Swap5)*(Peps_1p*Swap6)
 A=A*iden_h

 B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 B=(B*Swap12)*((Peps_2p*Swap11p)*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))
 B=((B*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))

 #print A.printDiagram(), B.printDiagram()
 val1=A*B

 A=E_left*mps_boundry_left[2*Location]
 A=A*mps_boundry_left[2*Location+1]
 A=A*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_right[2*Location+1]
 A=((A*Swap4p)*Swap5)*(Peps_1*Swap6)
 A=A*Ham

 B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right


 #print Swap11.printDiagram(), Swap11p.printDiagram()
 B=B*Swap12
 B1=(Peps_2*Swap11)
 B=B*B1
 B=B*(mps_boundry_left[2*Location+3]*mps_boundry_left[2*Location+2])
 B=((B*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))

 val2=A*B

 A=E_left*mps_boundry_left[2*Location]
 A=A*mps_boundry_left[2*Location+1]
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*(mps_boundry_right[2*Location]*mps_boundry_right[2*Location+1])

 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 A=A*Ham

 B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 B=(B*Swap12)*((Peps_2p*Swap11p)*(mps_boundry_left[2*Location+3]*mps_boundry_left[2*Location+2]))
 B=((B*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))

 val3=A*B
 #print val, val1, val2, val3
 #if val1[0]<0: print  "Oh, norm<0"
 #return val1[0]-2.0*val2[0]
 #return val[0]+val1[0]-val2[0]-val3[0]
 return val1[0]-val2[0]-val3[0], val1[0]





#@profile
def Obtain_grad_two_coulmn(Peps_1, Peps_1d, Peps_2, Peps_2d, Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p):


 Env_s=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_s=(Env_s*Swap12)*(Peps_2p*(Swap11p*(mps_boundry_left[2*Location+3]*mps_boundry_left[2*Location+2])))

 Env_s=((Env_s*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))
 Env_s=Env_s*iden_h
 Env_s=((Env_s))*(Swap6)

 A=E_left*((mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1)

 A=A*((Peps_1pd*Swap2)*Swap3p)
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_right[2*Location+1]
 A=(A*Swap4p)*Swap5
 Env_s=Env_s*A
 Env_s.permute([6,35,37,55,34],5)


 Env_s1=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_s1=(Env_s1*Swap12)*((Peps_2p)*((mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])*Swap11p))

 Env_s1=((Env_s1*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))
 Env_s1=Env_s1*iden_h
 Env_s1=((Env_s1)*(Swap4p*Swap5))

 A=E_left*((mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1)

 A=A*(Peps_1p*(Swap6))
 Env_s1=Env_s1*A
 Env_s1=Env_s1*mps_boundry_right[2*Location]
 Env_s1=Env_s1*mps_boundry_right[2*Location+1]

 Env_s1=Env_s1*((Swap2)*Swap3p)

 Env_s1.permute([53,56,39,54,57],0)

 Env_s1.transpose()
 Env_s=Env_s+Env_s1
 D_a=Env_s*1.0


 Env_s=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_s=(Env_s*Swap12)*(Peps_2p*(Swap11p*(mps_boundry_left[2*Location+3]*mps_boundry_left[2*Location+2])))

 Env_s=((Env_s*Swap10))*(Swap9*((Peps_2d*Swap7)*Swap8))
 Env_s=Env_s*Ham
 Env_s=((Env_s))*(Swap6)

 A=E_left*((mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1)

 A=A*((Peps_1d*Swap2)*Swap3)
 A=A*mps_boundry_right[2*Location]
 A=A*mps_boundry_right[2*Location+1]
 A=(A*Swap4)*Swap5
 Env_s=Env_s*A
 Env_s.permute([6,35,37,55,34],5)


 Env_s1=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right

 Env_s1=(Env_s1*Swap12)*((Peps_2)*((mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])*Swap11))

 Env_s1=((Env_s1*Swap10))*(Swap9*((Peps_2pd*Swap7)*Swap8p))
 Env_s1=Env_s1*Ham
 Env_s1=((Env_s1)*(Swap4p*Swap5))

 A=E_left*((mps_boundry_left[2*Location+1]*mps_boundry_left[2*Location])*Swap1)

 A=A*(Peps_1*(Swap6))
 Env_s1=Env_s1*A
 Env_s1=Env_s1*mps_boundry_right[2*Location]
 Env_s1=Env_s1*mps_boundry_right[2*Location+1]

 Env_s1=Env_s1*((Swap2)*Swap3p)

 Env_s1.permute([53,56,39,54,57],0)

 Env_s1.transpose()
 Env_s=Env_s+Env_s1
 D_a=D_a+Env_s*(-1.0)


 Env_s=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])

 Env_s=Env_s*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 Env_s=Env_s*(mps_boundry_right[2*Location]*mps_boundry_right[2*Location+1])

 Env_s=((Env_s*Swap4p)*Swap5)*(Peps_1p*Swap6)
 Env_s=Env_s*iden_h
 Env_s=((Env_s*Swap10)*Swap11p)
 Env_s=(Env_s*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))*Peps_2p

 Env_s=Env_s*E_right
 Env_s=Env_s*((mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12)

 Env_s=Env_s*Swap7
 Env_s=(Env_s)*(Swap9*Swap8p)
 Env_s.permute([48,58,46,49,57],5)

 Env_s1=E_left*mps_boundry_left[2*Location]
 Env_s1=Env_s1*mps_boundry_left[2*Location+1]
 Env_s1=Env_s1*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 Env_s1=Env_s1*mps_boundry_right[2*Location]
 Env_s1=Env_s1*mps_boundry_right[2*Location+1]

 Env_s1=((Env_s1*Swap4p)*Swap5)*(Peps_1p*Swap6)
 Env_s1=Env_s1*iden_h
 A=(((mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12)*(Swap7*(Swap8p*Peps_2pd)))*E_right

 Env_s1=A*(Env_s1*Swap9)
 Env_s1=((Env_s1*Swap10)*Swap11p)
 Env_s1=(Env_s1*mps_boundry_left[2*Location+2])
 Env_s1=(Env_s1*mps_boundry_left[2*Location+3])

 Env_s1.permute([11,31,32,33,15],0)

 Env_s1.transpose()
 Env_s=Env_s+Env_s1
 D_b=Env_s*1.0

 Env_s=E_left*(mps_boundry_left[2*Location]*mps_boundry_left[2*Location+1])

 Env_s=Env_s*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 Env_s=Env_s*(mps_boundry_right[2*Location]*mps_boundry_right[2*Location+1])

 Env_s=((Env_s*Swap4p)*Swap5)*(Peps_1*Swap6)
 Env_s=Env_s*Ham
 Env_s=((Env_s*Swap10)*Swap11)
 Env_s=(Env_s*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))*Peps_2

 Env_s=Env_s*E_right
 Env_s=Env_s*((mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12)

 Env_s=Env_s*Swap7
 Env_s=(Env_s)*(Swap9*Swap8p)
 Env_s.permute([48,58,46,49,57],5)


 Env_s1=E_left*mps_boundry_left[2*Location]
 Env_s1=Env_s1*mps_boundry_left[2*Location+1]
 Env_s1=Env_s1*(Swap1*((Peps_1d*Swap3)*Swap2))
 Env_s1=Env_s1*mps_boundry_right[2*Location]
 Env_s1=Env_s1*mps_boundry_right[2*Location+1]

 Env_s1=((Env_s1*Swap4)*Swap5)*(Peps_1p*Swap6)
 Env_s1=Env_s1*Ham
 A=(((mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*Swap12)*(Swap7*(Swap8*Peps_2d)))*E_right

 Env_s1=A*(Env_s1*Swap9)
 Env_s1=((Env_s1*Swap10)*Swap11p)
 Env_s1=(Env_s1*mps_boundry_left[2*Location+2])
 Env_s1=(Env_s1*mps_boundry_left[2*Location+3])

 Env_s1.permute([11,31,32,33,15],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s
 D_b=D_b+Env_s*(-1.0)


 D_a.transpose()
 #D_b.transpose()
 D_a.permute([6,35,37,55,34],5)
 #D_b.permute([48,58,46,49,57],5)
 
 return D_a, D_b


#@profile
def Do_optimization_Grad_col( Peps_1, Peps_1d, Peps_2, Peps_2d, Peps_1p, Peps_1pd, Peps_2p, Peps_2pd, Ham_ham, iden_h, H_ham, Swap1, Swap2, Swap3, Swap4, Swap5, Swap6, Swap7, Swap8, Swap9, Swap10, Swap11, Swap12, mps_boundry_right, mps_boundry_left, Location, E_right, E_left,Sys,H_orig,Swap3p,Swap8p,Swap11p,Swap4p):
  Opt_method="ST"
  #Opt_method="grad"

#    Peps_1p.permute([6,35,37,55,34],5)
#    Peps_2p.permute([11,31,32,33,15],5)
#    Peps_1pd.permute([53,56,39,54,57],5)
#    Peps_2pd.permute([48,58,46,49,57],5)


  Peps_1p_um=copy.copy(Peps_1p)
  Peps_2p_um=copy.copy(Peps_2p)
  Peps_1p_umd=Peps_1p_um*1.0
  Peps_2p_umd=Peps_2p_um*1.0
  Peps_1p_umd.transpose()
  Peps_2p_umd.transpose()
  Peps_1p_umd.setLabel([53,56,39,54,57])
  Peps_2p_umd.setLabel([48,58,46,49,57])

  time_val=0
  Es, E_positive=norm_f_val_inv_single( Peps_1, Peps_1d, Peps_2, Peps_2d, Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)

  Ef=0
  E2_val=0
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=1.e+18

  EE_previous=1.e+18
  EEf=0
  EE2_val=0

  count=0
  D_list=[0]*4
  H_list=[0]*4
  H_a=0; H_b=0;H_c=0;H_d=0;
  for i in xrange(Sys[7]):
   count+=1
   t0=time.time()

   E1_val,E_positive=norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)

#   E_energy, Norm_val=Energy_cal_grad(Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,iden_h,H_orig,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)



   #if i%200==0:
   #print 'E1=', E1_val, abs((E_previous-E1_val)/E1_val), i, count, time_val
   #print 'E_energy=', E_energy, abs((EE_previous-E_energy)/E_energy), i, count, time_val


   D_a, D_b=Obtain_grad_two_coulmn(Peps_1, Peps_1d, Peps_2, Peps_2d, Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)
   D_a=(-1.0)*D_a
   D_b=(-1.0)*D_b

   if i is 0:
    H_a=D_a
    H_b=D_b
   else:
    Z_a=D_a+(-1.0)*D_list[0]
    Z_b=D_b+(-1.0)*D_list[1]
    A=Z_a*D_a
    B=Z_b*D_b
    A1=D_list[0]*D_list[0]
    A2=D_list[1]*D_list[1]
    Gamma_grad=(A[0]+B[0]) / (A1[0]+A2[0])
    if Opt_method is 'ST':Gamma_grad=0;
    H_a=D_a+(Gamma_grad)*H_list[0]
    H_b=D_b+(Gamma_grad)*H_list[1]

#    A=D_a*D_list[0]
#    B=D_b*D_list[1]
#    C=D_c*D_list[2]
#    D=D_d*D_list[3]
#    check=A[0]+B[0]+C[0]+D[0] 
#    print "check", check 

   D_list[0]=copy.copy(D_a)
   D_list[1]=copy.copy(D_b)

   H_list[0]=copy.copy(H_a)
   H_list[1]=copy.copy(H_b)



   A=D_a*H_a
   B=D_b*H_b
   
   Norm_Z=A[0]+B[0]
   
   #print 'Norm', Norm_Z
   if (E1_val<E_previous) or (i is 0):
    if (abs(E1_val) > 1.0e-10):
     if abs((E_previous-E1_val)/E1_val) < 1.0e-10:
      print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)/E1_val), i
      break
     else: 
      if abs((E_previous-E1_val)) < 1.0e-15:
       print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)), i
       break


#   if (E_energy<EE_previous) or (i is 0):
#    if (abs(E_energy) > 1.0e-10):
#     if abs((EE_previous-E_energy)/E_energy) < 1.0e-10:
#      print 'Differnance Satisfied!', EE_previous, E_energy, abs((EE_previous-E_energy)/E_energy), i
#      break
#     else: 
#      if abs((EE_previous-E1_val)) < 1.0e-15:
#       print 'Differnance Satisfied!', EE_previous, E_energy, abs((EE_previous-E_energy)), i
#       break



   #print E1_val>E_previous  


   if  E_positive<0  :
    print "norm negetive", E_positive
    Peps_1p=copy.copy(Peps_1p_um) 
    Peps_2p=copy.copy(Peps_2p_um) 
    Peps_1pd=copy.copy(Peps_1p_umd) 
    Peps_2pd=copy.copy(Peps_2p_umd) 
    break



#   if E_energy>EE_previous  :
#    print "break, energy_c", E_energy, EE_previous, i
#    Peps_1p=copy.copy(Peps_1p_um) 
#    Peps_2p=copy.copy(Peps_2p_um) 
#    Peps_1pd=copy.copy(Peps_1p_umd) 
#    Peps_2pd=copy.copy(Peps_2p_umd) 
#    break
#   else:
#    Peps_1p_um=copy.copy(Peps_1p) 
#    Peps_2p_um=copy.copy(Peps_2p) 
#    Peps_1p_umd=Peps_1p_um*1.0
#    Peps_2p_umd=Peps_2p_um*1.0
#    Peps_1p_umd.setLabel([53,56,39,54,57])
#    Peps_2p_umd.setLabel([48,58,46,49,57])




   if E1_val>E_previous  :
    print "break, not satisfied", E1_val, E_previous
    Peps_1p=copy.copy(Peps_1p_um) 
    Peps_2p=copy.copy(Peps_2p_um) 
    Peps_1pd=copy.copy(Peps_1p_umd) 
    Peps_2pd=copy.copy(Peps_2p_umd) 
    break
   else:
    Peps_1p_um=copy.copy(Peps_1p) 
    Peps_2p_um=copy.copy(Peps_2p) 
    Peps_1p_umd=Peps_1p_um*1.0
    Peps_2p_umd=Peps_2p_um*1.0
    Peps_1p_umd.transpose()
    Peps_2p_umd.transpose()
    Peps_1p_umd.setLabel([53,56,39,54,57])
    Peps_2p_umd.setLabel([48,58,46,49,57])

   E_previous=E1_val
   #EE_previous=E_energy

   break_val="off"
   if abs(Norm_Z) < 1.0e-12:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   if (i%15)==0: 
    Gamma=1
   else: 
    if Gamma >= 1: Gamma*=(1.00/100)
    if Gamma < 1: Gamma*=100
   #Gamma=1
   #print "Gamma", Gamma
   while Break_loop is 1:
    count+=1
    #print Peps_1p.printDiagram(), H_a.printDiagram()
    #print Peps_2p.printDiagram(), H_b.printDiagram()

    a_ut=Peps_1p+(2.00)*Gamma*H_a
    b_ut=Peps_2p+(2.00)*Gamma*H_b
    a_utd=a_ut*1.0
    b_utd=b_ut*1.0
    a_utd.transpose()
    b_utd.transpose()
    a_utd.setLabel([53,56,39,54,57])
    b_utd.setLabel([48,58,46,49,57])
    
    #E2_val=cost_f(A_f, a_ut, b_ut, c_ut, d_ut,Swap1,Swap2,Swap3)
    E2_val,E_positive=norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,a_ut,a_utd,b_ut,b_utd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)
    if  E_positive<0: 
     break_val="on"
     break




    if abs((0.5)*Norm_Z*Gamma) > 1.0e+12 or  abs(Gamma)>1.0e+12 :
     print "break1", E1_val, abs((0.5)*Norm_Z*Gamma), E2_val, Gamma
     Gamma=1
     break
    if E1_val-E2_val >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0

   Break_loop=1
   while Break_loop is 1:
    count+=1
    a_ut=Peps_1p+(1.00)*Gamma*H_a
    b_ut=Peps_2p+(1.00)*Gamma*H_b

    a_utd=a_ut*1.0
    b_utd=b_ut*1.0
    a_utd.transpose()
    b_utd.transpose()

    a_utd.setLabel([53,56,39,54,57])
    b_utd.setLabel([48,58,46,49,57])



#    E2_val=cost_f(A_f, a_ut, b_ut, c_ut, d_ut,Swap1,Swap2,Swap3)
    E2_val,E_positive=norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,a_ut,a_utd,b_ut,b_utd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)
    if E_positive<0:
     break_val="on"
     break


    #print "Gamma", Gamma
    if abs((0.5)*Norm_Z*Gamma) <1.0e-16 or  (abs((E1_val-E2_val)/E2_val))<1.0e-16 or abs(Gamma)<1.0e-16 :
     print "break2", E1_val, E2_val, Gamma, abs((0.5)*Norm_Z*Gamma), (abs((E1_val-E2_val)/E2_val))
     break
     
    if E1_val-E2_val < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0

   if break_val=="on": 
    print "break_val", break_val
    break

   Peps_1p=Peps_1p+(1.00)*Gamma*H_a
   Peps_2p=Peps_2p+(1.00)*Gamma*H_b
   Peps_1p.setLabel([6,35,37,55,34])
   Peps_2p.setLabel([11,31,32,33,15])

   Peps_1pd=Peps_1p*1.0
   Peps_2pd=Peps_2p*1.0
   Peps_1pd.transpose()
   Peps_2pd.transpose()
   Peps_1pd.setLabel([53,56,39,54,57])
   Peps_2pd.setLabel([48,58,46,49,57])


   time_val=time.time() - t0

  Peps_1p.setLabel([6,35,37,55,34])
  Peps_2p.setLabel([11,31,32,33,15])
  Peps_1p.permute([6,35,37,55,34],3)
  Peps_2p.permute([11,31,32,33,15],3)


  return Peps_1p, Peps_1pd, Peps_2p, Peps_2pd


def Energy_cal_grad(Peps_1, Peps_1d,Peps_2, Peps_2d,iden_h,H_orig,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left):

 A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*mps_boundry_right[Location*2]
 A=A*mps_boundry_right[Location*2+1]

 A=((A*Swap4)*Swap5)*(Peps_1*Swap6)


 B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

 B=(B*Swap12)*((Peps_2*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2d*Swap7)*Swap8))


 N_ten=A*B
 #print N_ten
 #N_ten.permute([200,-400,-100,300],2)
 #N_ten.setLabel([-200,-400,-100,-300])

 H_orig.setLabel([44,45,42,43])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 return h_h[0]/Norm_h[0], Norm_h[0]


def  Update_twotensor_local_grad( PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, Ham, H_orig, D, threshold, interval, Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 Peps_1.setLabel([6,35,37,55,34])
 Peps_1.permute([6,35,37,55,34],3)
 Peps_2.setLabel([11,31,32,33,15])
 Peps_2.permute([11,31,32,33,15],3)

 Peps_1d=Peps_1*1.0
 Peps_1d.transpose()
 Peps_1d.setLabel([54,57,53,56,39])

 Peps_2d=Peps_2*1.0
 Peps_2d.transpose()
 Peps_2d.setLabel([49,57,48,58,46])



 Swap1=fermionicOPT(Sys, Peps_1.bond(1), Peps_1.bond(0))
 Swap1.setLabel([70,36,35,-60])
 Swap2=Swap1*1.0
 Swap2.setLabel([56,53,-7,36])
 Swap3=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap3.setLabel([40,-9,57,54])
 Swap4=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap4.setLabel([41,58,40,90])
 Swap6=fermionicOPT(Sys,Peps_1.bond(2), Peps_1d.bond(0))
 Swap6.setLabel([42,55,37,38])
 Swap5=fermionicOPT(Sys,Peps_1d.bond(0), Peps_1.bond(2))
 Swap5.setLabel([38,39,58,44])


#####################################################


 Swap7=fermionicOPT(Sys,Peps_2d.bond(0), Peps_2d.bond(1))
 Swap7.setLabel([-14,50,49,57])
 Swap8=fermionicOPT(Sys, Peps_2.bond(0), Peps_2.bond(1))
 Swap8.setLabel([48,58,47,41])
 Swap9=fermionicOPT(Sys,Peps_2.bond(0), Peps_2.bond(2))
 Swap9.setLabel([47,46,51,45])
 Swap10=fermionicOPT(Sys,Peps_2.bond(0), Peps_2.bond(2))
 Swap10.setLabel([51,43,52,32])
 Swap11=fermionicOPT(Sys,Peps_2.bond(0), Peps_2.bond(1))
 Swap11.setLabel([52,34,-110,31])
 Swap12=fermionicOPT(Sys,Peps_2d.bond(0), Peps_2d.bond(1))
 Swap12.setLabel([33,-150,140,50])

######################################################################
######################################################################

 mps_boundry_left[Location*2].setLabel([16,-60,18])
 mps_boundry_left[Location*2+1].setLabel([18,6,21])
 mps_boundry_left[Location*2+2].setLabel([21,-110,22])
 mps_boundry_left[Location*2+3].setLabel([22,11,25])

 mps_boundry_right[Location*2].setLabel([17,-9,19])
 mps_boundry_right[Location*2+1].setLabel([19,90,20])
 mps_boundry_right[Location*2+2].setLabel([20,-14,23])
 mps_boundry_right[Location*2+3].setLabel([23,140,24])

######################################################
 A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*mps_boundry_right[Location*2]
 A=A*mps_boundry_right[Location*2+1]

 A=((A*Swap4)*Swap5)*(Peps_1*Swap6)


 B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

 B=(B*Swap12)*((Peps_2*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2d*Swap7)*Swap8))


 N_ten=A*B
 #print N_ten
 #N_ten.permute([200,-400,-100,300],2)
 #N_ten.setLabel([-200,-400,-100,-300])

 H_orig.setLabel([44,45,42,43])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_1c", Norm_h[0], h_h[0]/Norm_h[0]





 if Sys[8]=="Previous":
  Peps_1p=Peps_1*1.0
  Peps_2p=Peps_2*1.0
  Peps_1pd=Peps_1d*1.0
  Peps_2pd=Peps_2d*1.0

  Swap8p=Swap8*1.0
  Swap11p=Swap11*1.0
  Swap3p=Swap3*1.0
  Swap4p=Swap4*1.0

  Peps_1p.permute([6,35,37,55,34],5)
  Peps_2p.permute([11,31,32,33,15],5)
  Peps_1pd.permute([53,56,39,54,57],0)
  Peps_2pd.permute([48,58,46,49,57],0)



 if  Sys[8]=="TEBD_SVD":
   A=Peps_1*1.0
   A.setLabel([1,2,3,4,5])
   A.setLabel([ 1, 2, 3, 4, 5])
   A.permute([ 1, 2, 3, 4, 5], 6)
   Swap=fermionicOPT(Sys,A.bond(2), A.bond(3))
   Swap.setLabel([-3,-4,3,4])
   A=Swap*A
   A.permute([ 1, 2, -3, -4, 5], 3)
   A.setLabel([ 1, 2, 3, 4, 5])
   PEPS_A=A*1.0

   PEPS_B=Peps_2*1.0
   PEPS_B.setLabel([6,5,8,9,10])
   Ham.setLabel([11,12,3,8])


   Theta=(PEPS_A*Ham)*PEPS_B
   Theta.permute([1,2,11,4,6,12,9,10],4)

   #print len(D)
   U, V, S=TU.setTruncation(Theta, len(D))
   U.setLabel([1,2,11,4,-100])
   V.setLabel([-100,6,12,9,10])
   S=Sqrt(S)
   S.setLabel([-100,100])
   U=U*S
   S.setLabel([100,-100])
   V=V*S
   V.permute([6,100,12,9,10],3)
   V.setLabel([11,31,32,33,15])

   U.permute([1,2,11,4,100],4)
   Swap=fermionicOPT(Sys,U.bond(2), U.bond(3))
   Swap.setLabel([-11,-4,11,4])
   U=Swap*U
   U.permute([ 1, 2, -11, -4, 100], 3)
   U.setLabel([6,35,37,55,34])

   Peps_1p=U*1.0
   Peps_2p=V*1.0

   Peps_1pd=Peps_1p*1.0
   Peps_1pd.transpose()
   Peps_1pd.setLabel([54,57,53,56,39])


   Peps_2pd=Peps_2p*1.0
   Peps_2pd.transpose()
   Peps_2pd.setLabel([49,57,48,58,46])

   Swap8p=fermionicOPT(Sys,Peps_2p.bond(0), Peps_2p.bond(1))
   Swap8p.setLabel([48,58,47,41])
   Swap11p=fermionicOPT(Sys,Peps_2p.bond(0), Peps_2p.bond(1))
   Swap11p.setLabel([52,34,-110,31])

   Swap3p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap3p.setLabel([40,-9,57,54])
   Swap4p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap4p.setLabel([41,58,40,90])
   
   
   Peps_1p.permute([6,35,37,55,34],5)
   Peps_2p.permute([11,31,32,33,15],5)
   Peps_1pd.permute([53,56,39,54,57],0)
   Peps_2pd.permute([48,58,46,49,57],0)



   
 #print "Swap", Swap11.printDiagram(), Swap11p.printDiagram()
 #print "PEPS", Peps_2.printDiagram(), Peps_2p.printDiagram() 
 Ham.setLabel([51,52,53,54])
 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([44,45,42,43])
 Ham.setLabel([44,45,42,43])
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])



#  valf, E_pos=norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)
#  print valf
# 
# 
# 
# 

 Peps_1p, Peps_1pd,Peps_2p,Peps_2pd=Do_optimization_Grad_col(Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Sys,H_orig,Swap3p,Swap8p,Swap11p,Swap4p)

#######################################################



 E_val, Norm_val=Energy_cal_grad(Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,iden_h,H_orig,Swap1,Swap2,Swap3p,Swap4p,Swap5,Swap6,Swap7,Swap8p,Swap9,Swap10,Swap11p,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left)

 #print "E_2c", E_val, Norm_val

 if Norm_val<0: 
  print "Norm is <0"
  mps_boundry_right[2*Location+2]=mps_boundry_right[2*Location+2]*(-1.0)

 if Norm_val> threshold[1] or Norm_val<(threshold[0]):
  Peps_1p, Peps_2p=renormalize_tensor_c_single(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3p,Swap4p,Swap5,Swap6,Swap7,Swap8p,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,threshold,interval)



 return Peps_1p, Peps_2p, E_val




def  Update_twotensor_local_inv(PEPS_listten, E_list_down, E_list_up, mps_boundry_left, mps_boundry_right, Location, Ham, H_orig, D, threshold, interval, Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=copy.copy(PEPS_listten[Location])
 Peps_2=copy.copy(PEPS_listten[Location+1])

 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([70])
  Resul=Resul*iden
  iden.setLabel([-7])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_down[Location-1])
  E_left.setLabel([16,70,-7,17])

 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([25])
  Resul=copy.copy(iden)
  iden.setLabel([15])
  Resul=Resul*iden
  iden.setLabel([-150])
  Resul=Resul*iden
  iden.setLabel([24])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_up[Location+2])
  E_right.setLabel([25,15,-150,24])

 Peps_1.setLabel([6,35,37,55,34])
 Peps_1.permute([6,35,37,55,34],3)
 Peps_2.setLabel([11,31,32,33,15])
 Peps_2.permute([11,31,32,33,15],3)

 Peps_1d=Peps_1*1.0
 Peps_1d.transpose()
 Peps_1d.setLabel([54,57,53,56,39])

 Peps_2d=Peps_2*1.0
 Peps_2d.transpose()
 Peps_2d.setLabel([49,57,48,58,46])



 Swap1=fermionicOPT(Sys, Peps_1.bond(1), Peps_1.bond(0))
 Swap1.setLabel([70,36,35,-60])
 Swap2=Swap1*1.0
 Swap2.setLabel([56,53,-7,36])
 Swap3=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap3.setLabel([40,-9,57,54])
 Swap4=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap4.setLabel([41,58,40,90])
 Swap6=fermionicOPT(Sys,Peps_1.bond(2), Peps_1d.bond(0))
 Swap6.setLabel([42,55,37,38])
 Swap5=fermionicOPT(Sys,Peps_1d.bond(0), Peps_1.bond(2))
 Swap5.setLabel([38,39,58,44])


#####################################################


 Swap7=fermionicOPT(Sys,Peps_2d.bond(0), Peps_2d.bond(1))
 Swap7.setLabel([-14,50,49,57])
 Swap8=fermionicOPT(Sys, Peps_2.bond(0), Peps_2.bond(1))
 Swap8.setLabel([48,58,47,41])
 Swap9=fermionicOPT(Sys,Peps_2.bond(0), Peps_2.bond(2))
 Swap9.setLabel([47,46,51,45])
 Swap10=fermionicOPT(Sys,Peps_2.bond(0), Peps_2.bond(2))
 Swap10.setLabel([51,43,52,32])
 Swap11=fermionicOPT(Sys,Peps_2.bond(0), Peps_2.bond(1))
 Swap11.setLabel([52,34,-110,31])
 Swap12=fermionicOPT(Sys,Peps_2d.bond(0), Peps_2d.bond(1))
 Swap12.setLabel([33,-150,140,50])

######################################################################
######################################################################

 mps_boundry_left[Location*2].setLabel([16,-60,18])
 mps_boundry_left[Location*2+1].setLabel([18,6,21])
 mps_boundry_left[Location*2+2].setLabel([21,-110,22])
 mps_boundry_left[Location*2+3].setLabel([22,11,25])

 mps_boundry_right[Location*2].setLabel([17,-9,19])
 mps_boundry_right[Location*2+1].setLabel([19,90,20])
 mps_boundry_right[Location*2+2].setLabel([20,-14,23])
 mps_boundry_right[Location*2+3].setLabel([23,140,24])

######################################################
 A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
 A=A*(Swap1*((Peps_1d*Swap3)*Swap2))
 A=A*mps_boundry_right[Location*2]
 A=A*mps_boundry_right[Location*2+1]

 A=((A*Swap4)*Swap5)*(Peps_1*Swap6)


 B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

 B=(B*Swap12)*((Peps_2*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2d*Swap7)*Swap8))


 N_ten=A*B
 #print N_ten
 #N_ten.permute([200,-400,-100,300],2)
 #N_ten.setLabel([-200,-400,-100,-300])

 H_orig.setLabel([44,45,42,43])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_1c", Norm_h[0], h_h[0]/Norm_h[0]





 if Sys[8]=="Previous":
  Peps_1p=Peps_1*1.0
  Peps_2p=Peps_2*1.0
  Peps_1pd=Peps_1d*1.0
  Peps_2pd=Peps_2d*1.0

  Swap8p=Swap8*1.0
  Swap11p=Swap11*1.0
  Swap3p=Swap3*1.0
  Swap4p=Swap4*1.0




 if  Sys[8]=="TEBD_SVD":
   A=Peps_1*1.0
   A.setLabel([1,2,3,4,5])
   A.setLabel([ 1, 2, 3, 4, 5])
   A.permute([ 1, 2, 3, 4, 5], 6)
   Swap=fermionicOPT(Sys,A.bond(2), A.bond(3))
   Swap.setLabel([-3,-4,3,4])
   A=Swap*A
   A.permute([ 1, 2, -3, -4, 5], 3)
   A.setLabel([ 1, 2, 3, 4, 5])
   PEPS_A=A*1.0

   PEPS_B=Peps_2*1.0
   PEPS_B.setLabel([6,5,8,9,10])
   Ham.setLabel([11,12,3,8])


   Theta=(PEPS_A*Ham)*PEPS_B
   Theta.permute([1,2,11,4,6,12,9,10],4)

   #print len(D)
   U, V, S=TU.setTruncation(Theta, len(D))
   U.setLabel([1,2,11,4,-100])
   V.setLabel([-100,6,12,9,10])
   S=Sqrt(S)
   S.setLabel([-100,100])
   U=U*S
   S.setLabel([100,-100])
   V=V*S
   V.permute([6,100,12,9,10],3)
   V.setLabel([11,31,32,33,15])

   U.permute([1,2,11,4,100],4)
   Swap=fermionicOPT(Sys,U.bond(2), U.bond(3))
   Swap.setLabel([-11,-4,11,4])
   U=Swap*U
   U.permute([ 1, 2, -11, -4, 100], 3)
   U.setLabel([6,35,37,55,34])

   Peps_1p=U*1.0
   Peps_2p=V*1.0

   Peps_1pd=Peps_1p*1.0
   Peps_1pd.transpose()
   Peps_1pd.setLabel([54,57,53,56,39])

   Peps_2pd=Peps_2p*1.0
   Peps_2pd.transpose()
   Peps_2pd.setLabel([49,57,48,58,46])

   Swap8p=fermionicOPT(Sys,Peps_2p.bond(0), Peps_2p.bond(1))
   Swap8p.setLabel([48,58,47,41])
   Swap11p=fermionicOPT(Sys,Peps_2p.bond(0), Peps_2p.bond(1))
   Swap11p.setLabel([52,34,-110,31])

   Swap3p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap3p.setLabel([40,-9,57,54])
   Swap4p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap4p.setLabel([41,58,40,90])







 Ham.setLabel([51,52,53,54])
 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([44,45,42,43])
 Ham.setLabel([44,45,42,43])
 iden_h=Ham*1.0
 iden_h.setLabel([-3,-13,3,13])
 iden_h.identity()
 iden_h.setLabel([44,45,42,43])

 Peps_1p_init=Peps_1p*1.0
 Peps_2p_init=Peps_2p*1.0
 Peps_2pd_init=Peps_2pd*1.0
 Peps_1pd_init=Peps_1pd*1.0

 valf, E_pos=norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(Sys[3]):
  count=count+1
  E_2=E_1*1.0
  val, E_pos=norm_f_val_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   Peps_1p=Peps_1p_init*1.0
   Peps_1pd=Peps_1pd_init*1.0
   Peps_2p=Peps_2p_init*1.0
   Peps_2pd=Peps_2pd_init*1.0
   break
  else:
   Peps_1p_init=Peps_1p*1.0
   Peps_1pd_init=Peps_1pd*1.0
   Peps_2p_init=Peps_2p*1.0
   Peps_2pd_init=Peps_2pd*1.0


  Peps_1p, Peps_1pd=optimum_0_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)

  Peps_2p, Peps_2pd=optimum_1_inv_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,Swap3p,Swap8p,Swap11p,Swap4p)



 #val=norm_f_val(N_ten, l_u, r_u, r_d, l_d, l_up, r_up,l_dp, r_dp, Ham,iden_h,H)
 #print "Tru_row",  abs(valf), abs(val), (abs(valf)-abs(val))/abs(valf), count
 #print "valf", valf
 if abs(valf)<1.0e-12 or abs(valf)>1.0e+12:
  print "warning_norm_in_optimization",  abs(valf), "count", count








######################################################

 A=E_left*mps_boundry_left[Location*2]
 A=mps_boundry_left[Location*2+1]*A
 A=A*(Swap1*((Peps_1pd*Swap3p)*Swap2))
 A=A*mps_boundry_right[Location*2]
 A=A*mps_boundry_right[Location*2+1]

 A=((A*Swap4p)*Swap5)*(Peps_1p*Swap6)


 B=(mps_boundry_right[2*Location+3]*mps_boundry_right[2*Location+2])*E_right


 B=(B*Swap12)*(Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3]))
 B=((B*Swap10)*Swap11p)*(Swap9*((Peps_2pd*Swap7)*Swap8p))


 N_ten=A*B

 H_orig.setLabel([44,45,42,43])
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_2c", Norm_h[0], h_h[0]/Norm_h[0]

 if Norm_h[0]> threshold[1] or Norm_h[0]<threshold[0]:
  Peps_1p, Peps_2p=renormalize_tensor_c_single(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3p,Swap4p,Swap5,Swap6,Swap7,Swap8p,Swap9,Swap10,Swap11p,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,threshold,interval)

 return Peps_1p, Peps_2p, h_h[0]/Norm_h[0]



def optimum_1_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p):

 Env_r=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_r=Env_r*(Peps_1p*Swap1)
 Env_r=Env_r*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 Env_r=Env_r*((Swap3p*Peps_1pd)*Swap2)
 Env_r=(Env_r*(Swap8*Swap9))*Swap10p
 Env_r=Env_r*iden_h

 Env_r=Env_r*Swap4p
 Env_r=Env_r*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_r=Env_r*Swap5p
 Env_r=Env_r*Swap12
 Env_r=Env_r*Swap11

 A=(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])*E_right

 A=Swap7*A
 A=Swap6*A

 Env_r=(Env_r*A)

 Env_r.permute([43,42,47,45,44,54,52,51,53,15],5)
 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1

 Env_s=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_s=Env_s*(Peps_1p*Swap1)
 Env_s=Env_s*mps_boundry_down[2*Location]
 Env_s=Env_s*mps_boundry_down[2*Location+1]
 Env_s=Env_s*((Swap3*Peps_1d)*Swap2)
 Env_s=(Env_s*(Swap8*Swap9))*Swap10p
 Env_s=Env_s*Ham
 Env_s=((Swap11*Env_s)*Swap12)*Swap4
 Env_s=Env_s*mps_boundry_down[2*Location+2]
 Env_s=Env_s*mps_boundry_down[2*Location+3]
 Env_s=Env_s*((Swap5*Peps_2d)*Swap6)
 Env_s=Env_s*E_right
 Env_s=Swap7*Env_s

 Env_s=Env_s*(mps_boundry_up[2*Location+3]*mps_boundry_up[2*Location+2])


 Env_s.permute([ 54,52,51,53,15],5)

 Env_s1=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_s1=Env_s1*(Peps_1*Swap1)
 Env_s1=Env_s1*mps_boundry_down[2*Location]
 Env_s1=Env_s1*mps_boundry_down[2*Location+1]

 Env_s1=Env_s1*((Swap3p*Peps_1pd)*Swap2)
 Env_s1=(Env_s1*(Swap8*Swap9))*Swap10
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*((Peps_2*Swap11))
 
 Env_s1=(Env_s1*Swap12)*Swap4p

 Env_s1=Env_s1*mps_boundry_up[2*Location+2]
 Env_s1=Env_s1*mps_boundry_up[2*Location+3]

 Env_s1=Env_s1*mps_boundry_down[2*Location+2]
 Env_s1=Env_s1*mps_boundry_down[2*Location+3]

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*E_right
 Env_s1=Env_s1*(Swap6*Swap5p)
 Env_s1.permute([ 43,42,47,45,44],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([54,52,51,53,15,43,42,47,45,44])
 
 
 Peps_2p=A2_inv*Env_s
 Peps_2p.permute([54,52,51,53,15],3)
 Peps_2pd=Peps_2p*1.0
 Peps_2pd.transpose()
 Peps_2pd.setLabel([45,44,43,42,47])

 return Peps_2p, Peps_2pd







######@profile
def renormalize_tensor_c_single(Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_right,mps_boundry_left,Location,E_right,E_left,threshold,interval):

 A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location*2]*mps_boundry_right[Location*2+1]

 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)


 B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

 B=(B*Swap12)*((Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))


 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "norm_val", norm_val
 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold[1]:
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
   A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
   A=A*mps_boundry_right[Location*2]*mps_boundry_right[Location*2+1]

   A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)


   B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

   B=(B*Swap12)*((Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
   B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))


   N_ten=A*B
   iden_h.setLabel([44,45,42,43])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<(threshold[0]):
   Peps_1p=Peps_1p*alpha
   Peps_1pd=Peps_1pd*alpha
   Peps_2p=Peps_2p*alpha
   Peps_2pd=Peps_2pd*alpha
   A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
   A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
   A=A*mps_boundry_right[Location*2]*mps_boundry_right[Location*2+1]

   A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)


   B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right

   B=(B*Swap12)*((Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
   B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))

   N_ten=A*B
   iden_h.setLabel([44,45,42,43])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])

 A=E_left*(mps_boundry_left[Location*2]*mps_boundry_left[Location*2+1])
 A=A*(Swap1*((Peps_1pd*Swap3)*Swap2))
 A=A*mps_boundry_right[Location*2]*mps_boundry_right[Location*2+1]
 A=((A*Swap4)*Swap5)*(Peps_1p*Swap6)
 B=(mps_boundry_right[2*Location+2]*mps_boundry_right[2*Location+3])*E_right
 B=(B*Swap12)*((Peps_2p*(mps_boundry_left[2*Location+2]*mps_boundry_left[2*Location+3])))
 B=((B*Swap10)*Swap11)*(Swap9*((Peps_2pd*Swap7)*Swap8))
 N_ten=A*B
 iden_h.setLabel([44,45,42,43])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "Fixednorm_local", abs(norm_val)
 return Peps_1p, Peps_2p



def renormalize_tensor_r_single(Peps_1, Peps_1d,Peps_2, Peps_2d,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,threshold,interval):

 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2*Swap11)

 N_ten=A*B
 iden_h.setLabel([40,49,56,55])
 #print N_ten.printDiagram()
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "norm_val", norm_val
 count=0

 while count<100:
  #print norm_val, count, 
  if abs(norm_val)<threshold[1] and abs(norm_val)>threshold[0]:break
  count=count+1
  alpha=1.0-interval
  if abs(norm_val)>threshold[1]:
   Peps_1=Peps_1*alpha
   Peps_1d=Peps_1d*alpha
   Peps_2=Peps_2*alpha
   Peps_2d=Peps_2d*alpha
   A=E_left*mps_boundry_up[2*Location]
   A=A*mps_boundry_up[2*Location+1]

   A=A*(Peps_1*Swap1)
   A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

   A=A*((Swap3*Peps_1d)*Swap2)
   A=(A*(Swap8*Swap9))*Swap10

   B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

   B=Swap7*B
   B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

   B=B*((Peps_2d*Swap6)*Swap5)
   B=(B*Swap4)*Swap12
   B=B*(Peps_2*Swap11)

   N_ten=A*B
   iden_h.setLabel([40,49,56,55])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])
   #norm_val=Cal_norm(PEPS_listten, N_x, D, chi_try, d)

  alpha=1.0+interval
  if abs(norm_val)<(threshold[0]):
   Peps_1=Peps_1*alpha
   Peps_1d=Peps_1d*alpha
   Peps_2=Peps_2*alpha
   Peps_2d=Peps_2d*alpha
   A=E_left*mps_boundry_up[2*Location]
   A=A*mps_boundry_up[2*Location+1]

   A=A*(Peps_1*Swap1)
   A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

   A=A*((Swap3*Peps_1d)*Swap2)
   A=(A*(Swap8*Swap9))*Swap10

   B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

   B=Swap7*B
   B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

   B=B*((Peps_2d*Swap6)*Swap5)
   B=(B*Swap4)*Swap12
   B=B*(Peps_2*Swap11)

   N_ten=A*B
   iden_h.setLabel([40,49,56,55])
   Norm_h=N_ten*iden_h
   norm_val=abs(Norm_h[0])

 A=E_left*mps_boundry_up[2*Location]
 A=A*mps_boundry_up[2*Location+1]

 A=A*(Peps_1*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2*Swap11)

 N_ten=A*B
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 norm_val=abs(Norm_h[0])
 #print "Fixednorm_local", abs(norm_val)
 return Peps_1, Peps_2



def optimum_0_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p):


 Env_r=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_r=Swap7*Env_r
 Env_r=Env_r*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_r=Env_r*((Peps_2pd*Swap6)*Swap5p)
 Env_r=(Env_r*Swap4p)*Swap12
 Env_r=Env_r*(Peps_2p*Swap11)
 Env_r=Env_r*iden_h

 Env_r=Env_r*(Swap10p*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1]))

 Env_r=(Env_r)*(Swap9)
 Env_r=(Env_r)*(Swap8)
 Env_r=Env_r*Swap3p


 A=E_left*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])
 A=A*Swap1
 A=A*Swap2


 Env_r=(Env_r)*(A)


 Env_r.permute([ 32,33,37, 35,34,6,31,58,59,10],5)

 Env_r1=Env_r*1.0
 Env_r1.transpose()
 Env_r=Env_r+Env_r1


 Env_s=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_s=Swap7*Env_s
 Env_s=Env_s*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_s=Env_s*((Peps_2d*Swap6)*Swap5)
 Env_s=(Env_s*Swap4)*Swap12
 Env_s=Env_s*(Peps_2p*Swap11)
 Env_s=Env_s*Ham


 Env_s=(Env_s*(Swap8*Swap9))*Swap10p


 Env_s=(Env_s*((Peps_1d*Swap3)*Swap2))*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 Env_s=(Env_s)*(E_left*Swap1)

 Env_s=Env_s*mps_boundry_up[2*Location]
 Env_s=Env_s*mps_boundry_up[2*Location+1]

 Env_s.permute([6,31,58,59,10],5)

 Env_s1=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_s1=Env_s1*((Peps_2pd*Swap6)*Swap5p)
 Env_s1=(Env_s1*Swap4p)*Swap12
 Env_s1=Env_s1*(Peps_2*Swap11)
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*((mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])*Swap10)

 Env_s1=(Env_s1*(Swap8*Swap9))
 Env_s1=(Env_s1*(Swap3p))
 Env_s1=Env_s1*Peps_1
 Env_s1=Env_s1*E_left
 Env_s1=Env_s1*(((mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])*Swap1)*Swap2)
 Env_s1.permute([32,33,37,35,34],0)
 Env_s1.transpose()
 Env_s=Env_s1+Env_s


 row, colm=cal_rowcol(Env_r)
 if (row<=colm):
  U,V,S=TU.setTruncation(Env_r,row)
 else:
  U,V,S=TU.setTruncation(Env_r,colm)


 U.transpose()
 V.transpose()
 S=inverse(S)
 
 U.setLabel([7,8,9,10,11,12])
 V.setLabel([1,2,3,4,5,6])
 S.setLabel([6,7])

 A2_inv=V*S*U
 A2_inv.permute([1,2,3,4,5,8,9,10,11,12],5)
 A2_inv.setLabel([6,31,58,59,10,32,33,37,35,34])

 Peps_1p=A2_inv*Env_s
 Peps_1p.permute([6,31,58,59,10],3)
 Peps_1pd=Peps_1p*1.0
 Peps_1pd.transpose()
 Peps_1pd.setLabel([35,34,32,33,37])

 return Peps_1p, Peps_1pd





#@profile
def norm_f_val_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p):

 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1p*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3p*Peps_1pd)*Swap2)

 A=(A*(Swap8*Swap9))*Swap10p
 A=A*iden_h

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2pd*Swap6)*Swap5p)
 B=(B*Swap4p)*Swap12
 B=B*(Peps_2p*Swap11)

 val1=A*B

 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1p*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10p
 A=A*Ham

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2p*Swap11)

 val2=A*B

 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1*Swap1)
 A=A*mps_boundry_down[2*Location]
 A=A*mps_boundry_down[2*Location+1]

 A=A*((Swap3p*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10
 A=A*Ham

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2pd*Swap6)*Swap5p)
 B=(B*Swap4p)*Swap12
 B=B*(Peps_2*Swap11)

 val3=A*B
 #print val1, val2, val3
 return val1[0]-val2[0]-val3[0], val1[0]





#@profile
def Obtain_grad_two_row(Peps_1, Peps_1d, Peps_2, Peps_2d, Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p):

 Env_s=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_s=Swap7*Env_s
 Env_s=Env_s*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_s=Env_s*((Peps_2pd*Swap6)*Swap5p)
 Env_s=(Env_s*Swap4p)*Swap12
 Env_s=Env_s*(Peps_2p*Swap11)
 Env_s=Env_s*iden_h

 Env_s=(Env_s*(Swap8*Swap9))*Swap10p
 Env_s=(Env_s*((Peps_1pd*Swap3p)*Swap2))*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 Env_s=(Env_s)*(E_left*Swap1)

 Env_s=Env_s*mps_boundry_up[2*Location]
 Env_s=Env_s*mps_boundry_up[2*Location+1]

 #print Env_s.printDiagram()
 Env_s.permute([6,31,58,59,10],5)

 Env_s1=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_s1=Env_s1*((Peps_2pd*Swap6)*Swap5p)
 Env_s1=(Env_s1*Swap4p)*Swap12
 Env_s1=Env_s1*(Peps_2p*Swap11)
 Env_s1=Env_s1*iden_h

 Env_s1=Env_s1*((mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])*Swap10p)

 Env_s1=(Env_s1*(Swap8*Swap9))
 Env_s1=(Env_s1*(Swap3p))
 Env_s1=Env_s1*Peps_1p
 Env_s1=Env_s1*E_left
 Env_s1=Env_s1*(((mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])*Swap1)*Swap2)
 Env_s1.permute([32,33,37,35,34],0)
 Env_s1.transpose()
 Env_s=Env_s+Env_s1
 D_a=Env_s*1.0


 Env_s=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_s=Swap7*Env_s
 Env_s=Env_s*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_s=Env_s*((Peps_2d*Swap6)*Swap5)
 Env_s=(Env_s*Swap4)*Swap12
 Env_s=Env_s*(Peps_2p*Swap11)
 Env_s=Env_s*Ham


 Env_s=(Env_s*(Swap8*Swap9))*Swap10p


 Env_s=(Env_s*((Peps_1d*Swap3)*Swap2))*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 Env_s=(Env_s)*(E_left*Swap1)

 Env_s=Env_s*mps_boundry_up[2*Location]
 Env_s=Env_s*mps_boundry_up[2*Location+1]

 Env_s.permute([6,31,58,59,10],5)

 Env_s1=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 Env_s1=Env_s1*((Peps_2pd*Swap6)*Swap5p)
 Env_s1=(Env_s1*Swap4p)*Swap12
 Env_s1=Env_s1*(Peps_2*Swap11)
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*((mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])*Swap10)

 Env_s1=(Env_s1*(Swap8*Swap9))
 Env_s1=(Env_s1*(Swap3p))
 Env_s1=Env_s1*Peps_1
 Env_s1=Env_s1*E_left
 Env_s1=Env_s1*(((mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])*Swap1)*Swap2)
 Env_s1.permute([32,33,37,35,34],0)
 Env_s1.transpose()
 Env_s=Env_s1+Env_s
 D_a=D_a+Env_s*(-1.0)


 Env_s=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_s=Env_s*(Peps_1p*Swap1)
 Env_s=Env_s*mps_boundry_down[2*Location]
 Env_s=Env_s*mps_boundry_down[2*Location+1]
 Env_s=Env_s*((Swap3p*Peps_1pd)*Swap2)
 Env_s=(Env_s*(Swap8*Swap9))*Swap10p
 Env_s=Env_s*iden_h
 Env_s=((Swap11*Env_s)*Swap12)*Swap4p
 Env_s=Env_s*mps_boundry_down[2*Location+2]
 Env_s=Env_s*mps_boundry_down[2*Location+3]
 Env_s=Env_s*((Swap5p*Peps_2pd)*Swap6)
 Env_s=Env_s*E_right
 Env_s=Swap7*Env_s

 Env_s=Env_s*(mps_boundry_up[2*Location+3]*mps_boundry_up[2*Location+2])


 Env_s.permute([ 54,52,51,53,15],5)

 Env_s1=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_s1=Env_s1*(Peps_1p*Swap1)
 Env_s1=Env_s1*mps_boundry_down[2*Location]
 Env_s1=Env_s1*mps_boundry_down[2*Location+1]

 Env_s1=Env_s1*((Swap3p*Peps_1pd)*Swap2)
 Env_s1=(Env_s1*(Swap8*Swap9))*Swap10p
 Env_s1=Env_s1*iden_h

 Env_s1=Env_s1*((Peps_2p*Swap11))
 
 Env_s1=(Env_s1*Swap12)*Swap4p

 Env_s1=Env_s1*mps_boundry_up[2*Location+2]
 Env_s1=Env_s1*mps_boundry_up[2*Location+3]

 Env_s1=Env_s1*mps_boundry_down[2*Location+2]
 Env_s1=Env_s1*mps_boundry_down[2*Location+3]

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*E_right
 Env_s1=Env_s1*(Swap6*Swap5p)
 Env_s1.permute([ 43,42,47,45,44],0)

 Env_s1.transpose()
 Env_s=Env_s+Env_s1
 D_b=Env_s*1.0

 Env_s=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_s=Env_s*(Peps_1p*Swap1)
 Env_s=Env_s*mps_boundry_down[2*Location]
 Env_s=Env_s*mps_boundry_down[2*Location+1]
 Env_s=Env_s*((Swap3*Peps_1d)*Swap2)
 Env_s=(Env_s*(Swap8*Swap9))*Swap10p
 Env_s=Env_s*Ham
 Env_s=((Swap11*Env_s)*Swap12)*Swap4
 Env_s=Env_s*mps_boundry_down[2*Location+2]
 Env_s=Env_s*mps_boundry_down[2*Location+3]
 Env_s=Env_s*((Swap5*Peps_2d)*Swap6)
 Env_s=Env_s*E_right
 Env_s=Swap7*Env_s

 Env_s=Env_s*(mps_boundry_up[2*Location+3]*mps_boundry_up[2*Location+2])


 Env_s.permute([ 54,52,51,53,15],5)

 Env_s1=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 Env_s1=Env_s1*(Peps_1*Swap1)
 Env_s1=Env_s1*mps_boundry_down[2*Location]
 Env_s1=Env_s1*mps_boundry_down[2*Location+1]

 Env_s1=Env_s1*((Swap3p*Peps_1pd)*Swap2)
 Env_s1=(Env_s1*(Swap8*Swap9))*Swap10
 Env_s1=Env_s1*Ham

 Env_s1=Env_s1*((Peps_2*Swap11))
 
 Env_s1=(Env_s1*Swap12)*Swap4p

 Env_s1=Env_s1*mps_boundry_up[2*Location+2]
 Env_s1=Env_s1*mps_boundry_up[2*Location+3]

 Env_s1=Env_s1*mps_boundry_down[2*Location+2]
 Env_s1=Env_s1*mps_boundry_down[2*Location+3]

 Env_s1=Swap7*Env_s1
 Env_s1=Env_s1*E_right
 Env_s1=Env_s1*(Swap6*Swap5p)
 Env_s1.permute([ 43,42,47,45,44],0)

 Env_s1.transpose()
 Env_s=Env_s1+Env_s
 D_b=D_b+Env_s*(-1.0)


 #print D_a.printDiagram()
 D_a.transpose()
 D_a.permute([6,31,58,59,10],5)

 D_b.transpose()
 D_b.permute([ 54,52,51,53,15],5)


 return D_a, D_b


#@profile
def  Do_optimization_Grad_row( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Sys,H_orig,Swap3p,Swap4p,Swap5p,Swap10p):
  Opt_method="ST"
  #Opt_method="grad"


  Peps_1p_um=Peps_1p*1.0
  Peps_2p_um=Peps_2p*1.0
  Peps_1pd_um=Peps_1p_um*1.0
  Peps_2pd_um=Peps_2p_um*1.0
  
  Peps_1pd_um.transpose()
  Peps_2pd_um.transpose()
  
  Peps_1pd_um.setLabel([32,33,37,35,34])
  Peps_2pd_um.setLabel([43,42,47,45,44])

  time_val=0
  Es, E_posi=norm_f_val_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)

  Ef=0
  E2_val=0
  EE2_val=0
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=1.e+18
  EE_previous=1.e+18

  count=0
  D_list=[0]*4
  H_list=[0]*4
  H_a=0; H_b=0;H_c=0;H_d=0;
  for i in xrange(Sys[7]):
   count+=1
   t0=time.time()

   E1_val,E_posi=norm_f_val_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)


#   E_energy, Norm_val=Energy_row_grad( Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,H_orig,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Sys)


   if (E_posi<0):
    print "Norm<0", E_posi
    Peps_1p=copy.copy(Peps_1p_um) 
    Peps_2p=copy.copy(Peps_2p_um) 
    Peps_1pd=copy.copy(Peps_1pd_um) 
    Peps_2pd=copy.copy(Peps_2pd_um) 
    break



   #if i%200==0:
   #print 'E1=', E1_val, abs((E_previous-E1_val)/E1_val), i, count, time_val
   #print 'E_energy=', E_energy, abs((EE_previous-E_energy)/E_energy), i, count, time_val


   D_a,D_b=Obtain_grad_two_row(Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)
   D_a=(-1.0)*D_a
   D_b=(-1.0)*D_b

   if i is 0:
    H_a=D_a
    H_b=D_b
   else:
    Z_a=D_a+(-1.0)*D_list[0]
    Z_b=D_b+(-1.0)*D_list[1]
    A=Z_a*D_a
    B=Z_b*D_b
    A1=D_list[0]*D_list[0]
    A2=D_list[1]*D_list[1]
    Gamma_grad=(A[0]+B[0]) / (A1[0]+A2[0])
    if Opt_method is 'ST':Gamma_grad=0;
    H_a=D_a+(Gamma_grad)*H_list[0]
    H_b=D_b+(Gamma_grad)*H_list[1]

#    A=D_a*D_list[0]
#    B=D_b*D_list[1]
#    C=D_c*D_list[2]
#    D=D_d*D_list[3]
#    check=A[0]+B[0]+C[0]+D[0] 
#    print "check", check 

   D_list[0]=copy.copy(D_a)
   D_list[1]=copy.copy(D_b)

   H_list[0]=copy.copy(H_a)
   H_list[1]=copy.copy(H_b)


   pos_chack="off"
   A=D_a*H_a
   B=D_b*H_b
   
   Norm_Z=A[0]+B[0]
   
   #print 'Norm', Norm_Z
   if (E1_val<E_previous) or (i is 0):
    if (abs(E1_val) > 1.0e-10):
     if abs((E_previous-E1_val)/E1_val) < 1.0e-10:
      print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)/E1_val), i
      break
     else: 
      if abs((E_previous-E1_val)) < 1.0e-15:
       print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)), i
       break

#   if (E_energy<EE_previous) or (i is 0):
#    if (abs(E_energy) > 1.0e-10):
#     if abs((EE_previous-E_energy)/E_energy) < 1.0e-10:
#      print 'Differnance Satisfied!', EE_previous, E_energy, abs((EE_previous-E_energy)/E_energy), i
#      break
#     else: 
#      if abs((EE_previous-E_energy)) < 1.0e-15:
#       print 'Differnance Satisfied!', EE_previous, E_energy, abs((EE_previous-E_energy)), i
#       break




#   if (E_energy>EE_previous):
#    print "break, energy_r", E_energy, EE_previous, i
#    Peps_1p=copy.copy(Peps_1p_um) 
#    Peps_2p=copy.copy(Peps_2p_um) 
#    Peps_1pd=copy.copy(Peps_1pd_um) 
#    Peps_2pd=copy.copy(Peps_2pd_um) 
#    break
#   else:
#    Peps_1p_um=copy.copy(Peps_1p) 
#    Peps_2p_um=copy.copy(Peps_2p) 
#    Peps_1pd_um=Peps_1p_um*1.0
#    Peps_2pd_um=Peps_2p_um*1.0
#    Peps_1pd_um.setLabel([32,33,37,35,34])
#    Peps_2pd_um.setLabel([43,42,47,45,44])




   if (E1_val>E_previous):
    print "break, not satisfied", E1_val, E_previous
    Peps_1p=copy.copy(Peps_1p_um) 
    Peps_2p=copy.copy(Peps_2p_um) 
    Peps_1pd=copy.copy(Peps_1pd_um) 
    Peps_2pd=copy.copy(Peps_2pd_um) 

    break
   else:
    Peps_1p_um=copy.copy(Peps_1p) 
    Peps_2p_um=copy.copy(Peps_2p) 
    Peps_1pd_um=Peps_1p_um*1.0
    Peps_2pd_um=Peps_2p_um*1.0

    Peps_1pd_um.transpose()
    Peps_2pd_um.transpose()

    Peps_1pd_um.setLabel([32,33,37,35,34])
    Peps_2pd_um.setLabel([43,42,47,45,44])




   E_previous=E1_val
   #EE_previous=E_energy
   
   if abs(Norm_Z) < 1.0e-12:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   if (i%15)==0: 
    Gamma=1
   else: 
    if Gamma >= 1: Gamma*=(1.00/100)
    if Gamma < 1: Gamma*=100
   #Gamma=1
   #print "Gamma", Gamma
   while Break_loop is 1:
    count+=1
    a_ut=Peps_1p+(2.00)*Gamma*H_a
    b_ut=Peps_2p+(2.00)*Gamma*H_b
    a_utd=a_ut*1.0
    b_utd=b_ut*1.0

    a_utd.transpose()
    b_utd.transpose()

    a_utd.setLabel([32,33,37,35,34])
    b_utd.setLabel([43,42,47,45,44])
    
    
    
    
    E2_val,E_posi=norm_f_val_inv_row_single(Peps_1,Peps_1d,Peps_2,Peps_2d,a_ut,a_utd,b_ut,b_utd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)
    #print "E_2_1", E2_val

    if E_posi<0:
     pos_chack="on"
     break

    if abs((0.5)*Norm_Z*Gamma) > 1.0e+12 or  abs(Gamma)>1.0e+12 :
     print "break1", E1_val, abs((0.5)*Norm_Z*Gamma), E2_val, Gamma
     Gamma=1
     break
    if E1_val-E2_val >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0

   Break_loop=1
   while Break_loop is 1:
    count+=1
    a_ut=Peps_1p+(1.00)*Gamma*H_a
    b_ut=Peps_2p+(1.00)*Gamma*H_b

    a_utd=a_ut*1.0
    b_utd=b_ut*1.0

    a_utd.transpose()
    b_utd.transpose()

    a_utd.setLabel([32,33,37,35,34])
    b_utd.setLabel([43,42,47,45,44])



    E2_val,E_posi=norm_f_val_inv_row_single(Peps_1,Peps_1d,Peps_2,Peps_2d,a_ut,a_utd,b_ut,b_utd,Ham_ham,iden_h,H_ham,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)
    #print "E_2_0",E2_val
    if E_posi<0:
     pos_chack="on"
     break

    #print "Gamma", Gamma
    if abs((0.5)*Norm_Z*Gamma) <1.0e-16 or  (abs((E1_val-E2_val)/E2_val))<1.0e-16 or abs(Gamma)<1.0e-16 :
     print "break2", E1_val, E2_val, Gamma, abs((0.5)*Norm_Z*Gamma), (abs((E1_val-E2_val)/E2_val))
     break
     
    if E1_val-E2_val < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0

   if pos_chack=="on": 
    print "norm<0"

   
   Peps_1p=Peps_1p+(1.00)*Gamma*H_a
   Peps_2p=Peps_2p+(1.00)*Gamma*H_b
   Peps_1p.setLabel([6,31,58,59,10])
   Peps_2p.setLabel([54,52,51,53,15])


   Peps_1pd=Peps_1p*1.0
   Peps_2pd=Peps_2p*1.0

   Peps_1pd.transpose()
   Peps_2pd.transpose()

   Peps_1pd.setLabel([32,33,37,35,34])
   Peps_2pd.setLabel([43,42,47,45,44])


   time_val=time.time() - t0



  Peps_1p.permute([6,31,58,59,10],3)
  Peps_2p.permute([54,52,51,53,15],3)

  return Peps_1p, Peps_1pd, Peps_2p, Peps_2pd



def  Energy_row_grad( Peps_1, Peps_1d,Peps_2, Peps_2d,H_orig,iden_h,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Sys):

 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1*Swap1)
 A=A*mps_boundry_down[2*Location]
 A=A*mps_boundry_down[2*Location+1]

 A=A*((Swap3*Peps_1d)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2d*Swap6)*Swap5)
 B=(B*Swap4)*Swap12
 B=B*(Peps_2*Swap11)

 N_ten=A*B

 N_ten.permute([40,49,56,55],2)

 H_orig.setLabel([40,49,56,55])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig

 return h_h[0]/Norm_h[0], Norm_h[0]



def  Update_twotensor_local_row_grad(PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, Ham, H_orig, D,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])



 Peps_1.setLabel([6,31,58,59,10])
 Peps_1d=Peps_1*1.0
 Peps_1d.transpose()
 Peps_1d.setLabel([35,34,32,33,37])
 Peps_2.setLabel([54,52,51,53,15])
 Peps_2d=Peps_2*1.0
 Peps_2d.transpose()
 Peps_2d.setLabel([45,44,43,42,47])


 Swap1=fermionicOPT(Sys,Peps_1.bond(1), Peps_1.bond(0))
 Swap1.setLabel([22,30,31,17])
 Swap2=fermionicOPT(Sys,Peps_1.bond(1), Peps_1.bond(0))
 Swap2.setLabel([33,32,-7,30])
 Swap3=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap3.setLabel([38,36,34,35])
 Swap4=fermionicOPT(Sys,Peps_2.bond(1), Peps_2.bond(0))
 Swap4.setLabel([23,41,48,36])
 Swap5=fermionicOPT(Sys,Peps_2.bond(1), Peps_2.bond(0))
 Swap5.setLabel([42,43,-10,41])
 Swap6=fermionicOPT(Sys,Peps_2d.bond(1), Peps_2d.bond(0))
 Swap6.setLabel([46,-14,44,45])

 Swap7=fermionicOPT(Sys,Peps_2d.bond(1), Peps_2d.bond(0))
 Swap7.setLabel([24,53,46,20])

 Swap8=fermionicOPT(Sys,Peps_1.bond(2), Peps_1d.bond(1))
 Swap8.setLabel([37,39,40,38])

 Swap9=fermionicOPT(Sys,Peps_1.bond(2), Peps_1d.bond(1))
 Swap9.setLabel([56,57,58,39])

 Swap10=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap10.setLabel([25,59,57,54])

 Swap11=fermionicOPT(Sys,Peps_2.bond(2), Peps_2.bond(1))
 Swap11.setLabel([55,50,51,52])

 Swap12=fermionicOPT(Sys,Peps_2.bond(2), Peps_2.bond(1))
 Swap12.setLabel([47,48,49,50])



######################################################################

 mps_boundry_up[Location*2].setLabel([16,10,266])
 mps_boundry_up[Location*2+1].setLabel([266,25,277])
 mps_boundry_up[Location*2+2].setLabel([277,15,288])
 mps_boundry_up[Location*2+3].setLabel([288,24,21])

 mps_boundry_down[Location*2].setLabel([18,22,311])
 mps_boundry_down[Location*2+1].setLabel([311,-7,322])
 mps_boundry_down[Location*2+2].setLabel([322,23,333])
 mps_boundry_down[Location*2+3].setLabel([333,-10,19])

######################################################

# A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

# A=A*(Peps_1*Swap1)
# A=A*mps_boundry_down[2*Location]
# A=A*mps_boundry_down[2*Location+1]

# A=A*((Swap3*Peps_1d)*Swap2)
# A=(A*(Swap8*Swap9))*Swap10

# B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

# B=Swap7*B
# B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

# B=B*((Peps_2d*Swap6)*Swap5)
# B=(B*Swap4)*Swap12
# B=B*(Peps_2*Swap11)

# N_ten=A*B
# 
# N_ten.permute([40,49,56,55],2)



# H_orig.setLabel([40,49,56,55])
# iden_h=H_orig*1.0
# iden_h.identity()
# iden_h.setLabel([40,49,56,55])
# Norm_h=N_ten*iden_h
# h_h=N_ten*H_orig
# print "E_1r", Norm_h[0], h_h[0]/Norm_h[0]



# Peps_1p=Peps_1*1.0
# Peps_2p=Peps_2*1.0
# 
# Peps_1pd=Peps_1p*1.0
# Peps_1pd.setLabel([32,33,37,35,34])
# Peps_2pd=Peps_2p*1.0
# Peps_2pd.setLabel([43,42,47,45,44])





 if  Sys[8]=="Previous":


  Peps_1p=Peps_1*1.0
  Peps_2p=Peps_2*1.0

  Peps_1pd=Peps_1p*1.0
  Peps_1pd.transpose()
  Peps_1pd.setLabel([35,34,32,33,37])

  Peps_2pd=Peps_2p*1.0
  Peps_2pd.transpose()
  Peps_2pd.setLabel([45,44,43,42,47])


  Swap3p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
  Swap3p.setLabel([38,36,34,35])
  Swap4p=fermionicOPT(Sys,Peps_2pd.bond(1), Peps_2pd.bond(0))
  Swap4p.setLabel([23,41,48,36])


  Swap5p=fermionicOPT(Sys,Peps_2p.bond(1), Peps_2p.bond(0))
  Swap5p.setLabel([42,43,-10,41])

  Swap10p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
  Swap10p.setLabel([25,59,57,54])





 if  Sys[8]=="TEBD_SVD":
   A=Peps_1*1.0
   A.setLabel([1,2,3,4,5])
   PEPS_A=A*1.0

   PEPS_B=Peps_2*1.0
   PEPS_B.setLabel([6,5,8,9,10])
   PEPS_B.permute([6,5,8,9,10],5)
   Swap=fermionicOPT(Sys,PEPS_B.bond(1), PEPS_B.bond(2))
   Swap.setLabel([-5,-8,5,8])
   PEPS_B=Swap*PEPS_B
   PEPS_B.permute([6,-5,-8,9,10],3)
   PEPS_B.setLabel([4,7,8,9,10])


   Ham.setLabel([11,12,3,8])


   Theta=(PEPS_A*Ham)*PEPS_B
   Theta.permute([1,2,11,5,7,12,9,10],4)
   #print len(D)
   U, V, S=TU.setTruncation(Theta, len(D))
   U.setLabel([1,2,11,5,-100])
   V.setLabel([-100,7,12,9,10])
   S=Sqrt(S)
   S.setLabel([-100,100])
   U=U*S
   U.permute([1,2,11,100,5],3)
   U.setLabel([6,31,58,59,10])

   S.setLabel([100,-100])
   V=V*S
   V.permute([100,7,12,9,10],5)
   Swap=fermionicOPT(Sys,V.bond(1), V.bond(2))
   Swap.setLabel([-7,-12,7,12])
   V=Swap*V
   V.permute([100,-7,-12,9,10],3)
   V.setLabel([54,52,51,53,15])


   Peps_1p=U*1.0
   Peps_2p=V*1.0
   Peps_1pd=Peps_1p*1.0
   Peps_2pd=Peps_2p*1.0

   Peps_1pd=Peps_1p*1.0
   Peps_1pd.transpose()
   Peps_1pd.setLabel([35,34,32,33,37])
   Peps_2pd=Peps_2p*1.0
   Peps_2pd.transpose()
   Peps_2pd.setLabel([45,44,43,42,47])

   Swap3p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap3p.setLabel([38,36,34,35])
   Swap4p=fermionicOPT(Sys,Peps_2p.bond(1), Peps_2p.bond(0))
   Swap4p.setLabel([23,41,48,36])


   Swap5p=fermionicOPT(Sys,Peps_2p.bond(1), Peps_2p.bond(0))
   Swap5p.setLabel([42,43,-10,41])

   Swap10p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap10p.setLabel([25,59,57,54])

   Peps_1p.permute([6,31,58,59,10],5)
   Peps_2p.permute([54,52,51,53,15],5)
   Peps_1pd.permute([32,33,37,35,34],0)
   Peps_2pd.permute([43,42,47,45,44],0)






 Ham.setLabel([51,52,53,54])

 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([40,49,56,55])

 Ham.setLabel([40,49,56,55])
 
 iden_h=Ham*1.0
 iden_h.setLabel([40,49,56,55])
 iden_h.identity()

 iden_h.setLabel([40,49,56,55])







 Peps_1p, Peps_1pd, Peps_2p, Peps_2pd=Do_optimization_Grad_row( Peps_1, Peps_1d, Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Sys,H_orig,Swap3p,Swap4p,Swap5p,Swap10p)




 #valf=norm_f_val_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left)
 

 E_energy, Norm_val=Energy_row_grad( Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,H_orig,iden_h,Swap1,Swap2,Swap3p,Swap4p,Swap5p,Swap6,Swap7,Swap8,Swap9,Swap10p,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Sys)



 #print "E_2r", Norm_val, E_energy

 if Norm_val<0:mps_boundry_down[2*Location+2]=mps_boundry_down[2*Location+2]*(-1.0)

 if Norm_val> threshold[1] or Norm_val<threshold[0]:
  Peps_1p, Peps_2p=renormalize_tensor_r_single(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3p,Swap4p,Swap5p,Swap6,Swap7,Swap8,Swap9,Swap10p,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,threshold,interval)


 return Peps_1p, Peps_2p, E_energy







def Update_twotensor_local_row_inv(PEPS_listten, E_list_right, E_list_left, mps_boundry_up, mps_boundry_down, Location, Ham, H_orig, D,threshold, interval,Sys):

 bdi_mid=uni10.Bond(uni10.BD_IN,1)
 iden=uni10.UniTensor([bdi_mid])
 iden.identity()

 Peps_1=PEPS_listten[Location]*1.0
 Peps_2=PEPS_listten[Location+1]*1.0


 E_left=0
 if Location==0:
  iden.setLabel([16])
  Resul=copy.copy(iden)
  iden.setLabel([6])
  Resul=Resul*iden
  iden.setLabel([17])
  Resul=Resul*iden
  iden.setLabel([18])
  Resul=Resul*iden
  E_left=Resul
 else:
  E_left= copy.copy(E_list_left[Location-1])
  E_left.setLabel([16,6,17,18])


 E_right=0
 if Location>=(len(PEPS_listten)-2):
  iden.setLabel([21])
  Resul=copy.copy(iden)
  iden.setLabel([20])
  Resul=Resul*iden
  iden.setLabel([-14])
  Resul=Resul*iden
  iden.setLabel([19])
  Resul=Resul*iden
  E_right=Resul
 else:
  E_right=copy.copy(E_list_right[Location+2])
  E_right.setLabel([21,20,-14,19])


 Peps_1.setLabel([6,31,58,59,10])
 Peps_1d=Peps_1*1.0
 Peps_1d.transpose()
 Peps_1d.setLabel([35,34,32,33,37])
 Peps_2.setLabel([54,52,51,53,15])
 Peps_2d=Peps_2*1.0
 Peps_2d.transpose()
 Peps_2d.setLabel([45,44,43,42,47])


 Swap1=fermionicOPT(Sys,Peps_1.bond(1), Peps_1.bond(0))
 Swap1.setLabel([22,30,31,17])
 Swap2=fermionicOPT(Sys,Peps_1.bond(1), Peps_1.bond(0))
 Swap2.setLabel([33,32,-7,30])
 Swap3=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap3.setLabel([38,36,34,35])
 Swap4=fermionicOPT(Sys,Peps_2.bond(1), Peps_2.bond(0))
 Swap4.setLabel([23,41,48,36])
 Swap5=fermionicOPT(Sys,Peps_2.bond(1), Peps_2.bond(0))
 Swap5.setLabel([42,43,-10,41])
 Swap6=fermionicOPT(Sys,Peps_2d.bond(1), Peps_2d.bond(0))
 Swap6.setLabel([46,-14,44,45])

 Swap7=fermionicOPT(Sys,Peps_2d.bond(1), Peps_2d.bond(0))
 Swap7.setLabel([24,53,46,20])

 Swap8=fermionicOPT(Sys,Peps_1.bond(2), Peps_1d.bond(1))
 Swap8.setLabel([37,39,40,38])

 Swap9=fermionicOPT(Sys,Peps_1.bond(2), Peps_1d.bond(1))
 Swap9.setLabel([56,57,58,39])

 Swap10=fermionicOPT(Sys,Peps_1d.bond(1), Peps_1d.bond(0))
 Swap10.setLabel([25,59,57,54])

 Swap11=fermionicOPT(Sys,Peps_2.bond(2), Peps_2.bond(1))
 Swap11.setLabel([55,50,51,52])

 Swap12=fermionicOPT(Sys,Peps_2.bond(2), Peps_2.bond(1))
 Swap12.setLabel([47,48,49,50])



######################################################################

 mps_boundry_up[Location*2].setLabel([16,10,266])
 mps_boundry_up[Location*2+1].setLabel([266,25,277])
 mps_boundry_up[Location*2+2].setLabel([277,15,288])
 mps_boundry_up[Location*2+3].setLabel([288,24,21])

 mps_boundry_down[Location*2].setLabel([18,22,311])
 mps_boundry_down[Location*2+1].setLabel([311,-7,322])
 mps_boundry_down[Location*2+2].setLabel([322,23,333])
 mps_boundry_down[Location*2+3].setLabel([333,-10,19])

######################################################

# A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

# A=A*(Peps_1*Swap1)
# A=A*mps_boundry_down[2*Location]
# A=A*mps_boundry_down[2*Location+1]

# A=A*((Swap3*Peps_1d)*Swap2)
# A=(A*(Swap8*Swap9))*Swap10

# B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

# B=Swap7*B
# B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

# B=B*((Peps_2d*Swap6)*Swap5)
# B=(B*Swap4)*Swap12
# B=B*(Peps_2*Swap11)

# N_ten=A*B
# 
# N_ten.permute([40,49,56,55],2)



# H_orig.setLabel([40,49,56,55])
# iden_h=H_orig*1.0
# iden_h.identity()
# iden_h.setLabel([40,49,56,55])
# Norm_h=N_ten*iden_h
# h_h=N_ten*H_orig
# print "E_1r", Norm_h[0], h_h[0]/Norm_h[0]



# Peps_1p=Peps_1*1.0
# Peps_2p=Peps_2*1.0
# 
# Peps_1pd=Peps_1p*1.0
# Peps_1pd.setLabel([32,33,37,35,34])
# Peps_2pd=Peps_2p*1.0
# Peps_2pd.setLabel([43,42,47,45,44])





 if  Sys[8]=="Previous":


  Peps_1p=Peps_1*1.0
  Peps_2p=Peps_2*1.0

  Peps_1pd=Peps_1p*1.0
  Peps_1pd.transpose()
  Peps_1pd.setLabel([35,34,32,33,37])

  Peps_2pd=Peps_2p*1.0
  Peps_2pd.transpose()
  Peps_2pd.setLabel([45,44,43,42,47])


  Swap3p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
  Swap3p.setLabel([38,36,34,35])
  Swap4p=fermionicOPT(Sys,Peps_2pd.bond(1), Peps_2pd.bond(0))
  Swap4p.setLabel([23,41,48,36])


  Swap5p=fermionicOPT(Sys,Peps_2p.bond(1), Peps_2p.bond(0))
  Swap5p.setLabel([42,43,-10,41])

  Swap10p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
  Swap10p.setLabel([25,59,57,54])





 if  Sys[8]=="TEBD_SVD":
   A=Peps_1*1.0
   A.setLabel([1,2,3,4,5])
   PEPS_A=A*1.0

   PEPS_B=Peps_2*1.0
   PEPS_B.setLabel([6,5,8,9,10])
   PEPS_B.permute([6,5,8,9,10],5)
   Swap=fermionicOPT(Sys,PEPS_B.bond(1), PEPS_B.bond(2))
   Swap.setLabel([-5,-8,5,8])
   PEPS_B=Swap*PEPS_B
   PEPS_B.permute([6,-5,-8,9,10],3)
   PEPS_B.setLabel([4,7,8,9,10])


   Ham.setLabel([11,12,3,8])


   Theta=(PEPS_A*Ham)*PEPS_B
   Theta.permute([1,2,11,5,7,12,9,10],4)
   #print len(D)
   U, V, S=TU.setTruncation(Theta, len(D))
   U.setLabel([1,2,11,5,-100])
   V.setLabel([-100,7,12,9,10])
   S=Sqrt(S)
   S.setLabel([-100,100])
   U=U*S
   U.permute([1,2,11,100,5],3)
   U.setLabel([6,31,58,59,10])

   S.setLabel([100,-100])
   V=V*S
   V.permute([100,7,12,9,10],5)
   Swap=fermionicOPT(Sys,V.bond(1), V.bond(2))
   Swap.setLabel([-7,-12,7,12])
   V=Swap*V
   V.permute([100,-7,-12,9,10],3)
   V.setLabel([54,52,51,53,15])


   Peps_1p=U*1.0
   Peps_2p=V*1.0
   Peps_1pd=Peps_1p*1.0
   Peps_2pd=Peps_2p*1.0

   Peps_1pd=Peps_1p*1.0
   Peps_1pd.transpose()
   Peps_1pd.setLabel([35,34,32,33,37])
   Peps_2pd=Peps_2p*1.0
   Peps_2pd.transpose()
   Peps_2pd.setLabel([45,44,43,42,47])

   Swap3p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap3p.setLabel([38,36,34,35])
   Swap4p=fermionicOPT(Sys,Peps_2p.bond(1), Peps_2p.bond(0))
   Swap4p.setLabel([23,41,48,36])


   Swap5p=fermionicOPT(Sys,Peps_2p.bond(1), Peps_2p.bond(0))
   Swap5p.setLabel([42,43,-10,41])

   Swap10p=fermionicOPT(Sys,Peps_1pd.bond(1), Peps_1pd.bond(0))
   Swap10p.setLabel([25,59,57,54])




 Ham.setLabel([51,52,53,54])

 H1=Ham*1.0
 H1.transpose()
 H1.setLabel([-20,-40,51,52])
 H=Ham*H1
 H.permute([-20,-40,53,54],2)
 H.setLabel([-2,-6,2,6])
 H.setLabel([40,49,56,55])

 Ham.setLabel([40,49,56,55])
 
 iden_h=Ham*1.0
 iden_h.setLabel([40,49,56,55])
 iden_h.identity()

 iden_h.setLabel([40,49,56,55])

 Peps_1p_init=Peps_1p*1.0
 Peps_1pd_init=Peps_1pd*1.0
 Peps_2p_init=Peps_2p*1.0
 Peps_2pd_init=Peps_2pd*1.0

 valf, E_pos=norm_f_val_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)
 
 E_1=200
 E_2=0
 E_0=0
 count=1
 for i in xrange(Sys[3]):
  count=count+1
  E_2=E_1*1.0
  val, E_pos=norm_f_val_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)
  if i==0: E_0=val
  E_1=val
  #print i, val, abs((E_1-E_2)/E_1)
  if E_1>E_2 or abs(E_1)<1.0e-12  or abs((E_1-E_2)/E_1)<+1.0e-12:
   #print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   Peps_1p=Peps_1p_init*1.0
   Peps_1pd=Peps_1pd_init*1.0
   Peps_2p=Peps_2p_init*1.0
   Peps_2pd=Peps_2pd_init*1.0
   break
  else:
   Peps_1p_init=Peps_1p*1.0
   Peps_1pd_init=Peps_1pd*1.0
   Peps_2p_init=Peps_2p*1.0
   Peps_2pd_init=Peps_2pd*1.0


  Peps_1p, Peps_1pd=optimum_0_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)

  Peps_2p, Peps_2pd=optimum_1_inv_row_single( Peps_1, Peps_1d,Peps_2, Peps_2d,Peps_1p, Peps_1pd,Peps_2p, Peps_2pd,Ham,iden_h,H,Swap1,Swap2,Swap3,Swap4,Swap5,Swap6,Swap7,Swap8,Swap9,Swap10,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,Swap3p,Swap4p,Swap5p,Swap10p)








 if abs(valf)<1.0e-11 or abs(valf)>1.0e+11:
  print "warning_norm_in_optimization",  abs(valf), "count", count



 A=E_left*(mps_boundry_up[2*Location]*mps_boundry_up[2*Location+1])

 A=A*(Peps_1p*Swap1)
 A=A*(mps_boundry_down[2*Location]*mps_boundry_down[2*Location+1])

 A=A*((Swap3p*Peps_1pd)*Swap2)
 A=(A*(Swap8*Swap9))*Swap10p

 B=E_right*(mps_boundry_up[2*Location+2]*mps_boundry_up[2*Location+3])

 B=Swap7*B
 B=B*(mps_boundry_down[2*Location+2]*mps_boundry_down[2*Location+3])

 B=B*((Peps_2pd*Swap6)*Swap5p)
 B=(B*Swap4p)*Swap12
 B=B*(Peps_2p*Swap11)

 N_ten=A*B


 H_orig.setLabel([40,49,56,55])
 iden_h=H_orig*1.0
 iden_h.identity()
 iden_h.setLabel([40,49,56,55])
 Norm_h=N_ten*iden_h
 h_h=N_ten*H_orig
 #print "E_2r", Norm_h[0], h_h[0]/Norm_h[0]


 if Norm_h[0]> threshold[1] or Norm_h[0]<threshold[0]:
  Peps_1p, Peps_2p=renormalize_tensor_r_single(Peps_1p,Peps_1pd,Peps_2p, Peps_2pd,iden_h,Swap1,Swap2,Swap3p,Swap4p,Swap5p,Swap6,Swap7,Swap8,Swap9,Swap10p,Swap11,Swap12,mps_boundry_up,mps_boundry_down,Location,E_right,E_left,threshold,interval)


 return Peps_1p, Peps_2p, h_h[0]/Norm_h[0]







####@profile
def solve_linear_eq(A,Ap):
 #print A
 Ap_h=copy.copy(Ap)
 Ap_h.transpose()
 Result=uni10.UniTensor(Ap.bond())
 blk_qnums = A.blockQnum()
 #print blk_qnums, Ap.printDiagram()
 blk_qnums1 = Ap.blockQnum()
 for qnum in blk_qnums:
   if qnum in blk_qnums1:
    A_mat=A.getBlock(qnum)
    Ap_mat=Ap.getBlock(qnum)
    A_np=Mat_uni_to_np(A_mat)
    b_np=Mat_uni_to_np(Ap_mat)
    #x_np=np.linalg.lstsq(A_np, b_np, rcond=-1)[0] 
    x_np=sp.linalg.lstsq(A_np, b_np,cond=1.e-12,overwrite_a=False,overwrite_b=False,check_finite=True, lapack_driver='gelsy')[0] 
    #x_np=sp.linalg.lstsq(A_np, b_np,cond=1.e-12,overwrite_a=False,overwrite_b=False,check_finite=True)[0] 

    x=Mat_np_to_Uni(x_np)
    Result.putBlock(qnum, x)
 return Result






def  inint_abcs(a_u, A):

 try:
  blk_qnums = a_u.blockQnum()
  blk_qnumsf = A.blockQnum()

  for qnum in blk_qnums:
   if qnum in blk_qnumsf:
    mat_t=a_u.getBlock(qnum)
    dx=int(mat_t.row())
    dy=int(mat_t.col())
    mat_f=A.getBlock(qnum)
    a_u.putBlock( qnum, mat_f.resize(dx,dy) )
  return a_u*1.0
 except:
   print "error, init_abcd"
   return a_u*1.0







##@profile
def optimiz_d( A_f, a, b, c, d, Sys):

 c=c*1.0
 b=b*1.0
 d=d*1.0
 a=a*1.0

 a.permute([5000,10,2,9000,1],3)
 c.permute([500,1,7,900,40],3)
 
 Swap=fermionicOPT(Sys, c.bond(0), c.bond(1))
 Swap.setLabel([-500,-1,500,1])
 c=c*Swap
 c.permute([-1,-500,7,900,40],3)
 c.setLabel([1,500,7,900,40])

 T_tem=a*c
 Swap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 Swap1.setLabel([-5000,18,5000,-500])
 T_temP=T_tem*Swap1 
 Swap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 Swap2.setLabel([-10,17,10,18])
 T_temP=T_temP*Swap2
 Swap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 Swap3.setLabel([-2,16,2,17])
 T_temP=T_temP*Swap3
 Swap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 Swap4.setLabel([19,500,9000,16])
 T_temP=T_temP*Swap4
 Swap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 Swap5.setLabel([-9000,7,19,-7])
 T_temP=T_temP*Swap5
 Swap11=fermionicOPT(Sys, T_tem.bond(7), T_tem.bond(6))
 Swap11.setLabel([40,900,21,-900])
 T_temP=T_temP*Swap11 
 Swap22=fermionicOPT(Sys, T_tem.bond(7), T_tem.bond(3))
 Swap22.setLabel([21,-90000,-40,-9000])
 T_temP=T_temP*Swap22
 T_temP.permute([-500,-5000,-10,-2,-7,-40,-90000,-900],4)
 T_temP.setLabel([500,5000,10,2,7,40,9000,900])
 T_tem_ac=T_temP*1.0



 b.permute([9000,100,5,800,1],3)
 d.permute([900,1,9,80,400],3)

 Swap=fermionicOPT(Sys, d.bond(0), d.bond(1))
 Swap.setLabel([-900,-1,900,1])
 d=d*Swap
 d.permute([-1,-900,9,80,400],3)
 d.setLabel([1,900,9,80,400])

 T_tem=d*b
 T_tem.permute([9000,100,5,800,900,9,80,400],4)
 OSwap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 OSwap1.setLabel([-9000,18,9000,-900])
 T_temP=T_tem*OSwap1 
 OSwap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 OSwap2.setLabel([-100,17,100,18])
 T_temP=T_temP*OSwap2
 OSwap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 OSwap3.setLabel([-5,16,5,17])
 T_temP=T_temP*OSwap3
 OSwap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 OSwap4.setLabel([15,900,800,16])
 T_temP=T_temP*OSwap4
 OSwap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 OSwap5.setLabel([-800,9,15,-9])
 T_temP=T_temP*OSwap5
 T_temP.permute([-900,-9000,-100,-5,-9,-800,80,400],4)
 T_tem_bd=T_temP*1.0
 T_tem_bd.setLabel([900,9000,100,5,9,800,80,400])


 Swap1T=OSwap1*1.0
 Swap2T=OSwap2*1.0
 Swap3T=OSwap3*1.0
 Swap4T=OSwap4*1.0
 Swap5T=OSwap5*1.0

 d.setLabel([1,2900,209,80,400])
 OSwap1.setLabel([9000,18,29000,900])
 OSwap2.setLabel([100,17,2100,18])
 OSwap3.setLabel([5,16,205,17])
 OSwap4.setLabel([15,2900,2800,16])
 OSwap5.setLabel([800,209,15,9])

 dT=d*1.0
 dT.transpose()
 Swap1T.transpose()
 Swap2T.transpose()
 Swap3T.transpose()
 Swap4T.transpose()
 Swap5T.transpose()


 dT.setLabel([80,400,-1,-2900,-209])
 Swap1T.setLabel([-29000,-900, -9000,-18])
 Swap2T.setLabel([-2100,-18, 100,-17])
 Swap3T.setLabel([-205,-17,5,-16])
 Swap4T.setLabel([-2800,-16,-15,-2900])
 Swap5T.setLabel([-15,9, 800,-209])


 b.setLabel([29000,2100,205,2800,1])
 bT=b*1.0
 bT.transpose()
 bT.setLabel([-2800,-1,-29000,-2100,-205])


 T_tem_ac.setLabel([500,5000,10,2,7,40,9000,900])
 T_tem_ac_conj=T_tem_ac*1.0
 T_tem_ac_conj.transpose()
 T_tem_ac_conj.setLabel([7,40,-9000,-900, 500,5000,10,2])
 double_ac=T_tem_ac_conj*T_tem_ac
 double_ac.permute([-900,-9000,900,9000],2)


 T=(OSwap1*double_ac)
 T=(Swap1T*T)
 T=(b*T)
 T=(bT*T)
 T=T*(OSwap2*Swap2T)
 T=T*(OSwap3*Swap3T)
 T=T*(OSwap4*Swap4T)
 T=OSwap5*T
 T=Swap5T*T

 bdi = uni10.Bond(uni10.BD_IN, d.bond(3).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, d.bond(3).Qlist())
 Iden1=uni10.UniTensor([bdi, bdo])
 Iden1.setLabel([280,-280])
 Iden1.identity()

 bdi = uni10.Bond(uni10.BD_IN, d.bond(4).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, d.bond(4).Qlist())
 Iden2=uni10.UniTensor([bdi, bdo])
 Iden2.setLabel([2400,-2400])
 Iden2.identity()

 T=T*Iden1
 T=T*Iden2

 T.permute([1,2900,209,280,2400, -1,-2900,-209,-280,-2400],5)
 #print T.printDiagram()
 
 A_f.permute([500,5000,10,2,7,40,100,5,9,800,80,400],6)
 A_f_t=A_f*1.0
 A_f_t.transpose()
 A_f_t.setLabel([100,5,9,800,80,400,500,5000,10,2,7,40])

 TT=A_f_t*T_tem_ac
 TT=TT*(OSwap5)
 Tp=b*(OSwap1)
 Tp=Tp*(OSwap2)
 Tp=Tp*(OSwap3)
 Tp=Tp*(OSwap4)
 Tp=TT*Tp

 Tp.permute([1,2900,209,80,400],5)
 Tp.setLabel([1,2900,209,280,2400])
 #print Tp.printDiagram()

 if Sys[4]=="SVD":
   row, colm=cal_rowcol(T)
   if (row<=colm):
    U,V,S=TU.setTruncation(T,row)
   else:
    U,V,S=TU.setTruncation(T,colm)

   U.transpose()
   V.transpose()
   S=inverse(S)

   U.setLabel([7,8,9,10,11,12])
   V.setLabel([1,2,3,4,5,6])
   S.setLabel([6,7])

   A2_inv=V*S*U
   A2_inv.permute([1,2,3,4,5, 8,9,10,11,12],5)
   A2_inv.setLabel([-1,-2900,-209,-280,-2400,1,2900,209,280,2400])

   dp=A2_inv*Tp
   dp.permute([-1,-2900,-209,-280,-2400],3)
   dp.setLabel([1,900,9,80,400])

   Swap=fermionicOPT(Sys, dp.bond(0), dp.bond(1))
   Swap.setLabel([-1,-900,1,900])
   dp=dp*Swap
   dp.permute([-900,-1,9,80,400],3)
   dp.setLabel([900,1,9,80,400])
   dp.transpose()
   dp.setLabel([80,400,900,1,9])
   dp.permute([900,1,9,80,400],3)

 elif  Sys[4]=="cg":
  try:
   dp=solve_linear_eq(T,Tp)
   dp.setLabel([1,900,9,80,400])
   Swap=fermionicOPT(Sys, dp.bond(0), dp.bond(1))
   Swap.setLabel([-1,-900,1,900])
   dp=dp*Swap
   dp.permute([-900,-1,9,80,400],3)
   dp.setLabel([900,1,9,80,400])
   dp.transpose()
   dp.setLabel([80,400,900,1,9])
   dp.permute([900,1,9,80,400],3)
  except:
   print "error, SVD"
   row, colm=cal_rowcol(T)
   if (row<=colm):
    U,V,S=TU.setTruncation(T,row)
   else:
    U,V,S=TU.setTruncation(T,colm)

   U.transpose()
   V.transpose()
   S=inverse(S)

   U.setLabel([7,8,9,10,11,12])
   V.setLabel([1,2,3,4,5,6])
   S.setLabel([6,7])

   A2_inv=V*S*U
   A2_inv.permute([1,2,3,4,5, 8,9,10,11,12],5)
   A2_inv.setLabel([-1,-2900,-209,-280,-2400,1,2900,209,280,2400])

   dp=A2_inv*Tp
   dp.permute([-1,-2900,-209,-280,-2400],3)
   dp.setLabel([1,900,9,80,400])

   Swap=fermionicOPT(Sys, dp.bond(0), dp.bond(1))
   Swap.setLabel([-1,-900,1,900])
   dp=dp*Swap
   dp.permute([-900,-1,9,80,400],3)
   dp.setLabel([900,1,9,80,400])
   dp.setLabel([900,1,9,80,400])
   dp.transpose()
   dp.setLabel([80,400,900,1,9])
   dp.permute([900,1,9,80,400],3)

 return dp




##@profile
def optimiz_b( A_f, a, b, c, d , Sys):

 c=c*1.0
 b=b*1.0
 d=d*1.0
 a=a*1.0

 a.permute([5000,10,2,9000,1],3)
 c.permute([500,1,7,900,40],3)
 
 Swap=fermionicOPT(Sys, c.bond(0), c.bond(1))
 Swap.setLabel([-500,-1,500,1])
 c=c*Swap
 c.permute([-1,-500,7,900,40],3)
 c.setLabel([1,500,7,900,40])

 T_tem=a*c
 OSwap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 OSwap1.setLabel([-5000,18,5000,-500])
 T_temP=T_tem*OSwap1 
 OSwap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 OSwap2.setLabel([-10,17,10,18])
 T_temP=T_temP*OSwap2
 OSwap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 OSwap3.setLabel([-2,16,2,17])
 T_temP=T_temP*OSwap3
 OSwap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 OSwap4.setLabel([19,500,9000,16])
 T_temP=T_temP*OSwap4
 OSwap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 OSwap5.setLabel([-9000,7,19,-7])
 T_temP=T_temP*OSwap5
 OSwap11=fermionicOPT(Sys, T_tem.bond(7), T_tem.bond(6))
 OSwap11.setLabel([40,900,21,-900])
 T_temP=T_temP*OSwap11 
 OSwap22=fermionicOPT(Sys, T_tem.bond(7), T_tem.bond(3))
 OSwap22.setLabel([21,-90000,-40,-9000])
 T_temP=T_temP*OSwap22
 T_temP.permute([-500,-5000,-10,-2,-7,-40,-90000,-900],4)
 T_temP.setLabel([500,5000,10,2,7,40,9000,900])
 T_tem_ac=T_temP*1.0


 b.permute([9000,100,5,800,1],3)
 d.permute([900,1,9,80,400],3)

 Swap=fermionicOPT(Sys, d.bond(0), d.bond(1))
 Swap.setLabel([-900,-1,900,1])
 d=d*Swap
 d.permute([-1,-900,9,80,400],3)
 d.setLabel([1,900,9,80,400])


 T_tem=d*b
 T_tem.permute([9000,100,5,800,900,9,80,400],4)
 OSwap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 OSwap1.setLabel([-9000,18,9000,-900])
 T_temP=T_tem*OSwap1 
 OSwap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 OSwap2.setLabel([-100,17,100,18])
 T_temP=T_temP*OSwap2
 OSwap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 OSwap3.setLabel([-5,16,5,17])
 T_temP=T_temP*OSwap3
 OSwap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 OSwap4.setLabel([15,900,800,16])
 T_temP=T_temP*OSwap4
 OSwap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 OSwap5.setLabel([-800,9,15,-9])
 T_temP=T_temP*OSwap5
 T_temP.permute([-900,-9000,-100,-5,-9,-800,80,400],4)
 T_tem_bd=T_temP*1.0
 T_tem_bd.setLabel([900,9000,100,5,9,800,80,400])


 Swap1T=OSwap1*1.0
 Swap2T=OSwap2*1.0
 Swap3T=OSwap3*1.0
 Swap4T=OSwap4*1.0
 Swap5T=OSwap5*1.0

 d.setLabel([1,2900,209,80,400])
 OSwap1.setLabel([9000,18,29000,900])
 OSwap2.setLabel([100,17,2100,18])
 OSwap3.setLabel([5,16,205,17])
 OSwap4.setLabel([15,2900,2800,16])
 OSwap5.setLabel([800,209,15,9])

 dT=d*1.0
 dT.transpose()
 Swap1T.transpose()
 Swap2T.transpose()
 Swap3T.transpose()
 Swap4T.transpose()
 Swap5T.transpose()


 dT.setLabel([80,400,-1,-2900,-209])
 Swap1T.setLabel([-29000,-900, -9000,-18])
 Swap2T.setLabel([-2100,-18, 100,-17])
 Swap3T.setLabel([-205,-17,5,-16])
 Swap4T.setLabel([-2800,-16,-15,-2900])
 Swap5T.setLabel([-15,9, 800,-209])


 T_tem_ac.setLabel([500,5000,10,2,7,40,9000,900])
 T_tem_ac_conj=T_tem_ac*1.0
 T_tem_ac_conj.transpose()
 T_tem_ac_conj.setLabel([7,40,-9000,-900, 500,5000,10,2])
 double_ac=T_tem_ac_conj*T_tem_ac
 double_ac.permute([-900,-9000,900,9000],2)



 T=OSwap5*d
 T1=Swap5T*dT
 T=T*T1
 T=T*(OSwap4*Swap4T)
 T=T*(OSwap3*Swap3T)
 T=T*(OSwap2*Swap2T)
 T=T*((OSwap1*Swap1T)*double_ac)
 
 T.permute([29000,2100,205,2800,1,-29000,-2100,-205,-2800,-1],5)



 A_f.permute([500,5000,10,2,7,40,100,5,9,800,80,400],6)
 A_f_t=A_f*1.0
 A_f_t.transpose()
 A_f_t.setLabel([100,5,9,800,80,400,500,5000,10,2,7,40])


 Tp=A_f_t*T_tem_ac
 TT=d*(OSwap5)
 Tp=TT*Tp
 Tp=Tp*(OSwap4)
 Tp=Tp*(OSwap3)
 Tp=Tp*(OSwap2)
 Tp=Tp*(OSwap1)

 Tp.permute([29000,2100,205,2800,1],5)


 if Sys[4]=="SVD":
  row, colm=cal_rowcol(T)
  if (row<=colm):
   U,V,S=TU.setTruncation(T,row)
  else:
   U,V,S=TU.setTruncation(T,colm)

  U.transpose()
  V.transpose()
  S=inverse(S)
  
  U.setLabel([7,8,9,10,11,12])
  V.setLabel([1,2,3,4,5,6])
  S.setLabel([6,7])

  A2_inv=V*S*U
  A2_inv.permute([1,2,3,4,5, 8,9,10,11,12],5)
  A2_inv.setLabel([-29000,-2100,-205,-2800,-1,29000,2100,205,2800,1])

  bp=A2_inv*Tp
  bp.permute([-29000,-2100,-205,-2800,-1],3)
  bp.setLabel([9000,100,5,800,1])
  bp.permute([9000,100,5,800,1],3)
  bp.transpose()
  bp.setLabel([800,1,9000,100,5])
  bp.permute([9000,100,5,800,1],3)



 elif  Sys[4]=="cg":
  try:
   bp=solve_linear_eq(T,Tp)
   bp.setLabel([9000,100,5,800,1])
   bp.permute([9000,100,5,800,1],3)
   bp.transpose()
   bp.setLabel([800,1,9000,100,5])
   bp.permute([9000,100,5,800,1],3)
  except:
   print "error, SVD"
   row, colm=cal_rowcol(T)
   if (row<=colm):
    U,V,S=TU.setTruncation(T,row)
   else:
    U,V,S=TU.setTruncation(T,colm)

   U.transpose()
   V.transpose()
   S=inverse(S)
  
   U.setLabel([7,8,9,10,11,12])
   V.setLabel([1,2,3,4,5,6])
   S.setLabel([6,7])

   A2_inv=V*S*U
   A2_inv.permute([1,2,3,4,5, 8,9,10,11,12],5)
   A2_inv.setLabel([-29000,-2100,-205,-2800,-1,29000,2100,205,2800,1])

   bp=A2_inv*Tp
   bp.permute([-29000,-2100,-205,-2800,-1],3)
   bp.setLabel([9000,100,5,800,1])
   bp.permute([9000,100,5,800,1],3)
   bp.transpose()
   bp.setLabel([800,1,9000,100,5])
   bp.permute([9000,100,5,800,1],3)

 return bp





def optimiz_c( A_f, a, b, c, d, Sys):

 c=c*1.0
 b=b*1.0
 d=d*1.0
 a=a*1.0

 a.permute([5000,10,2,9000,1],3)
 c.permute([500,1,7,900,40],3)
 
 Swap=fermionicOPT(Sys, c.bond(0), c.bond(1))
 Swap.setLabel([-500,-1,500,1])
 c=c*Swap
 c.permute([-1,-500,7,900,40],3)
 c.setLabel([1,500,7,900,40])

 T_tem=a*c

 OSwap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 OSwap1.setLabel([-5000,18,5000,-500])
 T_temP=T_tem*OSwap1 

 OSwap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 OSwap2.setLabel([-10,17,10,18])
 T_temP=T_temP*OSwap2

 OSwap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 OSwap3.setLabel([-2,16,2,17])
 T_temP=T_temP*OSwap3

 OSwap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 OSwap4.setLabel([19,500,9000,16])
 T_temP=T_temP*OSwap4

 OSwap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 OSwap5.setLabel([-9000,7,19,-7])
 T_temP=T_temP*OSwap5

 OSwap11=fermionicOPT(Sys, T_tem.bond(7), T_tem.bond(6))
 OSwap11.setLabel([40,900,21,-900])
 T_temP=T_temP*OSwap11 

 OSwap22=fermionicOPT(Sys, T_tem.bond(7), T_tem.bond(3))
 OSwap22.setLabel([21,-90000,-40,-9000])
 T_temP=T_temP*OSwap22

 T_temP.permute([-500,-5000,-10,-2,-7,-40,-90000,-900],4)
 T_temP.setLabel([500,5000,10,2,7,40,9000,900])
 T_tem_ac=T_temP*1.0


 c.setLabel([1,250,70,1900,1400])
 OSwap1.setLabel([5000,18,25000,500])
 OSwap2.setLabel([10,17,210,18])
 OSwap3.setLabel([2,16,202,17])
 OSwap4.setLabel([19,250,290,16])
 OSwap5.setLabel([29,70,19,7])
 OSwap11.setLabel([1400,1900,32,900])
 OSwap22.setLabel([32,9000,40,29])

 cT=c*1.0
 Swap1T=OSwap1*1.0
 Swap2T=OSwap2*1.0
 Swap3T=OSwap3*1.0
 Swap4T=OSwap4*1.0
 Swap5T=OSwap5*1.0
 Swap11T=OSwap11*1.0
 Swap22T=OSwap22*1.0

 cT.transpose()
 Swap1T.transpose()
 Swap2T.transpose()
 Swap3T.transpose()
 Swap4T.transpose()
 Swap5T.transpose()
 Swap11T.transpose()
 Swap22T.transpose()

 
 cT.setLabel([-1900,-1400,-1,-250,-70])
 Swap1T.setLabel([-25000,500,5000,-18])
 Swap2T.setLabel([-210,-18,10,-17])
 Swap3T.setLabel([-202,-17,2,-16])
 Swap4T.setLabel([-290,-16,-19,-250])
 Swap5T.setLabel([-19,7,-29,-70])
 Swap11T.setLabel([-32,-900,-1400,-1900])
 Swap22T.setLabel([40,-29,-32,-9000])



 a.setLabel([25000,210,202,290,1])
 aT=a*1.0
 aT.transpose()
 aT.setLabel([-290,-1,-25000,-210,-202])



 b.permute([9000,100,5,800,1],3)
 d.permute([900,1,9,80,400],3)

 Swap=fermionicOPT(Sys, d.bond(0), d.bond(1))
 Swap.setLabel([-900,-1,900,1])
 d=d*Swap
 d.permute([-900,-1,9,80,400],3)
 d.setLabel([900,1,9,80,400])

 T_tem=b*d

 T_tem.permute([9000,100,5,800,900,9,80,400],4)
 Swap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 Swap1.setLabel([-9000,18,9000,-900])
 T_temP=T_tem*Swap1 
 Swap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 Swap2.setLabel([-100,17,100,18])
 T_temP=T_temP*Swap2
 Swap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 Swap3.setLabel([-5,16,5,17])
 T_temP=T_temP*Swap3
 Swap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 Swap4.setLabel([15,900,800,16])
 T_temP=T_temP*Swap4
 Swap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 Swap5.setLabel([-800,9,15,-9])
 T_temP=T_temP*Swap5
 T_temP.permute([-900,-9000,-100,-5,-9,-800,80,400],4)
 T_tem_bd=T_temP*1.0
 T_tem_bd.setLabel([900,9000,100,5,9,800,80,400])


 T_tem_bd.setLabel([900,9000,100,5,9,800,80,400])
 T_tem_bd_conj=T_tem_bd*1.0
 T_tem_bd_conj.setLabel([-900,-9000,100,5,9,800,80,400])
 T_tem_bd_conj.transpose()
 T_tem_bd_conj.setLabel([9,800,80,400,-900,-9000,100,5])
 
 double_bd=T_tem_bd_conj*T_tem_bd
 double_bd.permute([-900,-9000,900,9000],2)
 
 
 T=(((a*OSwap1)*OSwap2)*OSwap3)*OSwap4
 T1=(((aT*Swap1T)*Swap2T)*Swap3T)*Swap4T
 T=T1*T
 T1=(double_bd*(OSwap11*OSwap22))
 T1=T1*(Swap11T*Swap22T) 
 T1=T1*(OSwap5)
 T1=T1*(Swap5T)
 T=T*T1
 T.permute([1,250,70,1900,1400,-1,-250,-70,-1900,-1400],5)


 A_f.permute([500,5000,10,2,7,40,100,5,9,800,80,400],6)
 A_f_t=A_f*1.0
 A_f_t.transpose()
 A_f_t.setLabel([100,5,9,800,80,400,500,5000,10,2,7,40])



 Tp=A_f_t*T_tem_bd

 TT=(((a*OSwap1)*OSwap2)*OSwap3)*OSwap4
 TT=TT*Tp
 TT=TT*(OSwap5)
 Tp=(TT*(OSwap11*OSwap22)) 

 Tp.permute([1,250,70,1900,1400],5)


 if Sys[4]=="SVD":
  row, colm=cal_rowcol(T)
  if (row<=colm):
   U,V,S=TU.setTruncation(T,row)
  else:
   U,V,S=TU.setTruncation(T,colm)

  U.transpose()
  V.transpose()
  S=inverse(S)
  
  U.setLabel([7,8,9,10,11,12])
  V.setLabel([1,2,3,4,5,6])
  S.setLabel([6,7])

  A2_inv=V*S*U
  A2_inv.permute([1,2,3,4,5, 8,9,10,11,12],5)
  A2_inv.setLabel([-1,-250,-70,-1900,-1400,1,250,70,1900,1400])
  cp=A2_inv*Tp
  cp.permute([-1,-250,-70,-1900,-1400],3)
  cp.setLabel([1,500,7,900,40])
  Swap=fermionicOPT(Sys, cp.bond(0), cp.bond(1))
  Swap.setLabel([-1,-500,1,500])
  cp=cp*Swap
  cp.permute([-500,-1,7,900,40],3)
  cp.setLabel([500,1,7,900,40])  
  cp.permute([500,1,7,900,40],3)
  cp.transpose()
  cp.setLabel([900,40,500,1,7])
  cp.permute([500,1,7,900,40],3)

 elif  Sys[4]=="cg":
  try:
   cp=solve_linear_eq(T,Tp)
   cp.setLabel([1,500,7,900,40])
   Swap=fermionicOPT(Sys, cp.bond(0), cp.bond(1))
   Swap.setLabel([-1,-500,1,500])
   cp=cp*Swap
   cp.permute([-500,-1,7,900,40],3)
   cp.setLabel([500,1,7,900,40])  
   cp.permute([500,1,7,900,40],3)
   cp.transpose()
   cp.setLabel([900,40,500,1,7])
   cp.permute([500,1,7,900,40],3)

  except:
   print "error, SVD"
   row, colm=cal_rowcol(T)
   if (row<=colm):
    U,V,S=TU.setTruncation(T,row)
   else:
    U,V,S=TU.setTruncation(T,colm)

   U.transpose()
   V.transpose()
   S=inverse(S)
  
   U.setLabel([7,8,9,10,11,12])
   V.setLabel([1,2,3,4,5,6])
   S.setLabel([6,7])

   A2_inv=V*S*U
   A2_inv.permute([1,2,3,4,5, 8,9,10,11,12],5)
   A2_inv.setLabel([-1,-250,-70,-1900,-1400,1,250,70,1900,1400])
   cp=A2_inv*Tp
   cp.permute([-1,-250,-70,-1900,-1400],3)
   cp.setLabel([1,500,7,900,40])
   Swap=fermionicOPT(Sys, cp.bond(0), cp.bond(1))
   Swap.setLabel([-1,-500,1,500])
   cp=cp*Swap
   cp.permute([-500,-1,7,900,40],3)
   cp.setLabel([500,1,7,900,40])  
   cp.permute([500,1,7,900,40],3)
   cp.transpose()
   cp.setLabel([900,40,500,1,7])
   cp.permute([500,1,7,900,40],3)

 return cp





def optimiz_a( A_f, a, b, c, d,Sys):

 c=c*1.0
 b=b*1.0
 d=d*1.0
 a=a*1.0

 a.permute([5000,10,2,9000,1],3)
 c.permute([500,1,7,900,40],3)
 
 Swap=fermionicOPT(Sys, c.bond(0), c.bond(1))
 Swap.setLabel([-500,-1,500,1])
 c=c*Swap
 c.permute([-1,-500,7,900,40],3)
 c.setLabel([1,500,7,900,40])

 T_tem=a*c
 T_tem.permute([5000,10,2,9000,500,7,900,40],4)

 OSwap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 OSwap1.setLabel([-5000,18,5000,-500])
 T_temP=T_tem*OSwap1 

 OSwap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 OSwap2.setLabel([-10,17,10,18])
 T_temP=T_temP*OSwap2

 OSwap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 OSwap3.setLabel([-2,16,2,17])
 T_temP=T_temP*OSwap3

 OSwap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 OSwap4.setLabel([19,500,9000,16])
 T_temP=T_temP*OSwap4

 OSwap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 OSwap5.setLabel([-9000,7,19,-7])
 T_temP=T_temP*OSwap5

 OSwap11=fermionicOPT(Sys, T_tem.bond(7), T_tem.bond(6))
 OSwap11.setLabel([40,900,21,-900])
 T_temP=T_temP*OSwap11 

 OSwap22=fermionicOPT(Sys, T_tem.bond(7), T_tem.bond(3))
 OSwap22.setLabel([21,-90000,-40,-9000])

 T_temP=T_temP*OSwap22

 T_temP.permute([-500,-5000,-10,-2,-7,-40,-90000,-900],4)
 T_temP.setLabel([500,5000,10,2,7,40,9000,900])
 T_tem_ac=T_temP*1.0



 c.setLabel([1,250,70,1900,1400])
 OSwap1.setLabel([5000,18,25000,500])
 OSwap2.setLabel([10,17,210,18])
 OSwap3.setLabel([2,16,202,17])
 OSwap4.setLabel([19,250,290,16])
 OSwap5.setLabel([29,70,19,7])
 OSwap11.setLabel([1400,1900,32,900])
 OSwap22.setLabel([32,9000,40,29])

 cT=c*1.0
 Swap1T=OSwap1*1.0
 Swap2T=OSwap2*1.0
 Swap3T=OSwap3*1.0
 Swap4T=OSwap4*1.0
 Swap5T=OSwap5*1.0
 Swap11T=OSwap11*1.0
 Swap22T=OSwap22*1.0

 cT.transpose()
 Swap1T.transpose()
 Swap2T.transpose()
 Swap3T.transpose()
 Swap4T.transpose()
 Swap5T.transpose()
 Swap11T.transpose()
 Swap22T.transpose()

 
 cT.setLabel([-1900,-1400,-1,-250,-70])
 Swap1T.setLabel([-25000,500,5000,-18])
 Swap2T.setLabel([-210,-18,10,-17])
 Swap3T.setLabel([-202,-17,2,-16])
 Swap4T.setLabel([-290,-16,-19,-250])
 Swap5T.setLabel([-19,7,-29,-70])
 Swap11T.setLabel([-32,-900,-1400,-1900])
 Swap22T.setLabel([40,-29,-32,-9000])



 b.permute([9000,100,5,800,1],3)
 d.permute([900,1,9,80,400],3)

 Swap=fermionicOPT(Sys, d.bond(0), d.bond(1))
 Swap.setLabel([-900,-1,900,1])
 d=d*Swap
 d.permute([-1,-900,9,80,400],3)
 d.setLabel([1,900,9,80,400])


 T_tem=b*d
 T_tem.permute([9000,100,5,800,900,9,80,400],4)
 Swap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 Swap1.setLabel([-9000,18,9000,-900])
 T_temP=T_tem*Swap1 
 Swap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 Swap2.setLabel([-100,17,100,18])
 T_temP=T_temP*Swap2
 Swap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 Swap3.setLabel([-5,16,5,17])
 T_temP=T_temP*Swap3
 Swap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 Swap4.setLabel([15,900,800,16])
 T_temP=T_temP*Swap4
 Swap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 Swap5.setLabel([-800,9,15,-9])
 T_temP=T_temP*Swap5
 T_temP.permute([-900,-9000,-100,-5,-9,-800,80,400],4)
 T_tem_bd=T_temP*1.0
 T_tem_bd.setLabel([900,9000,100,5,9,800,80,400])


 T_tem_bd.setLabel([900,9000,100,5,9,800,80,400])
 T_tem_bd_conj=T_tem_bd*1.0
 T_tem_bd_conj.setLabel([-900,-9000,100,5,9,800,80,400])
 T_tem_bd_conj.transpose()
 T_tem_bd_conj.setLabel([9,800,80,400,-900,-9000,100,5])
 
 double_bd=T_tem_bd_conj*T_tem_bd
 double_bd.permute([-900,-9000,900,9000],2)
 
 
 T=(double_bd*(OSwap11*OSwap22))
 T=T*(Swap11T*Swap22T) 
 T=T*(OSwap5)
 T=T*(Swap5T)
 T=T*c
 T=T*cT
 T=T*(OSwap4*Swap4T)
 T=T*(OSwap3*Swap3T)
 T=T*(OSwap2*Swap2T)
 T=T*(OSwap1*Swap1T)


 T.permute([25000,210,202,290,1,-25000,-210,-202,-290,-1],5)
 
 A_f.permute([500,5000,10,2,7,40,100,5,9,800,80,400],6)
 A_f_t=A_f*1.0
 A_f_t.transpose()
 A_f_t.setLabel([100,5,9,800,80,400,500,5000,10,2,7,40])

 Tp=A_f_t*T_tem_bd
 TT=(c*(OSwap11*OSwap22)) 
 TT=TT*(OSwap5)
 Tp=TT*Tp
 Tp=Tp*(OSwap4)
 Tp=Tp*(OSwap3)
 Tp=Tp*(OSwap2)
 Tp=Tp*(OSwap1)

 Tp.permute([25000,210,202,290,1],5)


 if Sys[4]=="SVD":
  row, colm=cal_rowcol(T)
  if (row<=colm):
   U,V,S=TU.setTruncation(T,row)
  else:
   U,V,S=TU.setTruncation(T,colm)

  U.transpose()
  V.transpose()
  S=inverse(S)
  
  U.setLabel([7,8,9,10,11,12])
  V.setLabel([1,2,3,4,5,6])
  S.setLabel([6,7])

  A2_inv=V*S*U
  A2_inv.permute([1,2,3,4,5, 8,9,10,11,12],3)
  A2_inv.setLabel([-25000,-210,-202,-290,-1,25000,210,202,290,1])
#   A2_inv.transpose()
#   A2_inv.setLabel([25000,210,202,290,1, -25000,-210,-202,-290,-1])

  ap=A2_inv*Tp
  ap.permute([-25000,-210,-202,-290,-1],3)
  ap.setLabel([5000,10,2,9000,1])
  ap.permute([5000,10,2,9000,1],3)
  #print a.printDiagram()
  ap.transpose()
  ap.setLabel([9000,1,5000,10,2])
  ap.permute([5000,10,2,9000,1],3)
  #print ap.printDiagram()

 elif  Sys[4]=="cg":
  try:
   ap=solve_linear_eq(T,Tp)
   ap.setLabel([5000,10,2,9000,1])
   ap.permute([5000,10,2,9000,1],3)
   ap.transpose()
   ap.setLabel([9000,1,5000,10,2])
   ap.permute([5000,10,2,9000,1],3)

  except:
   row, colm=cal_rowcol(T)
   if (row<=colm):
    U,V,S=TU.setTruncation(T,row)
   else:
    U,V,S=TU.setTruncation(T,colm)

   U.transpose()
   V.transpose()
   S=inverse(S)
  
   U.setLabel([7,8,9,10,11,12])
   V.setLabel([1,2,3,4,5,6])
   S.setLabel([6,7])

   A2_inv=V*S*U
   A2_inv.permute([1,2,3,4,5, 8,9,10,11,12],3)
   A2_inv.setLabel([-25000,-210,-202,-290,-1,25000,210,202,290,1])
 #   A2_inv.transpose()
 #   A2_inv.setLabel([25000,210,202,290,1, -25000,-210,-202,-290,-1])

   ap=A2_inv*Tp
   ap.permute([-25000,-210,-202,-290,-1],3)
   ap.setLabel([5000,10,2,9000,1])
   ap.permute([5000,10,2,9000,1],3)
   #print a.printDiagram()
   ap.transpose()
   ap.setLabel([9000,1,5000,10,2])
   ap.permute([5000,10,2,9000,1],3)
   #print ap.printDiagram()



 ap.setLabel([5000,10,2,9000,1])

 return ap





###@profile
def cost_f( A_f, a, b, c, d ,Sys):

 c1=c*1.0
 b1=b*1.0
 d1=d*1.0
 a1=a*1.0

 a1.permute([5000,10,2,9000,1],3)
 c1.permute([1,500,7,900,40],3)
 
 Swap=fermionicOPT(Sys, c1.bond(0), c1.bond(1))
 Swap.setLabel([-1,-500,1,500])
 c1=c1*Swap
 c1.permute([-1,-500,7,900,40],3)
 c1.setLabel([1,500,7,900,40])

 T_tem=a1*c1

 #T_tem.permute([5000,5000,10,2,7,9000,900,40],4)
 T_tem.permute([5000,10,2,9000,500,7,900,40],4)

 Swap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 Swap1.setLabel([-5000,18,+5000,-500])
 T_temP=T_tem*Swap1 

 Swap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 Swap2.setLabel([-10,17,+10,18])
 T_temP=T_temP*Swap2

 Swap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 Swap3.setLabel([-2,16,2,17])
 T_temP=T_temP*Swap3

 Swap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 Swap4.setLabel([19,500,9000,16])
 T_temP=T_temP*Swap4

 Swap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 Swap5.setLabel([-9000,7,19,-7])
 T_temP=T_temP*Swap5


 T_temP.permute([-500,-5000,-10,-2,-7,-9000,900,40],4)
 T_tem_ac=T_temP*1.0
 T_tem_ac.setLabel([500,5000,10,2,7,9000,900,40])


 Swap11=fermionicOPT(Sys, T_tem_ac.bond(7), T_tem_ac.bond(6))
 Swap11.setLabel([40,900,16,-900])
 T_temP=T_tem_ac*Swap11 

 Swap22=fermionicOPT(Sys, T_tem_ac.bond(7), T_tem_ac.bond(5))
 Swap22.setLabel([16,9000,-40,-9000])
 T_temP=T_temP*Swap22
 T_temP.permute([500,5000,10,2,7,-40,-9000,-900],4)
 T_temP.setLabel([500,5000,10,2,7,40,9000,900])
 T_tem_ac=T_temP*1.0
 T_tem_ac_t=T_tem_ac*1.0
 T_tem_ac_t.transpose()
 T_tem_ac_t.setLabel([7,40,9000,900, 500,5000,10,2])
 
 
 
 b1.permute([9000,100,5,800,1],3)
 d1.permute([900,1,9,80,400],3)

 Swap=fermionicOPT(Sys, d1.bond(0), d1.bond(1))
 Swap.setLabel([-900,-1,900,1])
 d1=d1*Swap
 d1.permute([-1,-900,9,80,400],3)
 d1.setLabel([1,900,9,80,400])
 d1.permute([900,1,9,80,400],3)
 #print b1.printDiagram(), d1.printDiagram()
 T_tem=b1*d1

 #print T_tem.printDiagram()
 T_tem.permute([9000,100,5,800,900,9,80,400],4)

 Swap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 Swap1.setLabel([-9000,18,9000,-900])
 T_temP=T_tem*Swap1 

 Swap2=fermionicOPT(Sys, T_tem.bond(1), T_tem.bond(4))
 Swap2.setLabel([-100,17,100,18])
 T_temP=T_temP*Swap2

 Swap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 Swap3.setLabel([-5,16,5,17])
 T_temP=T_temP*Swap3

 Swap4=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(4))
 Swap4.setLabel([15,900,800,16])
 T_temP=T_temP*Swap4

 Swap5=fermionicOPT(Sys, T_tem.bond(3), T_tem.bond(5))
 Swap5.setLabel([-800,9,15,-9])
 T_temP=T_temP*Swap5

 #print T_temP.printDiagram()
 T_temP.permute([-900,-9000,-100,-5,-9,-800,80,400],4)
 T_tem_bd=T_temP*1.0
 T_tem_bd.setLabel([900,9000,100,5,9,800,80,400])

 T_tem_bd_t=T_tem_bd*1.0
 T_tem_bd_t.transpose()
 T_tem_bd_t.setLabel([9,800,80,400,900,9000,100,5])



 A_f.permute([500,5000,10,2,7,40,100,5,9,800,80,400],6)
 A_f_t=A_f*1.0
 A_f_t.transpose()
 A_f_t.setLabel([100,5,9,800,80,400,500,5000,10,2,7,40])
 
 
 val0=A_f_t*A_f
 A_tem_t=T_tem_bd_t*T_tem_ac_t
 A_tem=T_tem_bd*T_tem_ac
 val1=A_tem_t*A_tem
 val2=(T_tem_bd*T_tem_ac)*A_f_t
 #print val0, val1, val2
 #return val1[0]-2.*val2[0] 
 return val0[0]+val1[0]-2.*val2[0] 















##@profile
def init_abcd_SVD( A_ten, Q, Sys, q_D, D, Rg_var):

 A_ten.setLabel([0,1,20,3,4])
 Q.setLabel([2,-5,-7,9,20])
 Swap1=fermionicOPT(Sys, Q.bond(1), Q.bond(2))
 Swap1.setLabel([5,7,-5,-7])

 A_q=A_ten*Q
 A_q=A_q*Swap1
 A_q.permute([0,1,3,4,2,7,5,9],8)
 A_q.setLabel([0,1,3,4,2,7,5,9])

 bdi1 = uni10.Bond(uni10.BD_IN, A_ten.bond(1).Qlist())
 bdo1 = uni10.Bond(uni10.BD_OUT, A_ten.bond(1).Qlist())
 T1=uni10.UniTensor([bdi1, bdi1, bdo1])
 T1.setLabel([10,100,1])

 bdi = uni10.Bond(uni10.BD_IN, A_ten.bond(4).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, A_ten.bond(4).Qlist())
 T3=uni10.UniTensor([bdi, bdo, bdo])
 T3.setLabel([4, 40,400])


 if Rg_var[4]=="iden":
  T3.identity() 
  T1.identity() 

 if Rg_var[4]=="part":
  T3=rand_SVD(T3,bdi,bdo,Sys)
  T1=rand_SVD(T1,bdi1,bdo1,Sys)
  T1.setLabel([10,100,1])
  T3.setLabel([4, 40,400])



 T_tem=(A_q*T1)*T3
 T_tem.permute([0,10,100,2,7,5,9,3,400,40],7)

 Swap1=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(3))
 Swap1.setLabel([16,-2,100,2])
 T_temP=T_tem*Swap1 
 
 Swap2=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(4))
 Swap2.setLabel([17,-7,16,7])
 T_temP=T_temP*Swap2
 
 Swap3=fermionicOPT(Sys, T_tem.bond(2), T_tem.bond(9))
 Swap3.setLabel([-100,21,17,4000])
 T_temP=T_temP*Swap3
 
 Swap4=fermionicOPT(Sys, T_tem.bond(9), T_tem.bond(5))
 Swap4.setLabel([20,-5,21,5])
 T_temP=T_temP*Swap4

 Swap5=fermionicOPT(Sys, T_tem.bond(9), T_tem.bond(6))
 Swap5.setLabel([19,-9,20,9])
 T_temP=T_temP*Swap5

 Swap6=fermionicOPT(Sys, T_tem.bond(9), T_tem.bond(7))
 Swap6.setLabel([18,3,19,-3])
 T_temP=T_temP*Swap6

 Swap7=fermionicOPT(Sys, T_tem.bond(8), T_tem.bond(9))
 Swap7.setLabel([400,40,-400,18])
 T_temP=T_temP*Swap7



 T_temP.permute([0,10,-2,-7,4000,-100,-5,-9,-3, -400],5)
 T_temP.setLabel([0,10,2,7,40,100,5,9,3, 400])
 T_tem=T_temP*1.0
 
 A_tensor_swap=T_tem*1.0
 
 row, colm=cal_rowcol(T_tem)

 if (row<=colm and row<=sum(D)*sum(D)):
  U,V,s=TU.setTruncation(T_tem,row)
 elif (row<=colm and row>sum(D)*sum(D)):
  U,V,s=TU.setTruncation(T_tem,sum(D)*sum(D))
 elif (row>colm and colm<=sum(D)*sum(D)):
  U,V,s=TU.setTruncation(T_tem,colm)
 elif (row>colm and colm>sum(D)*sum(D)):
  U,V,s=TU.setTruncation(T_tem,sum(D)*sum(D))

 U.setLabel([0,10,2,7,40,-30])
 V.setLabel([-30,100,5,9,3,400])
 s1=Sqrt(s)
 s1.setLabel([-30,30])
 U1=U*s1
 s1.setLabel([30,-30])
 V1=s1*V
 U1.permute([0,10,2,7,40,30],3)
 V1.permute([30,100,5,9,3,400],3)

####################
 bdi1 = uni10.Bond( uni10.BD_IN, A_ten.bond(0).Qlist())
 bdo1 = uni10.Bond( uni10.BD_OUT, A_ten.bond(0).Qlist())

 T0=uni10.UniTensor([bdi1, bdi1, bdo1])
 T0.setLabel([5000,500,0])

 bdo = uni10.Bond(uni10.BD_OUT, q_D)
 bdi = uni10.Bond(uni10.BD_IN, s1.bond(0).Qlist())
 #bdo = uni10.Bond(uni10.BD_OUT, s1.bond(0).Qlist())

 T2=uni10.UniTensor([bdi, bdo, bdo ])
 T2.setLabel([30, 9000,900])


 if Rg_var[5]=="on":
  T0.setLabel([0, 500,5000])
  T2.setLabel([30,900,9000])


 if Rg_var[4]=="iden":
  T0.identity()
  T2.identity() 


 if Rg_var[4]=="part":
  T0=rand_SVD(T0,bdi1,bdo1,Sys) 
  T2=rand_SVD(T2,bdi,bdo,Sys) 
  T0.setLabel([5000,500,0])
  T2.setLabel([30, 9000,900])


 A_tensor_swap=A_tensor_swap*T0
 T_tem=(U1*T0)*T2

 T_tem.permute([500,5000,10,2,7,40,9000,900],4)
 Swap1=fermionicOPT(Sys, T_tem.bond(5), T_tem.bond(6))
 Swap1.setLabel([40,9000,16,-9000])
 T_temP=T_tem*Swap1 

 Swap2=fermionicOPT(Sys, T_tem.bond(5), T_tem.bond(7))
 Swap2.setLabel([16,900,-40,-900])
 T_temP=T_temP*Swap2

 T_temP.permute([500,5000,10,2,7,-9000,-900,-40],4)
 T_temP.setLabel([500,5000,10,2,7,9000,900,40])
 T_tem=T_temP*1.0 
 T_tem.permute([500,5000,10,2,7,9000,900,40],4)

 Swap1=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(1))
 Swap1.setLabel([16,-5000,500,5000])
 T_temP=T_tem*Swap1 

 Swap2=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(2))
 Swap2.setLabel([17,-10,16,10])
 T_temP=T_temP*Swap2

 Swap3=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(3))
 Swap3.setLabel([18,-2,17,2])
 T_temP=T_temP*Swap3

 Swap4=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(5))
 Swap4.setLabel([-500,19,18,-9000])
 T_temP=T_temP*Swap4

 Swap5=fermionicOPT(Sys, T_tem.bond(4), T_tem.bond(5))
 Swap5.setLabel([7,9000,-7,19])
 T_temP=T_temP*Swap5

 T_temP.permute([-5000,-10,-2,-9000,-500,-7,900,40],4)
 T_tem=T_temP*1.0
 T_tem.setLabel([5000,10,2,9000,500,7,900,40])

 row, colm=cal_rowcol(T_tem)
 #print sum(D), colm, row
 if (row<=colm and row<=sum(D)):
  U,V,s=TU.setTruncation(T_tem,row)
 elif (row<=colm and row>sum(D)):
  U,V,s=TU.setTruncation(T_tem,sum(D))
 elif (row>colm and colm<=sum(D)):
  U,V,s=TU.setTruncation(T_tem,colm)
 elif (row>colm and colm>sum(D)):
  U,V,s=TU.setTruncation(T_tem,sum(D))

 U.setLabel([5000,10,2,9000,-1])
 V.setLabel([-1,500,7,900,40])
 s=Sqrt(s)
 s.setLabel([-1,1])
 a=U*s
 s.setLabel([1,-1])
 c=s*V
 a.permute([5000,10,2,9000,1],3)
 c.permute([1,500,7,900,40],3)

 Swap=fermionicOPT(Sys, c.bond(0), c.bond(1))
 Swap.setLabel([-1,-500,1,500])
 c=c*Swap
 c.permute([-500,-1,7,900,40],3)
 c.setLabel([500,1,7,900,40])

######################################################
 bdi1 = uni10.Bond( uni10.BD_IN, q_D)
 bdo1 = uni10.Bond( uni10.BD_OUT, s1.bond(1).Qlist())

 T0p=uni10.UniTensor([bdi1, bdi1, bdo1])
 T0p.setLabel([9000,900,30])

 bdi = uni10.Bond(uni10.BD_IN, A_ten.bond(3).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, A_ten.bond(3).Qlist())

 T2p=uni10.UniTensor([bdi, bdo, bdo])
 T2p.setLabel([3,800,80])

 if Rg_var[5]=="on":
  T2p.setLabel([3,800,80])
  T0p.setLabel([900,9000,30])


 if Rg_var[4]=="iden":
  T0p.identity()
  T2p.identity() 


 
 if Rg_var[4]=="part":
  T0p=rand_SVD(T0p, bdi1, bdo1,Sys) 
  T2p=rand_SVD(T2p, bdi, bdo,Sys) 
  T0p.setLabel([9000,900,30])
  T2p.setLabel([800,80,3])


 A_tensor_swap=A_tensor_swap*T2p
 A_tensor_swap.permute([500,5000,10,2,7,40,100,5,9,800,80,400],6)




 V1.permute([30,100,5,9,3,400],3)
 T_tem=(V1*T0p)*T2p
 T_tem.permute([900,9000,100,5,9,800,80,400],4)

 Swap1=fermionicOPT(Sys, T_tem.bond(4), T_tem.bond(5))
 Swap1.setLabel([9,800,-9,-800])
 T_temP=T_tem*Swap1

 T_tem=T_temP*1.0
 T_tem.permute([900,9000,100,5,-800,-9,80,400],4)


 Swap2=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(1))
 Swap2.setLabel([16,-9000,900,9000])
 T_temP=T_tem*Swap2

 Swap3=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(2))
 Swap3.setLabel([17,-100,16,100])
 T_temP=T_temP*Swap3

 Swap4=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(3))
 Swap4.setLabel([18,-5,17,5])
 T_temP=T_temP*Swap4

 Swap5=fermionicOPT(Sys, T_tem.bond(0), T_tem.bond(4))
 Swap5.setLabel([-900,-800,18,800])
 T_temP=T_temP*Swap5

 T_temP.permute([-9000,-100,-5,800,-900,-9,80,400],4)
 T_tem=T_temP*1.0
 T_tem.setLabel([9000,100,5,800,900,9,80,400])


 row, colm=cal_rowcol(T_tem)
 if (row<=colm and row<=sum(D)):
  U,V,s=TU.setTruncation(T_tem,row)
 elif (row<=colm and row>sum(D)):
  U,V,s=TU.setTruncation(T_tem,sum(D))
 elif (row>colm and colm<=sum(D)):
  U,V,s=TU.setTruncation(T_tem,colm)
 elif (row>colm and colm>sum(D)):
  U,V,s=TU.setTruncation(T_tem,sum(D))

 U.setLabel([9000,100,5,800,-1])
 V.setLabel([-1,900,9,80,400])
 s=Sqrt(s)
 s.setLabel([-1,1])
 b=U*s
 s.setLabel([1,-1])
 d=s*V
 b.permute([9000,100,5,800,1],3)
 d.permute([1,900,9,80,400],3)

 Swap=fermionicOPT(Sys, d.bond(0), d.bond(1))
 Swap.setLabel([-1,-900,1,900])
 d=d*Swap
 d.permute([-900,-1,9,80,400],3)
 d.setLabel([900,1,9,80,400])

 return a, b, c, d, A_tensor_swap


##@profile
def Var_PEPSQ( A_ten, Q, q_D, Sys, D, Rg_var):


 bdi = uni10.Bond( uni10.BD_IN, q_D)
 bdo = uni10.Bond( uni10.BD_OUT, q_D)

 a=uni10.UniTensor([ A_ten.bond(0), A_ten.bond(1), Q.bond(2), bdi, bdi ])
 a.setLabel([-2,-3,2,3,4])

 b=uni10.UniTensor([ bdi, A_ten.bond(1), Q.bond(2),A_ten.bond(3),bdo ])
 b.setLabel([3,-4,5,800,600])
 b.identity()


 c=uni10.UniTensor([ A_ten.bond(0), bdi, Q.bond(2), bdo, A_ten.bond(4) ])
 c.setLabel([-1,4,70,8,-7])
 c.identity()

 d=uni10.UniTensor([ bdi , bdi, Q.bond(2), A_ten.bond(3), A_ten.bond(4) ])
 d.setLabel([8,60,90,-5,-8])
 d.identity()


 a, b, c, d, A_f=init_abcd_SVD( A_ten, Q, Sys, q_D, D, Rg_var)



 if Rg_var[0]=="iden":
  a.identity()
  b.identity()
  c.identity()
  d.identity()


 if Rg_var[0]=="randm":
  a.randomize()
  b.randomize()
  c.randomize()
  d.randomize()

 if Rg_var[0]=="randQ":
  a.orthoRand()
  b.orthoRand()
  c.orthoRand()
  d.orthoRand()

 if Rg_var[0]=="previous":
  a=inint_abcs(a, A_ten)
  b=inint_abcs(b, A_ten)
  c=inint_abcs(c, A_ten)
  d=inint_abcs(d, A_ten)


 if Rg_var[0]=="SVD":
  a, b, c, d, A_f=init_abcd_SVD( A_ten, Q, Sys, q_D, D, Rg_var)



 a_init=a*1.0
 b_init=b*1.0
 c_init=c*1.0
 d_init=d*1.0
 valf=cost_f(A_f,a,b,c,d, Sys)
 print "valf", valf
 E_1=float("inf")
 E_2=0
 E_0=0
 count=1
 for i in xrange(Rg_var[2]):
  count=count+1
  E_2=E_1*1.0
  val=cost_f( A_f, a, b, c, d, Sys)
  if i==0: E_0=val
  E_1=val
  if abs(E_1)<1.0e-16: 
   print "break, E_1 is small", E_1
   break

  print i, val, abs((E_1-E_2)/E_1)

  if E_1>E_2  or abs((E_1-E_2)/E_1)<+1.0e-12:
   print "break"
   #print E_1, E_2, abs((E_1-E_2)/E_1)
   a=a_init*1.0
   b=b_init*1.0
   c=c_init*1.0
   d=d_init*1.0
   break
  else:
   a_init=a*1.0
   b_init=b*1.0
   c_init=c*1.0
   d_init=d*1.0
  
  a=optimiz_a( A_f, a, b, c, d, Sys)
  c=optimiz_c( A_f, a, b, c, d, Sys)
  b=optimiz_b( A_f, a, b, c, d, Sys)
  d=optimiz_d( A_f, a, b, c, d, Sys)


 if Rg_var[3]=="on":
  a, b, c, d=Do_optimization_Grad(A_f, a, b, c, d, Sys,Swap1,Swap2,Swap3)

 a.setLabel([-2,-3,2,3,4])
 b.setLabel([3,-4,5,-6,6])
 c.setLabel([-1,4,7,8,-8])
 d.setLabel([8,6,9,-5,-7])

 a.permute([-2,-3,2,3,4],3)
 b.permute([3,-4,5,-6,6],3)
 c.permute([-1,4,7,8,-8],3)
 d.permute([8,6,9,-5,-7],3)

 F=a*b*c*d

 #print Sys[0]

 return  a, b, c, d




def Do_optimization_Grad(A_f, a, b, c, d, Sys,Swap1,Swap2,Swap3):
  Opt_method="ST"
  #Opt_method="grad"

  a.setLabel([-2,-3,2,3,4])
  b.setLabel([3,-4,5,800,600])
  c.setLabel([-1,4,70,8,-8])
  d.setLabel([8,60,90,-5,-7])

  a.permute([-2,-3,2,3,4],3)
  b.permute([3,-4,5,800,600],3)
  c.permute([-1,4,70,8,-8],3)
  d.permute([8,60,90,-5,-7],3)

  a_um=copy.copy(a)
  b_um=copy.copy(b)
  c_um=copy.copy(c)
  d_um=copy.copy(d)

  time_val=0
  Es=cost_f( A_f, a, b, c, d, Swap1, Swap2, Swap3)
  Ef=0
  E2_val=0
  #print '\n', '\n', '\n', '\n'
  Gamma=1.0
  E_previous=1.e+18
  count=0
  D_list=[0]*4
  H_list=[0]*4
  H_a=0; H_b=0;H_c=0;H_d=0;
  for i in xrange(200):
   count+=1
   t0=time.time()

   E1_val=cost_f( A_f, a, b, c, d,Swap1,Swap2,Swap3)
   Ef=E1_val
   #if i%200==0:
   print 'E1=', E1_val, abs((E_previous-E1_val)/E1_val), i, count, time_val


   D_a,D_b,D_c,D_d=Obtain_grad_four(A_f, a, b, c, d,Swap1,Swap2,Swap3)
   D_a=(-1.0)*D_a
   D_b=(-1.0)*D_b
   D_c=(-1.0)*D_c
   D_d=(-1.0)*D_d

   if i is 0:
    H_a=D_a
    H_b=D_b
    H_c=D_c
    H_d=D_d
   else:
    Z_a=D_a+(-1.0)*D_list[0]
    Z_b=D_b+(-1.0)*D_list[1]
    Z_c=D_c+(-1.0)*D_list[2]
    Z_d=D_d+(-1.0)*D_list[3]
    A=Z_a*D_a
    B=Z_b*D_b
    C=Z_c*D_c
    D=Z_d*D_d
    A1=D_list[0]*D_list[0]
    A2=D_list[1]*D_list[1]
    A3=D_list[2]*D_list[2]
    A4=D_list[3]*D_list[3]
    Gamma_grad=(A[0]+B[0]+C[0]+D[0]) / (A1[0]+A2[0]+A3[0]+A4[0])
    if Opt_method is 'ST':Gamma_grad=0;
    H_a=D_a+(Gamma_grad)*H_list[0]
    H_b=D_b+(Gamma_grad)*H_list[1]
    H_c=D_c+(Gamma_grad)*H_list[2]
    H_d=D_d+(Gamma_grad)*H_list[3]

#    A=D_a*D_list[0]
#    B=D_b*D_list[1]
#    C=D_c*D_list[2]
#    D=D_d*D_list[3]
#    check=A[0]+B[0]+C[0]+D[0] 
#    print "check", check 

   D_list[0]=copy.copy(D_a)
   D_list[1]=copy.copy(D_b)
   D_list[2]=copy.copy(D_c)
   D_list[3]=copy.copy(D_d)

   H_list[0]=copy.copy(H_a)
   H_list[1]=copy.copy(H_b)
   H_list[2]=copy.copy(H_c)
   H_list[3]=copy.copy(H_d)



   A=D_a*H_a
   B=D_b*H_b
   C=D_c*H_c
   D=D_d*H_d
   
   Norm_Z=A[0]+B[0]+C[0]+D[0]
   
   #print 'Norm', Norm_Z
   if (E1_val<E_previous) or (i is 0):
    if (abs(E1_val) > 1.0e-10):
     if abs((E_previous-E1_val)/E1_val) < 1.0e-12:
      print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)/E1_val), i
      break
     else: 
      if abs((E_previous-E1_val)) < 1.0e-15:
       print 'Differnance Satisfied!', E_previous, E1_val, abs((E_previous-E1_val)), i
       break

   #print E1_val>E_previous  
   if (E1_val>E_previous):
    print "break, not satisfied", E1_val, E_previous
    a=copy.copy(a_um) 
    b=copy.copy(b_um) 
    c=copy.copy(c_um) 
    d=copy.copy(d_um)
    break
   else:
    a_um=copy.copy(a) 
    b_um=copy.copy(b) 
    c_um=copy.copy(c) 
    d_um=copy.copy(d) 

   E_previous=E1_val
   
   if abs(Norm_Z) < 1.0e-12:
    print 'Break Norm=', Norm_Z
    break
   Break_loop=1
   if (i%15)==0: 
    Gamma=1
   else: 
    if Gamma >= 1: Gamma*=(1.00/100)
    if Gamma < 1: Gamma*=100
   #Gamma=1
   #print "Gamma", Gamma
   while Break_loop is 1:
    count+=1
    #print a.printDiagram(), H_a.printDiagram()
    a_ut=a+(2.00)*Gamma*H_a
    b_ut=b+(2.00)*Gamma*H_b
    c_ut=c+(2.00)*Gamma*H_c
    d_ut=d+(2.00)*Gamma*H_d
    E2_val=cost_f(A_f, a_ut, b_ut, c_ut, d_ut,Swap1,Swap2,Swap3)
    if abs((0.5)*Norm_Z*Gamma) > 1.0e+12 or  abs(Gamma)>1.0e+12 :
     print "break1", E1_val, abs((0.5)*Norm_Z*Gamma), E2_val, Gamma
     Gamma=1
     break
    if E1_val-E2_val >=(Norm_Z*Gamma):
     Gamma*=2.00
    else:
     Break_loop=0

   Break_loop=1
   while Break_loop is 1:
    count+=1
    a_ut=a+(1.00)*Gamma*H_a
    b_ut=b+(1.00)*Gamma*H_b
    c_ut=c+(1.00)*Gamma*H_c
    d_ut=d+(1.00)*Gamma*H_d
    E2_val=cost_f(A_f, a_ut, b_ut, c_ut, d_ut,Swap1,Swap2,Swap3)
    #print "Gamma", Gamma
    if abs((0.5)*Norm_Z*Gamma) <1.0e-16 or  (abs((E1_val-E2_val)/E2_val))<1.0e-16 or abs(Gamma)<1.0e-16 :
     print "break2", E1_val, E2_val, Gamma, abs((0.5)*Norm_Z*Gamma), (abs((E1_val-E2_val)/E2_val))
     break
     
    if E1_val-E2_val < (0.50)*Norm_Z*Gamma:
     Gamma*=0.5
    else:
     Break_loop=0


   a=a+(1.00)*Gamma*H_a
   b=b+(1.00)*Gamma*H_b
   c=c+(1.00)*Gamma*H_c
   d=d+(1.00)*Gamma*H_d
   time_val=time.time() - t0

  return a, b, c, d





def Obtain_grad_four(A_f, a, b, c, d,Swap1,Swap2,Swap3):

 bdi = uni10.Bond(uni10.BD_IN, a.bond(0).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, a.bond(0).Qlist())

 T0=uni10.UniTensor([bdi, bdo])
 T0.setLabel([-2,-20])
 T0.identity() 

 bdi = uni10.Bond(uni10.BD_IN, a.bond(1).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, a.bond(1).Qlist())

 T1=uni10.UniTensor([bdi, bdo])
 T1.setLabel([-3,-30])
 T1.identity() 

 bdi = uni10.Bond(uni10.BD_IN, a.bond(2).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, a.bond(2).Qlist())

 T2=uni10.UniTensor([bdi, bdo])
 T2.setLabel([2,20])
 T2.identity() 

 A_tem=((c*d)*((b*Swap1)*Swap2))*Swap3
 A_tem.permute( [ 4, 3, -1,-7, -8, -5, -6, -4, 5, 7, 9] , 10)

 b_tem=A_tem*1.0
 b_tem.setLabel( [ 40, 30, -1,-7, -8, -5, -6, -4, 5, 7, 9])

 T=b_tem*A_tem
 T=((T*T2)*T1)*T0
 T.permute([-2,-3,2,3,4,-20,-30,20,30,40],5)
 a1=a*1.0
 a1.setLabel([-20,-30,20,30,40])
 D_a=T*a1
 D_a.permute([-2,-3,2,3,4],3)
 
 Tp=A_f*((c*d)*((b*Swap1)*Swap2))*Swap3
 Tp.permute([-2,-3,2,3,4],3)
 D_a=D_a+(-1.0)*Tp


 bdi = uni10.Bond(uni10.BD_IN, b.bond(1).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, b.bond(1).Qlist())

 T0=uni10.UniTensor([bdi, bdo])
 T0.setLabel([-4,-40])
 T0.identity()


 bdi = uni10.Bond(uni10.BD_IN, b.bond(2).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, b.bond(2).Qlist())

 T2=uni10.UniTensor([bdi, bdo])
 T2.setLabel([5,50])
 T2.identity()

 A_tem=(((c*a)*Swap1)*(d*Swap3)*Swap2)
 A_tem.permute( [ -1, -2, -3, 3, 600, -6,-5, -7, -8, 2, 7, 9,800] , 10)


 b_tem=A_tem*1.0
 b_tem.setLabel( [ -1, -2, -3, 30, 6000, -6,-5, -7, -8, 2, 7, 9,8000])

 T=b_tem*A_tem
 T=((T*T2))*T0
 T.permute([3,-4,5,800,600,30,-40,50,8000,6000],5)

 b1=b*1.0
 b1.setLabel([30,-40,50,8000,6000])
 D_b=T*b1
 D_b.permute([3,-4,5,800,600],3)

 Tp=A_f*((((c*a)*Swap1)*(d*Swap3)*Swap2))
 Tp.permute([3,-4,5,800,600],3)
 D_b=D_b+(-1.0)*Tp


 bdi = uni10.Bond(uni10.BD_IN, c.bond(0).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, c.bond(0).Qlist())

 T0=uni10.UniTensor([bdi, bdo])
 T0.setLabel([-1,-10])
 T0.identity()

 bdi = uni10.Bond(uni10.BD_IN, c.bond(4).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, c.bond(4).Qlist())

 T1=uni10.UniTensor([bdi, bdo])
 T1.setLabel([-8,-80])
 T1.identity()

 bdi = uni10.Bond(uni10.BD_IN, c.bond(2).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, c.bond(2).Qlist())


 A_tem=((a*((b*Swap1)*Swap2))*d)*Swap3
 A_tem.permute( [ 4, 8,-2, -3, -4, -6, -5, -7, 2, 5, 9,7,70] , 10)

 b_tem=A_tem*1.0
 b_tem.setLabel( [ 40, 80,-2, -3, -4, -6, -5, -7, 2, 5, 9,7,700])

 T=b_tem*A_tem
 T=((T)*T1)*T0
 T.permute([-1,4,70,8,-8,-10,40,700,80,-80],5)
 c1=c*1.0
 c1.setLabel([-10,40,700,80,-80])
 D_c=T*c1
 D_c.permute([-1,4,70,8,-8],3)


 Tp=A_f*(((a*((b*Swap1)*Swap2))*d)*Swap3)
 Tp.permute([-1,4,70,8,-8],3)
 D_c=D_c+(-1.0)*Tp


 bdi = uni10.Bond(uni10.BD_IN, d.bond(3).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, d.bond(3).Qlist())

 T0=uni10.UniTensor([bdi, bdo])
 T0.setLabel([-5,-50])
 T0.identity()

 bdi = uni10.Bond(uni10.BD_IN, d.bond(4).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, d.bond(4).Qlist())

 T1=uni10.UniTensor([bdi, bdo])
 T1.setLabel([-7,-70])
 T1.identity()

 bdi = uni10.Bond(uni10.BD_IN, d.bond(2).Qlist())
 bdo = uni10.Bond(uni10.BD_OUT, d.bond(2).Qlist())



 A_tem=((a*((b*Swap1)*Swap2))*c)*Swap3
 A_tem.permute([ 8, 60,-1,-2, -3, -4, -6, -8, 2, 5, 7,9,90], 10)


 b_tem=A_tem*1.0
 b_tem.setLabel([ 80, 600,-1,-2, -3, -4, -6, -8, 2, 5, 7,9,900])


 T=b_tem*A_tem
 T=((T)*T1)*T0
 T.permute([8,60,90,-5,-7,80,600,900,-50,-70],5)
 d1=d*1.0
 d1.setLabel([80,600,900,-50,-70])
 D_d=T*d1
 D_d.permute([8,60,90,-5,-7],3)


 Tp=A_f*(((a*((b*Swap1)*Swap2))*c)*Swap3)
 Tp.permute([8,60,90,-5,-7],3)

 D_d=D_d+(-1.0)*Tp



 return D_a, D_b, D_c, D_d





