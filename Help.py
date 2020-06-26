
# #T.orthoRand()
# blk_qnums = T.blockQnum()
# for qnum in blk_qnums:
#  T_mat=T.getBlock(qnum)
#  U=T_mat.svd()
#  if T_mat.row()>=T_mat.col():
#   U_mat=U[0]
#   T.putBlock(qnum,U_mat)
#  else:
#   U_mat=U[2]
#   T.putBlock(qnum,U_mat)




# a.setLabel([-2,-3,2,3,4])
# b.setLabel([3,-4,5,800,600])
# c.setLabel([-1,4,70,8,-7])
# d.setLabel([8,60,90,-5,-8])


# Swap1=fermionicOPT(Sys, c.bond(2), b.bond(4))
# Swap1.setLabel([70,600,900,60])

# Swap2=fermionicOPT(Sys, c.bond(2),b.bond(3))
# Swap2.setLabel([900,800,7,-600])

# Swap3=fermionicOPT(Sys, d.bond(2),b.bond(3))
# Swap3.setLabel([90,-600,9,-6])

# if Rg_var[1]=="iden":
#  Swap1.identity()
#  Swap2.identity()
#  Swap3.identity()
#  #print "Hi"

 A_new=((a*((b*Swap1)*Swap2))*(c*d))*Swap3

# val=A_new*A_new
# val1=A_f*A_f
# print val, val1
 #print A_new.printDiagram()

 A_new.permute([-1,-2,-3,-4,2,5,7,9,-6,-5,-8,-7],12)
 #print A_f.printDiagram()

 #val=A_new*A_f
 #print  val

 val=cost_f(A_f,a,b,c,d,Swap1,Swap2,Swap3)
 a_init=a*1.0
 b_init=b*1.0
 c_init=c*1.0
 d_init=d*1.0

 valf=cost_f(A_f,a,b,c,d,Swap1,Swap2,Swap3)
 
 E_1=float("inf")
 E_2=0
 E_0=0
 count=1
 for i in xrange(Rg_var[2]):
  count=count+1
  E_2=E_1*1.0
  val=cost_f( A_f, a, b, c, d,Swap1,Swap2,Swap3)
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

  a=optimiz_a( A_f, a, b, c, d, Sys, Swap1, Swap2, Swap3)
  b=optimiz_b( A_f, a, b, c, d, Sys, Swap1, Swap2, Swap3)
  c=optimiz_c( A_f, a, b, c, d, Sys, Swap1, Swap2, Swap3)
  d=optimiz_d( A_f, a, b, c, d, Sys, Swap1, Swap2, Swap3)


# a=inint_abcs(a, A_ten)
# b=inint_abcs(b, A_ten)
# c=inint_abcs(c, A_ten)
# d=inint_abcs(d, A_ten)

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

 return  a, b, c, d

















#def cost_T(A,T0,T1,T2,T3,a,b,c,d,Q):
#  A_f=((((A*T0)*T1)*T2)*T3)*Q
#  val0=A_f*A_f
#  A_tem=((a*b)*(c*d))
#  val1=A_tem*A_tem
#  val2=((a*b)*(c*d))*A_f
#  #print val0, val1, val2
#  #return val1[0]-2.*val2[0] 
#  return val0[0]+val1[0]-2.*val2[0] 




##@profile
#def optimiz_T0( A,T0,T1,T2,T3,a,b,c,d,Sys,Q):


# bdi = uni10.Bond(uni10.BD_IN, T0.bond(0).Qlist())
# bdo = uni10.Bond(uni10.BD_OUT, T0.bond(0).Qlist())
# T_h=uni10.UniTensor([bdi, bdo])
# T_h.setLabel([-1,-10])
# T_h.identity()

# bdi = uni10.Bond(uni10.BD_IN, T0.bond(1).Qlist())
# bdo = uni10.Bond(uni10.BD_OUT, T0.bond(1).Qlist())
# T_h1=uni10.UniTensor([bdi, bdo])
# T_h1.setLabel([-2,-20])
# T_h1.identity()




# A_tem=((((A)*T1)*T2)*T3)*Q
# A_tem.permute([0,-3,-4,2,5,7,9,-5,-6,-7,-8],12)

# b_tem=A_tem*1.0
# b_tem.setLabel([ 10,-3,-4,2,5,7,9,-5,-6,-7,-8])




# T=b_tem*A_tem
# T=((T*T_h1)*T_h)
# T.permute([-1,-2,0,-10,-20,10],3)


# A_f=((((A*T0)*T1)*T2)*T3)*Q

# Tp=A_f*( ((((A)*T1)*T2)*T3)*Q )
# Tp.permute([-1,-2,0],3)






# if Sys[4]=="SVD":
#   row, colm=cal_rowcol(T)
#   if (row<=colm):
#    U,V,S=TU.setTruncation(T,row)
#   else:
#    U,V,S=TU.setTruncation(T,colm)


#   U.transpose()
#   V.transpose()
#   S=inverse(S)
#   
#   U.setLabel([5,6,7,8])
#   V.setLabel([1,2,3,4])
#   S.setLabel([4,5])

#   A2_inv=V*S*U
#   A2_inv.permute([1,2,3,6,7,8],3)
#   A2_inv.setLabel([-10,-20,10,-1,-2,0])

#   T0=A2_inv*Tp
#   T0.permute([-10,-20,10],3)
#   T0.setLabel([-1,-2,0])
# elif  Sys[4]=="cg":
#  try:
#   T0=solve_linear_eq(T,Tp)
#   T0.setLabel([-1,-2,0])
#  except:
#   print "error, SVD"
#   row, colm=cal_rowcol(T)
#   if (row<=colm):
#    U,V,S=TU.setTruncation(T,row)
#   else:
#    U,V,S=TU.setTruncation(T,colm)


#   U.setLabel([5,6,7,8])
#   V.setLabel([1,2,3,4])
#   S.setLabel([4,5])

#   A2_inv=V*S*U
#   A2_inv.permute([1,2,3,6,7,8],3)
#   A2_inv.setLabel([-10,-20,10,-1,-2,0])

#   T0=A2_inv*Tp
#   T0.permute([-10,-20,10],3)
#   T0.setLabel([-1,-2,0])

# return T0



#############################################################################################
# val=cost_f( A_f, a, b, c, d)
# print "val", val



# val=cost_T(A,T0,T1,T2,T3,a,b,c,d, Q)
# print "val", val

# T0_init=T0*1.0
# T1_init=T1*1.0
# T2_init=T2*1.0
# T3_init=T3*1.0

# valf=cost_T(A,T0,T1,T2,T3,a,b,c,d, Q)

# E_1=float("inf")
# E_2=0
# E_0=0
# count=1
# for i in xrange(40):
#  count=count+1
#  E_2=E_1*1.0
#  val=cost_T(A,T0,T1,T2,T3,a,b,c,d, Q)
#  if i==0: E_0=val
#  E_1=val
#  if abs(E_1)<1.0e-16: 
#   print "break, E_1 is small", E_1
#   break

#  print i, val, abs((E_1-E_2)/E_1)


#  if E_1>E_2  or abs((E_1-E_2)/E_1)<+1.0e-12:
#   print "break"
#   #print E_1, E_2, abs((E_1-E_2)/E_1)
#   T0=T0_init*1.0
#   T1=T1_init*1.0
#   T2=T2_init*1.0
#   T3=T3_init*1.0
#   break
#  else:
#   T0_init=T0*1.0
#   T1_init=T1*1.0
#   T2_init=T2*1.0
#   T3_init=T3*1.0

#  T=optimiz_T0( A,T0,T1,T2,T3,a,b,c,d ,Sys, Q)
#  T1=optimiz_T1( A,T0,T1,T2,T3,a,b,c,d, Sys, Q)
#  T2=optimiz_T2( A,T0,T1,T2,T3,a,b,c,d, Sys, Q)
#  T3=optimiz_T3( A,T0,T1,T2,T3,a,b,c,d, Sys, Q)


def particle_val( H_col, Q_list, i, j, rho_col, N_x, Sys,P_f):

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
  Q_trans.setLabel([-1,2,-3,4,-5])


  Q_list[i][j-1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j-1]*1.0
  Q_transN.setLabel([6,7,8,9,-10])
  H_col[2*i][2*j].setLabel([-1,-3,1,3])
  rho_col[i][j-1].setLabel([-10,-5,10,5])

  result1=(((Q_list_ferm*H_col[2*i][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)
  P_f[2*i][2*j]=result1[0]
################################################

  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-3,-4,3,4])
  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,2,-3,-4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.setLabel([1,-2,3,-4,-5])


  Q_list[i][j-1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j-1]*1.0
  Q_transN.setLabel([6,7,8,9,-10])

  H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
  rho_col[i][j-1].setLabel([-10,-5,10,5])
  P_f[2*i+1][2*j]=(((Q_list_ferm*H_col[2*i+1][2*j])*Q_trans)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_transN)


#######################################################


  Swap.setLabel([12,4,6,-4])

  Q_list[i][j-1].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j-1]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])


  rho_col[i][j-1].setLabel([-5,-10,5,10])
  H_col[2*i][2*j-1].setLabel([-3,-12,3,12])

  P_f[2*i][2*j-1]=((Q_list[i][j]*H_col[2*i][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list_fermin*Q_trans)

  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j-1]*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.setLabel([6,-7,8,9,-10])

  H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  P_f[2*i+1][2*j-1]=((Q_list_ferminP*H_col[2*i+1][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_trans)





 if j<N_x-1:

  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-2,-3,2,3])

  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,-2,-3,4,5],5)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.setLabel([-1,2,-3,4,-5])

  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.setLabel([6,7,8,9,-10])
 
  H_col[2*i][2*j].setLabel([-1,-3,1,3])
  rho_col[i][j].setLabel([-5,-10,5,10])
  P_f[2*i][2*j]=((Q_list_ferm*H_col[2*i][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]

################################################


  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-3,-4,3,4])
  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,2,-3,-4,5],4)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.setLabel([1,-2,3,-4,-5])


  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.setLabel([6,7,8,9,-10])

  H_col[2*i+1][2*j].setLabel([-2,-4,2,4])
  rho_col[i][j].setLabel([-5,-10,5,10])
  P_f[2*i+1][2*j]=((Q_list_ferm*H_col[2*i+1][2*j]*Q_trans)*(Q_list[i][j+1]*Q_transN))*rho_col[i][j]
#######################################################

  Swap.setLabel([12,4,6,-4])

  Q_list[i][j].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])


  H_col[2*i][2*j+1].setLabel([-3,-12,3,12])
  P_f[2*i][2*j+1]=(Q_list_fermin*H_col[2*i][2*j+1]*Q_trans)*((Q_list[i][j+1]*Q_transN)*rho_col[i][j])


  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j+1].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j+1]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.setLabel([6,-7,8,9,-10])


  H_col[2*i+1][2*j+1].setLabel([-4,-7,4,7])
  rho_col[i][j].setLabel([-5,-10,5,10])

  P_f[2*i+1][2*j+1]=((Q_list[i][j]*H_col[2*i+1][2*j+1]*Q_trans)*rho_col[i][j])*(Q_list_ferminP*Q_transN)






###############################################################################
 if j-1>=0 and j!=N_x-1:

  Swap.setLabel([12,4,6,-4])

  Q_list[i][j-1].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j-1]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])


  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])


  rho_col[i][j-1].setLabel([-5,-10,5,10])

  H_col[2*i][2*j-1].setLabel([-3,-12,3,12])
  P_f[2*i][2*j-1]=((Q_list[i][j]*H_col[2*i][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list_fermin*Q_trans)

  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j-1]*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.setLabel([6,-7,8,9,-10])

  H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  P_f[2*i+1][2*j-1]=((Q_list_ferminP*H_col[2*i+1][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_trans)





































































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
    Q_trans.setLabel([-1,2,-3,4,-5])


    Q_tem=Q_list[i][j]*1.0
    Q_tem.setLabel([6,7,8,9,10])
    Q_temN=Q_tem*1.0
    Q_temN.setLabel([6,7,8,9,-10])

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
    Q_trans.setLabel([1,-2,3,-4,-5])


    Q_tem=Q_list[i][j]*1.0
    Q_tem.setLabel([6,7,8,9,10])
    Q_temN=Q_tem*1.0
    Q_temN.setLabel([6,7,8,9,-10])

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
    Q_trans.setLabel([-1,2,-3,4,-5])

    Q_list[i][j+1].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i][j+1]*1.0
    Q_transN.setLabel([6,7,8,9,-10])

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
    Q_trans.setLabel([1,-2,3,-4,-5])

    Q_list[i][j+1].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i][j+1]*1.0
    Q_transN.setLabel([6,7,8,9,-10])

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
    Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

    Q_list[i][j+1].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i][j+1]*1.0
    Q_transN.setLabel([-6,7,8,9,-10])

    H_list[2*i][2*j+1].setLabel([-3,-12,3,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*H_list[2*i][2*j+1])*Q_trans)*(Q_list[i][j+1]*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]

    Swap.setLabel([7,6,-7,-6])

    Q_list[i][j].setLabel([1,2,3,4,5])
    Q_trans=Q_list[i][j]*1.0
    Q_trans.setLabel([1,2,3,-4,-5])

    Q_list[i][j+1].setLabel([-6,-7,8,9,10])
    Q_list_ferminP=Q_list[i][j+1]*Swap
    Q_list_ferminP.permute([6,7,8,9,10],5)
    Q_transN=Q_list_ferminP*1.0
    Q_transN.setLabel([6,-7,8,9,-10])


    H_list[2*i+1][2*j+1].setLabel([-4,-7,4,7])
    result=((Q_list[i][j]*H_list[2*i+1][2*j+1])*Q_trans)*(Q_list_ferminP*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]
    #print 2*i+1, 2*j+1, H_list[2*i+1][2*j+1]

####################   Long  #################################

    Swap1=Swap*1.0
    Swap2=Swap*1.0
    Swap.setLabel([2,12,-2,13])
    Swap1.setLabel([3,13,-3,14])
    Swap2.setLabel([4,14,-4,6])

    Q_list[i][j].setLabel([1,-2,-3,-4,5])
    Q_list_fermin=((Q_list[i][j]*Swap)*Swap1)*Swap2
    Q_list_fermin.permute([1,2,3,4,5,12,6],7)
    Q_trans=Q_list_fermin*1.0
    Q_trans.setLabel([-1,2,3,4,-5,-12,-6])

    Q_list[i][j+1].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i][j+1]*1.0
    Q_transN.setLabel([-6,7,8,9,-10])

    h_long.setLabel([-1,-12,1,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*h_long)*Q_trans)*(Q_list[i][j+1]*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]
###################################################
    Swap1=Swap*1.0
    Swap2=Swap*1.0
    Swap.setLabel([3,12,-3,13])
    Swap1.setLabel([4,13,-4,14])
    Swap2.setLabel([6,14,-6,7])

    Q_list[i][j].setLabel([1,2,-3,-4,5])
    Q_list_fermin=((Q_list[i][j]*Swap)*Swap1)
    Q_list_fermin.permute([1,2,3,4,5,12,14],7)
    Q_trans=Q_list_fermin*1.0
    Q_trans.setLabel([1,-2,3,4,-5,-12,-14])

    Q_list[i][j+1].setLabel([-6,7,8,9,10])
    Q_list_fermin1=((Q_list[i][j+1]*Swap2))
    Q_list_fermin1.permute([6,8,9,10,14],6)
    Q_transN=Q_list_fermin1*1.0
    Q_transN.setLabel([6,8,9,-10,-14])

    h_long.setLabel([-2,-12,2,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*h_long)*Q_trans)*(Q_list_fermin1*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]

#################
    Swap1=Swap*1.0
    Swap2=Swap*1.0
    Swap.setLabel([4,12,-4,13])
    Swap1.setLabel([6,13,-6,14])
    Swap2.setLabel([7,14,-7,8])

    Q_list[i][j].setLabel([1,2,3,-4,5])
    Q_list_fermin=((Q_list[i][j]*Swap))
    Q_list_fermin.permute([1,2,3,4,5,12,13],7)
    Q_trans=Q_list_fermin*1.0
    Q_trans.setLabel([1,2,-3,4,-5,-12,-13])

    Q_list[i][j+1].setLabel([-6,-7,8,9,10])
    Q_list_fermin1=((Q_list[i][j+1]*Swap1))*Swap2
    Q_list_fermin1.permute([6,7,9,10,13],6)
    Q_transN=Q_list_fermin1*1.0
    Q_transN.setLabel([6,7,9,-10,-13])

    h_long.setLabel([-3,-12,3,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*h_long)*Q_trans)*(Q_list_fermin1*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]

####################
    Swap1=Swap*1.0
    Swap2=Swap*1.0
    Swap.setLabel([6,12,-6,13])
    Swap1.setLabel([7,13,-7,14])
    Swap2.setLabel([8,14,-8,9])

    Q_list[i][j].setLabel([1,2,3,4,5])
    Q_list_fermin=((Q_list[i][j]))
    Q_trans=Q_list_fermin*1.0
    Q_trans.setLabel([1,2,3,-4,-5])

    Q_list[i][j+1].setLabel([-6,-7,-8,9,10])
    Q_list_fermin1=((Q_list[i][j+1]*Swap)*Swap1)*Swap2
    Q_list_fermin1.permute([6,7,8,10,12],5)
    Q_transN=Q_list_fermin1*1.0
    Q_transN.setLabel([6,7,8,-10,-12])

    h_long.setLabel([-4,-12,4,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*h_long)*Q_trans)*(Q_list_fermin1*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]


 return HA_list





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



# c_iu, c_iu_dag, iden=C_i_spinUP()
# c_id, c_id_dag, iden=C_i_spinDOWN()

 for i in xrange(N_x):
  for j in xrange(N_x):

   if j==N_x-1:
     ham = uni10.otimes(iden,c_iu_dag*c_iu)
     ham =ham+ uni10.otimes(iden,c_id_dag*c_id)
   else:
     ham = uni10.otimes(c_iu_dag*c_iu,iden)
     ham = ham+uni10.otimes(c_id_dag*c_id,iden)

   H.setRawElem(ham)
   H_list[i][j]=H*1.0

 return H_list



































def Q_cost_val_middle( H_col, H_row, Q_list, i, j , rho_row, rho_col, N_x, Sys,HH_long):
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
  Q_trans.setLabel([-1,2,-3,4,-5])


  Q_list[i][j-1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j-1]*1.0
  Q_transN.setLabel([6,7,8,9,-10])
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
  Q_trans.setLabel([1,-2,3,-4,-5])


  Q_list[i][j-1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j-1]*1.0
  Q_transN.setLabel([6,7,8,9,-10])

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
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])


  rho_col[i][j-1].setLabel([-5,-10,5,10])
  H_col[2*i][2*j-1].setLabel([-3,-12,3,12])

  result1=((Q_list[i][j]*H_col[2*i][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list_fermin*Q_trans)
  result=result1[0]+result
  #print "10",result1

  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j-1]*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.setLabel([6,-7,8,9,-10])

  H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=((Q_list_ferminP*H_col[2*i+1][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_trans)
  result=result1[0]+result

  #print "9",result1





####################   Long  #################################

  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([2,12,-2,13])
  Swap1.setLabel([3,13,-3,14])
  Swap2.setLabel([4,14,-4,6])

  Q_list[i][j-1].setLabel([1,-2,-3,-4,5])
  Q_list_fermin=((Q_list[i][j-1]*Swap)*Swap1)*Swap2
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([-1,2,3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])

  h_long.setLabel([-1,-12,1,12])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j-1])*(Q_list[i][j]*Q_transN)
  result=result1[0]+result

  #print "10",result1
###################################################


  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([3,12,-3,13])
  Swap1.setLabel([4,13,-4,14])
  Swap2.setLabel([6,14,-6,7])

  Q_list[i][j-1].setLabel([1,2,-3,-4,5])
  Q_list_fermin=((Q_list[i][j-1]*Swap)*Swap1)
  Q_list_fermin.permute([1,2,3,4,5,12,14],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,-2,3,4,-5,-12,-14])

  Q_list[i][j].setLabel([-6,7,8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap2))
  Q_list_fermin1.permute([6,8,9,10,14],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,8,9,-10,-14])

  h_long.setLabel([-2,-12,2,12])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j-1])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "11",result1

#################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([4,12,-4,13])
  Swap1.setLabel([6,13,-6,14])
  Swap2.setLabel([7,14,-7,8])

  Q_list[i][j-1].setLabel([1,2,3,-4,5])
  Q_list_fermin=((Q_list[i][j-1]*Swap))
  Q_list_fermin.permute([1,2,3,4,5,12,13],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-13])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap1))*Swap2
  Q_list_fermin1.permute([6,7,9,10,13],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,9,-10,-13])

  h_long.setLabel([-3,-12,3,12])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j-1])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "12",result1

####################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([6,12,-6,13])
  Swap1.setLabel([7,13,-7,14])
  Swap2.setLabel([8,14,9,-8])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_list_fermin=((Q_list[i][j-1]))
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,-8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap)*Swap1)*Swap2
  Q_list_fermin1.permute([6,7,8,10,12],5)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,8,-10,-12])

  h_long.setLabel([-4,-12,4,12])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j-1])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result
  #print "13",result1











 if j<N_x-1:

  Q_list[i][j].setLabel([1,2,3,4,5])
  Swap.setLabel([-2,-3,2,3])

  Q_list_ferm=Swap*Q_list[i][j]
  Q_list_ferm.permute([1,-2,-3,4,5],5)

  Q_list_ferm.setLabel([1,2,3,4,5])
  Q_trans=Q_list_ferm*1.0
  Q_trans.setLabel([-1,2,-3,4,-5])

  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.setLabel([6,7,8,9,-10])
 
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
  Q_trans.setLabel([1,-2,3,-4,-5])


  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.setLabel([6,7,8,9,-10])

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
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])

  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])


  H_col[2*i][2*j+1].setLabel([-3,-12,3,12])
  result1=(Q_list_fermin*H_col[2*i][2*j+1]*Q_trans)*((Q_list[i][j+1]*Q_transN)*rho_col[i][j])
  result=result1[0]+result
  #print "6", result1


  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j+1].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j+1]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.setLabel([6,-7,8,9,-10])


  H_col[2*i+1][2*j+1].setLabel([-4,-7,4,7])
  rho_col[i][j].setLabel([-5,-10,5,10])

  result1=((Q_list[i][j]*H_col[2*i+1][2*j+1]*Q_trans)*rho_col[i][j])*(Q_list_ferminP*Q_transN)
  result=result1[0]+result
  #print "5",result1








####################   Long  #################################

  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([2,12,-2,13])
  Swap1.setLabel([3,13,-3,14])
  Swap2.setLabel([4,14,-4,6])

  Q_list[i][j].setLabel([1,-2,-3,-4,5])
  Q_list_fermin=((Q_list[i][j]*Swap)*Swap1)*Swap2
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([-1,2,3,4,-5,-12,-6])

  Q_list[i][j+1].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j+1]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])

  h_long.setLabel([-1,-12,1,12])
  rho_col[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j])*(Q_list[i][j]*Q_transN)
  result=result1[0]+result

  #print "10",result1
###################################################


  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([3,12,-3,13])
  Swap1.setLabel([4,13,-4,14])
  Swap2.setLabel([6,14,-6,7])

  Q_list[i][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=((Q_list[i][j]*Swap)*Swap1)
  Q_list_fermin.permute([1,2,3,4,5,12,14],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,-2,3,4,-5,-12,-14])

  Q_list[i][j+1].setLabel([-6,7,8,9,10])
  Q_list_fermin1=((Q_list[i][j+1]*Swap2))
  Q_list_fermin1.permute([6,8,9,10,14],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,8,9,-10,-14])

  h_long.setLabel([-2,-12,2,12])
  rho_col[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "11",result1

#################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([4,12,-4,13])
  Swap1.setLabel([6,13,-6,14])
  Swap2.setLabel([7,14,-7,8])

  Q_list[i][j].setLabel([1,2,3,-4,5])
  Q_list_fermin=((Q_list[i][j]*Swap))
  Q_list_fermin.permute([1,2,3,4,5,12,13],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-13])

  Q_list[i][j+1].setLabel([-6,-7,8,9,10])
  Q_list_fermin1=((Q_list[i][j+1]*Swap1))*Swap2
  Q_list_fermin1.permute([6,7,9,10,13],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,9,-10,-13])

  h_long.setLabel([-3,-12,3,12])
  rho_col[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "12",result1

####################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([6,12,-6,13])
  Swap1.setLabel([7,13,-7,14])
  Swap2.setLabel([8,14,9,-8])

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_list_fermin=((Q_list[i][j]))
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j+1].setLabel([-6,-7,-8,9,10])
  Q_list_fermin1=((Q_list[i][j+1]*Swap)*Swap1)*Swap2
  Q_list_fermin1.permute([6,7,8,10,12],5)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,8,-10,-12])

  h_long.setLabel([-4,-12,4,12])
  rho_col[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result
  #print "13",result1












###############################################################################
 if j-1>=0 and j!=N_x-1:

  Swap.setLabel([12,4,6,-4])

  Q_list[i][j-1].setLabel([1,2,3,-4,5])
  Q_list_fermin=Q_list[i][j-1]*Swap
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-6])


  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])


  rho_col[i][j-1].setLabel([-5,-10,5,10])

  H_col[2*i][2*j-1].setLabel([-3,-12,3,12])
  result1=((Q_list[i][j]*H_col[2*i][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list_fermin*Q_trans)
  result=result1[0]+result
  #print "outb0", result1[0]
  #print "4",result1

  Swap.setLabel([7,6,-7,-6])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j-1]*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=Q_list[i][j]*Swap
  Q_list_ferminP.permute([6,7,8,9,10],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.setLabel([6,-7,8,9,-10])

  H_col[2*i+1][2*j-1].setLabel([-4,-7,4,7])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=((Q_list_ferminP*H_col[2*i+1][2*j-1]*Q_transN)*rho_col[i][j-1])*(Q_list[i][j-1]*Q_trans)
  result=result1[0]+result
  #print "3", result1




####################   Long  #################################

  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([2,12,-2,13])
  Swap1.setLabel([3,13,-3,14])
  Swap2.setLabel([4,14,-4,6])

  Q_list[i][j-1].setLabel([1,-2,-3,-4,5])
  Q_list_fermin=((Q_list[i][j-1]*Swap)*Swap1)*Swap2
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([-1,2,3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])

  h_long.setLabel([-1,-12,1,12])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j-1])*(Q_list[i][j]*Q_transN)
  result=result1[0]+result

  #print "10",result1
###################################################


  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([3,12,-3,13])
  Swap1.setLabel([4,13,-4,14])
  Swap2.setLabel([6,14,-6,7])

  Q_list[i][j-1].setLabel([1,2,-3,-4,5])
  Q_list_fermin=((Q_list[i][j-1]*Swap)*Swap1)
  Q_list_fermin.permute([1,2,3,4,5,12,14],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,-2,3,4,-5,-12,-14])

  Q_list[i][j].setLabel([-6,7,8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap2))
  Q_list_fermin1.permute([6,8,9,10,14],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,8,9,-10,-14])

  h_long.setLabel([-2,-12,2,12])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j-1])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "11",result1

#################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([4,12,-4,13])
  Swap1.setLabel([6,13,-6,14])
  Swap2.setLabel([7,14,-7,8])

  Q_list[i][j-1].setLabel([1,2,3,-4,5])
  Q_list_fermin=((Q_list[i][j-1]*Swap))
  Q_list_fermin.permute([1,2,3,4,5,12,13],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-13])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap1))*Swap2
  Q_list_fermin1.permute([6,7,9,10,13],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,9,-10,-13])

  h_long.setLabel([-3,-12,3,12])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j-1])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "12",result1

####################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([6,12,-6,13])
  Swap1.setLabel([7,13,-7,14])
  Swap2.setLabel([8,14,9,-8])

  Q_list[i][j-1].setLabel([1,2,3,4,5])
  Q_list_fermin=((Q_list[i][j-1]))
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,-8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap)*Swap1)*Swap2
  Q_list_fermin1.permute([6,7,8,10,12],5)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,8,-10,-12])

  h_long.setLabel([-4,-12,4,12])
  rho_col[i][j-1].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_col[i][j-1])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result
  #print "13",result1





#################### Row ####################################


 if i==N_x-1:

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.setLabel([1,2,3,4,-5])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([6,7,-8,-9,-10])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  H_row[2*i][2*j+1].setLabel([-8,-9,8,9])
  result1=((Q_list[i][j]*H_row[2*i][2*j+1]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
  result=result1[0]+result
  #print "1", result1


  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.setLabel([1,2,3,4,-5])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,-7,8,9,-10])

  H_row[2*i][2*j].setLabel([-6,-7,6,7])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=((Q_list[i][j]*H_row[2*i][2*j]*Q_transN)*rho_row[i-1][j])*(Q_list[i-1][j]*Q_trans)
  result=result1[0]+result
  #print "2",result1





  Swap.setLabel([15,6,16,-6])
  Swap1.setLabel([16,7,8,-7])

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
  Q_list_ferminP.permute([6,7,9,10,15],5)
  Q_transN=Q_list_ferminP*1.0
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
  Q_trans.setLabel([1,-2,3,4,-5,-12,-6])


  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])


  H_row[2*i-1][2*j].setLabel([-2,-12,2,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=((Q_list[i][j]*H_row[2*i-1][2*j]*Q_transN)*rho_row[i-1][j])*(Q_list_fermin*Q_trans)

  result=result1[0]+result
  #print "17",result1






####################   Long  #################################

  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([2,12,-2,13])
  Swap1.setLabel([3,13,-3,14])
  Swap2.setLabel([4,14,-4,6])

  Q_list[i-1][j].setLabel([1,-2,-3,-4,5])
  Q_list_fermin=((Q_list[i-1][j]*Swap)*Swap1)*Swap2
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([-1,2,3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])

  h_long.setLabel([-1,-12,1,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i-1][j])*(Q_list[i][j]*Q_transN)
  result=result1[0]+result

  #print "10",result1
###################################################


  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([3,12,-3,13])
  Swap1.setLabel([4,13,-4,14])
  Swap2.setLabel([6,14,-6,7])

  Q_list[i-1][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=((Q_list[i-1][j]*Swap)*Swap1)
  Q_list_fermin.permute([1,2,3,4,5,12,14],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,-2,3,4,-5,-12,-14])

  Q_list[i][j].setLabel([-6,7,8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap2))
  Q_list_fermin1.permute([6,8,9,10,14],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,8,9,-10,-14])

  h_long.setLabel([-2,-12,2,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i-1][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "11",result1

#################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([4,12,-4,13])
  Swap1.setLabel([6,13,-6,14])
  Swap2.setLabel([7,14,-7,8])

  Q_list[i-1][j].setLabel([1,2,3,-4,5])
  Q_list_fermin=((Q_list[i-1][j]*Swap))
  Q_list_fermin.permute([1,2,3,4,5,12,13],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-13])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap1))*Swap2
  Q_list_fermin1.permute([6,7,9,10,13],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,9,-10,-13])

  h_long.setLabel([-3,-12,3,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i-1][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "12",result1

####################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([6,12,-6,13])
  Swap1.setLabel([7,13,-7,14])
  Swap2.setLabel([8,14,9,-8])

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_list_fermin=((Q_list[i-1][j]))
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,-8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap)*Swap1)*Swap2
  Q_list_fermin1.permute([6,7,8,10,12],5)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,8,-10,-12])

  h_long.setLabel([-4,-12,4,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i-1][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result
  #print "13",result1





 if i<N_x-1:
 #######################################################
  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.setLabel([1,2,-3,-4,-5])

  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.setLabel([6,7,8,9,-10])

  rho_row[i][j].setLabel([-5,-10,5,10])
  H_row[2*i][2*j+1].setLabel([-3,-4,3,4])
  result1=(Q_list[i][j]*H_row[2*i][2*j+1]*Q_trans)*((Q_list[i+1][j]*Q_transN)*rho_row[i][j])
  result=result1[0]+result
  #print "19", result1

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i][j]*1.0
  Q_trans.setLabel([-1,-2,3,4,-5])

  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.setLabel([6,7,8,9,-10])

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
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i+1][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=(Q_list[i+1][j]*Swap)*Swap1
  Q_list_ferminP.permute([6,7,9,10,15],5)
  Q_transN=Q_list_ferminP*1.0
  Q_transN.setLabel([6,7,9,-10,-15])

  rho_row[i][j].setLabel([-5,-10,5,10])
  H_row[2*i+1][2*j+1].setLabel([-4,-15,4,15])
  result1=(Q_list[i][j]*H_row[2*i+1][2*j+1]*Q_trans)*((Q_list_ferminP*Q_transN)*rho_row[i][j])
  result=result1[0]+result
  #print "20",result1


  Swap.setLabel([12,3,13,-3])
  Swap1.setLabel([13,4,6,-4])

  Q_list[i][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=(Q_list[i][j]*Swap)*Swap1
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,-2,3,4,-5,-12,-6])


  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])

  H_row[2*i+1][2*j].setLabel([-2,-12,2,12])
  rho_row[i][j].setLabel([-5,-10,5,10])

  result1=((Q_list_fermin*H_row[2*i+1][2*j]*Q_trans)*rho_row[i][j])*(Q_list[i+1][j]*Q_transN)
  result=result1[0]+result
  #print "21",result1


 ##############################################################################################




####################   Long  #################################

  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([2,12,-2,13])
  Swap1.setLabel([3,13,-3,14])
  Swap2.setLabel([4,14,-4,6])

  Q_list[i][j].setLabel([1,-2,-3,-4,5])
  Q_list_fermin=((Q_list[i][j]*Swap)*Swap1)*Swap2
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([-1,2,3,4,-5,-12,-6])

  Q_list[i+1][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i+1][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])

  h_long.setLabel([-1,-12,1,12])
  rho_row[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i][j])*(Q_list[i+1][j]*Q_transN)
  result=result1[0]+result

  #print "10",result1
###################################################


  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([3,12,-3,13])
  Swap1.setLabel([4,13,-4,14])
  Swap2.setLabel([6,14,-6,7])

  Q_list[i][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=((Q_list[i][j]*Swap)*Swap1)
  Q_list_fermin.permute([1,2,3,4,5,12,14],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,-2,3,4,-5,-12,-14])

  Q_list[i+1][j].setLabel([-6,7,8,9,10])
  Q_list_fermin1=((Q_list[i+1][j]*Swap2))
  Q_list_fermin1.permute([6,8,9,10,14],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,8,9,-10,-14])

  h_long.setLabel([-2,-12,2,12])
  rho_row[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "11",result1

#################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([4,12,-4,13])
  Swap1.setLabel([6,13,-6,14])
  Swap2.setLabel([7,14,-7,8])

  Q_list[i][j].setLabel([1,2,3,-4,5])
  Q_list_fermin=((Q_list[i][j]*Swap))
  Q_list_fermin.permute([1,2,3,4,5,12,13],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-13])

  Q_list[i+1][j].setLabel([-6,-7,8,9,10])
  Q_list_fermin1=((Q_list[i+1][j]*Swap1))*Swap2
  Q_list_fermin1.permute([6,7,9,10,13],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,9,-10,-13])

  h_long.setLabel([-3,-12,3,12])
  rho_row[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "12",result1

####################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([6,12,-6,13])
  Swap1.setLabel([7,13,-7,14])
  Swap2.setLabel([8,14,9,-8])

  Q_list[i][j].setLabel([1,2,3,4,5])
  Q_list_fermin=((Q_list[i][j]))
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i+1][j].setLabel([-6,-7,-8,9,10])
  Q_list_fermin1=((Q_list[i+1][j]*Swap)*Swap1)*Swap2
  Q_list_fermin1.permute([6,7,8,10,12],5)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,8,-10,-12])

  h_long.setLabel([-4,-12,4,12])
  rho_row[i][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result
  #print "13",result1














 if i-1>=0 and i!=N_x-1:

  Swap.setLabel([15,6,16,-6])
  Swap1.setLabel([16,7,8,-7])

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_trans=Q_list[i-1][j]*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_ferminP=(Q_list[i][j]*Swap)*Swap1
  Q_list_ferminP.permute([6,7,9,10,15],5)
  Q_transN=Q_list_ferminP*1.0
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
  Q_trans.setLabel([1,-2,3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])


  H_row[2*i-1][2*j].setLabel([-2,-12,2,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=((Q_list[i][j]*H_row[2*i-1][2*j]*Q_transN)*rho_row[i-1][j])*(Q_list_fermin*Q_trans)
  result=result1[0]+result
  #print "23",result1





####################   Long  #################################

  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([2,12,-2,13])
  Swap1.setLabel([3,13,-3,14])
  Swap2.setLabel([4,14,-4,6])

  Q_list[i-1][j].setLabel([1,-2,-3,-4,5])
  Q_list_fermin=((Q_list[i-1][j]*Swap)*Swap1)*Swap2
  Q_list_fermin.permute([1,2,3,4,5,12,6],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([-1,2,3,4,-5,-12,-6])

  Q_list[i][j].setLabel([6,7,8,9,10])
  Q_transN=Q_list[i][j]*1.0
  Q_transN.setLabel([-6,7,8,9,-10])

  h_long.setLabel([-1,-12,1,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i-1][j])*(Q_list[i][j]*Q_transN)
  result=result1[0]+result

  #print "10",result1
###################################################


  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([3,12,-3,13])
  Swap1.setLabel([4,13,-4,14])
  Swap2.setLabel([6,14,-6,7])

  Q_list[i-1][j].setLabel([1,2,-3,-4,5])
  Q_list_fermin=((Q_list[i-1][j]*Swap)*Swap1)
  Q_list_fermin.permute([1,2,3,4,5,12,14],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,-2,3,4,-5,-12,-14])

  Q_list[i][j].setLabel([-6,7,8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap2))
  Q_list_fermin1.permute([6,8,9,10,14],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,8,9,-10,-14])

  h_long.setLabel([-2,-12,2,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i-1][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "11",result1

#################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([4,12,-4,13])
  Swap1.setLabel([6,13,-6,14])
  Swap2.setLabel([7,14,-7,8])

  Q_list[i-1][j].setLabel([1,2,3,-4,5])
  Q_list_fermin=((Q_list[i-1][j]*Swap))
  Q_list_fermin.permute([1,2,3,4,5,12,13],7)
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,-3,4,-5,-12,-13])

  Q_list[i][j].setLabel([-6,-7,8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap1))*Swap2
  Q_list_fermin1.permute([6,7,9,10,13],6)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,9,-10,-13])

  h_long.setLabel([-3,-12,3,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i-1][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result

  #print "12",result1

####################
  Swap1=Swap*1.0
  Swap2=Swap*1.0
  Swap.setLabel([6,12,-6,13])
  Swap1.setLabel([7,13,-7,14])
  Swap2.setLabel([8,14,9,-8])

  Q_list[i-1][j].setLabel([1,2,3,4,5])
  Q_list_fermin=((Q_list[i-1][j]))
  Q_trans=Q_list_fermin*1.0
  Q_trans.setLabel([1,2,3,-4,-5])

  Q_list[i][j].setLabel([-6,-7,-8,9,10])
  Q_list_fermin1=((Q_list[i][j]*Swap)*Swap1)*Swap2
  Q_list_fermin1.permute([6,7,8,10,12],5)
  Q_transN=Q_list_fermin1*1.0
  Q_transN.setLabel([6,7,8,-10,-12])

  h_long.setLabel([-4,-12,4,12])
  rho_row[i-1][j].setLabel([-5,-10,5,10])

  result1=(((Q_list_fermin*h_long)*Q_trans)*rho_row[i-1][j])*(Q_list_fermin1*Q_transN)
  result=result1[0]+result
  #print "13",result1


 return result























































































####################   Long  #################################

    Swap1=Swap*1.0
    Swap2=Swap*1.0
    Swap.setLabel([2,12,-2,13])
    Swap1.setLabel([3,13,-3,14])
    Swap2.setLabel([4,14,-4,6])

    Q_list[i][j].setLabel([1,-2,-3,-4,5])
    Q_list_fermin=((Q_list[i][j]*Swap)*Swap1)*Swap2
    Q_list_fermin.permute([1,2,3,4,5,12,6],7)
    Q_trans=Q_list_fermin*1.0
    Q_trans.setLabel([-1,2,3,4,-5,-12,-6])

    Q_list[i][j+1].setLabel([6,7,8,9,10])
    Q_transN=Q_list[i][j+1]*1.0
    Q_transN.setLabel([-6,7,8,9,-10])

    h_long.setLabel([-1,-12,1,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*h_long)*Q_trans)*(Q_list[i][j+1]*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]
###################################################
    Swap1=Swap*1.0
    Swap2=Swap*1.0
    Swap.setLabel([3,12,-3,13])
    Swap1.setLabel([4,13,-4,14])
    Swap2.setLabel([6,14,-6,7])

    Q_list[i][j].setLabel([1,2,-3,-4,5])
    Q_list_fermin=((Q_list[i][j]*Swap)*Swap1)
    Q_list_fermin.permute([1,2,3,4,5,12,14],7)
    Q_trans=Q_list_fermin*1.0
    Q_trans.setLabel([1,-2,3,4,-5,-12,-14])

    Q_list[i][j+1].setLabel([-6,7,8,9,10])
    Q_list_fermin1=((Q_list[i][j+1]*Swap2))
    Q_list_fermin1.permute([6,8,9,10,14],6)
    Q_transN=Q_list_fermin1*1.0
    Q_transN.setLabel([6,8,9,-10,-14])

    h_long.setLabel([-2,-12,2,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*h_long)*Q_trans)*(Q_list_fermin1*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]

#################
    Swap1=Swap*1.0
    Swap2=Swap*1.0
    Swap.setLabel([4,12,-4,13])
    Swap1.setLabel([6,13,-6,14])
    Swap2.setLabel([7,14,-7,8])

    Q_list[i][j].setLabel([1,2,3,-4,5])
    Q_list_fermin=((Q_list[i][j]*Swap))
    Q_list_fermin.permute([1,2,3,4,5,12,13],7)
    Q_trans=Q_list_fermin*1.0
    Q_trans.setLabel([1,2,-3,4,-5,-12,-13])

    Q_list[i][j+1].setLabel([-6,-7,8,9,10])
    Q_list_fermin1=((Q_list[i][j+1]*Swap1))*Swap2
    Q_list_fermin1.permute([6,7,9,10,13],6)
    Q_transN=Q_list_fermin1*1.0
    Q_transN.setLabel([6,7,9,-10,-13])

    h_long.setLabel([-3,-12,3,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*h_long)*Q_trans)*(Q_list_fermin1*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]

####################
    Swap1=Swap*1.0
    Swap2=Swap*1.0
    Swap.setLabel([6,12,-6,13])
    Swap1.setLabel([7,13,-7,14])
    Swap2.setLabel([8,14,9,-8])

    Q_list[i][j].setLabel([1,2,3,4,5])
    Q_list_fermin=((Q_list[i][j]))
    Q_trans=Q_list_fermin*1.0
    Q_trans.setLabel([1,2,3,-4,-5])

    Q_list[i][j+1].setLabel([-6,-7,-8,9,10])
    Q_list_fermin1=((Q_list[i][j+1]*Swap)*Swap1)*Swap2
    Q_list_fermin1.permute([6,7,8,10,12],5)
    Q_transN=Q_list_fermin1*1.0
    Q_transN.setLabel([6,7,8,-10,-12])

    h_long.setLabel([-4,-12,4,12])
    #print 2*i, 2*j+1, H_list[2*i][2*j+1]

    result=((Q_list_fermin*h_long)*Q_trans)*(Q_list_fermin1*Q_transN)
    result.permute([-5,-10,5,10],2)
    HA_list[i][j]=result+HA_list[i][j]





























def make_H_col( N_x, d_phys, h, Model):

 H_list=[None]*(N_x)
 for i in xrange(N_x):
  H_list[i]=[None]*(N_x-1) 

 if Model[0]=="Fer_Z2" or Model[0]=="Fer":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
  c_i, c_i_dag, iden=C_i(len(d_phys))

  for i in xrange(N_x):
   for j in xrange(N_x-1):
    if i==0 or i==N_x-1:
     if j==0:
      ham =uni10.otimes(c_i_dag*c_i,iden)
     elif j==N_x-2:
      ham =(2.*uni10.otimes(c_i_dag*c_i,iden)+uni10.otimes(iden,c_i_dag*c_i))
     else:
      ham =(2.*uni10.otimes(c_i_dag*c_i,iden))
    else:
     if j==0:
      ham =2.*uni10.otimes(c_i_dag*c_i,iden)
     elif j==N_x-2:
      ham =(4.*uni10.otimes(c_i_dag*c_i,iden)+2.*uni10.otimes(iden,c_i_dag*c_i))
     else:
      ham =(4.*uni10.otimes(c_i_dag*c_i,iden))

    H.setRawElem(ham)
    H_list[i][j]=H*0.50*0.25*h[2]


  for i in xrange(N_x):
   for j in xrange(N_x-1):
    if i==0 or i==N_x-1:
     if j==0:
      ham =uni10.otimes(c_i_dag*c_i,iden)
     elif j==N_x-2:
      ham =(2.*uni10.otimes(c_i_dag*c_i,iden)+uni10.otimes(iden,c_i_dag*c_i))
     else:
      ham =(2.*uni10.otimes(c_i_dag*c_i,iden))
    else:
     if j==0:
      ham =2.*uni10.otimes(c_i_dag*c_i,iden)
     elif j==N_x-2:
      ham =(4.*uni10.otimes(c_i_dag*c_i,iden)+2.*uni10.otimes(iden,c_i_dag*c_i))
     else:
      ham =(4.*uni10.otimes(c_i_dag*c_i,iden))

    H.setRawElem(ham)
    H_list[i][j]=H*0.50*h[1]+H_list[i][j]

  for i in xrange(N_x):
   for j in xrange(N_x-1):
    if i==0 or i==N_x-1:
     if j==0 or j==N_x-2:
      ham=3.0*(uni10.otimes(c_i,c_i_dag)+uni10.otimes(c_i_dag,c_i))
     else:
      ham=4.0*(uni10.otimes(c_i,c_i_dag)+uni10.otimes(c_i_dag,c_i))
    else:
     if j==0 or j==N_x-2:
      ham=6.0*(uni10.otimes(c_i,c_i_dag)+uni10.otimes(c_i_dag,c_i))
     else:
      ham=8.0*(uni10.otimes(c_i,c_i_dag)+uni10.otimes(c_i_dag,c_i))

    H.setRawElem(ham)
    H_list[i][j]=H*h[0]*0.25+H_list[i][j]



 return H_list


def  make_H_row (N_x, d_phys,h, Model):

 bdi = uni10.Bond(uni10.BD_IN, d_phys)
 bdo = uni10.Bond(uni10.BD_OUT, d_phys)
 H = uni10.UniTensor([bdi, bdi, bdo, bdo], "Heis")
 #H_list=[None]*(N_x-1)

 
 H_list=[None]*(N_x-1)
 for i in xrange(N_x-1):
  H_list[i]=[None]*N_x 


 if Model[0]=="Fer_Z2" or Model[0]=="Fer":
  bdi = uni10.Bond( uni10.BD_IN, d_phys)
  bdo = uni10.Bond( uni10.BD_OUT, d_phys)
  H = uni10.UniTensor( [bdi, bdi, bdo, bdo], "Heis")
  c_i, c_i_dag, iden=C_i(len(d_phys))
  for i in xrange(N_x-1):
   for j in xrange(N_x):
    if j==0 or j==N_x-1:
     if i==0:
      ham =uni10.otimes(c_i_dag*c_i,iden)
     elif i==N_x-2:
      ham =(2.*uni10.otimes(c_i_dag*c_i,iden)+uni10.otimes(iden,c_i_dag*c_i))
     else:
      ham =(2.*uni10.otimes(c_i_dag*c_i,iden))
    else:
     if i==0:
      ham =2.*uni10.otimes(c_i_dag*c_i,iden)
     elif i== N_x-2:
      ham =(4.*uni10.otimes(c_i_dag*c_i,iden)+2.*uni10.otimes(iden,c_i_dag*c_i))
     else:
      ham =(4.*uni10.otimes(c_i_dag*c_i,iden))



    H.setRawElem(ham)
    H_list[i][j]=H*0.50*0.25*h[2]


  for i in xrange(N_x-1):
   for j in xrange(N_x):
    if j==0 or j==N_x-1:
     if i==0:
      ham =uni10.otimes(c_i_dag*c_i,iden)
     elif i==N_x-2:
      ham =(2.*uni10.otimes(c_i_dag*c_i,iden)+uni10.otimes(iden,c_i_dag*c_i))
     else:
      ham =(2.*uni10.otimes(c_i_dag*c_i,iden))
    else:
     if i==0:
      ham =2.*uni10.otimes(c_i_dag*c_i,iden)
     elif i== N_x-2:
      ham =(4.*uni10.otimes(c_i_dag*c_i,iden)+2.*uni10.otimes(iden,c_i_dag*c_i))
     else:
      ham =(4.*uni10.otimes(c_i_dag*c_i,iden))



    H.setRawElem(ham)
    #print h[1]
    H_list[i][j]=H*0.50*h[1]+H_list[i][j]
    #print H_list[i][j]


  for i in xrange(N_x-1):
   for j in xrange(N_x):
    if j==0 or j==N_x-1:
     if i==0 or i==N_x-2:
      ham=3.0*(uni10.otimes(c_i,c_i_dag)+uni10.otimes(c_i_dag,c_i))
     else:
      ham=4.0*(uni10.otimes(c_i,c_i_dag)+uni10.otimes(c_i_dag,c_i))
    else:
     if i==0 or i==N_x-2:
      ham=6.0*(uni10.otimes(c_i,c_i_dag)+uni10.otimes(c_i_dag,c_i))
     else:
      ham=8.0*(uni10.otimes(c_i,c_i_dag)+uni10.otimes(c_i_dag,c_i))



    H.setRawElem(ham)
    H_list[i][j]=H*h[0]*0.25+H_list[i][j]




 return H_list

