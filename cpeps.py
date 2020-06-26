import pyUni10 as uni10
import copy
import numpy as np
from   numpy import linalg as LA
import MPSclass
import cUD as UD
import time
import math

L=6.0
N_x=6
N_y=6
La_S=(L*L)/((N_x+1.0)*(N_x+1.0))

#No-symmetry
D=[4]
chi_boundry=[20]
chi_single=[25]
chi_try=[10]
d_in=[2]
d_out=[4]

#Z2-symmetry
D=[3,2]
chi_boundry=[20]
chi_single=[120]
chi_try=[80]
d_in=[2,2]
d_out=[4,4]

#U1-symmetry
D=[1,2,1]
chi_boundry=[30]
chi_single=[60]
chi_try=[20]
d_in=[2,2]
d_out=[1,2,3,4,3,2,1]
##########################################
interval=+1.0e-2
threshold=[1.0,1.0e+1]
accuracy=+1.0e-9
Model=["Fer_U1","on"]         #ITF,ITF_Z2, Heis, Heis_Z2, Heis_U1, Fer_Z2, Fer_U1, FFI_Z2, Fer_BOS, Fer_BOS_Z2#
N_iter_total=20
N_tebd=[10, "on"]

N_part=20.0
RG_Int_Coupl_final=-16.6720
La_S=1.0

h_coupling=[-1.0/La_S, +4.0/La_S, +13.0, 2.*RG_Int_Coupl_final, 0.0]

Sys=['Fer','single', "QR", 15, "cg", "iden", 20, 140, "TEBD_SVD", "U1" ]                 #"Swap:Fer, Bos", "Opt_energy:double, single, simple", "Opt_truncation=QR, Inv, Grad", QR_update_iteration,  "Inv=SVD, cg" , "Q_init=rand, iden, part", Q_update_iteration, Grad_update_iteration,  "Previous, TEBD_SVD"

#start_itebd=La_S*4.0
start_itebd=1.0
division_itebd=5.0
N_iter_SU=[40, "full"]     #"full" or "QR"


Rg_var=["SVD", "ferm", 10 , "off", "iden", "off", "on"]                             
#"Init_abcd:SVD,randm,randQ,iden,previous","Swap:iden,ferm","Opt_iter:, >0","on, off","Iso:iden, part","Switch_2nd="on", "Fine-graing Tranfomation=off, on"
N_tebd_RG=[10, "off"]
start_itebd_varPEPS=0.001
division_itebd_varPEPS=5.0
N_iter_SU_varPEPS=[0, "full", "off"]     #"full" or "QR",  "on, off"

#Z2
D_varPEPS=[3,3]
D_simple=[4,4]
D_TEBD=[3,3]
#D_varPEPS=[6]

#U1
D_varPEPS=[2,2,2]
D_simple=[1,2,1]
D_TEBD=[2,2,2]


print Model, Sys, "h=", h_coupling, "D=", D, "d_out", d_out, "chi_single", chi_single, "chi_boundry", chi_boundry

Mag_f_list=[]
E_f_list=[]
h_list=[]
E_iter_list=[]
E_iter_list1=[]
E_iter_listRG=[]
E_iter_listRG1=[]
E_iter_list1=[]
count_list=[]

#E_iter_listQ=[]


time_list=[]
Fidel_val=1
E_coulmn=[]
E_row=[]
E_mag_coulmn=[]
E_mag_row=[]
E_0=1.0
E_1=1.0
E_00=1.0
E_11=1.0


PEPS_mps=[None]*(N_x/2)
PEPS_mps_left=[None]*(N_x/2)
PEPS_mps_right=[None]*(N_x/2)
PEPS_listten=[None]*(N_x/2)
PEPS_listtenRG=[None]*(N_x)

PEPS_mps_leftU=[None]*(N_x/2)
PEPS_mps_rightU=[None]*(N_x/2)
PEPS_listtenU=[None]*(N_x/2)

mps_boundry_left=[None]*(N_x/2)
mps_boundry_right=[None]*(N_x/2)
mps_boundry_temp=[None]*(N_x/2)

mps_boundry_left1=[None]*(N_x)
mps_boundry_right1=[None]*(N_x)
mps_boundry_temp1=[None]*(N_x)


q_D, q_chi_boundry, q_chi_single, q_chi_try, q_d_in, q_d_out=UD.full_make_bond( Model, D, chi_boundry, chi_single, chi_try, d_in, d_out)

print "q_d_in", q_d_in
print "q_d_out", q_d_out
print "q_D", q_D


for i in xrange(N_x):
 PEPS_listtenRG[i] =UD.Init_PEPS( N_x, q_D, q_d_in, i)


for i in xrange(N_x/2):
 PEPS_listten[i] =UD.Init_PEPS( N_x/2, q_D, q_d_out, i)
 PEPS_listtenU[i]=UD.Init_PEPS( N_x/2, q_D, q_d_out, i)

#PEPS_listten, norm_val, count_val=UD.Normalize_PEPS( PEPS_listten, N_x/2, q_D, q_chi_try, q_d_out, threshold, interval,Sys)




Q_list=UD.Init_Q_list( N_x/2, q_d_in, q_d_out,Sys)
H_col=UD.make_H_col( N_x, q_d_in, h_coupling, Model, L)
H_row=UD.make_H_row( N_x, q_d_in, h_coupling, Model, L)
H_long=UD.make_H_long( N_x, q_d_in, h_coupling, Model)

#print H_col[1][1]


Landa_col=UD.Landa_f_col( q_D, N_x/2)
Landa_row=UD.Landa_f_row( q_D, N_x/2)
UD.Store_Gamma( PEPS_listten, N_x/2)
UD.Store_Landa_row( Landa_row, N_x/2)
UD.Store_Landa_col( Landa_col, N_x/2)

N_col=UD.make_N_col( N_x, q_d_in, h_coupling, Model)
N_row=UD.make_N_row( N_x, q_d_in, h_coupling, Model)
N_long=UD.make_N_long( N_x, q_d_in, h_coupling, Model)

Sz_col=UD.make_Sz_col( N_x, q_d_in, h_coupling, Model)
Sz_row=UD.make_Sz_row( N_x, q_d_in, h_coupling, Model)


particle_col=UD.make_particle_col( N_x, q_d_in, h_coupling, Model)
Magz_col=UD.make_Magz_col( N_x, q_d_in, h_coupling, Model)



Magz_col_direct=UD.make_Magz_col_direct( N_x, q_d_in, h_coupling, Model)
particle_col_direct=UD.make_particle_col_direct( N_x, q_d_in, h_coupling, Model)


#for i in xrange(N_x/2):
# for j in xrange(N_y/2):
#   print "peps", i, j , PEPS_listten[i][j].printDiagram()




#UD.Store_f(PEPS_listten, N_x)
#HA_col=UD.Ascend_f_col( Q_list, H_col, N_x/2, Sys,H_long)
#HA_row=UD.Ascend_f_row( Q_list, H_row, N_x/2, Sys,H_long)
#PEPS_listten, Landa_col, Landa_row=UD.simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iter_SU, HA_col, HA_row, q_d_out, Model, N_x/2, q_D, q_chi_try, threshold, interval,Sys)

############# Start: init guess ##################
Q_list=UD.Q_update_init( H_col, H_row, N_x/2, q_D, accuracy, q_d_out, Q_list, Sys, H_long)
HA_col=UD.Ascend_f_col( Q_list, H_col, N_x/2, Sys, H_long)
HA_row=UD.Ascend_f_row( Q_list, H_row, N_x/2, Sys, H_long)


PEPS_listten, Landa_col, Landa_row=UD.simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iter_SU, HA_col, HA_row, q_d_out, Model, N_x/2, q_D, q_chi_try, threshold, interval, Sys)
PEPS_listten=UD.make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x/2)
# 
# #print Q_list[0][0].printDiagram()
# #print PEPS_listten[0][0].printDiagram()
# 
# #UD.Store_Gamma( PEPS_listten, N_x/2)
# #UD.Store_Landa_row( Landa_row, N_x/2)
# #UD.Store_Landa_col( Landa_col, N_x/2)
############# End: init guess ##################



#UD.Reload_Q_list( Q_list, N_x/2)
############# Init guess 2 ################## 
#UD.Reload_f( PEPS_listten, N_x/2)
#Q_list=UD.Init_Q_list( N_x/2, q_d_in, PEPS_listten[0][0].bond(2).Qlist(),Sys)
#UD.increase_physicalbond( Q_list, PEPS_listten, N_x/2, q_d_out)
##print PEPS_listten[2][2].printDiagram()
##print Q_list[2][2].printDiagram(), q_d_out
##Q_list=UD.Q_update_init( H_col, H_row, N_x/2, q_D, accuracy, q_d_out, Q_list, Sys, H_long)
#for i in xrange(4):

#       UD.Q_update( H_col, H_row, N_x/2, PEPS_listten, E_iter_list1, q_D, accuracy, q_d_out, q_chi_single, q_chi_try, Q_list,Sys, H_long)
#       HA_col=UD.Ascend_f_col( Q_list, H_col, N_x/2,Sys,H_long)
#       HA_row=UD.Ascend_f_row( Q_list, H_row, N_x/2,Sys,H_long)
#       NA_col=UD.Ascend_f_col( Q_list, N_col, N_x/2, Sys,N_long)
#       NA_row=UD.Ascend_f_row( Q_list, N_row, N_x/2, Sys,N_long)
#       rho_row, rho_col=UD.make_density_matrix_sinlgeLayer( PEPS_listten, N_x/2, q_chi_single, q_d_out, q_D,Sys)
#       E_0=UD.Energy_from_Density(rho_row, rho_col, HA_col, HA_row, N_x/2)
#       N_0=UD.Energy_from_Density(rho_row, rho_col, NA_col, NA_row, N_x/2)
#       print "N_0", N_0
#       print "E_final",(E_0-h_coupling[2]*N_0)
#       print "E_final",(E_0-h_coupling[2]*N_part)
#       UD.Store_Q_list(Q_list, N_x/2)
############# End: init guess ##################


#Q_list, PEPS_listten=UD.Q_feed_p( 4, PEPS_listten, Q_list)

# UD.Reload_f( PEPS_listten, N_x/2)
# UD.Reload_Q_list( Q_list, N_x/2)
# UD.increase_physicalbond( Q_list, PEPS_listten, N_x/2, q_d_out)

# PEPS_listten, norm_val, count_val=UD.Normalize_PEPS( PEPS_listten, N_x/2, q_D, q_chi_try, q_d_out, threshold, interval,Sys)
# print 'norm_val', norm_val

#UD.Reload_fRG( PEPS_listten, N_x/2)


#print Q_list[0][0].printDiagram()



for iter in xrange(N_iter_total):
 print iter
 E_min=1
 E_0=1
 E_1=100
 count_iter=0
 t0=time.time()


 UD.Q_update( H_col, H_row, N_x/2, PEPS_listten, E_iter_list1, q_D, accuracy, q_d_out, q_chi_single, q_chi_try, Q_list,Sys, H_long)

 HA_col=UD.Ascend_f_col( Q_list, H_col, N_x/2,Sys,H_long)
 HA_row=UD.Ascend_f_row( Q_list, H_row, N_x/2,Sys,H_long)

 NA_col=UD.Ascend_f_col( Q_list, N_col, N_x/2, Sys,N_long)
 NA_row=UD.Ascend_f_row( Q_list, N_row, N_x/2, Sys,N_long)

 SzA_col=UD.Ascend_f_col( Q_list, Sz_col, N_x/2, Sys,N_long)
 SzA_row=UD.Ascend_f_row( Q_list, Sz_row, N_x/2, Sys,N_long)


 if Sys[1] is "double" or Sys[1] is "simple":
  rho_row, rho_col=UD.make_density_matrix_double( PEPS_listten, N_x/2, q_chi_single, q_d_out, q_D,Sys)
 if Sys[1] is "single" or Sys[1] is "simple":
  rho_row, rho_col=UD.make_density_matrix_sinlgeLayer( PEPS_listten, N_x/2, q_chi_single, q_d_out, q_D,Sys)

 E_0=UD.Energy_from_Density(rho_row, rho_col, HA_col, HA_row, N_x/2)
 print "E_f",E_0

 N_0=UD.Energy_from_Density(rho_row, rho_col, NA_col, NA_row, N_x/2)
 print "N_f",N_0
 count_list.append(N_0)

 Sz_0=UD.Energy_from_Density(rho_row, rho_col, SzA_col, SzA_row, N_x/2)
 print "Sz_f", Sz_0



 particle_list=UD.Particle_from_Density( rho_col, particle_col, N_x, Q_list, Sys)

 Magz_list=UD.Particle_from_Density( rho_col, Magz_col, N_x, Q_list, Sys)


 print "E_final",(E_0-h_coupling[2]*N_0)*0.5
 print "E_final",(E_0-h_coupling[2]*N_part)*0.5

 E_peps=(E_0-h_coupling[2]*N_part)*0.5
 E_fermi=math.pi*(N_part/(N_x*N_x))
 E_unit=0.5*N_part*E_fermi
 print "E_total",  E_unit, E_peps/E_unit




#  print "E_U",(E_0-h_coupling[2]*N_0)*(1.0/(N_part*math.pi*(1.0/(N_x*N_x))))
#  print "E_U",(E_0-h_coupling[2]*N_part)*(1.0/(N_part*math.pi*(1.0/(N_x*N_x))))


 #print particle_list

 sum=0
 file = open("particle.txt", "w")
 for index in range(N_x):
  for index1 in range(N_x):
   file.write(str((index+1.0)*(L/(N_x+1))) + " "+str((index1+1.0)*(L/(N_x+1)))+" " + str(particle_list[index][index1])+" "+ "\n")
   sum=particle_list[index][index1]+sum
 file.close()
 print "sumN", sum



 sumZ=0
 file = open("Magz.txt", "w")
 for index in range(N_x):
  for index1 in range(N_x):
   file.write(str((index+1.0)*(L/(N_x+1))) + " "+str((index1+1.0)*(L/(N_x+1)))+" " + str(Magz_list[index][index1])+" "+ "\n")
   sumZ=Magz_list[index][index1]+sumZ
 file.close()
 print "sumZ", sumZ



 E_iter_list.append(   (E_0-h_coupling[2]*N_0)*0.5  )
 E_iter_list1.append(  (E_0-h_coupling[2]*N_part)*0.5  )
 file = open("Energy.txt", "w")
 for index in range(len(E_iter_list)):
   file.write(str(index) + " " + str(E_iter_list[index])+ " " + str(E_iter_list1[index])+" "+ "\n")
 file.close()


 count_list.append(N_0)
 file = open("NumberQ.txt", "w")
 for index in range(len(count_list)):
   file.write(str(index) + " " + str(count_list[index])+" "+ "\n")
 file.close()


 time_list.append( time.time() - t0 )
 file = open("time.txt", "w")
 for index in range(len(time_list)):
   file.write(str(index) + " " + str(time_list[index])+" "+ "\n")
 file.close()


 if Sys[1] is "single":
  UD.TEBD_Full( HA_col, HA_row, N_x/2, PEPS_listten, E_iter_list, q_D, accuracy, N_tebd, i, q_d_out, q_chi_single, q_chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model, count_list, E_iter_list1, Sys)

 if Sys[1] is "double":
  UD.TEBD_Full_double( HA_col, HA_row, N_x/2, PEPS_listten, E_iter_list, q_D, accuracy, N_tebd, i, q_d_out, q_chi_boundry, q_chi_try, mps_boundry_temp, mps_boundry_left, mps_boundry_right, threshold, interval, Model, count_list, E_iter_list1, Sys)



 if Sys[1] is "simple":

#  UD.Reload_Landa_row(Landa_row, N_x/2)
#  UD.Reload_Landa_col(Landa_col, N_x/2)
#  UD.Reload_Gamma(PEPS_listten, N_x/2)
  
  start_itebd=La_S*2.0
  division_itebd=5.0

  Landa_col=UD.Landa_f_col_iden( Landa_col,N_x/2)
  Landa_row=UD.Landa_f_row_iden( Landa_row,N_x/2)
  PEPS_listten, Landa_col, Landa_row=UD.simple_update( PEPS_listten, Landa_col, Landa_row, start_itebd, division_itebd, N_iter_SU, HA_col, HA_row, q_d_out, Model, N_x/2, q_D,q_chi_try, threshold, interval,Sys)

  UD.Store_Gamma( PEPS_listten, N_x/2)
  UD.Store_Landa_row( Landa_row, N_x/2)
  UD.Store_Landa_col( Landa_col, N_x/2)
  PEPS_listten=UD.make_PEPS_tensors(PEPS_listten, Landa_row, Landa_col,N_x/2)
  UD.Store_f(PEPS_listten, N_x/2)



 file = open("PEPSEnergy.txt", "w")
 for index in range(len(E_iter_list1)):
   file.write(str(index) + " " + str(E_iter_list1[index])+" "+ "\n")
 file.close()



NA_col=UD.Ascend_f_col( Q_list, N_col, N_x/2, Sys,N_long)
NA_row=UD.Ascend_f_row( Q_list, N_row, N_x/2, Sys,N_long)
HA_col=UD.Ascend_f_col( Q_list, H_col, N_x/2, Sys, H_long)
HA_row=UD.Ascend_f_row( Q_list, H_row, N_x/2, Sys, H_long)

rho_row, rho_col=UD.make_density_matrix_sinlgeLayer( PEPS_listten, N_x/2, q_chi_single, q_d_out, q_D, Sys)
N_0=UD.Energy_from_Density(rho_row, rho_col, NA_col, NA_row, N_x/2)
Energy_val=UD.Energy_cal( PEPS_listten, q_d_out, q_chi_single, N_x/2, q_D, HA_col, HA_row, Sys)
E_0=Energy_val
print "E_peps", E_0
print "N_f",N_0
print "E_final",(E_0-h_coupling[2]*N_0)*0.5
print "E_final",(E_0-h_coupling[2]*N_part)*0.5

q_D1, q_chi_boundry1, q_chi_single1, q_chi_try1, q_d_in1, q_d_out1=UD.full_make_bond( Model, D_varPEPS, chi_boundry, chi_single, chi_try, d_in, d_out)
q_D_simple, q_chi_boundry1, q_chi_single1, q_chi_try1, q_d_in1, q_d_out1=UD.full_make_bond( Model, D_simple, chi_boundry, chi_single, chi_try, d_in, d_out)
q_D_TEBD, q_chi_boundry1, q_chi_single1, q_chi_try1, q_d_in1, q_d_out1=UD.full_make_bond( Model, D_TEBD, chi_boundry, chi_single, chi_try, d_in, d_out)



print q_D1
if Rg_var[6]=="on":
 for i in xrange(N_x/2):
  for j in xrange(N_x/2):
   print "i,j", i, j
   PEPS_listtenRG[2*i][2*j], PEPS_listtenRG[2*i+1][2*j], PEPS_listtenRG[2*i][2*j+1] , PEPS_listtenRG[2*i+1][2*j+1]=UD.Var_PEPSQ( PEPS_listten[i][j], Q_list[i][j], q_D1, Sys, D_varPEPS, Rg_var)



#PEPS_listtenRG, norm_val, count_val=UD.Normalize_PEPS( PEPS_listtenRG, N_x, q_D, q_chi_try, q_d_in, threshold, interval,Sys)




if N_iter_SU_varPEPS[2]=="on":
 Landa_col=UD.Landa_f_col( q_D, N_x)
 Landa_row=UD.Landa_f_row( q_D, N_x)

 UD.Landa_f_row_rebonding(PEPS_listtenRG, Landa_row, N_x)
 UD.Landa_f_col_rebonding(PEPS_listtenRG, Landa_col, N_x)

 #UD.Reload_Landa_rowRG(Landa_row, N_x)
 #UD.Reload_Landa_colRG(Landa_col, N_x)
 #UD.Reload_GammaRG(PEPS_listtenRG, N_x)


 PEPS_listtenRG, Landa_col, Landa_row=UD.simple_update( PEPS_listtenRG, Landa_col, Landa_row, start_itebd_varPEPS, division_itebd_varPEPS, N_iter_SU_varPEPS, H_col, H_row, q_d_in, Model, N_x, q_D_simple, q_chi_try, threshold, interval, Sys)
 PEPS_listtenRG=UD.make_PEPS_tensors(PEPS_listtenRG, Landa_row, Landa_col,N_x)

 UD.Store_GammaRG( PEPS_listtenRG, N_x)
 UD.Store_Landa_rowRG( Landa_row, N_x)
 UD.Store_Landa_colRG( Landa_col, N_x)


#UD.Reload_fRG(PEPS_listtenRG, N_x)
#print PEPS_listtenRG[0][0].printDiagram()


rho_row, rho_col=UD.make_density_matrix_sinlgeLayer( PEPS_listtenRG, N_x, q_chi_single, q_d_in, q_D1, Sys)
Magz_col_direct=UD.Sz_from_Density_direct( rho_col, Magz_col_direct, N_x, Sys)
particle_col_direct=UD.Sz_from_Density_direct( rho_col, particle_col_direct, N_x, Sys)

sumZ=0
file = open("Zdirect.txt", "w")
for index in range(N_x):
 for index1 in range(N_x):
  file.write(str(index) + " "+str(index1)+" " + str(Magz_col_direct[index][index1])+" "+ "\n")
  sumZ=Magz_col_direct[index][index1]+sumZ
file.close()
print "sumZ", sumZ


sum=0
file = open("particle.txt", "w")
for index in range(N_x):
 for index1 in range(N_x):
  file.write(str((index+1.0)*(L/(N_x+1))) + " "+str((index1+1.0)*(L/(N_x+1)))+" " + str(particle_col_direct[index][index1])+" "+ "\n")
  sum=particle_col_direct[index][index1]+sum
file.close()
print "sumN", sum



N_0=UD.Energy_from_Density(rho_row, rho_col, N_col, N_row, N_x)
Sz_0=UD.Energy_from_Density(rho_row, rho_col, Sz_col, Sz_row, N_x)

Energy_val=UD.Energy_cal( PEPS_listtenRG, q_d_in, q_chi_single, N_x, q_D1, H_col, H_row, Sys)
E_0=Energy_val
print "E_peps", E_0
print "Sz_0",Sz_0
print "N_f",N_0

print "E_final",(E_0-h_coupling[2]*N_0)*0.5
print "E_final",(E_0-h_coupling[2]*N_part)*0.5

UD.TEBD_Full_RG( H_col, H_row, N_x, PEPS_listtenRG, E_iter_listRG, q_D_TEBD, accuracy, N_tebd_RG, i, q_d_in, q_chi_single, q_chi_try, mps_boundry_temp1, mps_boundry_left1, mps_boundry_right1, threshold, interval, Model, count_list, E_iter_listRG1, Sys,N_part,h_coupling,N_col,N_row, Sz_col, Sz_row)

UD.Store_fRG(PEPS_listtenRG, N_x)






