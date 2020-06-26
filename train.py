import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
import MPSclass
import cUD as UD
import time
import math

for i in xrange(100):
 L=1.0
 N_x=i+1
 N_y=i+1
 #La_S=1.0/((N_x-1)*(N_x-1))
 La_S=(L*L)/((N_x+1.0)*(N_x+1.0))
 #La_S=1.0/((N_x)*(N_x))

 N_part=2.
 RG_Int_Coupl_N=(math.pi/La_S)

 RG_Int_Coupl_N=(4.0*math.pi)

 RG_Int_Coupl_D=.50-math.log(0.49758*(N_part/(((N_x)*(N_x))))**0.5)
 RG_Int_Coupl=(-1.0*RG_Int_Coupl_N)/RG_Int_Coupl_D

 print N_x, RG_Int_Coupl, math.log(0.80261*(N_part/(((N_x+1)*(N_x+1))**0.5)))
