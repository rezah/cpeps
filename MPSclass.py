import pyUni10 as uni10
import copy
import numpy as np
from numpy import linalg as LA
#########  prerequisite functions  #############
def setTruncation(theta, chi):
 #print theta
 LA=uni10.UniTensor(theta.bond())
 GA=uni10.UniTensor(theta.bond())
 GB=uni10.UniTensor(theta.bond())

 bond_In_list,bond_OUT_list =bond_list_dist(theta)
 #print bond_In_list,bond_OUT_list
 svds = {}
 blk_qnums = theta.blockQnum()
 dim_svd=[]
 for qnum in blk_qnums:
     M_tem=theta.getBlock(qnum)
     svds[qnum] = M_tem.svd()
     #print qnum, svds[qnum][1]
     dim_svd.append(int(svds[qnum][1].col()))
 svs = []
 bidxs = []
 for bidx in xrange(len(blk_qnums)):
     svs, bidxs = sv_merge(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi, len(blk_qnums))
     #print "lol", svs, bidxs
 dims = [0] * len(blk_qnums)
 for bidx in bidxs:
     dims[bidx] += 1  
 qnums = []
 for bidx in xrange(len(blk_qnums)):
     qnums += [blk_qnums[bidx]] * dims[bidx]
 bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
 #print "Hi", chi, svs, qnums, bdi_mid
 #print bdi_mid
 bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
 bond_In_list.append(bdo_mid)
 GA.assign(bond_In_list)
 bond_OUT_list.insert(0,bdi_mid)
 GB.assign(bond_OUT_list)
 LA.assign([bdi_mid, bdo_mid])
 degs = bdi_mid.degeneracy()
 for qnum, dim in degs.iteritems():
     if qnum not in svds:
         raise Exception("In setTruncaton(): Fatal error.")
     svd = svds[qnum]
     GA.putBlock(qnum, svd[0].resize(svd[0].row(), dim))
     GB.putBlock(qnum, svd[2].resize(dim, svd[2].col()))
     LA.putBlock(qnum, svd[1].resize(dim, dim)  )
 #print LA
 return GA, GB, LA
def sv_merge(svs, bidxs, bidx, sv_mat, chi, len_qn):
    if(len(svs)):
        length = len(svs) + sv_mat.elemNum()
        length = length if length < chi else chi
        #print "length", length
        ori_svs = svs
        ori_bidxs = bidxs
        svs = [0] * length
        bidxs = [0] * length
        svs = []
        bidxs = []
        cnt  = 0
        cur1 = 0
        cur2 = 0
        while cnt <  length:
            if (cur1 < len(ori_svs)) and cur2 < sv_mat.elemNum():
                if ori_svs[cur1] >= sv_mat[cur2]:
                    if (ori_svs[cur1] > -1.0e-12):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                     cur1 += 1
                else:
                    if (sv_mat[cur2] > -1.0e-12):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                     cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > -1.0e-12):
                     if len(svs)==length:break
                     svs.append(sv_mat[cur2]) 
                     bidxs.append(bidx) 
                     cur2 += 1
                break
            else:
                for i in xrange(cur1, len(ori_svs)):
                 if len(svs)==length:break
                 svs.append(ori_svs[i])
                 bidxs.append(ori_bidxs[i]) 
                break
            cnt += 1
    else:
       if (len_qn is 1):
        bidxs = [bidx] * chi  
        svs = [sv_mat[i] for i in xrange(chi)]
       elif (sv_mat[0] > -1.0e-12):
        if sv_mat.elemNum()<chi:
         bidxs = [bidx] * sv_mat.elemNum()
         svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
        else:
         bidxs = [bidx] * chi
         svs = [sv_mat[i] for i in xrange(chi)]
       else: bidxs = [bidx];  svs = [sv_mat[0]];  
    return svs, bidxs




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







def  cal_rowcol(T):
 blk_qnums = T.blockQnum()
 row_list=[]
 col_list=[]
 for qnum in blk_qnums:
  M_tem=T.getBlock(qnum)
  col_list.append(M_tem.col())
  row_list.append(M_tem.row())
 return  sum(row_list),  sum(col_list)


def bond_list_dist(T):
 bond_list=list(T.bond())
 bond_Inlist=[]
 bond_OUTlist=[]
 for i in xrange(len(bond_list)):
  if bond_list[i].type()== 1:
   bond_Inlist.append(bond_list[i])
  else: 
   bond_OUTlist.append(bond_list[i])
 return bond_Inlist, bond_OUTlist


def   svd_parity(theta):
    #print theta,theta.getBlock().svd()
    bdo=uni10.Bond(uni10.BD_OUT,theta.bond(1).Qlist())
    bdo1=uni10.Bond(uni10.BD_OUT,theta.bond(2).Qlist())

    #A_{m<n}=U_{mm}S_{mm}V_{mn}

    LA=uni10.UniTensor([theta.bond(1), bdo])
    GA=uni10.UniTensor([theta.bond(0), theta.bond(1),bdo])
    GB=uni10.UniTensor([theta.bond(1), bdo1])
    svds = {}
    blk_qnums = theta.blockQnum()
    dim_svd=[]
    for qnum in blk_qnums:
        svds[qnum] = theta.getBlock(qnum).svd()
    for qnum in blk_qnums:
        svd = svds[qnum]
        GA.putBlock(qnum, svd[0])
        GB.putBlock(qnum, svd[2])
        LA.putBlock(qnum, svd[1])
#    print LA
    return GA, LA,GB

def   qr_parity(theta):

        #bd1=copy.copy(theta.bond(3))
        #bd1.change(uni10.BD_IN)
        bd1=uni10.Bond(uni10.BD_IN,theta.bond(1).Qlist())

        GA=uni10.UniTensor(uni10.CTYPE,[theta.bond(0),theta.bond(1)])
        LA=uni10.UniTensor(uni10.CTYPE,[bd1, theta.bond(1)])

        svds = {}
        blk_qnums = theta.blockQnum()
        dim_svd=[]
        for qnum in blk_qnums:
                svds[qnum] = theta.getBlock(qnum).qr()
                GA.putBlock(qnum, svds[qnum][0])
                LA.putBlock(qnum, svds[qnum][1])

        #    print LA
        return GA, LA
##########################################################
def    qr_parity1(theta):

        #bd1=copy.copy(theta.bond(3))
        #bd1.change(uni10.BD_IN)
        bd1=uni10.Bond(uni10.BD_IN,theta.bond(2).Qlist())

        GA=uni10.UniTensor(uni10.CTYPE,[theta.bond(0),theta.bond(1),theta.bond(2)])
        LA=uni10.UniTensor(uni10.CTYPE,[bd1, theta.bond(2)])

        svds = {}
        blk_qnums = theta.blockQnum()
        dim_svd=[]
        for qnum in blk_qnums:
                svds[qnum] = theta.getBlock(qnum).qr()
                GA.putBlock(qnum, svds[qnum][0])
                LA.putBlock(qnum, svds[qnum][1])

        #    print LA
        return GA, LA






############### MPS---Square-root #############################
class    MPS:

#################################################################################
#Use the __init__() function to assign values to object properties
#The self parameter is a reference to the class instance itself, and is used to access variables that belongs to the class.
 #@profile
 def __init__(self, physical=2, Dimension=2, Number=2, rand_fuc='rand'):
  self.N = Number
  self.D = Dimension
  self.d = physical
  self.tensor=[None]*Number
  bdi = uni10.Bond(uni10.BD_IN, self.D)
  bdo = uni10.Bond(uni10.BD_OUT, self.D)
  bdi1 = uni10.Bond(uni10.BD_IN, 1)
  bdo1 = uni10.Bond(uni10.BD_OUT, 1)
  bdi_pys = uni10.Bond(uni10.BD_IN, self.d)
  #print "Hi", bdi_pys
  A_fixed=uni10.UniTensor([bdi,bdi_pys,bdo], "A_middle")
  A_fixed.randomize()

  for i in xrange(self.N):
   if i == 0:
    A=uni10.UniTensor([bdi1,bdi_pys,bdo], "A_0")
    if rand_fuc is 'rand':
     A.randomize()
     self.tensor[i]=A
     #print "Hi0A"
    elif rand_fuc is 'ortho':
     A.orthoRand()
     #print "Hi10"
     self.tensor[i]=A
    elif rand_fuc is 'iden':
     A.identity()
     self.tensor[i]=A
     #print "Hi20"

   elif i ==((self.N)-1):
    A=uni10.UniTensor([bdi,bdi_pys,bdo1], "A_N")
    if rand_fuc is 'rand' :
     A.randomize()
     self.tensor[i]=A
    elif rand_fuc is 'ortho':
     A.orthoRand()
     self.tensor[i]=A
    elif rand_fuc is 'iden':
     #print "HIIIIIIIIIIIIIII"
     A.identity()
     self.tensor[i]=A
   else:
    A=uni10.UniTensor([bdi,bdi_pys,bdo], "A_middle")
    if rand_fuc is 'rand':
     A.randomize()
     self.tensor[i]=A
    elif rand_fuc is 'ortho':
     A.orthoRand()
     self.tensor[i]=A
    elif rand_fuc is 'randuni':
     self.tensor[i]=copy.copy(A_fixed)
    elif rand_fuc is 'iden':
     A.identity()
     self.tensor[i]=A


#################################################################################
 def __getitem__(self,i):
  return self.tensor[i]
#################################################################################
 def __setitem__(self, i, A):
  self.tensor[i]=A

#################################################################################
 def norm(self):
   A_l=[  copy.copy(self.tensor[i])  for i in xrange(len(self.tensor))  ]
   E_a=0
   for i in xrange(len(A_l)):
     #print "i", i
     if i == 0:
       bdi, bdo=bond_list_dist(A_l[i])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdi))]
       indexOUT=[ i_iter+len(bdi)+1  for i_iter in xrange(len(bdo))]
       A_l[i].setLabel(indexIN+indexOUT)
       A_l_dag=A_l[i]*1.0
       A_l_dag.transpose()
       indexOUT_minus=[   -1*indexOUT[i_iter]   for   i_iter   in   xrange(len(indexOUT))   ]
       A_l_dag.setLabel(indexOUT_minus+indexIN)
       E_a=A_l_dag*A_l[i]
       E_a.permute(indexOUT+indexOUT_minus,0)
     elif i == (len(A_l)-1):
       bdi, bdo=bond_list_dist(A_l[i-1])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdo))]
       indexIN_minus=[ -1*(i_iter+1) for i_iter in xrange(len(bdo))]
       bdi1, bdo1=bond_list_dist(A_l[i])
       index_in=indexIN+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_in_minus=indexIN_minus+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_out=[ i_iter+len(bdi1)+1  for i_iter in xrange(len(bdo1))]
       index_out_minus=[ -i_iter-len(bdi1)-1  for i_iter in xrange(len(bdo1))]
       A_l[i].setLabel(index_in+index_out)
       A_l_dag=A_l[i]*1.0
       A_l_dag.transpose()
       A_l_dag.setLabel(index_out+index_in_minus)
       E_a.setLabel(indexIN+indexIN_minus)
       E_a=A_l_dag*(E_a*A_l[i])
       #E_a.permute(index_out+index_out_minus,len(index_out))
     else:
       bdi, bdo=bond_list_dist(A_l[i-1])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdo))]
       indexIN_minus=[ -1*(i_iter+1) for i_iter in xrange(len(bdo))]
       bdi1, bdo1=bond_list_dist(A_l[i])
       index_in=indexIN+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_in_minus=indexIN_minus+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_out=[ i_iter+len(bdi1)+1  for i_iter in xrange(len(bdo1))]
       index_out_minus=[ -i_iter-len(bdi1)-1  for i_iter in xrange(len(bdo1))]
       A_l[i].setLabel(index_in+index_out)
       A_l_dag=A_l[i]*1.0
       A_l_dag.transpose()
       A_l_dag.setLabel(index_out_minus+index_in_minus)
       E_a.setLabel(indexIN+indexIN_minus)
       #print A_l_dag.printDiagram(), A_l[i].printDiagram(),  E_a.printDiagram()
       E_a=A_l_dag*(E_a*A_l[i])
       E_a.permute(index_out+index_out_minus,len(index_out))

   return E_a[0]


##################################################################################
 def __mul__(self, m):
  mps_t=MPS(self.d,self.D,self.N)
  for i in xrange(self.N):
    mps_t[i]=copy.copy(self[i])
  for i in xrange(mps_t.N):
   if m<0:
    mps_t[i]=mps_t[i]*(abs(m)**(1.0/mps_t.N))
    if i==0:
     mps_t[0]=mps_t[0]*(-1.0)
   else:
    mps_t[i]=mps_t[i]*(abs(m)**(1.0/mps_t.N))
  return mps_t

#################################################################################
 def __add__(self, mps):
  list_bond=[]
  list_bond1=[]

  for q in xrange(self.N):
   list_bond.append(self[q].bond(2).dim())

  for q in xrange(self.N):
   list_bond1.append(mps[q].bond(2).dim())

  #print "Max", max(list_bond), max(list_bond1)
  mps_f=MPS(self.d,max(list_bond)+max(list_bond1),self.N)

  for q in xrange(self.N):
   if q ==0:
    D=self[q].bond(2).dim()
    D1=mps[q].bond(2).dim()
    self_mat=self[q].getBlock()
    mps_mat=mps[q].getBlock()
    mps_f_mat=uni10.Matrix(mps_f.d,D+D1)
    mps_f_mat.set_zero()
    for i in xrange(mps_f.d):
     for j in xrange(D+D1):
       if j<D:
        mps_f_mat[i*(D+D1)+j]=self_mat[i*D+j]
       elif j>=D and (j-D)<D1:
        mps_f_mat[i*(D+D1)+j]=mps_mat[i*D1+(j-D)]

    bdi1=uni10.Bond(uni10.BD_IN, 1)
    bdi=uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo=uni10.Bond(uni10.BD_OUT, D+D1)

    mps_f[q]=uni10.UniTensor([bdi1,bdi,bdo], "A_0")
    mps_f[q].putBlock(mps_f_mat)

   elif q ==(self.N-1):
    D=self[q].bond(0).dim()
    D1=mps[q].bond(0).dim()
    self_mat=self[q].getBlock()
    mps_mat=mps[q].getBlock()
    mps_f_mat=uni10.Matrix(mps_f.d*(D+D1),1)
    mps_f_mat.set_zero()
    for i in xrange(mps_f.d):
     for j in xrange(D+D1):
       if j<D:
        mps_f_mat[j*mps_f.d+i]=self_mat[j*self.d+i]
       elif j>=D and (j-D)<D1:
        mps_f_mat[j*mps_f.d+i]=mps_mat[(j-D)*mps.d+i]

    bdi1=uni10.Bond(uni10.BD_IN, D+D1)
    bdi=uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo=uni10.Bond(uni10.BD_OUT, 1)

    mps_f[q]=uni10.UniTensor([bdi1,bdi,bdo], "A_N")
    mps_f[q].putBlock(mps_f_mat)

   else:
    D=self[q].bond(0).dim()
    Dy=self[q].bond(2).dim()
    D1=mps[q].bond(0).dim()
    Dy1=mps[q].bond(2).dim()
    self_mat=self[q].getBlock()
    mps_mat=mps[q].getBlock()
    mps_f_mat=uni10.Matrix(mps_f.d*(D+D1),Dy1+Dy)
    mps_f_mat.set_zero()
    for i in xrange(mps_f.d):
     for j in xrange(Dy1+Dy):
      for m in xrange(D+D1):
       if j<Dy and m<D:
        mps_f_mat[m*(Dy1+Dy)*mps_f.d+i*(Dy1+Dy)+j]=self_mat[m*Dy*self.d+i*Dy+j]
       elif j>=Dy and m>=D and (j-Dy)<Dy1 and (m-D)<D1:
        mps_f_mat[m*(Dy1+Dy)*mps_f.d+i*(Dy1+Dy)+j]=mps_mat[(m-D)*Dy1*mps.d+i*Dy1+(j-Dy)]

    bdi1=uni10.Bond(uni10.BD_IN, D+D1)
    bdi=uni10.Bond(uni10.BD_IN, mps_f.d)
    bdo=uni10.Bond(uni10.BD_OUT, Dy1+Dy)

    mps_f[q]=uni10.UniTensor([bdi1,bdi,bdo], "A_midle")
    mps_f[q].putBlock(mps_f_mat)

  return mps_f
###################################################################################
 def __sub__(self, other):
  return self+other*(-1.0)
#################################################################################


#########################################################################
 def fidel(self,mps):
   mps_1=MPS(mps.d,mps.D,mps.N,'randuni')
   for i in xrange(self.N):
    mps_1[i]=copy.copy(mps[i])
   self_1=MPS(self.d,self.D,self.N,'randuni')
   for i in xrange(self.N):
    self_1[i]=copy.copy(self[i])

   mps_1=mps_1.normalize()
   self_1=self_1.normalize()

   A_l=[  copy.copy(self_1.tensor[i])  for i in xrange(len(self_1.tensor))  ]
   B_l=[  copy.copy(mps_1.tensor[i])  for i in xrange(len(mps_1.tensor))  ]
   E_a=0
   for i in xrange(len(A_l)):
     if i == 0:
       bdi, bdo=bond_list_dist(A_l[i])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdi))]
       indexOUT=[ i_iter+len(bdi)+1  for i_iter in xrange(len(bdo))]
       A_l[i].setLabel(indexIN+indexOUT)
       A_l_dag=B_l[i]*1.0
       bdi, bdo=bond_list_dist(A_l_dag)
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdi))]
       indexOUT_minus=[ -1*(i_iter+len(bdi)+1)  for i_iter in xrange(len(bdo))]
       A_l_dag.transpose()
       A_l_dag.setLabel(indexOUT_minus+indexIN)

       E_a=A_l_dag*A_l[i]
       E_a.permute(indexOUT+indexOUT_minus,0)
     elif i == (len(A_l)-1):
       bdi, bdo=bond_list_dist(A_l[i-1])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdo))]
       indexIN_minus=[ -1*(i_iter+1) for i_iter in xrange(len(bdo))]
       bdi1, bdo1=bond_list_dist(A_l[i])
       index_in=indexIN+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_in_minus=indexIN_minus+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_out=[ i_iter+len(bdi1)+1  for i_iter in xrange(len(bdo1))]
       index_out_minus=[ -i_iter-len(bdi1)-1  for i_iter in xrange(len(bdo1))]
       A_l[i].setLabel(index_in+index_out)
       A_l_dag=B_l[i]*1.0
       #A_l_dag.transpose()
       #A_l_dag.setLabel(index_out+index_in_minus)
       A_l_dag.setLabel(index_in_minus+index_out)

       E_a.setLabel(indexIN+indexIN_minus)
       E_a=A_l_dag*(E_a*A_l[i])
       #E_a.permute(index_out+index_out_minus,len(index_out))
     else:
       bdi, bdo=bond_list_dist(A_l[i-1])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdo))]
       indexIN_minus=[ -1*(i_iter+1) for i_iter in xrange(len(bdo))]
       bdi1, bdo1=bond_list_dist(A_l[i])
       index_in=indexIN+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_in_minus=indexIN_minus+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_out=[ i_iter+len(bdi1)+1  for i_iter in xrange(len(bdo1))]
       index_out_minus=[ -i_iter-len(bdi1)-1  for i_iter in xrange(len(bdo1))]
       A_l[i].setLabel(index_in+index_out)

       bdi, bdo=bond_list_dist(B_l[i-1])
       indexINB=[ i_iter+1  for i_iter in xrange(len(bdo))]
       indexIN_minusB=[ -1*(i_iter+1) for i_iter in xrange(len(bdo))]
       bdi1, bdo1=bond_list_dist(B_l[i])
       index_inB=indexIN+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_in_minusB=indexIN_minus+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_outB=[ i_iter+len(bdi1)+1  for i_iter in xrange(len(bdo1))]
       index_out_minusB=[ -i_iter-len(bdi1)-1  for i_iter in xrange(len(bdo1))]

       A_l_dag=B_l[i]*1.0


       A_l_dag.transpose()
       A_l_dag.setLabel(index_out_minusB+index_in_minusB)

       E_a.setLabel(indexIN+indexIN_minusB)
       E_a=A_l_dag*(E_a*A_l[i])
       E_a.permute(index_out+index_out_minusB,len(index_out))
   
   return E_a[0]
##################################################################################


##############################Product###########################################
 def product(self, mps):
   mps_1=MPS(mps.d,mps.D,mps.N,'randuni')
   for i in xrange(self.N):
    mps_1[i]=copy.copy(mps[i])
   self_1=MPS(self.d,self.D,self.N,'randuni')
   for i in xrange(self.N):
    self_1[i]=copy.copy(self[i])
   
   A_l=[  copy.copy(self_1.tensor[i])  for i in xrange(len(self_1.tensor))  ]
   B_l=[  copy.copy(mps_1.tensor[i])  for i in xrange(len(mps_1.tensor))  ]


   for i in xrange(len(A_l)):
     if i == 0:
       bdi, bdo=bond_list_dist(A_l[i])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdi))]
       indexOUT=[ i_iter+len(bdi)+1  for i_iter in xrange(len(bdo))]
       A_l[i].setLabel(indexIN+indexOUT)
       A_l_dag=B_l[i]*1.0
       A_l_dag.transpose()
       indexOUT_minus=[   -1*indexOUT[i_iter]   for   i_iter   in   xrange(len(indexOUT))   ]
       A_l_dag.setLabel(indexOUT_minus+indexIN)
       #A_l_dag.setLabel(indexIN+indexOUT_minus)

       E_a=A_l_dag*A_l[i]
       E_a.permute(indexOUT+indexOUT_minus,len(indexIN))
     elif i == (len(A_l)-1):
       bdi, bdo=bond_list_dist(A_l[i-1])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdo))]
       indexIN_minus=[ -1*(i_iter+1) for i_iter in xrange(len(bdo))]
       bdi1, bdo1=bond_list_dist(A_l[i])
       index_in=indexIN+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_in_minus=indexIN_minus+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_out=[ i_iter+len(bdi1)+1  for i_iter in xrange(len(bdo1))]
       index_out_minus=[ -i_iter-len(bdi1)-1  for i_iter in xrange(len(bdo1))]
       A_l[i].setLabel(index_in+index_out)
       A_l_dag=B_l[i]*1.0
       A_l_dag.transpose()
       A_l_dag.setLabel(index_out+index_in_minus)
       #A_l_dag.setLabel(index_in_minus+index_out)

       E_a.setLabel(indexIN+indexIN_minus)
       E_a=A_l_dag*(E_a*A_l[i])
       #E_a.permute(index_out+index_out_minus,len(index_out))
     else:
       bdi, bdo=bond_list_dist(A_l[i-1])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdo))]
       indexIN_minus=[ -1*(i_iter+1) for i_iter in xrange(len(bdo))]
       bdi1, bdo1=bond_list_dist(A_l[i])
       index_in=indexIN+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_in_minus=indexIN_minus+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_out=[ i_iter+len(bdi1)+1  for i_iter in xrange(len(bdo1))]
       index_out_minus=[ -i_iter-len(bdi1)-1  for i_iter in xrange(len(bdo1))]
       A_l[i].setLabel(index_in+index_out)
       A_l_dag=B_l[i]*1.0
       A_l_dag.transpose()
       A_l_dag.setLabel(index_out_minus+index_in_minus)
       #A_l_dag.setLabel(index_in_minus+index_out_minus)
       E_a.setLabel(indexIN+indexIN_minus)
       E_a=A_l_dag*(E_a*A_l[i])
       E_a.permute(index_out+index_out_minus,len(index_out))

   #print E_a
   return E_a[0]









##############################Product###########################################
 def product_nonsymm(self, mps):
   mps_1=MPS(mps.d,mps.D,mps.N,'randuni')
   for i in xrange(self.N):
    mps_1[i]=copy.copy(mps[i])
   self_1=MPS(self.d,self.D,self.N,'randuni')
   for i in xrange(self.N):
    self_1[i]=copy.copy(self[i])
   
   A_l=[  copy.copy(self_1.tensor[i])  for i in xrange(len(self_1.tensor))  ]
   B_l=[  copy.copy(mps_1.tensor[i])  for i in xrange(len(mps_1.tensor))  ]


   for i in xrange(len(A_l)):
     if i == 0:
       bdi, bdo=bond_list_dist(A_l[i])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdi))]
       indexOUT=[ i_iter+len(bdi)+1  for i_iter in xrange(len(bdo))]
       A_l[i].setLabel(indexIN+indexOUT)
       A_l_dag=B_l[i]*1.0
       #A_l_dag.transpose()
       indexOUT_minus=[   -1*indexOUT[i_iter]   for   i_iter   in   xrange(len(indexOUT))   ]
       #A_l_dag.setLabel(indexOUT_minus+indexIN)
       A_l_dag.setLabel(indexIN+indexOUT_minus)

       E_a=A_l_dag*A_l[i]
       E_a.permute(indexOUT+indexOUT_minus,len(indexIN))
     elif i == (len(A_l)-1):
       bdi, bdo=bond_list_dist(A_l[i-1])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdo))]
       indexIN_minus=[ -1*(i_iter+1) for i_iter in xrange(len(bdo))]
       bdi1, bdo1=bond_list_dist(A_l[i])
       index_in=indexIN+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_in_minus=indexIN_minus+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_out=[ i_iter+len(bdi1)+1  for i_iter in xrange(len(bdo1))]
       index_out_minus=[ -i_iter-len(bdi1)-1  for i_iter in xrange(len(bdo1))]
       A_l[i].setLabel(index_in+index_out)
       A_l_dag=B_l[i]*1.0
       #A_l_dag.transpose()
       #A_l_dag.setLabel(index_out+index_in_minus)
       A_l_dag.setLabel(index_in_minus+index_out)

       E_a.setLabel(indexIN+indexIN_minus)
       E_a=A_l_dag*(E_a*A_l[i])
       #E_a.permute(index_out+index_out_minus,len(index_out))
     else:
       bdi, bdo=bond_list_dist(A_l[i-1])
       indexIN=[ i_iter+1  for i_iter in xrange(len(bdo))]
       indexIN_minus=[ -1*(i_iter+1) for i_iter in xrange(len(bdo))]
       bdi1, bdo1=bond_list_dist(A_l[i])
       index_in=indexIN+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_in_minus=indexIN_minus+[ len(indexIN)+i_iter+1  for i_iter in xrange(len(bdi1)-len(bdo))]
       index_out=[ i_iter+len(bdi1)+1  for i_iter in xrange(len(bdo1))]
       index_out_minus=[ -i_iter-len(bdi1)-1  for i_iter in xrange(len(bdo1))]
       A_l[i].setLabel(index_in+index_out)
       A_l_dag=B_l[i]*1.0
       #A_l_dag.transpose()
       #A_l_dag.setLabel(index_out_minus+index_in_minus)
       A_l_dag.setLabel(index_in_minus+index_out_minus)
       E_a.setLabel(indexIN+indexIN_minus)
       E_a=A_l_dag*(E_a*A_l[i])
       E_a.permute(index_out+index_out_minus,len(index_out))

   #print E_a
   return E_a[0]















 def distance(self,mps):
   return self.product(self)+mps.product(mps)-2.0*(mps.product(self))

####################################################################################
 def normalize(self):
  mps_t=MPS(self.d,self.D,self.N,'randuni')
  for i in xrange(self.N):
    mps_t[i]=copy.copy(self[i])
  mps_t=mps_t*(1/(mps_t.norm()**(0.5)))
  return mps_t




#####################################################################################
 #@profile
 def appSVD(self,chi):

   chi_dim=[]
   bdi = uni10.Bond(uni10.BD_IN, chi)
   degs = bdi.degeneracy()
   for qnum, dim in degs.iteritems():
    chi_dim.append(dim)

   D_dim=[]
   bdi = uni10.Bond(uni10.BD_IN, self.D)
   degs = bdi.degeneracy()
   for qnum, dim in degs.iteritems():
    D_dim.append(dim)

   #assert (sum(D_dim)>=sum(chi_dim))
   A_l=[  copy.copy(self.tensor[i])  for  i  in  xrange(len(self.tensor))  ]
   A_l_can=[  copy.copy(self.tensor[i])  for  i  in  xrange(len(self.tensor))  ]

   for i in xrange(self.N-1):
           row, colm=cal_rowcol(A_l[i])
           if (row<=colm):
            U,V,s=setTruncation(A_l[i],row)
           else:
            U,V,s=setTruncation(A_l[i],colm)
           A_l_can[i]=U*1.0
           bdi, bdo=bond_list_dist(V)
           index=[ i_iter+2  for i_iter in xrange(len(bdo))]
           index.insert(0,1)
           V.setLabel(index)
           V=s*V
           index[0]=0
           V.permute(index,1)
           bdi, bdo=bond_list_dist(A_l[i+1])
           index=[ i_iter+2  for i_iter in xrange(len(bdi)+len(bdo))]
           A_l[i+1].setLabel(index)
           A_l[i+1]=V*A_l[i+1]
           n_bonds=list(A_l[i+1].bond())
           index=[ i_iter for i_iter in xrange(len(n_bonds))  ]
           A_l[i+1].setLabel(index)
           A_l[i+1].permute( index, len(n_bonds)-len(bdo))
           A_l_can[i+1]=A_l[i+1]*1.0


   bdi, bdo=bond_list_dist(A_l_can[0])
   bdi_tem1=[]
   indexIN=[ i_iter+1  for i_iter in xrange(len(bdi))]
   indexOUT=[ i_iter+len(bdi)+1  for i_iter in xrange(len(bdo))]
   A_l_can[0].setLabel(indexIN+indexOUT)
   bdi_tem=[]
   for i in xrange(len(bdi)):
    if bdi[i].dim()== 1:
     bdi_tem.append(bdi[i])
   if len(bdi_tem)>1:
    bdi_tem1=[ bdi_tem[i] for i in xrange(len(bdi_tem)-1) ]
    Iden=uni10.UniTensor(bdi_tem1)
    Iden.identity()
    index=[ i_iter+1  for i_iter in xrange(len(bdi_tem1))]
    Iden.setLabel(index)
    A_l_can[0]=Iden*A_l_can[0]
    index=[ i_iter+1+len(bdi_tem1)  for i_iter in xrange(len(bdi)-len(bdi_tem1))]
    A_l_can[0].permute(index+indexOUT,len(index))


   bdi, bdo=bond_list_dist(A_l_can[self.N-1])
   bdo_tem1=[]
   indexIN=[ i_iter+1  for i_iter in xrange(len(bdi))]
   indexOUT=[ i_iter+len(bdi)+1  for i_iter in xrange(len(bdo))]
   A_l_can[self.N-1].setLabel(indexIN+indexOUT)
   bdo_tem=[]
   for i in xrange(len(bdo)):
    if bdo[i].dim()== 1:
     bdo_tem.append(bdo[i])
   if len(bdo_tem)>1:
    bdo_tem1=[ bdo_tem[i] for i in xrange(len(bdo_tem)-1) ]
    Iden=uni10.UniTensor(bdo_tem1)
    Iden.identity()
    index=[ i_iter+1+len(bdi)  for i_iter in xrange(len(bdo_tem1))]
    Iden.setLabel(index)
    A_l_can[self.N-1]=Iden*A_l_can[self.N-1]
    index=[ i_iter+1+len(bdo_tem1)+len(bdi)  for i_iter in xrange(len(bdo)-len(bdo_tem1))]
    A_l_can[self.N-1].permute(indexIN+index,len(indexIN))

   mps_t=MPS( self.d, self.D, self.N)
   for i in xrange(self.N):
     mps_t[i]=copy.copy(A_l_can[i])
     #print  mps_t[i].printDiagram()
   #print "inside", mps_t.norm()
   for i in xrange(self.N-1):           #self.N-1
      #print "i" , i, A_l_can[self.N-1-i].printDiagram()
      bdi, bdo=bond_list_dist(A_l_can[self.N-1-i])
      index0=[0]
      indexIN=[ i_iter+1  for i_iter in xrange(len(bdi)-1)]
      indexOUT=[ i_iter+len(bdi)  for i_iter in xrange(len(bdo))]
      A_l_can[self.N-1-i].setLabel(index0+indexIN+indexOUT)
      A_l_can[self.N-1-i].permute(indexIN+indexOUT+index0,len(indexIN)+len(indexOUT))
      row, colm=cal_rowcol(A_l_can[self.N-1-i])

      #print "i" , i, A_l_can[self.N-1-i].printDiagram()


      #print "1",row, colm, sum(chi_dim)
      if (row<=colm and row<=sum(chi_dim)):
       U,V,s=setTruncation(A_l_can[self.N-1-i],row)
      elif (row<=colm and row>sum(chi_dim)):
       #print "Y"
       U,V,s=setTruncation(A_l_can[self.N-1-i],sum(chi_dim))
      elif (row>colm and colm<=sum(chi_dim)):
       U,V,s=setTruncation(A_l_can[self.N-1-i],colm)
      elif (row>colm and colm>sum(chi_dim)):
       #print "Y"
       U,V,s=setTruncation(A_l_can[self.N-1-i],sum(chi_dim))

      row, colm=cal_rowcol(U)
      #print "2", row, colm
      if (row<=colm):
       U1,V1,s1=setTruncation(U,row)
      else:
       #print "YY"
       U1,V1,s1=setTruncation(U,colm)

      U1.setLabel(indexIN+indexOUT+index0)
      U1.permute(index0+indexIN+indexOUT,len(index0)+len(indexIN))
      A_l_can[self.N-1-i]=U1*1.0

      s1.setLabel([1,2])
      V1.setLabel([2,3])
      s.setLabel([3,4])
      V.setLabel([4,5])
      V_f=(s*V)*(s1*V1)
      V_f.permute([1,5],1)
      V_f.setLabel([-1,-5])
      
      bdi, bdo=bond_list_dist(A_l_can[self.N-i-2])
      indexIN=[ i_iter+1  for i_iter in xrange(len(bdi))]
      indexOUT=[ i_iter+len(bdi)  for i_iter in xrange(len(bdo))]

      A_l_can[self.N-i-2].setLabel(indexIN+[-5])
      A_l_can[self.N-i-2]=V_f*A_l_can[self.N-i-2]
      A_l_can[self.N-i-2].permute(indexIN+[-1],len(indexIN))




   mps_app=MPS( 2, 2, self.N)
   for i in xrange(self.N):
     mps_app[i]=A_l_can[i]*1.0

   return mps_app
###################################################################
 def appINV(self,chi):
   assert (self.D>=chi)
   A_l=[  copy.copy(self.tensor[i])  for i in xrange(len(self.tensor))  ]
   A_l_can=[  copy.copy(self.tensor[i])  for i in xrange(len(self.tensor))  ]
   return mps_app
#####################################################################################

