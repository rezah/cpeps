import pyUni10 as uni10
#import sys
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
#import pylab
#import random
import copy
#import time
#import line_profiler
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



#@profile
def setTruncation(theta, chi):
 
 LA=uni10.UniTensor(theta.bond())
 GA=uni10.UniTensor(theta.bond())
 GB=uni10.UniTensor(theta.bond())

 bond_In_list,bond_OUT_list =bond_list_dist(theta)
 #print "chi", chi
 svds = {}
 blk_qnums = theta.blockQnum()
 dim_svd=[]
 for qnum in blk_qnums:
     M_tem=theta.getBlock(qnum)
     svds[qnum] = M_tem.svd()
     dim_svd.append(int(svds[qnum][1].col()))
 svs = []
 bidxs = []
 for bidx in xrange(len(blk_qnums)):
     svs, bidxs = sv_merge1(svs, bidxs, bidx, svds[blk_qnums[bidx]][1], chi, len(blk_qnums))
 dims = [0] * len(blk_qnums)
 for bidx in bidxs:
     dims[bidx] += 1  
 qnums = []
 for bidx in xrange(len(blk_qnums)):
     qnums += [blk_qnums[bidx]] * dims[bidx]
 bdi_mid = uni10.Bond(uni10.BD_IN, qnums)
 #print bdi_mid, "chi", chi
 bdo_mid = uni10.Bond(uni10.BD_OUT, qnums)
 bond_In_list.append(bdo_mid)
 GA.assign(bond_In_list)
 bond_OUT_list.insert(0, bdi_mid)
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
 return GA, GB, LA











def sv_merge(svs, bidxs, bidx, sv_mat, chi, len_qn):
    #print "chi_in", chi
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
                    if (ori_svs[cur1] > 1.0e-12):
                     svs.append(ori_svs[cur1]) 
                     bidxs.append(ori_bidxs[cur1])
                    cur1 += 1
                else:
                    if (sv_mat[cur2] > 1.0e-12):
                     svs.append( sv_mat[cur2])
                     bidxs.append(bidx) 
                    cur2 += 1
            elif cur2 < sv_mat.elemNum() :
                for i in xrange(cnt, length):
                    if (sv_mat[cur2] > 1.0e-12):
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
       elif (sv_mat[0] > 1.0e-12):
        bidxs = [bidx] * sv_mat.elemNum()
        svs = [sv_mat[i] for i in xrange(sv_mat.elemNum())]
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

