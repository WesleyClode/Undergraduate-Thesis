# -*- coding: utf-8 -*-
"""
@author: tianyun
"""
import math
import numpy as np
from coefficient import *
from gauss import *
from numpy import pi

WORDLIST_FILENAME = "cross_out_with_err&spe.txt"
WORDLIST_FILENAME1 = "12.txt"
rank23=2
rank=9
dis_min=0.2
dis_max=0.8
vel_limit=20000.0

def pm_trans(ra,de,b,pmra,pmde):
    #galactic north pole
    ra_gp=192.85948/180.0*pi
    de_gp=27.12825/180.0*pi
    
    sin_phi=math.sin(ra-ra_gp)*math.cos(de_gp)/math.cos(b)
    cos_phi=(math.cos(de)*math.sin(de_gp)-math.sin(de)*math.cos(ra-ra_gp)*math.cos(de_gp))/math.cos(b)

    pm_l=cos_phi*pmra+sin_phi*pmde
    pm_b=-sin_phi*pmra+cos_phi*pmde
    return pm_l,pm_b

def rade2lb(ra,de):
    b=math.asin(math.cos(ra)*math.cos(de)*(-0.867666149)+math.sin(ra)*math.cos(de)*(-0.1980763734)+math.sin(de)*0.4559837762)
    #math.asin=double sin
    l2=(math.cos(ra)*math.cos(de)*(0.4941094279)+math.sin(ra)*math.cos(de)*(-0.44482963)+math.sin(de)*0.7439822445)/math.cos(b)
    l1=(math.cos(ra)*math.cos(de)*(-0.0548755604)+math.sin(ra)*math.cos(de)*(-0.8734370902)+math.sin(de)*(-0.4838350155))/math.cos(b)
    if l2 >= 0:
        l=math.acos(l1)
    else:
        l=-math.acos(l1)+2*pi
    return l,b

def mat_dot(A ,B):
    N = len(A[0])
    x = np.zeros(N)
    for i in range(0,N):
        for j in range(0,N):
            x[i] = x[i]+A[j][i]*B[j]
    return x
    
def write_data(list):
    fl=open('list.txt', 'w')
    for i in list:
        fl.write(i)
        fl.write("\n")
    fl.close()
    
def load_data():
    """
    Returns a list of data.
    """
    print ("Loading data from file...")
    # inFile: file
    inFile = open(WORDLIST_FILENAME, 'r')
    # wordlist: list of strings
    datalist = []
    for line in inFile:
        datalist.append(line.strip().split())
    print ("  ", len(datalist), "lines loaded.")
    return datalist
    

def get_parameter(datalist):
    """
    Return a list of data. 
    """
    mat=np.zeros((rank,rank))
    mat13=np.zeros(rank)
    count=0
    for i in range(0,len(datalist)-1):
        radeg=float(datalist[i][1])
        dedeg=float(datalist[i][2])
        plx=float(datalist[i][3])
        pmra=float(datalist[i][5])
        pmde=float(datalist[i][6])
        ra=radeg/180*pi
        de=dedeg/180*pi
        dis=1/plx
#        print(ra,de,dis)
        rv=0
        if (rank23 == 3):
            vel=math.sqrt((4.74*dis*math.sqrt(pmra**2+pmde**2))**2+rv**2)
        elif (rank23 == 2):
            vel=4.74*dis*math.sqrt(pmra**2+pmde**2)
        else:
            print('undefined rank23')
            
        global a
        
        if dis>dis_min and dis<dis_max and vel<vel_limit:
            
            count+=1
            l , b = rade2lb(ra, de)
            #定义变量a不可以同时符合12，9两个不同阶数
            if ((rank == 12) and (rank23 == 3)):
                a=value_a12_3(l, b, dis) #l,b & dis to matrix a
#                print(a)
            elif ((rank == 9) and (rank23 == 2)):
                a=value_a9_2(l, b, dis) 
#                print(a)#'9参数结果：剔除了膨胀项的9参数，无视向速度'
            else: 
                print('undefined rank')
            
            pm_l, pm_b = pm_trans(ra,de,b,pmra,pmde)
#    for i in range(0,len(datalist)-1):
#            l,b,pm_l,pm_b,ra,de,dis,rv=chang_line(datalist,i)
            
            
            for m in range(0,rank):
                for n in range(m,rank):
                    for k in range(0,rank23):
                        mat[m][n]=mat[m][n]+a[m][k]*a[n][k]
                    mat[n][m]=mat[m][n]
            
            ulbv=np.zeros(3)
            ulbv[0]=4.74*pm_l
            ulbv[1]=4.74*pm_b
            if (rank23 == 3):
                ulbv[2]=rv/dis
                
            
            for m in range(0,rank):
                for k in range(0,rank23):
                    mat13[m]=mat13[m]+ulbv[k]*a[m][k]
    print('符合距离条件的恒星数目:',count,len(datalist))
    return mat13,mat
    
    
if __name__=="__main__":
    datalist=load_data()
    mat13, mat = get_parameter(datalist)
#-----------------------------------------
    mat_inv =np.linalg.inv(mat)
    x = np.dot(mat_inv, mat13)
#    mat_inv=inv(mat) 
    print(x)           
        
    