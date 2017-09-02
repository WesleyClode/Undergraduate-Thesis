# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 23:50:18 2017

@author: tianyun
"""

import math
import numpy as np
from numpy import pi
from coefficient import *
import matplotlib.pyplot as plt

dis_min=0.2#是距离还是在xy平面画一个圈呢？
dis_max=0.8
xydis_min=0.2
xydis_max=0.8
beishu=7
vel_limit=20000
rank=9
rank23=2
WORDLIST_FILENAME = "cross_out_with_err&spe.txt"

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


class star_list(object):
    """
    A data list.
    """
    def __init__(self):
        """
        Initializes a datalist object.

        postcondition: words are read in from a file (WORDLIST_FILENAME, a
        global constant) and stored as a list.
        """
        self.residual=np.zeros(rank23)
        self.filter=np.zeros(rank23)
        self.solution=np.zeros(rank)
        inputFile = open(WORDLIST_FILENAME)
        try:
            self.datalist = []
            for line in inputFile:
                self.datalist.append(line.strip().split())
        finally:
            inputFile.close()
            
    def length(self):
        return len(self.datalist)
    
    def update_for_xyRange(self):
        """
        平面距离约束
        
        Return a list of data.
        Ra(deg) -> Ra(rad)
        dec(deg) -> dec(rad)
        plx(degree) -> dis
        
        print:(     ) stars have been rejected for the speed limit
        
        cost:152.6
        """

        
        datalist=self.datalist
        count=0
        rejected=[]              
        for line in datalist:
            ra=float(line[1])/180*pi
            de=float(line[2])/180*pi
            dis=1/float(line[3])
            l , b = rade2lb(ra, de)
            x=math.cos(l)*math.cos(b)*dis
            y=math.sin(l)*math.cos(b)*dis
            
            r=math.sqrt(x**2+y**2)
            
            if not (r>xydis_min and r<xydis_max):
                rejected.append(line)
                count += 1
            
        for i in rejected:
            datalist.remove(i)
        print("The datalist have been updata for XYplane range limit,and %i lines have been regected"%count)
        return rejected        
        
    
    def update_for_Rdis_speed(self):
        """
        球面距离约束
        
        Return a list of data.
        Ra(deg) -> Ra(rad)
        dec(deg) -> dec(rad)
        plx(degree) -> dis
        
        print:(     ) stars have been rejected for the speed limit
        
        cost:152.6
        """

        
        datalist=self.datalist
        count=0
        rejected=[]              
        
        for line in datalist:
            plx=float(line[3])
            pmra=float(line[5])
            pmde=float(line[6])
            dis=1/plx
            rv=0
            if (rank23 == 3):
                vel=math.sqrt((4.74*dis*math.sqrt(pmra**2+pmde**2))**2+rv**2)
            elif (rank23 == 2):
                vel=4.74*dis*math.sqrt(pmra**2+pmde**2)
            else:
                print('undefined rank23')
            
            if not (dis>dis_min and dis<dis_max and vel<vel_limit):
                rejected.append(line)
                count += 1
            
        for i in rejected:
            datalist.remove(i)
        print("The datalist have been updata for speed,and %i lines have been regected"%count)
#        return rejected
    
    def write_data(self):
        """
        tyc ra de plx plx_err pm_ra pm_de V B rv
        """
        fl=open('list.txt', 'w')
        count=0
        for line in self.datalist:
            count+=1
            fl.writelines(format("  "+line[0],">14")\
              +"  "+format(line[1],">25")+"  "+format(line[2],">25")+"  "+format(line[3],">25")\
              +"  "+format(line[4],">4")+"  "+format(line[5],">9")+"  "+format(line[6],">9")\
              +"  "+format(line[7],">6")+"  "+format(line[8],">6")\
              +"  "+format(0,">1")+"\n")
        fl.close()
        print("The datalist for %i lines have been written"%count)
        
    def get_parameter(self):
        """
        After update_for_speed and update_to_galaxy.
        
        return:
        the cofficient matix is the last element of each line.
        """
        mat=np.zeros((rank,rank))
        mat13=np.zeros(rank)
        datalist=self.datalist
        for line in datalist:
            ra=float(line[1])/180*pi 
            de=float(line[2])/180*pi 
            pmra=float(line[5])
            pmde=float(line[6])
            dis=1/float(line[3])
            rv=0
            l , b = rade2lb(ra, de)
            pm_l, pm_b = pm_trans(ra,de,b,pmra,pmde)
            
            if ((rank == 12) and (rank23 == 3)):
                a=value_a12_3(l, b, dis) #l,b & dis to matrix a
#                line.append(a)
#                print(a)
            elif ((rank == 9) and (rank23 == 2)):
                a=value_a9_2(l, b, dis) 
#                line.append(a)
#                print(a)#'9参数结果：剔除了膨胀项的9参数，无视向速度'
            else: 
                print('undefined rank')
                
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
        
        mat_inv =np.linalg.inv(mat)
        x = np.dot(mat_inv, mat13)
        self.solution=x
        return mat,mat13,x
    
    def xyplot(self):
        X=list()
        Y=list()
        datalist=self.datalist
        for line in datalist:
            ra=float(line[1])/180*pi
            de=float(line[2])/180*pi
            dis=1/float(line[3])
            l , b = rade2lb(ra, de)
            x=math.cos(l)*math.cos(b)*dis
            y=math.sin(l)*math.cos(b)*dis
            X.append(x)
            Y.append(y) 
        plt.figure(0,figsize=(6,6),dpi=100) 
        plt.scatter(X,Y,s=0.1)
        plt.xlim(-1.5,1.5)
        plt.ylim(-1.5,1.5)
        plt.xlabel('X(kpc)')
        plt.ylabel('Y(Kpc)')
        plt.savefig("plane figure for %.1f-%.1f kpc.eps"%(dis_min,dis_max),dpi = 1000,bbox_inches='tight')
                
    def Residual(self):
        datalist=self.datalist
        dres=np.zeros(rank23)
        res=np.zeros(rank23)
        x=self.solution

        for line in datalist:
            ra=float(line[1])/180*pi 
            de=float(line[2])/180*pi 
            pmra=float(line[5])
            pmde=float(line[6])
            dis=1/float(line[3])
            l , b = rade2lb(ra, de)
            pm_l, pm_b = pm_trans(ra,de,b,pmra,pmde)
            rv=0
            if ((rank == 12) and (rank23 == 3)):
                a=value_a12_3(l, b, dis) #l,b & dis to matrix a
#                line.append(a)
#                print(a)
            elif ((rank == 9) and (rank23 == 2)):
                a=value_a9_2(l, b, dis) 
#                line.append(a)
#                print(a)#'9参数结果：剔除了膨胀项的9参数，无视向速度'
            else: 
                print('undefined rank')   
                
            ulbv=np.zeros(3)
            ulbv[0]=4.74*pm_l
            ulbv[1]=4.74*pm_b
            if (rank23 == 3):
                ulbv[2]=rv/dis   
            
            for m in range(0,rank23):
                for n in range(0,rank):
                    dres[m]=dres[m]+x[n]*a[n][m]
                res[m]=res[m]+(ulbv[m]-dres[m])**2
                
            
        self.residual=res
        return res
    
    def filtering(self):
        res=self.residual
        n=len(self.datalist)
        filt=np.zeros(rank23)
        for m in range(0,rank23):
            filt[m]=res[m]/n*beishu
        self.filter=filt
        
    def update_for_filter(self):
        """
        After filtering,the datalist will be updated for the digital filtering
        
        return the rejected data lines
        
        """
        filt=self.filter
        datalist=self.datalist
        x=self.solution
        dres=np.zeros(rank23)
        res=np.zeros(rank23)
        rejected=[]
        count=0
        for line in datalist:
            ra=float(line[1])/180*pi 
            de=float(line[2])/180*pi 
            pmra=float(line[5])
            pmde=float(line[6])
            dis=1/float(line[3])
            l , b = rade2lb(ra, de)
            pm_l, pm_b = pm_trans(ra,de,b,pmra,pmde)
            rv=0
            if ((rank == 12) and (rank23 == 3)):
                a=value_a12_3(l, b, dis) #l,b & dis to matrix a
#                line.append(a)
#                print(a)
            elif ((rank == 9) and (rank23 == 2)):
                a=value_a9_2(l, b, dis) 
#                line.append(a)
#                print(a)#'9参数结果：剔除了膨胀项的9参数，无视向速度'
            else: 
                print('undefined rank')   
                
            ulbv=np.zeros(3)
            ulbv[0]=4.74*pm_l
            ulbv[1]=4.74*pm_b
            if (rank23 == 3):
                ulbv[2]=rv/dis   
            
            for m in range(0,rank23):
                for n in range(0,rank):
                    dres[m]=dres[m]+x[n]*a[n][m]
                res[m]=(ulbv[m]-dres[m])**2     
            
            if not (res[1]<filt[1] and res[2]<filt[2]):
                rejected.append(line)
                count += 1
            
        for i in rejected:
            datalist.remove(i)
        print("The datalist have been updata for digital filtering,and %i lines have been regected"%count)
#        return rejected                
        
        
        
        
        
        
        
        
        
        
            
if __name__=="__main__":
    star=star_list()
    star.update_for_Rdis_speed()
    temp=0
#    star.update_to_galaxy()
    for i in range(0,10):
        if star.residual[0]>=temp:
            star.get_parameter()
            
            star.Residual()
            temp=star.residual[0]
            star.filtering()
            star.update_for_filter()
        else:
            break
        
    x=star.solution
    print(x)