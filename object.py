# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 10:05:41 2017

@author: tianyun
"""
import math
import numpy as np
from numpy import pi
from coefficient import *
from gauss import *


rank=9
rank23=2
WORDLIST_FILENAME = "cross_out_with_err&spe.txt"

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
        inputFile = open(WORDLIST_FILENAME)
        try:
            self.datalist = []
            for line in inputFile:
#                self.wordlist.append((line.strip().split()[:][0],line.strip().split()[:][1:-1]))
                self.datalist.append(line.strip().split())
        finally:
            inputFile.close()
#        self.dic=dict(wordlist)
    def length(self):
        return len(self.datalist)

    def update_for_speed(self):
        """
        Return a list of data.
        Ra(deg) -> Ra(rad)
        dec(deg) -> dec(rad)
        plx(degree) -> dis
        
        print:(     ) stars have been rejected for the speed limit
        """
        dis_min=0
        dis_max=1
        vel_limit=2000
        datalist=self.datalist
        count=0
        for line in datalist:
            radeg=float(line[1])
            dedeg=float(line[2])
            plx=float(line[3])
            pmra=float(line[5])
            pmde=float(line[6])
            
            line[1]=radeg/180*pi 
            line[2]=dedeg/180*pi
            line[3]=1/plx
        
            dis=1/plx
#        
#        for i in range(0,len(datalist)-1):
#            radeg=float(datalist[i][1])
#            dedeg=float(datalist[i][2])
#            plx=float(datalist[i][3])
#            pmra=float(datalist[i][5])
#            pmde=float(datalist[i][6])
#            #Ra de dis 改写
#            datalist[i][1]=radeg/180*pi 
#            datalist[i][2]=dedeg/180*pi
#            datalist[i][3]=1/plx
#            
#            dis=1/plx
    #        print(ra,de,dis)
            rv=0
            if (rank23 == 3):
                vel=math.sqrt((4.74*dis*math.sqrt(pmra**2+pmde**2))**2+rv**2)
            elif (rank23 == 2):
                vel=4.74*dis*math.sqrt(pmra**2+pmde**2)
            else:
                print('undefined rank23')
            
            if not (dis>dis_min and dis<dis_max and vel<vel_limit):
                datalist.remove(line)
                count += 1
        print("The datalist have been updata for speed,and %i lines have been regected"%count)
        
    def update_to_galaxy(self):
        """
        Return a list of data.
        Ra(rad) -> l(rad)
        dec(rad) -> b(rad)
        pm_ra -> pm_l
        pm_dec -> pm_b
        """
        datalist=self.datalist
        for i in range(0,len(datalist)-1):
            ra=float(datalist[i][1])
            de=float(datalist[i][2])
            pmra=float(datalist[i][5])
            pmde=float(datalist[i][6])
            
            l , b = rade2lb(ra, de)
            #定义变量a不可以同时符合12，9两个不同阶数
            pm_l, pm_b = pm_trans(ra,de,b,pmra,pmde)
            datalist[i][1]=l
            datalist[i][2]=b
            datalist[i][5]=pmra
            datalist[i][6]=pmde
        print("The datalist for %i lines have been updated to galaxy"%(len(datalist)-1))

    def write_data(self):
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
            l=float(line[1])
            b=float(line[2])
            dis=float(line[3])
            pm_l=float(line[5])
            pm_b=float(line[6])
            
            if ((rank == 12) and (rank23 == 3)):
                a=value_a12_3(l, b, dis) #l,b & dis to matrix a
                line.append(a)
#                print(a)
            elif ((rank == 9) and (rank23 == 2)):
                a=value_a9_2(l, b, dis) 
                line.append(a)
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
        
        return mat,mat13,x
    
    
if __name__=="__main__":
    star=star_list()
    star.update_for_speed()
    star.update_to_galaxy()
    mat,mat13,x=star.get_parameter()
    print(x)
#-----------------------------------------

