# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 14:52:00 2017

@author: tianyun
"""
import math
import numpy as np

def value_a12_3(l,b,dis):
    a=np.zeros((12,3))
    a[0][0]=math.sin(l)/dis
    a[0][1]=math.cos(l)*math.sin(b)/dis
    a[0][2]=-math.cos(l)*math.cos(b)/dis
    a[1][0]=-math.cos(l)/dis
    a[1][1]=math.sin(l)*math.sin(b)/dis
    a[1][2]=-math.sin(l)*math.cos(b)/dis
    a[2][0]=0
    a[2][1]=-math.cos(b)/dis
    a[2][2]=-math.sin(b)/dis
    a[3][0]=-math.cos(l)*math.sin(b)
    a[3][1]=math.sin(l)
    a[3][2]=0
    a[4][0]=-math.sin(l)*math.sin(b)
    a[4][1]=-math.cos(l)
    a[4][2]=0
    a[5][0]=math.cos(b)
    a[5][1]=0
    a[5][2]=0
    a[6][0]=math.cos(2*l)*math.cos(b)
    a[6][2]=-math.sin(2*l)*math.sin(2*b)/2.0
    a[6][2]=math.sin(2*l)*((math.cos(b))**2)
    a[7][0]=-math.sin(l)*math.sin(b)
    a[7][1]=math.cos(l)*math.cos(2*b)
    a[7][2]=math.cos(l)*math.sin(2*b)
    a[8][0]=math.cos(l)*math.sin(b)
    a[8][1]=math.sin(l)*math.cos(2*b)
    a[8][2]=math.sin(l)*math.sin(2*b)
    a[9][0]=-math.sin(2*l)*math.cos(b)/2.0
    a[9][1]=-((math.cos(l))**2)*math.sin(2*b)/2.0
    a[9][2]=(math.cos(l)*math.cos(b))**2
    a[10][0]=math.sin(2*l)*math.cos(b)/2.0
    a[10][1]=-((math.sin(l))**2)*math.sin(2*b)/2.0
    a[10][2]=(math.sin(l)*math.cos(b))**2
    a[11][0]=0
    a[11][1]=math.sin(2*b)/2
    a[11][2]=(math.sin(b))**2.0
    return a

def value_a9_2(l,b,dis):
    a=np.zeros((9,2))
    a[0][0]=math.sin(l)/dis
    a[0][1]=math.cos(l)*math.sin(b)/dis
#    a[0][2]=-math.cos(l)*math.cos(b)/dis
    a[1][0]=-math.cos(l)/dis
    a[1][1]=math.sin(l)*math.sin(b)/dis
#    a[1][2]=-math.sin(l)*math.cos(b)/dis
    a[2][0]=0
    a[2][1]=-math.cos(b)/dis
#    a[2][2]=-math.sin(b)/dis
    a[3][0]=-math.cos(l)*math.sin(b)
    a[3][1]=math.sin(l)
#    a[3][2]=0
    a[4][0]=-math.sin(l)*math.sin(b)
    a[4][1]=-math.cos(l)
#    a[4][2]=0
    a[5][0]=math.cos(b)
    a[5][1]=0
#    a[5][2]=0
    a[6][0]=math.cos(2*l)*math.cos(b)
    a[6][1]=-math.sin(2*l)*math.sin(2*b)/2.0
#    a[6][2]=math.sin(2*l)*((math.cos(b))**2)
    a[7][0]=-math.sin(l)*math.sin(b)
    a[7][1]=math.cos(l)*math.cos(2*b)
#    a[7][2]=math.cos(l)*math.sin(2*b)
    a[8][0]=math.cos(l)*math.sin(b)
    a[8][1]=math.sin(l)*math.cos(2*b)
#    a[8][2]=math.sin(l)*math.sin(2*b)
#    a[9][0]=-math.sin(2*l)*math.cos(b)/2.0
#    a[9][1]=-((math.cos(l))**2)*math.sin(2*b)/2.0
#    a[9][2]=(math.cos(l)*math.cos(b))**2
#    a[10][0]=math.sin(2*l)*math.cos(b)/2.0
#    a[10][1]=-((math.sin(l))**2)*math.sin(2*b)/2.0
#    a[10][2]=(math.sin(l)*math.cos(b))**2
#    a[11][0]=0
#    a[11][1]=math.sin(2*b)/2
#    a[11][2]=(math.sin(b))**2.0
    return a

    