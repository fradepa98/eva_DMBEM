# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 10:08:25 2021

@author: Admin
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem
import math


no_t = no_q = sum(wall['Slices'][2:4]) + 1                      
                        
# Conductance matrix
R = np.zeros([no_q])
R[0] = R_cv['out'] + R_cd['Brick'] / 4
R[1] = R_cd['Brick'] / 2
R[2] = R_cd['Brick'] / 4 + R_cd['Insulation'] / 4
R[3] = R_cd['Insulation'] / 2
R[4] = R_cd['Insulation'] / 4 + R_cd['Oak']/4
R[5] = R_cd['Oak'] / 2
R[6] = R_cd['Oak'] / 4 + R_cv['in']
G_roof = np.diag(np.reciprocal(R)                        
                        
C_roof = np.zeros(no_t)
C_roof = np.diag([0,C_wall['Brick'],0,C_wall['Insulation'],0,C_wall['Oak'],0])

A_roof = np.eye(no_q, no_t + 1)
A_roof = -np.diff(A_roof, n=1, axis=1)

b_roof = np.zeros(no_q)
f_roof = np.zeros(no_t)

no_t = no_q = sum(wall['Slices'][5]) + 1                      
                        
# Conductance matrix
R = np.zeros([no_q])
R[0] = R_cv['out'] + R_cd['Glass'] / 2
R[1] = R_cd['Glass'] / 2 + R_cv['in']
G_roof = np.diag(np.reciprocal(R)                        
                        
C_glass = np.zeros(no_t)
C_glass = np.diag([0,0])

A_glass = np.eye(no_q, no_t + 1)
A_glass = -np.diff(A_glass, n=1, axis=1)

b_glass = np.zeros(no_q)
f_glass = np.zeros(no_t)


            