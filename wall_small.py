# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 10:09:07 2021

@author: Francesco
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem
import math

# Physical properties
# ===================
wall = {'Conductivity': [0.18,0.040,0.12,0.08,1.4],
        'Density': [800, 16, 600,2000,2500],
        'Specific heat': [2385, 1210, 2385,920,750],
        'Width': [0.1, 0.2, 0.1, 0.2,0.004],
        'Slices': [ 2, 2, 2, 2, 1]}

wall = pd.DataFrame(wall, index=['Beech', 'Insulation', 'Oak','Brick','Glass'])

air = {'Density': 1.007,
       'Specific heat': 1150}

# altitude: 1608 m

# convection coefficients, W/m² K

# estimation of the outer convection coefficient by using the following empirical formula:

nu = 12.94e-6
Too = 268.15
beta = 1/Too
g = 9.81
DT = 10
L = 10
Gr = (g*beta*DT*(L**3))/(nu**2)
k_air = 0.046
mu = nu*1.007
Pr = (mu*1150)/k_air
Ra = Gr*Pr
Nu = 0.1*(Ra**(1/3))
h_out = Nu*k_air/L 

h = pd.DataFrame([{'in': 4., 'out': h_out}])

S_wall = 10*2.5
V_air = 9.6 * 9.6 * 3 + ((10/math.sqrt(2))*9.6*0.5)  # m³

# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_wall)

# convection
R_cv = 1 / (h * S_wall)

# thermal capacity big wall
C_wall = wall['Density'] * wall['Specific heat'] * wall['Width'] * S_wall

# number of temperature nodes and flow branches
no_t = no_q = sum(wall['Slices'][1:3]) + 1

# Conductance matrix
R = np.zeros([no_q])
R[0] = R_cv['out'] + R_cd['Beech'] / 4
R[1] = R_cd['Beech'] / 2
R[2] = R_cd['Beech'] / 4 + R_cd['Insulation'] / 4
R[3] = R_cd['Insulation'] / 2
R[4] = R_cd['Insulation'] / 4 + R_cd['Oak']/4
R[5] = R_cd['Oak'] / 2
R[6] = R_cd['Oak'] / 4 + R_cv['in']
G_small = np.diag(np.reciprocal(R))

# Capacity matrix
C_small = np.zeros(no_t)
C_small = np.diag([0,C_wall['Beech'],0,C_wall['Insulation'],0,C_wall['Oak'],0])

# Arc-node matrix A
A_small = np.eye(no_q, no_t + 1)
A_small = -np.diff(A_small, n=1, axis=1)

# Source vectors
b_small = np.zeros(no_q)     # temperature sources
f_small = np.zeros(no_t)     # heat flow sources


