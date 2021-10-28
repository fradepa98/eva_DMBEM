# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 09:56:05 2021

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

S_wall_big = 10*2.5 + 25
V_air = 9.6 * 9.6 * 3 + ((10/math.sqrt(2))*9.6*0.5)  # m³

# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_wall_big)

# convection
R_cv = 1 / (h * S_wall_big)

# thermal capacity big wall
C_wall = wall['Density'] * wall['Specific heat'] * wall['Width'] * S_wall_big

# number of temperature nodes and flow branches
no_t = no_q = sum(wall['Slices']) + 1

# Conductance matrix
R = np.zeros([no_q])
R[0] = R_cv['out'] + R_cd['Beech'] / 4
R[1] = R_cd['Beech'] / 2
R[2] = R_cd['Beech'] / 4 + R_cd['Insulation'] / 4
R[3] = R_cd['Insulation'] / 2
R[4] = R_cd['Insulation'] / 4 + R_cd['Oak']/4
R[5] = R_cd['Oak'] / 2
R[6] = R_cd['Oak'] / 4 + R_cv['in']
G = np.diag(np.reciprocal(R))

# Capacity matrix
C = np.zeros(no_t)
C = np.diag([0,C_wall['Beech'],0,C_wall['Insulation'],0,C_wall['Oak'],0])
C_det = np.linalg.det(C)

# Arc-node matrix A
A = np.eye(no_q, no_t + 1)
A = -np.diff(A, n=1, axis=1)

# Source vectors
b = np.zeros(no_q)     # temperature sources
f = np.zeros(no_t)     # heat flow sources

# Steady state solution with fixed outdoor temperature
b[0] = -4
temp_steady_To = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
np.set_printoptions(precision=3)
print('When To = -4°C, the temperatures in steady-state are:', temp_steady_To, '°C')
print(f'The indoor temperature is: {temp_steady_To[-1]:.3f} °C')

# Steady state solution with additional indoor heat flux
b[0] = -4
f[-1] = 100
temp_steady_Qh = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
print('When Qh = 10, the temperatures in steady-state are:', temp_steady_Qh, '°C')
print(f'The indoor temperature is: {temp_steady_Qh[-1]:.3f} °C')
