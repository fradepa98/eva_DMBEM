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
S_wall = 10*2.5
S_roof = S_wall + 25

V_air = 9.6 * 9.6 * 3 + ((10/math.sqrt(2))*9.6*0.5)  # m³

# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_wall_big)

# convection
R_cv = 1 / (h * S_wall_big)

# thermal capacity big wall
C_wall_big = wall['Density'] * wall['Specific heat'] * wall['Width'] * S_wall_big
C_wall = wall['Density'] * wall['Specific heat'] * wall['Width'] * S_wall
C_wall_roof = wall['Density'] * wall['Specific heat'] * wall['Width'] * S_roof

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
G_big = np.diag(np.reciprocal(R))

# Capacity matrix
C_big = np.zeros(no_t)
C_big = np.diag([0,C_wall_big['Beech'],0,C_wall_big['Insulation'],0,C_wall_big['Oak'],0])

# Arc-node matrix A
A_big = np.eye(no_q, no_t + 1)
A_big = -np.diff(A_big, n=1, axis=1)

# Source vectors
b = np.zeros(no_q)     # temperature sources
f = np.zeros(no_t)     # heat flow sources

# Small wall

# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_wall)

# convection
R_cv = 1 / (h * S_wall)

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

# Roof

# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_roof)

# convection
R_cv = 1 / (h * S_roof)

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
G_roof = np.diag(np.reciprocal(R) )

C_roof = np.zeros(no_t)                      
C_roof = np.diag([0,C_wall_roof['Brick'],0,C_wall_roof['Insulation'],0,C_wall_roof['Oak'],0])

A_roof = np.eye(no_q, no_t + 1)
A_roof = -np.diff(A_roof, n=1, axis=1)

b_roof = np.zeros(no_q)
f_roof = np.zeros(no_t)

# Glass

# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_wall_big)

# convection
R_cv = 1 / (h * S_wall_big)

R = np.zeros([no_q])
R[0] = R_cv['out'] + R_cd['Glass'] / 2
R[1] = R_cd['Glass'] / 2 + R_cv['in']
G_glass = np.diag(np.reciprocal(R))                        

no_t = no_q = sum(wall['Slices'][5]) + 1                      
                        
C_glass = np.zeros(no_t)

A_glass = np.eye(no_q, no_t + 1)
A_glass = -np.diff(A_glass, n=1, axis=1)

b_glass = np.zeros(no_q)
f_glass = np.zeros(no_t)

