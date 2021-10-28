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
from dm4bem import read_epw, sol_rad_tilt_surf
import tuto

rad_E = pd.read_csv(r'C:\Users\Mario\Desktop\ZHAW_modules\HS2021\EVA_DMBEM\depalfra_project\tot_rad_E.csv')

rad_ER = pd.read_csv(r'C:\Users\Mario\Desktop\ZHAW_modules\HS2021\EVA_DMBEM\depalfra_project\tot_rad_ER.csv')

rad_G = pd.read_csv(r'C:\Users\Mario\Desktop\ZHAW_modules\HS2021\EVA_DMBEM\depalfra_project\tot_rad_G.csv')

rad_N = pd.read_csv(r'C:\Users\Mario\Desktop\ZHAW_modules\HS2021\EVA_DMBEM\depalfra_project\tot_rad_N.csv')

rad_S = pd.read_csv(r'C:\Users\Mario\Desktop\ZHAW_modules\HS2021\EVA_DMBEM\depalfra_project\tot_rad_S.csv')

rad_W = pd.read_csv(r'C:\Users\Mario\Desktop\ZHAW_modules\HS2021\EVA_DMBEM\depalfra_project\tot_rad_W.csv')

rad_WR = pd.read_csv(r'C:\Users\Mario\Desktop\ZHAW_modules\HS2021\EVA_DMBEM\depalfra_project\tot_rad_WR.csv')

T_out = pd.read_csv(r'C:\Users\Mario\Desktop\ZHAW_modules\HS2021\EVA_DMBEM\depalfra_project\T_out.csv')


#%% Physical properties
# ===================
wall = {'Conductivity': [0.18,0.040,0.12,0.08,1.4],
        'Density': [800, 16, 600,2000,2500],
        'Specific heat': [2385, 1210, 2385,920,750],
        'Width': [0.1, 0.2, 0.1, 0.2,0.02],
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

#%% Big Walls

S_wall_big = 10*2.5 + 25
S_wall = 10*2.5
S_roof = S_wall + 25

V_air = 9.6 * 9.6 * 3 + ((10/math.sqrt(2))*9.6*0.5)  # m³
m_dot = 0.3*V_air*air['Density']/3600

# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_wall_big)

# convection
R_cv = 1 / (h * S_wall_big)

# thermal capacity big wall
C_wall_big = wall['Density'] * wall['Specific heat'] * wall['Width'] * S_wall_big
C_wall = wall['Density'] * wall['Specific heat'] * wall['Width'] * S_wall
C_wall_roof = wall['Density'] * wall['Specific heat'] * wall['Width'] * S_roof

# number of temperature nodes and flow branches
no_t = no_q = 8

# Conductance matrix
R = np.zeros([no_q])
R[0] = R_cv['out'] + R_cd['Beech'] / 4
R[1] = R_cd['Beech'] / 2
R[2] = R_cd['Beech'] / 4 + R_cd['Insulation'] / 4
R[3] = R_cd['Insulation'] / 2
R[4] = R_cd['Insulation'] / 4 + R_cd['Oak']/4
R[5] = R_cd['Oak'] / 2
R[6] = R_cd['Oak'] / 4 
R[7] = R_cv['in']
G_big = np.diag(np.reciprocal(R))

# Capacity matrix
C_big = np.zeros(no_t)
C_big = np.diag([0,C_wall_big['Beech'],0,C_wall_big['Insulation'],0,C_wall_big['Oak'],0,0])

# Arc-node matrix A
A_big = np.eye(no_q, no_t + 1)
A_big = -np.diff(A_big, n=1, axis=1)

# Source vectors
b = np.zeros(no_q)     # temperature sources
f = np.zeros(no_t)     # heat flow sources

# f_N = np.zeros(no_t)
# f_N[0] = rad_N
# f_N[-1] = 0.85*rad_S/6

# f_E = np.zeros(no_t)
# f_E[0] = rad_E
# f_E[-1] = 0.85*rad_S/6

# f_W = np.zeros(no_t)
# f_W[0] = rad_W
# f_W[-1] = 0.85*rad_S/6

f[0] = f[6] = 1

b[0] = 1

y = np.zeros(no_t)
y[-1] = 1

#%% Small wall

print(no_q)

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
R[6] = R_cd['Oak'] / 4
R[7] = R_cv['in']
G_small = np.diag(np.reciprocal(R))

# Capacity matrix
C_small = np.zeros(no_t)
C_small = np.diag([0,C_wall['Beech'],0,C_wall['Insulation'],0,C_wall['Oak'],0,0])

# Arc-node matrix A
A_small = np.eye(no_q, no_t + 1)
A_small = -np.diff(A_small, n=1, axis=1)

#%% Roof

# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_roof)

# convection
R_cv = 1 / (h * S_roof)
                        
# Conductance matrix
R = np.zeros([no_q])
R[0] = R_cv['out'] + R_cd['Brick'] / 4
R[1] = R_cd['Brick'] / 2
R[2] = R_cd['Brick'] / 4 + R_cd['Insulation'] / 4
R[3] = R_cd['Insulation'] / 2
R[4] = R_cd['Insulation'] / 4 + R_cd['Oak']/4
R[5] = R_cd['Oak'] / 2
R[6] = R_cd['Oak'] / 4
R[7] = R_cv['in']
G_roof = np.diag(np.reciprocal(R) )

C_roof = np.zeros(no_t)                      
C_roof = np.diag([0,C_wall_roof['Brick'],0,C_wall_roof['Insulation'],0,C_wall_roof['Oak'],0,0])

A_roof = np.eye(no_q, no_t + 1)
A_roof = -np.diff(A_roof, n=1, axis=1)

b_roof = np.zeros(no_q)
f_roof = np.zeros(no_t)

# f_ER = np.zeros(no_t)
# f_ER[0] = rad_ER
# f_ER[-1] = 0.85*rad_S/6

# f_WR = np.zeros(no_t)
# f_WR[0] = rad_WR
# f_WR[-1] = 0.85*rad_S/6

b_roof[0] = 1
f_roof[0] = f_roof[6] = 1

y_roof = np.zeros(no_t)
y_roof[-1] = 1

#%% Glass

no_t = no_q = 3

# conduction
R_cd = wall['Width'] / (wall['Conductivity'] * S_wall_big)
# G_cd = (wall['Conductivity'] * S_wall_big)/wall['Width']

# convection
R_cv = 1/(h * S_wall_big)

# G_glass = np.zeros([no_q])
# G_glass[0] = 2*h['out']*S_wall_big*wall['Conductivity']['Glass']/(2*wall['Conductivity']['Glass'] + wall['Width']['Glass']*h['out'] )
# G_glass[1] = 2*h['in']*S_wall_big*wall['Conductivity']['Glass']/(2*wall['Conductivity']['Glass'] + wall['Width']['Glass']*h['in'] )
# G_glass = np.diag(G_glass)                        

R = np.zeros([no_q])
R[0] = R_cv['out'] + R_cd['Glass'] / 2
R[1] = R_cd['Glass'] / 2
R[2] = R_cv['in']
G_glass = np.diag(np.reciprocal(R))                      
                        
C_glass = np.zeros([no_t,no_q])

A_glass = np.eye(no_q, no_t + 1)
A_glass = -np.diff(A_glass, n=1, axis=1)

b_glass = np.zeros(no_q)
f_glass = np.zeros(no_t)

b_glass[0] = 1

y_glass = np.zeros(no_t)
y_glass[-1] = 1

#%% Ventilation and Control

Gv = air['Specific heat']*m_dot
Kp = 1e-3
Qa = 2000
Ca = air['Density'] * air['Specific heat'] * V_air

A_vc = np.array([[1],
                  [1]])
G_vc = np.diag(np.array([Gv, Kp]))
b_vc = np.array([1, 1])
C_vc = np.array([Ca])
f_vc = 1
y_vc = 1

#%% Assembling circuits

TC_North_wall = {'A': A_big, 'G': G_big, 'b': b, 'C': C_big, 'f': f, 'y': y}

TC_East_wall = {'A': A_small, 'G': G_small, 'b': b, 'C': C_small, 'f': f, 'y': y}

TC_West_wall = {'A': A_small, 'G': G_small, 'b': b, 'C': C_small, 'f': f, 'y': y}

TC_South_wall = {'A': A_glass, 'G': G_glass, 'b': b_glass, 'C': C_glass, 'f': f_glass, 'y': y_glass}

TC_East_roof = {'A': A_roof, 'G': G_roof, 'b': b_roof, 'C': C_roof, 'f': f_roof, 'y': y_roof}

TC_West_roof = {'A': A_roof, 'G': G_roof, 'b': b_roof, 'C': C_roof, 'f': f_roof, 'y': y_roof}

TC_VC = {'A': A_vc, 'G': G_vc, 'b': b_vc, 'C': C_vc, 'f': f_vc, 'y': y_vc}

TCd = {'0': TC_North_wall,
       '1': TC_East_wall,
       '2': TC_West_wall,
       '3': TC_South_wall,
       '4': TC_East_roof,
       '5': TC_West_roof,
       '6': TC_VC}

# Assembling matrix

AssX = np.array([[0, 6, 1, 6],
                   [0, 6, 2, 6],
                   [0, 6, 3, 2],
                   [0, 6, 4, 6],
                   [0, 6, 5, 6],
                   [0, 6, 6, 0]])


TCa = dm4bem.TCAss(TCd, AssX)
print(TCa)

#%% Thermal circuit -> state-space
[As, Bs, Cs, Ds] = dm4bem.tc2ss(
    TCa['A'], TCa['G'], TCa['b'], TCa['C'], TCa['f'], TCa['y'])

# Maximum time-step
dtmax = min(-2. / np.linalg.eig(As)[0])
print(f'Maximum time step: {dtmax:.2f} s')

#%% Step Response
dt = 500
duration = 24*3600       # [s]

days = T_out.shape[0] / 24
# n = int(np.floor(3600 / dt * 24 * days))
n = int(np.ceil(3600 / dt * 24 * days))
no_q = As.shape[0]
no_t = As.shape[1]

# time
t = np.arange(0, n * dt, dt)

fig, axs = plt.subplots(2, 2)

u = np.block([[np.ones([1, n])],
              [np.zeros([1, n])]])

# initial values for temperatures obtained by explicit and implicit Euler
temp_exp = np.zeros([no_t, t.shape[0]])
temp_imp = np.zeros([no_t, t.shape[0]])
for k in range(t.shape[0] - 1):
    temp_exp[:, k + 1] = (np.eye(no_t) + dt * As) @\
        temp_exp[:, k] + dt * Bs @ u[:, k]
    temp_imp[:, k + 1] = np.linalg.inv(np.eye(no_t) - dt * As) @\
        (temp_imp[:, k] + dt * Bs @ u[:, k])

# axs[0, 0].plot(t / 3600, temp_exp[-1, :], t / 3600, temp_imp[-1, :])
axs[0, 0].plot(t / 3600, temp_exp[-1, :], t / 3600, temp_imp[-1, :])
# axs[0, 0].set_ylabel('Air temperature [°C]')
# axs[0, 0].set_title('Step input: To')
axs[0, 0].set(ylabel='Air temperature [°C]', title='Step input: To')              


# tuto.step_response(duration, dt, As, Bs, Cs, Ds)

# #%% Simulation with weather Data

# # Interpolate weather data for time step dt
#     data = pd.concat([T_out, rad_E, axis=1)
#     data = data.resample(str(dt) + 'S').interpolate(method='linear')
#     data = data.rename(columns={'temp_air': 'To'})

#     # Indoor temperature set-point
#     data['Ti'] = 20 * np.ones(data.shape[0])

#     # Indoor auxiliary heat flow rate
#     data['Qa'] = 0 * np.ones(data.shape[0])

#     # time
#     t = dt * np.arange(data.shape[0])

#     u = pd.concat([data['To'], data['To'], data['To'], data['Ti'],
#                    α_wSW * wall['Surface']['Concrete'] * data['Φt1'],
#                    τ_gSW * α_wSW * wall['Surface']['Glass'] * data['Φt1'],
#                    data['Qa'],
#                    α_gSW * wall['Surface']['Glass'] * data['Φt1']], axis=1)

#     # initial values for temperatures
#     temp_exp = 20 * np.ones([As.shape[0], u.shape[0]])

#     # integration in time
#     I = np.eye(As.shape[0])
#     for k in range(u.shape[0] - 1):
#         temp_exp[:, k + 1] = (I + dt * As) @ temp_exp[:, k]\
#             + dt * Bs @ u.iloc[k, :]
#     # Indoor temperature
#     y_exp = Cs @ temp_exp + Ds @ u.to_numpy().T
#     # HVAC heat flow
#     q_HVAC = Kp * (data['Ti'] - y_exp[0, :])

#     # plot indoor and outdoor temperature
#     fig, axs = plt.subplots(2, 1)
#     axs[0].plot(t / 3600, y_exp[0, :], label='$T_{indoor}$')
#     axs[0].plot(t / 3600, data['To'], label='$T_{outdoor}$')
#     axs[0].set(xlabel='Time [h]',
#                ylabel='Temperatures [°C]',
#                title='Simulation for weather')
#     axs[0].legend(loc='upper right')

#     # plot total solar radiation and HVAC heat flow
#     axs[1].plot(t / 3600,  q_HVAC, label='$q_{HVAC}$')
#     axs[1].plot(t / 3600, data['Φt1'], label='$Φ_{total}$')
#     axs[1].set(xlabel='Time [h]',
#                ylabel='Heat flows [W]')
#     axs[1].legend(loc='upper right')

#     fig.tight_layout()
