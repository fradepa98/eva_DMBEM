# -*- coding: utf-8 -*-
"""
EVA_project_WeatherData
"""
from dm4bem import read_epw, sol_rad_tilt_surf
import numpy as np
import pandas as pd

filename = 'CHE_VS_Zermatt.067480_TMYx.2004-2018.epw'

# Data Extraction

[data, meta] = read_epw(filename, coerce_year=None)
month_year = data['month'].astype(str) + '-' + data['year'].astype(str)
print(f"Months - years in the dataset: {sorted(set(month_year))}")

weather_data = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
weather_data.index = weather_data.index.map(lambda t: t.replace(year=2012))

start_date = '2012-1-1'
end_date = '2012-2-1'

weather_data = weather_data[(weather_data.index >= start_date) & (
    weather_data.index < end_date)]
del data
weather_data.plot()

# Surface Radiation definition

albedo_ground = 0.2
albedo_snow = 0.9
albedo = 0.9*albedo_snow + 0.1*albedo_ground

# North wall

surface_orientation_N = {'slope': 90,
                       'azimuth': 180,
                       'latitude': 46}

rad_surf_N = sol_rad_tilt_surf(
    weather_data, surface_orientation_N, albedo)

rad_surf_N.plot()


B = surface_orientation_N['slope']
Z = surface_orientation_N['azimuth']
L = surface_orientation_N['latitude']

# Transform degrees in radians
B = B * np.pi / 180
Z = Z * np.pi / 180
L = L * np.pi / 180

n = weather_data.index.dayofyear

declination_angle = 23.45 * np.sin(360 * (284 + n) / 365 * np.pi / 180)
d = declination_angle * np.pi / 180

hour = weather_data.index.hour
minute = weather_data.index.minute + 60
hour_angle = 15 * ((hour + minute / 60) - 12)
h = hour_angle * np.pi / 180

theta = np.sin(d) * np.sin(L) * np.cos(B)
theta -= np.sin(d) * np.cos(L) * np.sin(B) * np.cos(Z)
theta += np.cos(d) * np.cos(L) * np.cos(B) * np.cos(h)
theta += np.cos(d) * np.sin(L) * np.sin(B) * np.cos(Z) * np.cos(h)
theta += np.cos(d) * np.sin(B) * np.sin(Z) * np.sin(h)
theta = np.array(np.arccos(theta))
theta[theta > (np.pi / 2)] = np.pi / 2

dir_rad = weather_data["dir_n_rad"] * np.cos(theta)
dir_rad[dir_rad < 0] = 0

dif_rad = weather_data["dif_h_rad"] * (1 + np.cos(B)) / 2

gamma = np.cos(d) * np.cos(L) * np.cos(h)
gamma += np.sin(d) * np.sin(L)
gamma = np.array(np.arcsin(gamma))
gamma[gamma < 1e-5] = 1e-5

dir_h_rad = weather_data["dir_n_rad"] * np.sin(gamma)

ref_rad = (dir_h_rad + weather_data["dif_h_rad"]) * albedo
ref_rad *= (1 - np.cos(B) / 2)

tot_rad_N = ref_rad + dif_rad + dir_rad
tot_rad_N.plot()

tot_rad_N.to_csv(r'C:\Users\Admin\OneDrive - Politecnico di Milano\Documenti\SCooL_Pro\ZHAW\DM4BEM\dm4bem-main\dm4bem-main\Tot_Rad\tot_rad_N.csv')

# East wall

surface_orientation_E = {'slope': 90,
                       'azimuth': -90,
                       'latitude': 46}

rad_surf_E = sol_rad_tilt_surf(
    weather_data, surface_orientation_E, albedo)

rad_surf_E.plot()


B = surface_orientation_E['slope']
Z = surface_orientation_E['azimuth']
L = surface_orientation_E['latitude']

# Transform degrees in radians
B = B * np.pi / 180
Z = Z * np.pi / 180
L = L * np.pi / 180

n = weather_data.index.dayofyear

declination_angle = 23.45 * np.sin(360 * (284 + n) / 365 * np.pi / 180)
d = declination_angle * np.pi / 180

hour = weather_data.index.hour
minute = weather_data.index.minute + 60
hour_angle = 15 * ((hour + minute / 60) - 12)
h = hour_angle * np.pi / 180

theta = np.sin(d) * np.sin(L) * np.cos(B)
theta -= np.sin(d) * np.cos(L) * np.sin(B) * np.cos(Z)
theta += np.cos(d) * np.cos(L) * np.cos(B) * np.cos(h)
theta += np.cos(d) * np.sin(L) * np.sin(B) * np.cos(Z) * np.cos(h)
theta += np.cos(d) * np.sin(B) * np.sin(Z) * np.sin(h)
theta = np.array(np.arccos(theta))
theta[theta > (np.pi / 2)] = np.pi / 2

dir_rad = weather_data["dir_n_rad"] * np.cos(theta)
dir_rad[dir_rad < 0] = 0

dif_rad = weather_data["dif_h_rad"] * (1 + np.cos(B)) / 2

gamma = np.cos(d) * np.cos(L) * np.cos(h)
gamma += np.sin(d) * np.sin(L)
gamma = np.array(np.arcsin(gamma))
gamma[gamma < 1e-5] = 1e-5

dir_h_rad = weather_data["dir_n_rad"] * np.sin(gamma)

ref_rad = (dir_h_rad + weather_data["dif_h_rad"]) * albedo
ref_rad *= (1 - np.cos(B) / 2)

tot_rad_E = ref_rad + dif_rad + dir_rad
tot_rad_E.plot()

tot_rad_E.to_csv(r'C:\Users\Admin\OneDrive - Politecnico di Milano\Documenti\SCooL_Pro\ZHAW\DM4BEM\dm4bem-main\dm4bem-main\Tot_Rad\tot_rad_E.csv')

# South wall

surface_orientation_S = {'slope': 90,
                       'azimuth': 0,
                       'latitude': 46}

rad_surf_S = sol_rad_tilt_surf(
    weather_data, surface_orientation_S, albedo)

rad_surf_S.plot()


B = surface_orientation_S['slope']
Z = surface_orientation_S['azimuth']
L = surface_orientation_S['latitude']

# Transform degrees in radians
B = B * np.pi / 180
Z = Z * np.pi / 180
L = L * np.pi / 180

n = weather_data.index.dayofyear

declination_angle = 23.45 * np.sin(360 * (284 + n) / 365 * np.pi / 180)
d = declination_angle * np.pi / 180

hour = weather_data.index.hour
minute = weather_data.index.minute + 60
hour_angle = 15 * ((hour + minute / 60) - 12)
h = hour_angle * np.pi / 180

theta = np.sin(d) * np.sin(L) * np.cos(B)
theta -= np.sin(d) * np.cos(L) * np.sin(B) * np.cos(Z)
theta += np.cos(d) * np.cos(L) * np.cos(B) * np.cos(h)
theta += np.cos(d) * np.sin(L) * np.sin(B) * np.cos(Z) * np.cos(h)
theta += np.cos(d) * np.sin(B) * np.sin(Z) * np.sin(h)
theta = np.array(np.arccos(theta))
theta[theta > (np.pi / 2)] = np.pi / 2

dir_rad = weather_data["dir_n_rad"] * np.cos(theta)
dir_rad[dir_rad < 0] = 0

dif_rad = weather_data["dif_h_rad"] * (1 + np.cos(B)) / 2

gamma = np.cos(d) * np.cos(L) * np.cos(h)
gamma += np.sin(d) * np.sin(L)
gamma = np.array(np.arcsin(gamma))
gamma[gamma < 1e-5] = 1e-5

dir_h_rad = weather_data["dir_n_rad"] * np.sin(gamma)

ref_rad = (dir_h_rad + weather_data["dif_h_rad"]) * albedo
ref_rad *= (1 - np.cos(B) / 2)

tot_rad_S = ref_rad + dif_rad + dir_rad
tot_rad_S.plot()

tot_rad_S.to_csv(r'C:\Users\Admin\OneDrive - Politecnico di Milano\Documenti\SCooL_Pro\ZHAW\DM4BEM\dm4bem-main\dm4bem-main\Tot_Rad\tot_rad_S.csv')

# West wall

surface_orientation_W = {'slope': 90,
                       'azimuth': 90,
                       'latitude': 46}

rad_surf_W = sol_rad_tilt_surf(
    weather_data, surface_orientation_W, albedo)

rad_surf_W.plot()


B = surface_orientation_W['slope']
Z = surface_orientation_W['azimuth']
L = surface_orientation_W['latitude']

# Transform degrees in radians
B = B * np.pi / 180
Z = Z * np.pi / 180
L = L * np.pi / 180

n = weather_data.index.dayofyear

declination_angle = 23.45 * np.sin(360 * (284 + n) / 365 * np.pi / 180)
d = declination_angle * np.pi / 180

hour = weather_data.index.hour
minute = weather_data.index.minute + 60
hour_angle = 15 * ((hour + minute / 60) - 12)
h = hour_angle * np.pi / 180

theta = np.sin(d) * np.sin(L) * np.cos(B)
theta -= np.sin(d) * np.cos(L) * np.sin(B) * np.cos(Z)
theta += np.cos(d) * np.cos(L) * np.cos(B) * np.cos(h)
theta += np.cos(d) * np.sin(L) * np.sin(B) * np.cos(Z) * np.cos(h)
theta += np.cos(d) * np.sin(B) * np.sin(Z) * np.sin(h)
theta = np.array(np.arccos(theta))
theta[theta > (np.pi / 2)] = np.pi / 2

dir_rad = weather_data["dir_n_rad"] * np.cos(theta)
dir_rad[dir_rad < 0] = 0

dif_rad = weather_data["dif_h_rad"] * (1 + np.cos(B)) / 2

gamma = np.cos(d) * np.cos(L) * np.cos(h)
gamma += np.sin(d) * np.sin(L)
gamma = np.array(np.arcsin(gamma))
gamma[gamma < 1e-5] = 1e-5

dir_h_rad = weather_data["dir_n_rad"] * np.sin(gamma)

ref_rad = (dir_h_rad + weather_data["dif_h_rad"]) * albedo
ref_rad *= (1 - np.cos(B) / 2)

tot_rad_W = ref_rad + dif_rad + dir_rad
tot_rad_W.plot()

tot_rad_W.to_csv(r'C:\Users\Admin\OneDrive - Politecnico di Milano\Documenti\SCooL_Pro\ZHAW\DM4BEM\dm4bem-main\dm4bem-main\Tot_Rad\tot_rad_W.csv')

# East Tilted wall

surface_orientation_ER = {'slope': 45,
                       'azimuth': -90,
                       'latitude': 46}

rad_surf_ER = sol_rad_tilt_surf(
    weather_data, surface_orientation_ER, albedo)

rad_surf_ER.plot()


B = surface_orientation_ER['slope']
Z = surface_orientation_ER['azimuth']
L = surface_orientation_ER['latitude']

# Transform degrees in radians
B = B * np.pi / 180
Z = Z * np.pi / 180
L = L * np.pi / 180

n = weather_data.index.dayofyear

declination_angle = 23.45 * np.sin(360 * (284 + n) / 365 * np.pi / 180)
d = declination_angle * np.pi / 180

hour = weather_data.index.hour
minute = weather_data.index.minute + 60
hour_angle = 15 * ((hour + minute / 60) - 12)
h = hour_angle * np.pi / 180

theta = np.sin(d) * np.sin(L) * np.cos(B)
theta -= np.sin(d) * np.cos(L) * np.sin(B) * np.cos(Z)
theta += np.cos(d) * np.cos(L) * np.cos(B) * np.cos(h)
theta += np.cos(d) * np.sin(L) * np.sin(B) * np.cos(Z) * np.cos(h)
theta += np.cos(d) * np.sin(B) * np.sin(Z) * np.sin(h)
theta = np.array(np.arccos(theta))
theta[theta > (np.pi / 2)] = np.pi / 2

dir_rad = weather_data["dir_n_rad"] * np.cos(theta)
dir_rad[dir_rad < 0] = 0

dif_rad = weather_data["dif_h_rad"] * (1 + np.cos(B)) / 2

gamma = np.cos(d) * np.cos(L) * np.cos(h)
gamma += np.sin(d) * np.sin(L)
gamma = np.array(np.arcsin(gamma))
gamma[gamma < 1e-5] = 1e-5

dir_h_rad = weather_data["dir_n_rad"] * np.sin(gamma)

ref_rad = (dir_h_rad + weather_data["dif_h_rad"]) * albedo
ref_rad *= (1 - np.cos(B) / 2)

tot_rad_ER = ref_rad + dif_rad + dir_rad
tot_rad_ER.plot()

tot_rad_ER.to_csv(r'C:\Users\Admin\OneDrive - Politecnico di Milano\Documenti\SCooL_Pro\ZHAW\DM4BEM\dm4bem-main\dm4bem-main\Tot_Rad\tot_rad_ER.csv')

# West Tilted wall

surface_orientation_WR = {'slope': 45,
                       'azimuth': 90,
                       'latitude': 46}

rad_surf_WR = sol_rad_tilt_surf(
    weather_data, surface_orientation_WR, albedo)

rad_surf_WR.plot()


B = surface_orientation_WR['slope']
Z = surface_orientation_WR['azimuth']
L = surface_orientation_WR['latitude']

# Transform degrees in radians
B = B * np.pi / 180
Z = Z * np.pi / 180
L = L * np.pi / 180

n = weather_data.index.dayofyear

declination_angle = 23.45 * np.sin(360 * (284 + n) / 365 * np.pi / 180)
d = declination_angle * np.pi / 180

hour = weather_data.index.hour
minute = weather_data.index.minute + 60
hour_angle = 15 * ((hour + minute / 60) - 12)
h = hour_angle * np.pi / 180

theta = np.sin(d) * np.sin(L) * np.cos(B)
theta -= np.sin(d) * np.cos(L) * np.sin(B) * np.cos(Z)
theta += np.cos(d) * np.cos(L) * np.cos(B) * np.cos(h)
theta += np.cos(d) * np.sin(L) * np.sin(B) * np.cos(Z) * np.cos(h)
theta += np.cos(d) * np.sin(B) * np.sin(Z) * np.sin(h)
theta = np.array(np.arccos(theta))
theta[theta > (np.pi / 2)] = np.pi / 2

dir_rad = weather_data["dir_n_rad"] * np.cos(theta)
dir_rad[dir_rad < 0] = 0

dif_rad = weather_data["dif_h_rad"] * (1 + np.cos(B)) / 2

gamma = np.cos(d) * np.cos(L) * np.cos(h)
gamma += np.sin(d) * np.sin(L)
gamma = np.array(np.arcsin(gamma))
gamma[gamma < 1e-5] = 1e-5

dir_h_rad = weather_data["dir_n_rad"] * np.sin(gamma)

ref_rad = (dir_h_rad + weather_data["dif_h_rad"]) * albedo
ref_rad *= (1 - np.cos(B) / 2)

tot_rad_WR = ref_rad + dif_rad + dir_rad
tot_rad_WR.plot()

tot_rad_WR.to_csv(r'C:\Users\Admin\OneDrive - Politecnico di Milano\Documenti\SCooL_Pro\ZHAW\DM4BEM\dm4bem-main\dm4bem-main\Tot_Rad\tot_rad_WR.csv')

# Ground Wall

surface_orientation_G = {'slope': 180,
                       'azimuth': 90,
                       'latitude': 46}

rad_surf_G = sol_rad_tilt_surf(
    weather_data, surface_orientation_G, albedo)

B = surface_orientation_G['slope']
Z = surface_orientation_G['azimuth']
L = surface_orientation_G['latitude']

# Transform degrees in radians
B = B * np.pi / 180
Z = Z * np.pi / 180
L = L * np.pi / 180

n = weather_data.index.dayofyear

declination_angle = 23.45 * np.sin(360 * (284 + n) / 365 * np.pi / 180)
d = declination_angle * np.pi / 180

hour = weather_data.index.hour
minute = weather_data.index.minute + 60
hour_angle = 15 * ((hour + minute / 60) - 12)
h = hour_angle * np.pi / 180

theta = np.sin(d) * np.sin(L) * np.cos(B)
theta -= np.sin(d) * np.cos(L) * np.sin(B) * np.cos(Z)
theta += np.cos(d) * np.cos(L) * np.cos(B) * np.cos(h)
theta += np.cos(d) * np.sin(L) * np.sin(B) * np.cos(Z) * np.cos(h)
theta += np.cos(d) * np.sin(B) * np.sin(Z) * np.sin(h)
theta = np.array(np.arccos(theta))
theta[theta > (np.pi / 2)] = np.pi / 2

dir_rad = weather_data["dir_n_rad"] * np.cos(theta)
dir_rad[dir_rad < 0] = 0

dif_rad = weather_data["dif_h_rad"] * (1 + np.cos(B)) / 2

gamma = np.cos(d) * np.cos(L) * np.cos(h)
gamma += np.sin(d) * np.sin(L)
gamma = np.array(np.arcsin(gamma))
gamma[gamma < 1e-5] = 1e-5

dir_h_rad = weather_data["dir_n_rad"] * np.sin(gamma)

ref_rad = (dir_h_rad + weather_data["dif_h_rad"]) * albedo
ref_rad *= (1 - np.cos(B) / 2)

tot_rad_G = ref_rad + dif_rad + dir_rad
tot_rad_G.plot()

tot_rad_WR.to_csv(r'C:\Users\Admin\OneDrive - Politecnico di Milano\Documenti\SCooL_Pro\ZHAW\DM4BEM\dm4bem-main\dm4bem-main\Tot_Rad\tot_rad_G.csv')