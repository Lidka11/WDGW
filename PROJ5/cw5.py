# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 16:32:41 2023

@author: mgrzy
"""

import numpy as np
import matplotlib.pyplot as plt

import gpxpy
import plotly.graph_objects as go
import plotly.io as pio
import gpxpy.gpx
import pandas as pd
import tkinter as tk
from tkinter import ttk



def gpx2df(gpx_file_name, new_file_name = 'route_df.csv'):
    with open(gpx_file_name, 'r') as gpx_file:
        gpx = gpxpy.parse(gpx_file)
    route_info = []
    
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                route_info.append({
                    'latitude': point.latitude,
                    'longitude': point.longitude,
                    'elevation': point.elevation,
                    'time':point.time
                })
    
    route_df = pd.DataFrame(route_info)
    
    save_df(route_df, new_file_name)
    

def save_df(route_df, new_file_name= 'route_df.csv'):
    route_df.to_csv(new_file_name, index=False)
    print(f'file saved to {new_file_name}')
def read_route_df(file):
    route_df = pd.read_csv(file)
    return route_df

def blh2xyz(phi, lam, h):
    phi_rad = np.deg2rad(phi)
    lam_rad = np.deg2rad(lam)
    a = 6378137
    e2 = 0.00669438002290
    N = a / (np.sqrt(1 - e2 * np.sin(phi_rad) * np.sin(phi_rad)))

    X = (N + h) * np.cos(phi_rad) * np.cos(lam_rad)
    Y = (N + h) * np.cos(phi_rad) * np.sin(lam_rad)
    Z = (N * (1 - e2) + h) * np.sin(phi_rad)
    return np.array([X, Y, Z]) #w radianach




df = read_route_df('poznaj-historyczna-stolice-podlasia.csv')

pio.renderers.default = 'browser'

fig = go.Figure(data=go.Scattermapbox(lat = df['latitude'],
    lon=df['longitude'], mode='lines'))
fig.update_layout(mapbox_style='open-street-map')
fig.show()

phis = df['latitude']
lams = df['longitude']

dists = []
for i, _ in enumerate(phis):
    if i < len(phis)-1:
        phi1 = phis[i]
        lam1 = lams[i]
        phi2 = phis[i+1]
        lam2 = lams[i+1]
        
        xyz1 = blh2xyz(phi1, lam1, 200)
        xyz2 = blh2xyz(phi2, lam2, 200)
        dxyz = np.array(xyz2) - np.array(xyz1)
        
        d = np.linalg.norm(dxyz)
        dists.append(d)
dists.append(np.nan)

df['dists'] = dists
df['distance'] = df['dists'].cumsum()

# dlugosc trasy
suma_dists = df['dists'].sum()
print(f"Suma dla 'dists': {suma_dists}")


model = np.genfromtxt('siatka_gugik-geoid2011-PL-KRON86-NH.txt')

def zmniejszenieSiatki(model, phis, lams, grid_step=0.01):
    max_phi = np.max(phis)
    min_phi = np.min(phis)
    max_lam = np.max(lams)
    min_lam = np.min(lams)
    
    ind_phi = np.logical_and(model[:,0]<(max_phi+grid_step), model[:,0]>(min_phi-grid_step))
    ind_lam = np.logical_and(model[:,1]<(max_lam+grid_step), model[:,1]>(min_lam-grid_step))
    
    indeksy = np.logical_and(ind_phi, ind_lam)
    model2 = model[indeksy, :]
    return model2
model2 = zmniejszenieSiatki(model, phis, lams)
#%%


def interpolacja(model2, phi, lam, grid_step=0.01):
    ind_phi = np.logical_and(model2[:,0]<(phi+grid_step), 
                             model2[:,0]>(phi-grid_step))
    ind_lam = np.logical_and(model2[:,1]<(lam+grid_step), 
                             model2[:,1]>(lam-grid_step))
    
    indeksy = np.logical_and(ind_phi, ind_lam)
    model3 = model2[indeksy, :]

    (y1,x1,Q11),(_y1,x2,Q21),(y2,_x1,Q12),(_y2,_x2,Q22) = model3
    x = lam
    y = phi
    R1 = Q11 + (Q21-Q11)/(x2-x1) * (x-x1)
    R2 = Q12 + (Q22-Q12)/(x2-x1) * (x-x1)
    
    P = R1 + (R2-R1)/(y2-y1) * (y-y1)
    return P


dzety = []
for phi, lam in zip(phis, lams):
    dzeta = interpolacja(model2, phi, lam)
    dzety.append(dzeta)


df['hel'] = df['elevation'] + np.array(dzety)

modelA = np.genfromtxt('Model_quasi-geoidy-PL-geoid2021-PL-EVRF2007-NH.txt')
modelA2 = zmniejszenieSiatki(modelA, phis, lams)

dzetyA = []
for phi, lam in zip(phis, lams):
    dzeta = interpolacja(modelA2, phi, lam)
    dzetyA.append(dzeta)


df['h_norm_Amsterdam'] = df['hel'] - np.array(dzetyA)

df['roznica']= df['h_norm_Amsterdam'] - df['elevation']

df['hel_A']= df['h_norm_Amsterdam']+np.array(dzetyA)

#wykonanie profilu trasy do wysokosci normalnej dla Amsterdama i Kronsztad
plt.figure(figsize=(12, 6))
plt.plot(df['distance'], df['h_norm_Amsterdam'], label='Wysokość normalna(Amsterdam)')
plt.plot(df['distance'], df['elevation'], label='Wysokość normalna (Kronsztad)', linestyle='--')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.title('Profil Trasy')
plt.xlabel('Dystans (km)')
plt.ylabel('Wysokość (m)')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black', alpha=0.2)
plt.show()


plt.figure(figsize=(12, 6))
plt.plot(df['distance'], df['hel'], label='Wysokość elipsoidalna (Kronsztad)')
plt.plot(df['distance'], df['hel_A'], label='Wysokość elipsoidalna (Amsterdam)', linestyle='--')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.title('Profil Trasy')
plt.xlabel('Dystans (km)')
plt.ylabel('Wysokość (m)')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black', alpha=0.2)
plt.show()

plt.figure(figsize=(12, 6))
plt.plot(df['distance'], df['roznica'], label='różnica')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.title('Wykres różnicy między wysokościami normalnymi w układzie Kronsztad i Amsterdam')
plt.xlabel('Różnica (m)')
plt.ylabel('Wysokość (m)')
plt.minorticks_on()
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black', alpha=0.2)
plt.show()

def show_table():
    window = tk.Tk()
    window.title("Różnica między wysokościami normalnymi w układzie Kronsztad i Amsterdam")
    
    tree = ttk.Treeview(window, show="headings")  # Set show="headings" to hide the leftmost column
    tree["columns"] = ("index", "roznica")  # Update columns
    tree.column("index", width=50, anchor="center")  # Width of the index column
    tree.column("roznica", width=100, anchor="center")  # Width of the "roznica" column
    
    tree.heading("index", text="Np.")
    tree.heading("roznica", text="Różnica")

    for i, (index, row) in enumerate(df.iterrows(), start=1):
        rounded_value = round(row['roznica'], 6)
        tree.insert("", tk.END, values=(i, rounded_value))
        
    tree.pack(expand=True, fill=tk.BOTH)  # Expand the table to fill the entire window
    window.geometry("450x200")  # Set the window size, adjust as needed
    window.mainloop()

show_table()
