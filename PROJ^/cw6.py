import numpy as np 
import matplotlib.pyplot as plt
import geopandas as gpd

e2 = 0.00669438002290
a = 6378137
def Np(B, a=a, e2=e2):
    N = a/(1-e2*(np.sin(B)**2))**0.5
    return N

def hirv(x,y,z):
    p = np.sqrt(x**2+y**2)

    phi = np.arctan(z/(p*(1-e2)))
    
    while True:
        N = Np(phi)
        h = p/np.cos(phi) - N
        phi_poprzednie = phi
        phi = np.arctan((z/p) * (1-(N*e2)/(N+h))**(-1))
        if abs(phi_poprzednie-phi)<(0.000001/60/60/360):
            break
    lam = np.arctan2(y,x)
    return phi, lam, h

dane = np.genfromtxt('C:\sem3\geodezja wyzsza\cw6\AJAC.txyz2.txt')


t = dane[:, 2]
xyz = dane[:, 3:6]
dxyz = xyz - xyz[0, :]
fig, ax = plt.subplots(3, 1)
ax[0].plot(dxyz[:,0])
ax[1].plot(dxyz[:,1])
ax[2].plot(dxyz[:,2])

phi, lam, h = hirv(xyz[0, 0], xyz[0, 1], xyz[0, 2])
Rneu = np.array([[-np.sin(phi)*np.cos(lam), -np.sin(lam), np.cos(phi) * np.cos(lam)], #uzupelnic tak jak w instrukcji 
                        [-np.sin(phi) * np.sin(lam), np.cos(lam), np.cos(phi) * np.sin(lam)],
                        [np.cos(phi), 0, np.sin(phi)]])

dneu = []
for dx in dxyz:
    dneu.append(Rneu.T@dx)

dneu= np.array(dneu)

# dla wspolrzednej x
x = dxyz[:,0]
jedynki = np.ones(len(t))
A = np.column_stack((t, jedynki))
L = x
xx = np.linalg.inv(A.T@A)@(A.T@L)
model = A@xx
v = model - L


fig1, ax = plt.subplots()
ax.plot(t, dxyz[:,0])
ax.plot(t, model)
ax.plot(t, v)

#dla wspolrzednej y
y = dxyz[:,1]
Ly= y
yy= np.linalg.inv(A.T@A)@(A.T@Ly)
modely= A@yy
vy= modely-Ly

#dla wspolrzednej z
z = dxyz[:,2]
Lz= z
zz= np.linalg.inv(A.T@A)@(A.T@Lz)
modelz= A@zz
vz= modelz-Lz
"""
fig2, ax= plt.subplots()
ax.plot(t, dxyz[:,1])
ax.plot(t, modely)
ax.plot(t, vy)

fig3, ax= plt.subplots()
ax.plot(t, dxyz[:,2])
ax.plot(t, modelz)
ax.plot(t, vz)


x = dxyz[:, 0]

plt.show()
"""


y_linia_lista=[]
for i in range(3):
    y = dneu[:,i]
    xx= np.linalg.inv(A.T@A)@(A.T@y)
    y_linia = A@xx
    y_linia_lista.append(y_linia)
y_linia_lista= np.array(y_linia_lista).T


y_model= y_linia_lista

#--------------------------
x_linia_lista=[]
for i in range(3):
    x = dneu[:,i]
    xx= np.linalg.inv(A.T@A)@(A.T@x)
    x_linia = A@xx
    x_linia_lista.append(x_linia)
x_linia_lista= np.array(x_linia_lista).T
x_model= x_linia_lista

z_linia_lista=[]
for i in range(3):
    z = dneu[:,i]
    xx= np.linalg.inv(A.T@A)@(A.T@z)
    z_linia = A@xx
    z_linia_lista.append(z_linia)
z_linia_lista= np.array(z_linia_lista).T
z_model= z_linia_lista
#---------------------------------------------

#y_model= np.array(y_linia_lista)
#strzalka
fig, ax = plt.subplots(3,1)
for i in range(len(ax)):
    ax[i].plot(t, dneu[:,i])
    ax[i].plot(t, y_linia_lista[:,i])


filee = r"C:\sem3\geodezja wyzsza\cw6\CNTR_BN_03M_2020_4326.shp"
shpfile = gpd.read_file(filee)
fig, ax = plt.subplots()
shpfile.plot(ax=ax, zorder=0)
vn = (y_model[-1, 0] - y_model[0, 0]) / (t[-1] - t[0]) * 100
ve = (y_model[-1, 1] - y_model[0, 1]) / (t[-1] - t[0]) * 100
vnx = (x_model[-1, 0] - x_model[0, 0]) / (t[-1] - t[0]) * 100
vex = (x_model[-1, 1] - x_model[0, 1]) / (t[-1] - t[0]) * 100
vnz = (z_model[-1, 0] - z_model[0, 0]) / (t[-1] - t[0]) * 100
vez = (z_model[-1, 1] - z_model[0, 1]) / (t[-1] - t[0]) * 100
Q = ax.quiver(np.rad2deg(lam), np.rad2deg(phi), vn, ve, color='red', zorder=1)
Q1 = ax.quiver(np.rad2deg(lam), np.rad2deg(phi), vnx, vex, color='green', zorder=1)
Q2 = ax.quiver(np.rad2deg(lam), np.rad2deg(phi), vnz, vez, color='blue', zorder=1)
ax.quiverkey(Q, 0.9, 0.9, 1, r'$1\frac{cm}{year}$', labelpos='E', coordinates='figure', zorder=2)
ax.quiverkey(Q1, 0.9, 0.8, 1, r'$1\frac{cm}{year}$', labelpos='E', coordinates='figure', zorder=2)
ax.quiverkey(Q2, 0.9, 0.7, 1, r'$1\frac{cm}{year}$', labelpos='E', coordinates='figure', zorder=2)
plt.show()
