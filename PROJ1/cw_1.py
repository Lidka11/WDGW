import numpy as np
import matplotlib.pyplot as plt
from math import radians, degrees

def hms2rad(h, m, s):
    return (h + m / 60 + s / 3600) * (np.pi / 12)

def dms2rad(d, m, s):
    return (d + m / 60 + s / 3600) * (np.pi / 180)

def julday(y, m, d, h):
    '''
    Simplified Julian Date generator, valid only between
    1 March 1900 to 28 February 2100
    '''
    if m <= 2:
        y = y - 1
        m = m + 12

    jd = np.floor(365.25 * (y + 4716)) + np.floor(30.6001 * (m + 1)) + d + h / 24 - 1537.5
    return jd

def GMST(jd):
    '''
    Calculation of Greenwich Mean Sidereal Time - GMST in hours
    ----------
    jd : julian date
    '''
    T = (jd - 2451545) / 36525
    GMST = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T**2 - T**3 / 38710000
    GMST = (GMST % 360) / 15  # GMST w godzinach
    return GMST

def plot_skyplot_one_star(RA, Dec, hours):
    latitude=52.0  # szerokość geograficzna Warszawy
    longitude=21.0 # długość geograficzna Warszawy
    jd=julday(2023, 7, 1, hours) # data i godzina
    time_offset=2
    jd_utc=jd-time_offset/24
    GMST0 = GMST(jd_utc)
    LST = GMST0 * 15 + longitude
    
    t = np.deg2rad(LST) - hms2rad(*RA)
    h = np.arcsin(np.sin(np.deg2rad(latitude)) * np.sin(dms2rad(*Dec)) +
                  np.cos(np.deg2rad(latitude)) * np.cos(dms2rad(*Dec)) * np.cos(t))
    Az = np.arctan2(-np.cos(dms2rad(*Dec)) * np.sin(t),
                     np.cos(np.deg2rad(latitude)) * np.sin(dms2rad(*Dec)) -
                     np.sin(np.deg2rad(latitude)) * np.cos(dms2rad(*Dec)) * np.cos(t))
    
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_yticks(range(0, 90+10, 10))
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    ax.set_rlim(0, 90)
    ax.scatter(Az, 90 - np.rad2deg(h))
    ax.set_title("Wykres skyplot położenia gwiazdy dla Warszawy")
    plt.show()

def plot_skyplot_stars(RA_list, Dec_list, hours):
    latitude=52.0  # szerokość geograficzna Warszawy
    longitude=21.0 # długość geograficzna Warszawy
    jd=julday(2023, 7, 1, hours) # data i godzina
    time_offset=2
    jd_utc=jd-time_offset/24
    GMST0 = GMST(jd_utc)
    LST = GMST0 * 15 + longitude
    fig=plt.figure(figsize=(8,8))
    ax = fig.add_subplot(polar=True)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_yticks(range(0, 90+10, 10))
    yLabel = ['90', '', '', '60', '', '', '30', '', '', '']
    ax.set_yticklabels(yLabel)
    ax.set_rlim(0, 90)
    for i in range(len(RA_list)):
        RA, Dec = RA_list[i], Dec_list[i]
        t = np.deg2rad(LST) - hms2rad(*RA)
        h = np.arcsin(np.sin(np.deg2rad(latitude)) * np.sin(dms2rad(*Dec)) +
                      np.cos(np.deg2rad(latitude)) * np.cos(dms2rad(*Dec)) * np.cos(t))
        Az = np.arctan2(-np.cos(dms2rad(*Dec)) * np.sin(t),
                         np.cos(np.deg2rad(latitude)) * np.sin(dms2rad(*Dec)) -
                         np.sin(np.deg2rad(latitude)) * np.cos(dms2rad(*Dec)) * np.cos(t))
        ax.scatter(Az, 90 - np.rad2deg(h), label=f'Gwiazda {i+1}')
    ax.set_title("Wykres skyplot położenia gwiazd dla Warszawy")
    ax.legend(loc='upper right')
    plt.show()



def plot_3d_stars(RA_list, Dec_list, hours):
    latitude=52.0
    longitude=21.0
    r= 1
    # Siatka współrzędnych
    u, v = np.mgrid[0:(2 * np.pi + 0.1):0.1, 0:np.pi:0.1]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)

    z[z < 0] = 0  # Rysujemy tylko półkulę
    fig=plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(x, y, z, alpha=0.1)
    jd=julday(2023, 7, 1, hours)
    jd-=2/24
    G0=GMST(jd)
    LST=G0*15+longitude
    for RA, Dec in zip(RA_list, Dec_list):
        t = np.deg2rad(LST) - hms2rad(*RA)
        h = np.arcsin(np.sin(np.deg2rad(latitude)) * np.sin(dms2rad(*Dec)) +
                      np.cos(np.deg2rad(latitude)) * np.cos(dms2rad(*Dec)) * np.cos(t))
        Az = np.arctan2(-np.cos(dms2rad(*Dec)) * np.sin(t),
                         np.cos(np.deg2rad(latitude)) * np.sin(dms2rad(*Dec)) -
                         np.sin(np.deg2rad(latitude)) * np.cos(dms2rad(*Dec)) * np.cos(t))
        gx = r * np.sin(Az) * np.cos(h)
        gy = r * np.cos(Az) * np.cos(h)
        gz = r * np.sin(h)
        ax.plot3D(gx, gy, gz, 'o')
    ax.set_title("Położenie gwiazd na niebie dla Warszawy")
    plt.show()

def plot_3d_one_star(RA, Dec, hours):
    latitude=52.0
    longitude=21.0  
    r= 1
    # Siatka współrzędnych
    u, v = np.mgrid[0:(2 * np.pi + 0.1):0.1, 0:np.pi:0.1]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)

    z[z < 0] = 0  # Rysujemy tylko półkulę
    fig=plt.figure(figsize=(10,10))    
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(x, y, z, alpha=0.1)
    jd=julday(2023, 7, 1, hours)
    jd-=2/24
    G0=GMST(jd)
    LST=G0*15+longitude
    t = np.deg2rad(LST) - hms2rad(*RA)
    h = np.arcsin(np.sin(np.deg2rad(latitude)) * np.sin(dms2rad(*Dec)) +
                  np.cos(np.deg2rad(latitude)) * np.cos(dms2rad(*Dec)) * np.cos(t))    
    Az = np.arctan2(-np.cos(dms2rad(*Dec)) * np.sin(t),
                     np.cos(np.deg2rad(latitude)) * np.sin(dms2rad(*Dec)) - 
                     np.sin(np.deg2rad(latitude)) * np.cos(dms2rad(*Dec)) * np.cos(t))    
    gx = r * np.sin(Az) * np.cos(h)    
    gy = r * np.cos(Az) * np.cos(h)    
    gz = r * np.sin(h)                      
    ax.plot3D(gx, gy, gz, 'o')
    ax.set_title("Położenie gwiazd na niebie dla Warszawy") 
    plt.show()

if __name__ == '__main__':
    hours=np.arange(1,25)
    RA_list = [(13, 24,52.075), (11,3,14.669), (11,5,9.530),(12,55,3.395), (13,48,27.861), (11,55,3.388),(12,16,34.755)]
    Dec_list = [(54,48,11.50), (56,15,21.18), (61,37,24.44), (55,49,57.62), (49,11,48.04), (53,33,50.55), (56,54,7.80)]

    RA=[13, 26, 26.067]
    Dec=[-11, 16, 59.73]

    plot_skyplot_one_star(RA, Dec, hours)
    plot_3d_one_star(RA, Dec, hours)
    plot_skyplot_stars(RA_list, Dec_list, hours)
    plot_3d_stars(RA_list, Dec_list, hours)
    
