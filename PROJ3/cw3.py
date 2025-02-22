#numer:10
"""
Dane:
phi = 51,30
lam = 17,30

Azymuty i linie 3 linii geodezyjnych
    Długości   Azymuty
1-2: 40 km      0
2-3: 100 km     90
3-4: 40 km      180
4-1: 100 km     270
"""

#algorytm kivioji
import numpy as np
from pyproj import Geod
import folium
import webbrowser
from folium import plugins
from pyproj import Geod
from shapely.geometry import Polygon
from shapely.geometry import Point
a = 6378137
e2 = 0.00669438002290

Az_array = np.deg2rad([0, 90, 180, 270])
s_array = [40, 100, 40, 100]

phi_1= np.deg2rad(51.30)
lam_1 = np.deg2rad(17.30)
Az_1= np.deg2rad(0.0)
print(phi_1, lam_1, Az_1)
s=40000

def Np(B, a=a, e2=e2):
    N=a/(1-e2*(np.sin(B)**2))**0.5
    return N
def Mp(B, a=a, e2=e2):
    M= a*(1-e2)/((1-e2*np.sin(B)**2)**3)**(0.5)
    return M


#podzial linii
def Kivioj(phi_1, lam_1, Az_1, s):
    n = round(s/1000)
    ds = s/n

    for i in range(n):
        # 2. Obliczamy główne promienie krzywizny M i N w punkcie wyjściowym P1 oraz stałą c linii geodezyjnej
        M_i = Mp(phi_1)
        N_i = Np(phi_1)

        # 3. Pierwsze przybliżenie przyrostu szerokości i azymutu
        dphi_i = ds * np.cos(Az_1) / M_i
        dAz_i = (ds * np.sin(Az_1) * np.tan(phi_1)) / N_i

        # 4. Obliczenie szerokości i azymutu w punkcie środkowym (m) odcinka, na podstawie przyrostów
        phi_m = phi_1 + dphi_i/2
        Az_m = Az_1 + dAz_i/2

        # 5. Obliczenie promieni krzywizn w kierunkach głównych w punkcie m
        M_m = Mp(phi_m)
        N_m = Np(phi_m)

        # 6. Ostateczne przyrosty szerokości, długości i azymutu
        dphi_m = (ds * np.cos(Az_m)) / M_m
        dAz_m = (ds * np.sin(Az_m) * np.tan(phi_m)) / N_m
        dlam_m = (ds * np.sin(Az_m)) / (N_m * np.cos(phi_m))

        phi_1 = phi_1 + dphi_m
        #Az_1 = Az_1 + dAz_m
        lam_1 = lam_1 + dlam_m
    Az_1 = Az_1 + dAz_m
    # Return the final values after iterating through all segments
    return phi_1, lam_1, Az_1

def vincenty(BA,LA,BB,LB):
    '''
    Parameters
    ----------
    BA : szerokosc geodezyjna punktu A [RADIAN]
    LA : dlugosc geodezyjna punktu A [RADIAN]
    BB : szerokosc geodezyjna punktu B [RADIAN]
    LB : dlugosc geodezyjna punktu B [RADIAN]

    Returns
    -------
    sAB : dlugosc linii geodezyjnej AB [METR]
    A_AB : azymut linii geodezyjnej AB [RADIAN]
    A_BA : azymut odwrotny linii geodezyjne [RADIAN]
    '''
    b = a * np.sqrt(1-e2)
    f = 1-b/a
    dL = LB - LA
    UA = np.arctan((1-f)*np.tan(BA))
    UB = np.arctan((1-f)*np.tan(BB))
    L = dL
    while True:
        sin_sig = np.sqrt((np.cos(UB)*np.sin(L))**2 +\
                          (np.cos(UA)*np.sin(UB) - np.sin(UA)*np.cos(UB)*np.cos(L))**2)
        cos_sig = np.sin(UA)*np.sin(UB) + np.cos(UA) * np.cos(UB) * np.cos(L)
        sig = np.arctan2(sin_sig,cos_sig)
        sin_al = (np.cos(UA)*np.cos(UB)*np.sin(L))/sin_sig
        cos2_al = 1 - sin_al**2
        cos2_sigm = cos_sig - (2 * np.sin(UA) * np.sin(UB))/cos2_al
        C = (f/16) * cos2_al * (4 + f*(4 - 3 * cos2_al))
        Lst = L
        L = dL + (1-C)*f*sin_al*(sig+C*sin_sig*(cos2_sigm+C*cos_sig*(-1 + 2*cos2_sigm**2)))
        if abs(L-Lst)<(0.000001/206265):
            break
    
    u2 = (a**2 - b**2)/(b**2) * cos2_al
    A = 1 + (u2/16384) * (4096 + u2*(-768 + u2 * (320 - 175 * u2)))
    B = u2/1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    d_sig = B*sin_sig * (cos2_sigm + 1/4*B*(cos_sig*(-1+2*cos2_sigm**2)\
            - 1/6 *B*cos2_sigm * (-3 + 4*sin_sig**2)*(-3+4*cos2_sigm**2)))
    sAB = b*A*(sig-d_sig)
    A_AB = np.arctan2((np.cos(UB) * np.sin(L)),(np.cos(UA)*np.sin(UB) - np.sin(UA)*np.cos(UB)*np.cos(L)))
    A_BA = np.arctan2((np.cos(UA) * np.sin(L)),(-np.sin(UA)*np.cos(UB) + np.cos(UA)*np.sin(UB)*np.cos(L))) + np.pi
    return sAB, A_AB, A_BA


g= Geod(ellps="WGS84")
lam_2P, phi_2P, Az_odwP = g.fwd(lam_1, phi_1, Az_1, s)

pi_pol=np.pi/2
trzy_2_pi=np.pi*3/2
phi_2K, lam_2K, Az_odwK = Kivioj(phi_1, lam_1, Az_1, s)
phi_3K, lam_3K, Az_3K = Kivioj(phi_2K, lam_2K, pi_pol, 100000)
phi_4K, lam_4K, Az_4K = Kivioj(phi_3K, lam_3K, np.pi, 40000)
phi_5K, lam_5K, Az_5K = Kivioj(phi_4K, lam_4K, trzy_2_pi, 100000)

print("Results from Kivioj:")
print("phi_2K:", np.rad2deg(phi_2K))
print("lam_2K:", np.rad2deg(lam_2K))
print("Az_odwK:", np.rad2deg(Az_odwK))

print("phi_3K:", np.rad2deg(phi_3K))
print("lam_3K:", np.rad2deg(lam_3K))
print("Az_3K:", np.rad2deg(Az_3K))

print("phi_4K:", np.rad2deg(phi_4K))
print("lam_4K:", np.rad2deg(lam_4K))
print("Az_4K:", np.rad2deg(Az_4K))

print("phi_5K:", np.rad2deg(phi_5K))
print("lam_5K:", np.rad2deg(lam_5K))
print("Az_5K:", np.rad2deg(Az_5K))


print("\nResults from Geod:")
print("phi_2P:", phi_2P)
print("lam_2P:", lam_2P)
print("Az_odwP:", np.rad2deg(Az_odwP))


# Utwórz obiekt mapy
mymap = folium.Map(location=[51.30, 17.30], zoom_start=12)  # Tutaj podane są współrzędne i poziom przybliżenia

# Dodaj znacznik na mapie
folium.Marker(location=[np.rad2deg(phi_1), np.rad2deg(lam_1)], popup='Point 1').add_to(mymap)
folium.Marker(location=[np.rad2deg(phi_2K), np.rad2deg(lam_2K)], popup='Point 2K').add_to(mymap)
folium.Marker(location=[np.rad2deg(phi_3K), np.rad2deg(lam_3K)], popup='Point 3K').add_to(mymap)
folium.Marker(location=[np.rad2deg(phi_4K), np.rad2deg(lam_4K)], popup='Point 4K').add_to(mymap)
folium.Marker(location=[np.rad2deg(phi_5K), np.rad2deg(lam_5K)], popup='Point 5K').add_to(mymap)

line = folium.PolyLine(
    locations=[
        [np.rad2deg(phi_1), np.rad2deg(lam_1)],
        [np.rad2deg(phi_2K), np.rad2deg(lam_2K)],
        [np.rad2deg(phi_3K), np.rad2deg(lam_3K)],
        [np.rad2deg(phi_4K), np.rad2deg(lam_4K)],
        [np.rad2deg(phi_5K), np.rad2deg(lam_5K)],
    ],
    color='blue',
).add_to(mymap)


mymap.save('mapa.html')

mymap_path = 'mapa.html'

# Otwórz mapę w domyślnej przeglądarce
webbrowser.open('file://' + mymap_path, new=2)

#file://C:\sem3\geodezja wyzsza\cw3\mapa.html

s12,A12, A21 = vincenty(phi_1, lam_1, phi_2K, lam_2K)
print("s12:", s12)
print("A12:", np.rad2deg(A12))
print("A21:", np.rad2deg(A21))

s23,A23, A32 = vincenty(phi_2K, lam_2K, phi_3K, lam_3K)
print("s23:", s23)
print("A23:", np.rad2deg(A23))
print("A32:", np.rad2deg(A32))

s34,A34, A43 = vincenty(phi_3K, lam_3K, phi_4K, lam_4K)
print("s34:", s34)
print("A34:", np.rad2deg(A34))
print("A43:", np.rad2deg(A43))

s45,A45, A54 = vincenty(phi_4K, lam_4K, phi_5K, lam_5K)
print("s45:", s45)
print("A45:", np.rad2deg(A45))
print("A54:", np.rad2deg(A54))

s51,A51, A15 = vincenty(phi_5K, lam_5K, phi_1, lam_1)
print("s51:", s51)
print("A51:", np.rad2deg(A51))
print("A15:", np.rad2deg(A15))
#geod = Geod(ellps="WGS84")

points = [
    (np.rad2deg(phi_5K), np.rad2deg(lam_5K)),
    (np.rad2deg(phi_4K), np.rad2deg(lam_4K)),
    (np.rad2deg(phi_3K), np.rad2deg(lam_3K)),
    (np.rad2deg(phi_2K), np.rad2deg(lam_2K)),
    (np.rad2deg(phi_1), np.rad2deg(lam_1))
]
print('Punkty: ',points)
point_objects = [Point(lon, lat) for lat, lon in points]
print('Punkty1: ',point_objects)
polygon = Polygon(point_objects)
geod = Geod(ellps="WGS84")
area, perimeter = geod.geometry_area_perimeter(polygon)
area = area/1000000
print(f"Obszar: {area:.2f} km^2")
print(f"Obwód: {perimeter:.2f} m")
r = geod.inv(np.rad2deg(lam_1), np.rad2deg(phi_1), np.rad2deg(lam_5K), np.rad2deg(phi_5K))

# Print the result
print('Różnica położenia punktów 1 i 1* wynosi: {:.3f} m'.format(r[2]))

