import numpy as np
from pyproj import Transformer
import math
from pyproj import Geod, CRS
import pyproj
from shapely.geometry import Polygon
from shapely.geometry import Point
from geopy.distance import distance

lon_array = np.array([17.5, 17.5, 18.951338613315183, 18.951338613315183])
lat_array = np.array([51.5, 51.85951417510828, 51.85056129662844, 51.4910465703782])
p1=(17.5, 51.5)
p2=(17.5, 51.859514175)
p3=(18.951338613, 51.850561296)
p4=(18.951338613, 51.491046570)
punkty=[(51.5,17.5 ), (51.859514175, 17.5), (51.850561296, 18.951338613), (51.491046570, 18.951338613)]
dist1 = distance(p1, p2).kilometers
dist2= distance(p2, p3).kilometers
dist3= distance(p3, p4).kilometers
dist4= distance(p1, p4).kilometers

print(f"Distance between points: {dist1, dist2, dist3, dist4} kilometers")
def calculate_distance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

def calculate_punkty_posrednie(lista_punktow):
    punkty_posrednie = []
    
    for i in range(len(lista_punktow)):
        punkt1 = lista_punktow[i]
        punkt2 = lista_punktow[(i + 1) % len(lista_punktow)]
        
        punkt_posredni = ((punkt1[0] + punkt2[0]) / 2, (punkt1[1] + punkt2[1]) / 2)
        
        punkty_posrednie.append(punkt_posredni)
    
    return punkty_posrednie

#PLLAEA
input_proj2 = CRS.from_epsg(4326)
output_proj2 = CRS.from_epsg(3035)
philam_laea_transformer = Transformer.from_proj(input_proj2, output_proj2)
pllaea_1 = philam_laea_transformer.transform(p1[1], p1[0])
pllaea_2 = philam_laea_transformer.transform(p2[1], p2[0])
pllaea_3 = philam_laea_transformer.transform(p3[1], p3[0])
pllaea_4 = philam_laea_transformer.transform(p4[1], p4[0])

pllaea_points = [pllaea_1, pllaea_2, pllaea_3, pllaea_4]
print("pllaea_points",pllaea_points)
#UTM

input_proj3= CRS.from_epsg(4326)
output_proj3 = CRS.from_epsg(32633)
philam_utm_transformer = Transformer.from_proj(input_proj3, output_proj3)
plutm_1 = philam_utm_transformer.transform(p1[1], p1[0])
plutm_2 = philam_utm_transformer.transform(p2[1], p2[0])
plutm_3 = philam_utm_transformer.transform(p3[1], p3[0])
plutm_4 = philam_utm_transformer.transform(p4[1], p4[0])

plutm_points = [plutm_1, plutm_2, plutm_3, plutm_4]

#LLC
input_proj4 = CRS.from_epsg(4326)
output_proj4 = CRS.from_epsg(3034)
philam_llc_transformer = Transformer.from_proj(input_proj4, output_proj4)
plllc_1 = philam_llc_transformer.transform(p1[1], p1[0])
plllc_2 = philam_llc_transformer.transform(p2[1], p2[0])
plllc_3 = philam_llc_transformer.transform(p3[1], p3[0])
plllc_4 = philam_llc_transformer.transform(p4[1], p4[0])
pllc_points = [plllc_1, plllc_2, plllc_3, plllc_4]
print("pllc_points",pllc_points)
#PL-1992
input_proj1 = CRS.from_epsg(4326)
output_proj1 = CRS.from_epsg(2180)
philam_1992_transformer = Transformer.from_proj(input_proj1, output_proj1)
pl1992_1 = philam_1992_transformer.transform(p1[1], p1[0])
pl1992_2 = philam_1992_transformer.transform(p2[1], p2[0])
pl1992_3 = philam_1992_transformer.transform(p3[1], p3[0])
pl1992_4 = philam_1992_transformer.transform(p4[1], p4[0])
pl1992_points = [pl1992_1, pl1992_2, pl1992_3, pl1992_4]

x_gk_array_z_92=[]
y_gk__array_z_92=[]
for i in range(len(pl1992_points)):
    x_gk= (pl1992_points[i][0]+5300000)/0.9993
    y_gk= (pl1992_points[i][1]-500000)/0.9993
    x_gk_array_z_92.append(x_gk)
    y_gk__array_z_92.append(y_gk)

coordinates_pairs_gk_z_92 = list(zip(x_gk_array_z_92, y_gk__array_z_92))

#Odległości PL-1992

distances92 = [calculate_distance(coordinates_pairs_gk_z_92[i], coordinates_pairs_gk_z_92[(i+1) % len(coordinates_pairs_gk_z_92)]) for i in range(len(coordinates_pairs_gk_z_92))]
for i, distance in enumerate(distances92, start=1):
    print(f"Odległość w 92 w gk między punktem {i} a punktem {i % len(coordinates_pairs_gk_z_92) + 1}: {distance/1000:.3f} km")

distancespl1992 = [calculate_distance(pl1992_points[i], pl1992_points[(i+1) % len(pl1992_points)]) for i in range(len(pl1992_points))]
for i, distance in enumerate(distancespl1992, start=1):
    print(f"Odległość w 92 między punktem {i} a punktem {i % len(pl1992_points) + 1}: {distance/1000:.3f} km")
#PL-2000
pl_2000 = Transformer.from_crs('epsg:4326', 'epsg:2178')
input_proj = CRS.from_epsg(4326)
output_proj = CRS.from_epsg(2178)
philam_2000_transformer = Transformer.from_proj(input_proj, output_proj)
pl2000_1 = philam_2000_transformer.transform(p1[1], p1[0])
pl2000_2 = philam_2000_transformer.transform(p2[1], p2[0])
pl2000_3 = philam_2000_transformer.transform(p3[1], p3[0])
pl2000_4 = philam_2000_transformer.transform(p4[1], p4[0])
pl2000_points = [pl2000_1, pl2000_2, pl2000_3, pl2000_4]

x_gk_array_z_00=[]
y_gk__array_z_00=[]
nr= 6 

for i in range(len(pl2000_points)):
    x_gk= pl2000_points[i][0]/0.999923
    y_gk= (pl2000_points[i][1]-nr*1000000-500000)/0.999923
    x_gk_array_z_00.append(x_gk)
    y_gk__array_z_00.append(y_gk)

coordinates_pairs_gk_z_00 = list(zip(x_gk_array_z_00, y_gk__array_z_00))


#Odległości PL-2000
distances00 = [calculate_distance(coordinates_pairs_gk_z_00[i], coordinates_pairs_gk_z_00[(i+1) % len(coordinates_pairs_gk_z_00)]) for i in range(len(coordinates_pairs_gk_z_00))]
for i, distance in enumerate(distances00, start=1):
    print(f"Odległość w 00 gk między punktem {i} a punktem {i % len(coordinates_pairs_gk_z_00) + 1}: {distance/1000:.3f} km")

distancespl2000 = [calculate_distance(pl2000_points[i], pl2000_points[(i+1) % len(pl2000_points)]) for i in range(len(pl2000_points))]
for i, distance in enumerate(distancespl2000, start=1):
    print(f"Odległość w 00 między punktem {i} a punktem {i % len(pl2000_points) + 1}: {distance/1000:.3f} km")

#Obliczenie redukcji

a = 6378137
e2 = 0.00669438002290

def calculate_RM(a, e2, points):
    RM= []
    for i in range(len(points)):
        N= a/(1-e2*(np.sin(np.radians(points[i]))**2))**0.5
        M= a*(1-e2)/((1-e2*(np.sin(np.radians(points[i]))**2))**1.5)
        R= (N*M)**0.5
        RM.append(R) 
    return RM

phi_12_00=(pl2000_points[1][0]+pl2000_points[0][0])/2.0
phi_23_00=(pl2000_points[2][0]+pl2000_points[1][0])/2.0
phi_34_00=(pl2000_points[3][0]+pl2000_points[2][0])/2.0
phi_41_00=(pl2000_points[0][0]+pl2000_points[3][0])/2.0
posrednie_00= [phi_12_00, phi_23_00, phi_34_00, phi_41_00]

RM= calculate_RM(a, e2, posrednie_00)
print("RM", RM)

r12_00= distancespl2000[0]*(coordinates_pairs_gk_z_00[0][1]**2 +coordinates_pairs_gk_z_00[0][1]*coordinates_pairs_gk_z_00[1][1]+ coordinates_pairs_gk_z_00[1][1]**2)
r12_00= r12_00/(6*(RM[0]**2))
r23_00= distancespl2000[1]*(coordinates_pairs_gk_z_00[1][1]**2 +coordinates_pairs_gk_z_00[1][1]*coordinates_pairs_gk_z_00[2][1]+ coordinates_pairs_gk_z_00[2][1]**2)
r23_00= r23_00/(6*(RM[1]**2))
r34_00= distancespl2000[2]*(coordinates_pairs_gk_z_00[2][1]**2 +coordinates_pairs_gk_z_00[2][1]*coordinates_pairs_gk_z_00[3][1]+ coordinates_pairs_gk_z_00[3][1]**2)
r34_00= r34_00/(6*(RM[2]**2))
r41_00= distancespl2000[3]*(coordinates_pairs_gk_z_00[3][1]**2 +coordinates_pairs_gk_z_00[3][1]*coordinates_pairs_gk_z_00[0][1]+ coordinates_pairs_gk_z_00[0][1]**2)
r41_00= r41_00/(6*(RM[3]**2))
print("r12_00", r12_00, "r23_00", r23_00, "r34_00", r34_00, "r41_00", r41_00)

selip12_00= distances00[0]- r12_00
selip23_00= distances00[1]- r23_00
selip34_00= distances00[2]- r34_00
selip41_00= distances00[3]- r41_00

print("selip12_00", selip12_00/1000, "selip23_00", selip23_00/1000, "selip34_00", selip34_00/1000, "selip41_00", selip41_00/1000)

phi_12_92= (pl1992_points[1][0]+pl1992_points[0][0])/2.0
phi_23_92= (pl1992_points[2][0]+pl1992_points[1][0])/2.0
phi_34_92= (pl1992_points[3][0]+pl1992_points[2][0])/2.0
phi_41_92= (pl1992_points[0][0]+pl1992_points[3][0])/2.0
posrednie_92= [phi_12_92, phi_23_92, phi_34_92, phi_41_92]

RM1= calculate_RM(a, e2, posrednie_92)
print("RM", RM1)

r12_92= distancespl1992[0]*(coordinates_pairs_gk_z_92[0][1]**2 +coordinates_pairs_gk_z_92[0][1]*coordinates_pairs_gk_z_92[1][1]+ coordinates_pairs_gk_z_92[1][1]**2)
r12_92= r12_92/(6*(RM1[0]**2))
r23_92= distancespl1992[1]*(coordinates_pairs_gk_z_92[1][1]**2 +coordinates_pairs_gk_z_92[1][1]*coordinates_pairs_gk_z_92[2][1]+ coordinates_pairs_gk_z_92[2][1]**2)
r23_92= r23_92/(6*(RM1[1]**2))
r34_92= distancespl1992[2]*(coordinates_pairs_gk_z_92[2][1]**2 +coordinates_pairs_gk_z_92[2][1]*coordinates_pairs_gk_z_92[3][1]+ coordinates_pairs_gk_z_92[3][1]**2)
r34_92= r34_92/(6*(RM1[2]**2))
r41_92= distancespl1992[3]*(coordinates_pairs_gk_z_92[3][1]**2 +coordinates_pairs_gk_z_92[3][1]*coordinates_pairs_gk_z_92[0][1]+ coordinates_pairs_gk_z_92[0][1]**2)
r41_92= r41_92/(6*(RM1[3]**2))
print("r12_92", r12_92, "r23_92", r23_92, "r34_92", r34_92, "r41_92", r41_92)
print("distances92", distances92)
selip12_92= distances92[0]- r12_92
selip23_92= distances92[1]- r23_92
selip34_92= distances92[2]- r34_92
selip41_92= distances92[3]- r41_92

print("selip12_92", selip12_92/1000, "selip23_92", selip23_92/1000, "selip34_92", selip34_92/1000, "selip41_92", selip41_92/1000)

gamma = []
for i in range(len(punkty)):
    punkty = np.radians(punkty)
    lambda_0 = np.radians(18)
    b2 = a**2 * (1 - e2)
    e2_prime = (a**2 - b2)/ b2
    delta_lambda = punkty[i][1] - lambda_0
    t = math.tan(punkty[i][0])
    eta2 = e2_prime * math.cos(punkty[i][0])**2
    eta = math.sqrt(eta2)
    N = a / math.sqrt(1 - e2 * math.sin(punkty[i][0])**2)
    gamma_i = delta_lambda * math.sin(punkty[i][0]) + (delta_lambda**3 / 3) * math.sin(punkty[i][0]) * math.cos(punkty[i][0])**2 * (1 + 3 * eta**2 + 2 * eta**4) + (delta_lambda**5 / 15) * math.sin(punkty[i][0]) * math.cos(punkty[i][0])**4 * (2 - t**2)
    gamma.append(gamma_i)

print("gamma", gamma)

#redukcja azymutów
alfa_kier_12_00= np.arctan2(coordinates_pairs_gk_z_00[1][1]- coordinates_pairs_gk_z_00[0][1],coordinates_pairs_gk_z_00[1][0]- coordinates_pairs_gk_z_00[0][0])
alfa_kier_23_00= np.arctan2(coordinates_pairs_gk_z_00[2][1]- coordinates_pairs_gk_z_00[1][1],coordinates_pairs_gk_z_00[2][0]- coordinates_pairs_gk_z_00[1][0])
alfa_kier_34_00= np.arctan2(coordinates_pairs_gk_z_00[3][1]- coordinates_pairs_gk_z_00[2][1],coordinates_pairs_gk_z_00[3][0]- coordinates_pairs_gk_z_00[2][0])
alfa_kier_41_00= np.arctan2(coordinates_pairs_gk_z_00[0][1]- coordinates_pairs_gk_z_00[3][1],coordinates_pairs_gk_z_00[0][0]- coordinates_pairs_gk_z_00[3][0])
print("alfa_kier_ab ",np.rad2deg(alfa_kier_12_00))
print("alfa_kier_bc ",np.rad2deg(alfa_kier_23_00))
print("alfa_kier_cd ",np.rad2deg(2*np.pi+alfa_kier_34_00))
print("alfa_kier_da ",np.rad2deg(2*np.pi+alfa_kier_41_00))

alfa_kier_12_92= np.arctan2(coordinates_pairs_gk_z_92[1][1]- coordinates_pairs_gk_z_92[0][1],coordinates_pairs_gk_z_92[1][0]- coordinates_pairs_gk_z_92[0][0])
alfa_kier_23_92= np.arctan2(coordinates_pairs_gk_z_92[2][1]- coordinates_pairs_gk_z_92[1][1],coordinates_pairs_gk_z_92[2][0]- coordinates_pairs_gk_z_92[1][0])
alfa_kier_34_92= np.arctan2(coordinates_pairs_gk_z_92[3][1]- coordinates_pairs_gk_z_92[2][1],coordinates_pairs_gk_z_92[3][0]- coordinates_pairs_gk_z_92[2][0])
alfa_kier_41_92= np.arctan2(coordinates_pairs_gk_z_92[0][1]- coordinates_pairs_gk_z_92[3][1],coordinates_pairs_gk_z_92[0][0]- coordinates_pairs_gk_z_92[3][0])
print("alfa_kier_ab ",np.rad2deg(alfa_kier_12_92))
print("alfa_kier_bc ",np.rad2deg(alfa_kier_23_92))
print("alfa_kier_cd ",np.rad2deg(2*np.pi+alfa_kier_34_92))
print("alfa_kier_da ",np.rad2deg(2*np.pi+alfa_kier_41_92))

#obliczenie redukcji kierunków

teta12_00= ((coordinates_pairs_gk_z_00[1][0]- coordinates_pairs_gk_z_00[0][0])*(2*coordinates_pairs_gk_z_00[0][1]+ coordinates_pairs_gk_z_00[1][1]))/(6*RM[0]**2)
teta23_00= ((coordinates_pairs_gk_z_00[2][0]- coordinates_pairs_gk_z_00[1][0])*(2*coordinates_pairs_gk_z_00[1][1]+ coordinates_pairs_gk_z_00[2][1]))/(6*RM[1]**2)
teta34_00= ((coordinates_pairs_gk_z_00[3][0]- coordinates_pairs_gk_z_00[2][0])*(2*coordinates_pairs_gk_z_00[2][1]+ coordinates_pairs_gk_z_00[3][1]))/(6*RM[2]**2)
teta41_00= ((coordinates_pairs_gk_z_00[0][0]- coordinates_pairs_gk_z_00[3][0])*(2*coordinates_pairs_gk_z_00[3][1]+ coordinates_pairs_gk_z_00[0][1]))/(6*RM[3]**2)
print("teta12 ",np.rad2deg(teta12_00), "teta23 ",np.rad2deg(teta23_00), "teta34 ",np.rad2deg(teta34_00), "teta41 ",np.rad2deg(teta41_00))
teta12_92= ((coordinates_pairs_gk_z_92[1][0]- coordinates_pairs_gk_z_92[0][0])*(2*coordinates_pairs_gk_z_92[0][1]+ coordinates_pairs_gk_z_92[1][1]))/(6*RM[0]**2)
teta23_92= ((coordinates_pairs_gk_z_92[2][0]- coordinates_pairs_gk_z_92[1][0])*(2*coordinates_pairs_gk_z_92[1][1]+ coordinates_pairs_gk_z_92[2][1]))/(6*RM[1]**2)
teta34_92= ((coordinates_pairs_gk_z_92[3][0]- coordinates_pairs_gk_z_92[2][0])*(2*coordinates_pairs_gk_z_92[2][1]+ coordinates_pairs_gk_z_92[3][1]))/(6*RM[2]**2)
teta41_92= ((coordinates_pairs_gk_z_92[0][0]- coordinates_pairs_gk_z_92[3][0])*(2*coordinates_pairs_gk_z_92[3][1]+ coordinates_pairs_gk_z_92[0][1]))/(6*RM[3]**2)
print("teta12 ",np.rad2deg(teta12_92), "teta23 ",np.rad2deg(teta23_92), "teta34 ",np.rad2deg(teta34_92), "teta41 ",np.rad2deg(teta41_92))
A12_00= alfa_kier_12_00+teta12_00+gamma[0]
A23_00= alfa_kier_23_00+teta23_00+gamma[1]
A34_00= alfa_kier_34_00+teta34_00+2*np.pi+gamma[2]
A41_00= alfa_kier_41_00+teta41_00+2*np.pi+gamma[3]

A12_92= alfa_kier_12_92+teta12_92+gamma[0]
A23_92= alfa_kier_23_92+teta23_92+gamma[1]
A34_92= alfa_kier_34_92+teta34_92+2*np.pi+gamma[2]
A41_92= alfa_kier_41_92+teta41_92+2*np.pi+gamma[3]

print("A12_00 ",np.rad2deg(A12_00), "A23_00 ",np.rad2deg(A23_00), "A34_00 ",np.rad2deg(A34_00), "A41_00 ",np.rad2deg(A41_00))
print("A12_92 ",np.rad2deg(A12_92), "A23_92 ",np.rad2deg(A23_92), "A34_92 ",np.rad2deg(A34_92), "A41_92 ",np.rad2deg(A41_92))
#stopnie minuty sekundy
def radian_to_degrees_minutes_seconds(angle_in_radians):
    degrees = np.degrees(angle_in_radians)
    degrees_int = int(degrees)
    minutes = (degrees - degrees_int) * 60
    minutes_int = int(minutes)
    seconds = (minutes - minutes_int) * 60

    return degrees_int, minutes_int, seconds

# Konwersja kątów na stopnie, minuty i sekundy
A12_00_deg = radian_to_degrees_minutes_seconds(A12_00)
A23_00_deg = radian_to_degrees_minutes_seconds(A23_00)
A34_00_deg = radian_to_degrees_minutes_seconds(A34_00)
A41_00_deg = radian_to_degrees_minutes_seconds(A41_00)

A12_92_deg = radian_to_degrees_minutes_seconds(A12_92)
A23_92_deg = radian_to_degrees_minutes_seconds(A23_92)
A34_92_deg = radian_to_degrees_minutes_seconds(A34_92)
A41_92_deg = radian_to_degrees_minutes_seconds(A41_92)

# Wyświetlanie wyników
print("A12_00:", A12_00_deg)
print("A23_00:", A23_00_deg)
print("A34_00:", A34_00_deg)
print("A41_00:", A41_00_deg)

print("A12_92:", A12_92_deg)
print("A23_92:", A23_92_deg)
print("A34_92:", A34_92_deg)
print("A41_92:", A41_92_deg)

def oblicz_pole_powierzchni(dane):
    n = len(dane)
    
    if n < 3:
        raise ValueError("Za mało punktów do obliczenia pola powierzchni.")

    pole = 0.0

    for i in range(n):
        xi, yi = dane[i]
        xi_plus_1, yi_plus_1 = dane[(i + 1) % n]
        xi_minus_1, yi_minus_1 = dane[(i - 1) % n]

        pole += xi * (yi_plus_1 - yi_minus_1)

    pole /= 2.0

    return abs(pole)

wynik = oblicz_pole_powierzchni(coordinates_pairs_gk_z_00)
print("wynik 00 W KM2",wynik/1000000)
wynik1= oblicz_pole_powierzchni(coordinates_pairs_gk_z_92)
print("wynik 92 W KM2",wynik1/1000000)
wynik2= oblicz_pole_powierzchni(pllaea_points)
print("wynik pllaea W KM2",wynik2/1000000)






 
