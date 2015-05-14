'''
Created on Nov 25, 2014
Computes the single frequency GPS pseudorange correction for the ionospheric 
delay using the Klobuchar model. The main method is Klo_Iono, the remaining 
sub-methods are used in support of it. 
@author: Sonal R.
'''
import math
import numpy as np

'''
Computes the ionospheric error for a single frequency pseudorange (Ex. C1)
using the Klobuchar model.
latlong - Array that contains [Latitude,Longitude,Height] of the station
Azel - Azimuth and elevation Array to the satellite of interest
tsecs - GPS seconds of the current Epoch
alpha_i - Ionospheric parameters A0-A3 of the almanac transmitted in the Navigation data
beta_i - Ionospheric parameters B0-B3 of the almanac transmitted in the Navigation data
returns - Iono_delay - the Pseudorange correction for Ionospheric Delay
'''   
def Klo_Iono(latlong,Azel,tsecs,alpha_i,beta_i):
    Az = Azel[0]*(math.pi/180)
    Ele = Azel[1]*(math.pi/180)
    lat = latlong[0]*(math.pi/180)
    longi = latlong[1]*(math.pi/180)
    h = latlong[2]
    E_semi = Ele/math.pi
    #Earth Centered angle (Elevation E in semicircles)    
    ECA = (0.0137/(E_semi+0.11))-0.022
    #Latitude of the Ionospheric Pierce point
    phi_I = lat+(ECA*math.cos(Az))
    if phi_I > 0.416:
        phi_I = 0.416
    elif phi_I < -0.416:
        phi_I = -0.416
    else:
        phi_I = phi_I
        
    #longitude of the IPP
    lambda_I = longi+((ECA*math.sin(Az))/math.cos(phi_I))

    #geomagnetic latitude of the IPP
    phi_m = phi_I+0.064*(math.cos(lambda_I-1.617))

    #local time at IPP
    t = (43200*lambda_I)+tsecs
    if t >= 86400:
        t = t-86400
    elif t < 0:
        t = t+86400
    else:
        t = t
        
    #Amplitude of the Ionospheric Delay
    A_i = alpha_i[0]*math.pow(phi_m,0)+alpha_i[1]*math.pow(phi_m,1)+alpha_i[2]*math.pow(phi_m,2)+alpha_i*math.pow(phi_m,3)
    if A_i < 0:
        A_i = 0
    else:
        A_i = A_i
        
    #Period of ionospheric delay
    P_i = beta_i[0]*math.pow(phi_m,0)+beta_i[1]*math.pow(phi_m,1)+beta_i[2]*math.pow(phi_m,2)+beta_i*math.pow(phi_m,3)
    if P_i < 72000:
        P_i = 72000
        
    #Phase of ionospheric delay
    X_i = (2*math.pi*(t-50400))/P_i
    #Compute the slant factor
    F = 1.0+(16.0*math.pow((0.53-ECA),3))
    #The Ionospheric Delay
    amp = alpha_i[0]*math.pow(phi_m,0)+alpha_i[1]*math.pow(phi_m,1)+alpha_i[2]*math.pow(phi_m,2)+alpha_i*math.pow(phi_m,3)
    
    if X_i <= 1.5:
        Iono_delay = 5e-9+(amp*(1-(math.pow(X_i,2)/2)+(math.pow(X_i,4)/24)))
    else:
        Iono_delay = 5e-9*F

    return Iono_delay
'''
This function converts ECEF (XYZ) coordinate into Geodetic coordinates 
(Latitude,Longitude,Ellipsoidal Height)
x - ECEF x coordinate
y - ECEF y coordinate
z - ECEF z coordinate
returns - Array that contains the Geodetic coordinates of the station in the 
          form [Latitude,Longitude,Height] 
'''   
def getLatLong(x,y,z):
    #semi major axis of the WGS84 ellipsoid
    a = 6378137.0
    # semi minor axis of the WGS84 ellipsoid
    b = 6356752.314245
    #reciprocal of flattening
    reciprocal = 298.257223563
    #flattening f
    f = 1/reciprocal
    #eccentricity e
    e = math.sqrt((2*f-math.pow(f,2)))
    
    #distance of the point from the Z axis for height 0
    p = math.sqrt(math.pow(x,2)+math.pow(y,2))
    #the latitude of the point for height zero
    phi0 = math.atan(z/((1-math.pow(e,2))*p))
    #Radius of the curvature of the prime vertical section for height zero
    N0 = math.pow(a,2)/math.sqrt((math.pow(a,2)*math.pow(math.cos(phi0),2))+(math.pow(b,2)*math.pow(math.sin(phi0),2)))
    #calculate the height using the calculated values of p, phi0, and N0
    h0 = (p/math.cos(phi0))-N0
    #calculate the phi with the new h value
    phi = math.atan(z/(p*(1-(math.pow(e,2)*N0))/(N0-h0)))
    
    while abs(phi0-phi)>0.0000001:
        phi0 = phi
        N = math.pow(a,2)/math.sqrt((math.pow(a,2)*math.pow(math.cos(phi),2))+(math.pow(b,2)*math.pow(math.sin(phi),2)))
        h = (p/math.cos(phi))-N
        zp = z/p
        eN = (1-((math.pow(e,2)*N)/(N+h)))
        phi = math.atan2(zp,eN)
        
    #Compute Latitude and Longitude
    lat = phi*(180/math.pi)
    lon = math.acos(x/(N*math.cos(phi)))*(180/math.pi)

    return (lat,lon,h)
'''
Computes the topocentric coordinates (Azimuth,Elevation,Distance) of a satellite
given the receiver position in geodetic coordinates 
(Latitude,Longitude,Ellipsoidal Height) and the position vector from the 
receiver to the satellite ECEF (x,y,z)
u_latlong - Array that contains [Latitude,Longitude,Height] of the station
sat_vec - Array that contains the position vector from the receiver to the 
          satellite in the form [dx,dy,dz]
returns - Array that contains the topocentric coordinates of the satellite
          in the form [Azimuth,Elevation,Distance] 
''' 
def getElevAz(u_latlong,sat_vec):

    lat = u_latlong[0]
    longi = u_latlong[1]
    h = u_latlong[2]
    clong = math.cos(longi*(math.pi/180))
    slong = math.sin(longi*(math.pi/180))
    clat = math.sin(lat*(math.pi/180))
    slat = math.sin(lat*(math.pi/180))
    R = np.array([[-slong,(-slat*clong),(clat*clong)],[clong,(-slat*slong),(clat*slong)],[0,clat,slat]])
    vec = np.dot(R.transpose(),sat_vec)
    E = vec[0]
    N = vec[1]
    U = vec[2]
    hor_dis = math.sqrt(math.pow(N,2)+math.pow(E,2))

    if hor_dis < 1e-20:
        Az = 0
        El = 90;
    else:
        Az = math.atan2(E,N)*(180/math.pi)
        El = math.atan2(U,hor_dis)*(180/math.pi)

    if Az < 0:
        Az = Az+360
    D = math.sqrt(math.pow(sat_vec[0],2)+math.pow(sat_vec[1],2)+math.pow(sat_vec[2],2))

    return(Az,El,D)

'''
Determines topocentric Elevation and Azimuth angles to the satellite from the 
geodetic receiver coordinates (Latitude,Longitude,Ellipsoidal Height) and 
the position vector from the receiver to the satellite ECEF (x,y,z)

Reference: Coordinate Systems in Geodesy, E.J. Krakiwsky and D.E. 
            Wells May 1971, UNB. Page 101
u_latlong - Array that contains [Latitude,Longitude,Height] of the station
sat_vec - Array that contains the position vector from the receiver to the 
          satellite in the form [dx,dy,dz]
returns - Array that contains the topocentric coordinates of the satellite
          in the form [Azimuth,Elevation,Range] 
'''
def altAz(u_latlong,sat_vec):
    lat = u_latlong[0]
    longi = u_latlong[1]
    h = u_latlong[2]
    sat_range = math.sqrt(math.pow(sat_vec[0],2)+math.pow(sat_vec[1],2)+math.pow(sat_vec[2],2))
    #transform the range vector to Local geodetic frame
    xy_LG = np.dot(np.dot(np.dot(P2(),rot2(lat-90)),rot3(longi-180)),sat_vec)
    #Elevation and Azimuth angles
    alt = math.asin(xy_LG[2]/sat_range)*(180/math.pi)
    Az = math.atan2(xy_LG[1],xy_LG[0])*(180/math.pi)
    return(Az,alt,sat_range)
'''
Returns a rotation matrix about the X axis for the input angle
x - rotation angle in degrees
returns - A 3x3 rotation matrix in the X axis
'''    
def rot1(x):
    ang = x*(math.pi/180)
    return np.array([[1,0,0],[0,math.cos(ang),math.sin(ang)],[0,-math.sin(ang),math.cos(ang)]])
'''
Returns a rotation matrix about the Y axis for the input angle
x - rotation angle in degrees
returns - A 3x3 rotation matrix in the Y axis
'''  
def rot2(x):
    ang = x*(math.pi/180)
    y = np.array([[math.cos(ang),0,-math.sin(ang)],[0,1,0],[math.sin(ang),0,math.cos(ang)]])
    return y
'''
Returns a rotation matrix about the Z axis for the input angle
x - rotation angle in degrees
returns - A 3x3 rotation matrix in the Z axis
'''  
def rot3(x):
    ang = x*(math.pi/180)
    y = np.array([[math.cos(ang),math.sin(ang),0],[-math.sin(ang),math.cos(ang),0],[0,0,1]])
    return y
'''
Returns a reflection matrix about the Y axis
returns - A 3x3 reflection matrix in the Y axis
'''  
def P2():
    y = np.array([[1,0,0],[0,-1,0],[0,0,1]])
    return y

    
                  
