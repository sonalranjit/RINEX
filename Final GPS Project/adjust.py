'''
Created on Nov 13, 2014

@author: Patrick
'''
import numpy as np
from math import *
from SatPos import SatPos
from SatPos import getSecs
from dtropo import dTropo
from numpy import shape
import datetime,time

firstepoch =(('10:1:1:2:0:30.0000000'))
c = 2.99792458e8
''' Perform coarse adjustment using CA (on S1) pseudoranges'''
def singleFrequencyPosition(obs_data, nav_data):
    e_list = obs_data['OBS']['LIST']
    x0 = obs_data['HEAD']['POS']['X']
    y0 = obs_data['HEAD']['POS']['Y']
    z0 = obs_data['HEAD']['POS']['Z']
    #SatPos(nav_data,'G01','07:1:1:1:0:0.0',22581697.779)

''' Perform adjustment using P1 (S1) and P2 (S2) pseudoranges'''
def dualFrequencyPosition(obs_data, nav_data):
    e_list = obs_data['OBS']['LIST']
    #initial Estimates (From HEADER), Time offset is estimated to be zero
    x0=np.array([[obs_data['HEAD']['POS']['X']],[obs_data['HEAD']['POS']['Y']],[obs_data['HEAD']['POS']['Z']],[0.0]])
    #print x0
    results = []
    cx = []
    res = []
    sat_xyz = []
    GDOP = []
    PDOP = []
    stds = []
    of = open('results.txt',"w")
    of1 = open('cx.txt',"w")
    of2 = open('residuals.txt',"w")
    #iterate through the list of epochs
    for i in xrange(0,len(e_list)):
        epoch = e_list[i] #a list of all epochs strings of the format 07:1:1:1:0:0.0
        #x0=np.array([[obs_data['HEAD']['POS']['X']],[obs_data['HEAD']['POS']['Y']],[obs_data['HEAD']['POS']['Z']],[0.0]])
        
        #print epoch
        #get Observation vector
        sat_list = CUTLIST(obs_data['OBS'][epoch],'Dfr')
        #sat_list = obs_data['OBS'][epoch]['LIST']
        #print sat_list
        PR = IONOFREEOBS(READOBS(sat_list,obs_data['OBS'][epoch],'P1'),READOBS(sat_list,obs_data['OBS'][epoch],'P2'))
        SPOS = SATCONSTANTS(epoch,sat_list,PR,nav_data)
        #print SPOS
        sat_list=SPOS[2]
        
        #print sat_list
        if (len(sat_list)<4):
            #results.append([0,0,0,0])
            #cx.append([0])
            #res.append([0])
            #APPLY THINGS
            continue
        #W matrix
        #WTF is goin on here??
        PR = SPOS[1]
        SPOS=SPOS[0]
        
        #A Matrix and Tropospheric Corrections
        ATropo = Amat(SPOS, x0)
        
        A = ATropo[0]
        Tropo = ATropo[1]
        #print epoch
        #l matrix
        #print SPOS
        #print SPOS[:,3:]

        #WTF is going on here???????
        #l = np.subtract(np.add(PR,np.multiply(SPOS[:,3:],c)),Tropo)
        l = np.add(PR,np.multiply(SPOS[:,3:],c))
        fx = Fx(x0, SPOS)
        w = l-fx-Tropo
        #print fx.shape
        #print w
        
        #delm = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(A),A)),np.transpose(A)),w)
        #x0 = x0+delm
        #print x0.shape
#        while SMALLENOUGH(delm, 0.0001): 
        for i in range(0,5):
            ATropo = Amat(SPOS, x0)
            A = ATropo[0]
            Tropo = ATropo[1]
            fx = Fx(x0, SPOS)
            w = l-fx
            delm = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(A),A)),np.transpose(A)),w)
            x0 = x0+delm
            #print delm
        cxx = np.linalg.inv(np.dot(np.transpose(A),A))    
        cx.append(cxx)
        x00 = np.ndarray.flatten(x0)
        print x00
        results.append([x00[0],x00[1],x00[2],x00[3]])
        ress = np.subtract(w,np.dot(A,delm))
        res.append(ress)
        x00 = np.ndarray.flatten(x0)
  

        for pos in SPOS:
            
            sat_xyz.append(pos)

        gdop = sqrt(cxx[0][0]+cxx[1][1]+cxx[2][2]+c*(cxx[3][3]))
        pdop = sqrt(cxx[0][0]+cxx[1][1]+cxx[2][2])
        stds.append([sqrt(cxx[0][0]),sqrt(cxx[1][1]),sqrt(cxx[2][2]),sqrt(cxx[3][3]),])
        GDOP.append(gdop)
        PDOP.append(pdop)
        of.write(str(getSecs(epoch))+' '+str(x00[0])+' '+str(x00[1])+' '+str(x00[2])+' '+str(x00[3])+'\n')
        of1.write(str(getSecs(epoch))+' '+str(cxx[0][0])+' '+str(cxx[1][1])+' '+str(cxx[2][2])+' '+str(cxx[3][3])+'\n')
        of2.write(str(getSecs(epoch))+' '+str(ress[0])+' '+str(ress[1])+' '+str(ress[2])+' '+str(ress[3])+'\n')
    of.close()
    of1.close()
    of2.close()
    return [results,res,cx,sat_xyz,GDOP,PDOP,stds]      
         
''' obs_data =  obs_data['OBS'][epoch]'''
def READOBS(sat_list,obs_data,obstype):
    obs = []  
    #Iterate Through Satellites
    for sat in sat_list:
        obs.append([obs_data[sat][obstype]])
    return np.array(obs)
##wtf is this
def Amat(SPOS,x0):
    x0 = np.ndarray.flatten(x0)
    geodeticCoords = getLatLong(x0[0],x0[1],x0[2])
    A = []
    Tropo = []
    for pos in SPOS:
        p0 = sqrt(pow(pos[0]-x0[0],2)+pow(pos[1]-x0[1],2)+pow(pos[2]-x0[2],2))
        a1 = (-(pos[0]-x0[0])/p0)
        a2 = (-(pos[1]-x0[1])/p0)
        a3 = (-(pos[2]-x0[2])/p0)
        A.append([a1,a2,a3,c])
        elevation = 90 - getElevAz(geodeticCoords, [pos[0]-x0[0],pos[1]-x0[1],pos[2]-x0[2]])[1]
        Tropo.append([dTropo(elevation)])
    A=np.array(A)
    Tropo= np.array(Tropo)
    #return A
    return [A,Tropo]   
''' Sfr (single frequency (c1))
    Dfr (dual frequency (p1,p2))
    obs_data['OBS'][epoch]'''
#WTFF is this
def CUTLIST(obs_data,obstype):
    sat_list = obs_data['LIST']
    append_list = []
    for sat in sat_list:
        if obstype == 'Sfr':
            if not (isnan(obs_data[sat]['C1'])):
                append_list.append(sat) 
        elif obstype == 'Dfr':
            if not (isnan(obs_data[sat]['P1'])):
                if not (isnan(obs_data[sat]['P2'])):
                    append_list.append(sat) 
        else:
            if not (isnan(obs_data[sat][obstype])):
                append_list.append(sat) 
    return append_list
#evaluates function at x0 for every value in spos
#does not include -c(dtrs) which is applied to the observations
def Fx(x0, spos):
    x0 = np.ndarray.flatten(x0)
    #print x0
    fx = []
    for pos in spos:
        fx.append([sqrt(pow(pos[0]-x0[0],2)+pow(pos[1]-x0[1],2)+pow(pos[2]-x0[2],2))+c*x0[3]])
    return np.array(fx)

def IONOFREEOBS(P1OBS, P2OBS):
    gamma = pow(77.0/60.0,2)
    PR = ((P2OBS-gamma*P1OBS)/(1-gamma))
    return PR
    
    
def SMALLENOUGH(delm,pre):
    for val in delm:
        if val < pre:
            return True
    return False
    
def ParseTime(epoch):
    t = epoch.split(":")
    return datetime.datetime((2000+int(t[0])),int(t[1]),int(t[2]),int(t[3]),int(t[4]),int(floor(float(t[5]))),int(1000000*((float(t[5]))-int(floor((float(t[5])))))))
# how the hell is this working??????
def SATCONSTANTS(epoch,sat_list,rangeOBS,nav_data):
    i = 0
    satcon = []
    rOBS = []
    slist=[]
    '''FIX DIS'''
    
    for sat in sat_list:
        #print sat
        sp = SatPos(nav_data,sat,epoch,rangeOBS[i][0])
        #print sp[0]
        if sp[0]!=0.0:
            if sat == 'G01':
                dif = ParseTime(epoch)-ParseTime(firstepoch)
                #print str(dif.total_seconds())+' '+str(sp[0]) + ' ' + str(sp[1]) + ' ' + str(sp[2]) + ' ' + str(sp[3])
            slist.append(sat)
            satcon.append(sp)
            rOBS.append([rangeOBS[i][0]])
            
        i=i+1
    #print slist
    
    return [np.array(satcon),np.array(rOBS),slist]
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
    e = sqrt((2*f-pow(f,2)))
    #distance of the point from the Z axis for height 0
    p = sqrt(pow(x,2)+pow(y,2))
    #the latitude of the point for height zero
    phi0 = atan(z/((1-pow(e,2))*p))
    #Radius of the curvature of the prime vertical section for height zero
    N0 = pow(a,2)/sqrt((pow(a,2)*pow(cos(phi0),2))+(pow(b,2)*pow(sin(phi0),2)))
    #calculate the height using the calculated values of p, phi0, and N0
    h0 = (p/cos(phi0))-N0
    #calculate the phi with the new h value
    phi = atan(z/(p*(1-(pow(e,2)*N0))/(N0-h0)))
    
    while abs(phi0-phi)>0.0000001:
        phi0 = phi
        N = pow(a,2)/sqrt((pow(a,2)*pow(cos(phi),2))+(pow(b,2)*pow(sin(phi),2)))
        h = (p/cos(phi))-N
        zp = z/p
        eN = (1-((pow(e,2)*N)/(N+h)))
        phi = atan2(zp,eN)
        
    #Compute Latitude and Longitude
    lat = phi*(180/pi)
    lon = acos(x/(N*cos(phi)))*(180/pi)
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
    clong = cos(longi*(pi/180))
    slong = sin(longi*(pi/180))
    clat = sin(lat*(pi/180))
    slat = sin(lat*(pi/180))
    R = np.array([[-slong,(-slat*clong),(clat*clong)],[clong,(-slat*slong),(clat*slong)],[0,clat,slat]])
    vec = np.dot(R.transpose(),sat_vec)
    E = vec[0]
    N = vec[1]
    U = vec[2]
    hor_dis = sqrt(pow(N,2)+pow(E,2))

    if hor_dis < 1e-20:
        Az = 0
        El = 90;
    else:
        Az = atan2(E,N)*(180/pi)
        El = atan2(U,hor_dis)*(180/pi)

    if Az < 0:
        Az = Az+360
    D = sqrt(pow(sat_vec[0],2)+pow(sat_vec[1],2)+pow(sat_vec[2],2))

    return(Az,El,D)
