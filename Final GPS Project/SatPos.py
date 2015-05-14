""" This function reads in the RINEX navigation parameters and outputs the ECEF position of the satellite"""
from math import *
import datetime,time
import numpy as np

secsInWeek = 604800
secsInDay = 86400
gpsEpoch = (1980,1,6,0,0,0)
mu = 3.986005e+14
OHMDOTe = 7.2921151467e-5
c = 2.99792458E+8
F = -4.442807633e-10
PI = 3.1415926535898

def SatPos(nav_data,prn,epoch,pseudo):
    #Constants and the Keplerian elements
    
    #sat Position
    sat_epoch=cTime(epoch, nav_data[prn])
    if sat_epoch == 'none':
        return [0,0,0,0]
    #print "sat epoch " + str(getSecs(sat_epoch))
    #Declare Variables
    ecc = nav_data[prn][sat_epoch]['Ecc']
    A=pow(nav_data[prn][sat_epoch]['SqrtA'],2)
    deln=nav_data[prn][sat_epoch]['Deln']
    TOE=nav_data[prn][sat_epoch]['TOE']
    M0=nav_data[prn][sat_epoch]['Mo']
    a0=nav_data[prn][sat_epoch]['SV_CLB']
    a1=nav_data[prn][sat_epoch]['SV_CLD']
    a2=nav_data[prn][sat_epoch]['SV_CLR']
    tgd=nav_data[prn][sat_epoch]['TGD']
    #Argument of the Perigee
    ohm=nav_data[prn][sat_epoch]['Omega']
    
    cus=nav_data[prn][sat_epoch]['Cus']
    crs=nav_data[prn][sat_epoch]['Crs']
    cis=nav_data[prn][sat_epoch]['Cis']
    cuc=nav_data[prn][sat_epoch]['Cuc']
    crc=nav_data[prn][sat_epoch]['Crc']
    cic=nav_data[prn][sat_epoch]['Cic']
    
    i0=nav_data[prn][sat_epoch]['Io']
    IDOT=nav_data[prn][sat_epoch]['IDOT']
    OHM=nav_data[prn][sat_epoch]['OMEGA']
    OHMDOT=nav_data[prn][sat_epoch]['OMEGA_DOT']
    
    #Mean Motion
    n0 = sqrt(mu/(pow(A,3)))
    #Corrected Mean Motion
    n = n0 + deln
    traw = getSecs(epoch)-pseudo/c
    #print traw
    tk=check_t(traw-TOE)
    #Mean Anomoly
    Mk=M0+n*tk
    #Mk=fmod(Mk+2*PI,2*PI)
    Ek0=Mk
    Ek=Mk+ecc*sin(Ek0)
    #Ek=fmod(Ek+2*PI,2*PI)
    while abs(Ek-Ek0)>0.000000000001:
        Ek0=Ek
        Ek=Mk+ecc*sin(Ek0)
        #Ek=fmod(Ek+2*PI,2*PI)
        
    dtr = F*ecc*sqrt(A)*sin(Ek)
    dtcorr= a0+a1*tk+a2*tk*tk+dtr-tgd
    tsv = traw-dtcorr
    
    #Compute SatPOS
    #Compute True Anomaly Components
    ta_num = sin(Ek)*sqrt(1-pow(ecc,2))
    ta_den = (cos(Ek)-ecc)
    #Compute True Anomaly
    v_k = atan2(ta_num,ta_den)
    #Argument of the Latitude
    phi_k=v_k+ohm
    #phi_k=fmod(phi_k+2*PI,2*PI)
    #Correction to Argument of Latitude
    del_uk=cus*sin(2*phi_k)+cuc*cos(2*phi_k)
    #Correction to Radius
    del_rk=crs*sin(2*phi_k)+crc*cos(2*phi_k)
    #Correction to Inclination
    del_ik=cis*sin(2*phi_k)+cic*cos(2*phi_k)
    #Corrected Argument of Latitude
    uk = phi_k+del_uk
    #uk = fmod(uk+2*PI,2*PI)
    #Corrected Radius
    rk=A*(1-ecc*cos(Ek))+del_rk
    #Corrected Inclination
    ik=i0+del_ik+IDOT*tk
    #ik=fmod(ik+2*PI,2*PI)
    #Orbital Plan Positions
    x0 = rk*cos(uk)
    y0 = rk*sin(uk)
    #Corrected Longitude of ascending node
    OMEGAK=OHM+tk*(OHMDOT-OHMDOTe)-OHMDOTe*TOE
    #OMEGAK=fmod(OMEGAK+2*PI,2*PI)
    #Earthfixed Coordinates
    xk=x0*cos(OMEGAK)-y0*cos(ik)*sin(OMEGAK)
    yk=x0*sin(OMEGAK)+y0*cos(ik)*cos(OMEGAK)
    zk=y0*sin(ik)
    #Rotate
    rot=pseudo*OHMDOTe/c
    #xtrans=(np.dot(rot3(rot),np.array([[xk],[yk],[zk]])))
    #return [xtrans[0][0],xtrans[1][0],xtrans[2][0],tk]
    return [xk,yk,zk,dtcorr]
"""looks for the largest epoch that is less than the epoch we are searching for"""
def cTime(epoch, dic):
    sat_epoch = 'none'
    t  = ParseTime(epoch)
    #print t
    for key in dic:
        dif = ParseTime(key)-t
        if ParseTime(key)<t and dif.total_seconds()<=7200:
            if sat_epoch == 'none':
                sat_epoch = key
            elif ParseTime(sat_epoch) < ParseTime(key):
                sat_epoch = key
    #if sat_epoch != 'none':
        #print 'e '+epoch+ ' se '+sat_epoch
    return sat_epoch

def ParseTime(epoch):
    t = epoch.split(":")
    return datetime.datetime((2000+int(t[0])),int(t[1]),int(t[2]),int(t[3]),int(t[4]),int(floor(float(t[5]))),int(1000000*((float(t[5]))-int(floor((float(t[5])))))))
'''returns GPS seconds of the week with the inputted date
    kai borre julianday and gpsseconds function'''
def getSecs(epoch):
    t = epoch.split(":")
    yy = 2000+float(t[0])
    mm = float(t[1])
    dd = float(t[2])
    hh = float(t[3])
    mins = float(t[4])
    ss = float(t[5])
    hh = hh+(mins/60)+(ss/3600)
    if mm <= 2:
        yy = yy-1
        mm = mm+12
    jd = floor(365.25*(yy+4716))+floor(30.6001*(mm+1))+dd+(hh/24)-1537.5
    a = floor(jd+0.5)
    b = a+1537
    c = floor((b-122.1)/365.25)
    e = floor(365.25*c)
    f = floor((b-e)/30.6001)
    d = b-e-floor(30.6001*f)+(fmod((jd+0.5),1))
    day_of_week = fmod(floor(jd+0.5),7)
    secs = (fmod(d,1)+day_of_week+1)*86400
    return secs

def check_t(t):
    half_week = 302400
    if t > half_week:
        tt = t-(2*half_week)
    elif t < -1*half_week:
        tt = t+(2*half_week)
    else:
        tt=t
    return tt

'''Rotation matrix about the Z axis for a given angle'''
def rot3(x):
    ang=x
    y = np.array([[cos(ang),sin(ang),0],[-sin(ang),cos(ang),0],[0,0,1]])
    return y
    
