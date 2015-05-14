'''This function calculates the tropospheric delay using the
    Saastamoinen Closed Form Model. Inputs are the zenith angle to the satellite of interest'''
from math import *
'''Saastamoinen model is defined as:
    drop: 0.002277/cos(z)[pd+(1255/T + 0.05)pwv - (tan(z))^2]
    where:
    dtrop is the tropospheric delay
    z is the zenith angle of the satellite from the receiver
    pd is the partial pressure of dry air in mbar
    pwv is the partial pressure of water vapour in mbar
    T is the absolute temperature in Kelvins'''
def dTropo(z):
    T = 273.25
    pd = 1013.24
    pwv = 23.39

    dtrop = (0.002277/(cos(z*(pi/180)))*(pd+((1255/T)+0.005)*pwv - pow(tan(z),2)))

    return dtrop
