'''
Created on Oct 30, 2014

@author: Patrick and Sonal
'''
from Reader import RINEXREADER,NAVREADER
from adjust import dualFrequencyPosition
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import *

#Read Observation and Navigation File
obs_data= RINEXREADER('ALBH0010.10O')
x0 = np.array([[obs_data['HEAD']['POS']['X']],[obs_data['HEAD']['POS']['Y']],[obs_data['HEAD']['POS']['Z']]])
x0 = np.array([[-2341333.003],[-3539049.514],[4745791.300]])
nav_data = NAVREADER('ALBH0010.10N')
#Perform Adjustment using Dual Frequency Observations
adj = dualFrequencyPosition(obs_data,nav_data)
print 'DONE'
results = np.array(adj[0])
residuals = adj[1]
cov = np.array(adj[2])
sat_xyz = np.array(adj[3])
std = np.array(adj[6])
gdop = np.array(adj[4])
pdop = np.array(adj[5])
res1 = []
gdop1 = []
pdop1 = []

'''for result in results:
    dist = sqrt(pow(result[0]-x0[0],2)+pow(result[1]-x0[1],2)+pow(result[2]-x0[2],2))
    if dist < 25:
        res1.append(result)
    #else:
        #print 'Bigger'


results = np.array(res1)
print results'''
#gdop = np.array(gdop1)
#pdop = np.array(pdop1)

plt.figure(1)
ax = plt.subplot(111, projection='3d')
ax.plot(sat_xyz[:,0],sat_xyz[:,1],sat_xyz[:,2],'.',color = 'b')
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')
ax.set_title('Satellite Position for each Epoch of observation')

plt.figure(2)
plt.plot(gdop[:],label='GDOP')
plt.plot(pdop[:],label = 'PDOP')
plt.xlabel('Epoch')
plt.legend()
plt.title('DOP Values for each Epoch')

plt.figure(4)
#ax1 = plt.subplot(111,projection='3d')
#ax1.plot(results[:,0],results[:,1],results[:,2],'.',color = 'b')
#plt.scatter(results[:,0],results[:,1])
plt.plot(x0[0]-results[:,0],'b',label = 'X')
plt.plot(x0[1]-results[:,1],'r',label = 'Y')
plt.plot(x0[2]-results[:,2],'g',label = 'Z')
plt.legend()
plt.xlabel('Epoch')
plt.ylabel('Difference [m]')
plt.title('Difference in Solution Coordinates from Intial')
#ax1.set_xlabel('X [m]')
#ax1.set_ylabel('Y [m]')
#ax1.set_zlabel('Z [m]')
#ax1.set_title('Receiver ECEF Coordinates [m]')

plt.figure(5)
'''ax2 = plt.subplot(111, projection = '3d')
ax2.plot(results[:,0],results[:,1],results[:,2],'.',color='b')
ax2.set_xlabel('X [m]')
ax2.set_ylabel('Y [m]')
ax2.set_zlabel('Z [m]')
ax2.set_title('Receiver Coordinates')'''
plt.plot(results[:,0],results[:,1],'.')
plt.xlabel('X [m]')
plt.ylabel('Y [m]')
plt.title('Receiver ECEF Coordiantes')

plt.figure(6)
plt.plot(std[:,0],'r',label = 'X')
plt.plot(std[:,1],'g',label = 'Y')
plt.plot(std[:,2],'b',label = 'Z')
plt.plot(std[:,3],'k',label = 'Time')
plt.legend()
plt.xlabel('Epoch')
plt.ylabel('Standard Deviation')
plt.title('Standard deviation of each element of the solution')
plt.show()

'''print 'X: '+str(np.average(results[:,0]))
print 'Y: '+str(np.average(results[:,1]))
print 'Z: '+str(np.average(results[:,2]))
print 't: '+str(np.average(results[:,3]))

print 'dX: '+str(np.average(std[:,0]))
print 'dY: '+str(np.average(std[:,1]))
print 'dZ: '+str(np.average(std[:,2]))
print 'dt: '+str(np.average(std[:,3]))

print 'gdop: '+str(np.average(gdop[:]))
print 'pdop: '+str(np.average(pdop[:]))

print 'stdX: '+str(np.std(results[:,0]))
print 'stdY: '+str(np.std(results[:,1]))
print 'stdZ: '+str(np.std(results[:,2]))
print 'stdt: '+str(np.std(results[:,3]))

print 'difX: '+str(np.average(results[:,0])-x0[0])
print 'difY: '+str(np.average(results[:,1])-x0[1])
print 'difZ: '+str(np.average(results[:,2])-x0[2])'''
