#BRANDON NOBLE
#3Y03 Stellar Structure: Astrometry Project
#February 12, 2021

import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import csv
from math import sqrt

fileM35 = 'E:\\Documents\\School\\3Y03\\1611340333922O-result.csv' #M35
fileNGC = 'E:\\Documents\\School\\3Y03\\1611340488926O-result.csv' #NGC 6397
M35 = np.genfromtxt(fileM35, dtype=('float'), delimiter=',') #break up the file into each row
NGC = np.genfromtxt(fileNGC, dtype=('float'), delimiter=',') #break up the file into each row

#define variables for easier indexing
src = 0 #source id
ra = 1 #right ascension
raE = 2 #right ascension error
dec = 3 #declination
decE = 4 #declination error
par = 5 #parallax
parE = 6 #parallax error
pmA = 7 #proper motion in right ascension
pmAE = 8 #error in proper motion, RA
pmD = 9 #proper motion in declination
pmDE = 10 #error in peroper motion, DEC
phot = 11 #brightness in G magnitude
bprp = 12 #colour in BR filters
rV = 13 #radial velocity
rVE = 14 #radial velocity error


################################################################
################################################################
################################################################
#           PLOTTING OF M35 PROPER MOTION AND HR PLOT
################################################################
################################################################
################################################################
print("---M35---")

xVals = []
yVals = []

#points that have been established within the M35 cluster
M35vals = []
rad = 0.4

for i in range(1,len(M35)):
    x = M35[i][pmA]
    y = M35[i][pmD]

    if sqrt( (x-2.31)**2 + (y+2.99)**2 ) <= rad: #selects stars within the cluster, using the 'centre' and approximate radius
        if sqrt( (M35[i][ra]-92.25)**2 + (M35[i][dec]-24.33)**2 ) <= 0.5: #restricting to stars within a 0.5 degree radius following research paper (https://ui.adsabs.harvard.edu/abs/2009ApJ...695..679M/abstract)
            M35vals.append(M35[i])

    xVals.append(x)
    yVals.append(y)

#print(len(M35vals)) #1049 selected stars for M35

#~~~CLUSTER PLOT FOR M35~~~: Centre Around (2.31, -2.99)
f1 = plt.figure(1)
plt.subplot(121)
plt.scatter(xVals,yVals,s=0.001,c='b',marker='o')
circle1 = plt.Circle((2.31, -2.99), rad, color='r', fill=False)
fig = plt.gcf()
ax = fig.gca()
ax.add_patch(circle1)
plt.axis([-3,4,-8,4])
plt.xlabel('Proper Motion - Right Ascension')
plt.ylabel('Proper Motion - Declination')
plt.title('Proper Motion of M35 Cluster')

xVals = []
yVals = []
MxVals = []
MyVals = []
avgP = 0.0

for i in range(1,len(M35)):
    x = M35[i][bprp]
    y = M35[i][phot]
    xVals.append(x)
    yVals.append(y)

for i in range(1,len(M35vals)):
    x = M35vals[i][bprp]
    y = M35vals[i][phot]
    MxVals.append(x)
    MyVals.append(y)
    avgP += M35vals[i][pmA]

avgP /= len(M35vals)
M35pmA = avgP
print("Proper Motion RA: ", M35pmA)

#~~~HR DIAGRAM FOR M35~~~
plt.subplot(122)
plt.scatter(xVals,yVals,s=0.05,c='b',marker='o',label='All Data') #plot all the data collected
plt.scatter(MxVals,MyVals,s=0.5,c='r',marker='o',label='Located Cluster') #overplot of the actual cluster
plt.axis([-0.5,3,21,9])
plt.legend(loc='upper right', markerscale=10.0)
plt.legend.framealpha = 0.8
plt.xlabel('Colour in B-R Filters')
plt.ylabel('Brightness in G Magnitude')
plt.title('Colour-Magnitude Diagram for M35 Cluster')

#######################################################################################################

#Was used to plot the RA vs DEC of the data to help find the correct stars

xVals = []
yVals = []
MxVals = []
MyVals = []
avgP = 0.0

for i in range(1,len(M35)):
    x = M35[i][ra]
    y = M35[i][dec]
    xVals.append(x)
    yVals.append(y)

for i in range(1,len(M35vals)):
    x = M35vals[i][ra]
    y = M35vals[i][dec]
    MxVals.append(x)
    MyVals.append(y)
    avgP += M35vals[i][pmD]

avgP /= len(M35vals)
M35pmD = avgP
print("Proper Motion DEC: ", M35pmD)

#M35 cluster plot
f2 = plt.figure(2)
plt.scatter(xVals,yVals,s=0.001,c='b',marker='o')
plt.scatter(MxVals,MyVals,s=0.1,c='r',marker='o')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('RA vs DEC for M35 Cluster')

#Calcualte Radial Velocity Average
radV = []

for i in range(1,len(M35vals)):
    if np.isnan(M35vals[i][rV]):
        continue
    if M35vals[i][rVE] <= 6:
        radV.append(M35vals[i])

r = 0.0
re = 0.0

for i in range(0,len(radV)):
    r += radV[i][rV]
    re += radV[i][rVE]
r /= len(radV)
re /= len(radV)
print("Radial Velocity: ", r," +/- ",re)

#Tangential Velocity
Vt = []
VtE = []
mu = []
parallax = []

for i in range(1,len(M35vals)):
    a = M35vals[i][pmA]
    d = M35vals[i][pmD]
    p = M35vals[i][par]
    ea = M35vals[i][pmAE]
    ed = M35vals[i][pmDE]
    ep = M35vals[i][parE]

    u = sqrt((a/1000)**2 + (d/1000)**2)
    Vtan = 4.74*u/(p/1000)
    VtanE = 4.74*sqrt((ea/1000)**2 + (ed/1000)**2)/(ep/1000)
    if abs(Vtan) > 150:
        continue
    Vt.append(Vtan)
    VtE.append(VtanE)
    mu.append(u)
    parallax.append(p/1000)

Vtan = np.sum(Vt)/len(Vt)
VtanE = np.sum(VtE)/len(VtE)
muVal = np.sum(mu)/len(mu)
plax = np.sum(parallax)/len(parallax)

#Find the average velocity of stars in the cluster
VtVals = []
for i in range(0,len(Vt)):
    if abs(Vt[i] - Vtan) < 10:
        VtVals.append(Vt[i] - Vtan)

print("Star Velocity, 1STD: ",np.std(VtVals))

f12 = plt.figure(12)
plt.hist(VtVals, bins=100)
plt.xlabel('Velocity')
plt.ylabel('Amount')
plt.title('M35 Cluster: Star Velocity - Cluster Tangential Velocity')

#Assuming (quite inaccurately) that V_R = V_T = 4.74*mu*d, we can get a rough estimate for d
print("Tangential Velocity: ",Vtan," +/-",VtanE)
d = abs(r/(4.74*muVal))
print("Distance Using Radial: ", d)
d = 1/plax
print("Distance using Parallax: ", d)
d = Vtan/(4.74*muVal)
print("Distance Using Tangential: ", d)
print("")


################################################################
################################################################
################################################################
#       PLOTTING OF NGC 6397 PROPER MOTION AND HR PLOT
################################################################
################################################################
################################################################
print("---NGC 6397---")

xVals = []
yVals = []

#points that have been established within the NGC cluster
NGCvals = []
rad = 0.5 #radius of the cluster on plot, started with 2 and determined experimentally that 0.5 was a good fit for the data

for i in range(1,len(NGC)):
    x = NGC[i][pmA]
    y = NGC[i][pmD]

    if sqrt( (x-3.24)**2 + (y+17.58)**2 ) <= rad: #selects stars within the cluster, using the 'centre' and approximate radius
        if sqrt( (NGC[i][ra]-265.2)**2 + (NGC[i][dec]+53.67)**2 ) <= 0.5: #using 0.5 degree radius as suggested by (https://ui.adsabs.harvard.edu/abs/1997A%26AS..122....1K/abstract)
            NGCvals.append(NGC[i])

    xVals.append(x)
    yVals.append(y)

#print(len(NGCvals)) #21 718 selected stars for NGC

#~~~CLUSTER PLOT FOR NGC~~~: centre around (3.24, -17.58)
f3 = plt.figure(3)
plt.subplot(121)
plt.scatter(xVals,yVals,s=0.001,c='b',marker='o')
circle1 = plt.Circle((3.24, -17.58), rad, color='r', fill=False)
fig = plt.gcf()
ax = fig.gca()
ax.add_patch(circle1)
plt.axis([-11,8,-22,7])
plt.xlabel('Proper Motion - Right Ascension')
plt.ylabel('Proper Motion - Declination')
plt.title('Proper Motion of NGC 6397 Cluster')

xVals = []
yVals = []
NxVals = []
NyVals = []
avgP = 0.0

for i in range(1,len(NGC)):
    x = NGC[i][bprp]
    y = NGC[i][phot]
    xVals.append(x)
    yVals.append(y)

for i in range(1,len(NGCvals)):
    x = NGCvals[i][bprp]
    y = NGCvals[i][phot]
    NxVals.append(x)
    NyVals.append(y)
    avgP += NGCvals[i][pmA]

avgP /= len(NGCvals)
NGCpmA = avgP
print("Proper Motion RA: ", NGCpmA)

#~~~HR DIAGRAM FOR NGC~~~
plt.subplot(122)
plt.scatter(xVals,yVals,s=0.05,c='b',marker='o', label='All Data') #plot of all points that were downloaded
plt.scatter(NxVals,NyVals,s=0.5,c='r',marker='o', label='Located Cluster') #overplot of the actual cluster
plt.axis([-1,3,21,9])
plt.xlabel('Colour in B-R Filters')
plt.ylabel('Brightness in G Magnitude')
plt.title('Colour-Magnitude Diagram for NGC 6397 Cluster')
plt.legend(loc='upper left', markerscale=10.0)
plt.legend.framealpha = 0.8

#######################################################################################################

#Was used to plot the RA vs DEC of the data to help find the correct stars

xVals = []
yVals = []
NxVals = []
NyVals = []
avgP = 0.0

for i in range(1,len(NGC)):
    x = NGC[i][ra]
    y = NGC[i][dec]
    xVals.append(x)
    yVals.append(y)

for i in range(1,len(NGCvals)):
    x = NGCvals[i][ra]
    y = NGCvals[i][dec]
    NxVals.append(x)
    NyVals.append(y)
    avgP += NGCvals[i][pmD]

avgP /= len(NGCvals)
NGCpmD = avgP
print("Proper Motion DEC: ", NGCpmD)

#NGC cluster plot
f4 = plt.figure(4)
plt.scatter(xVals,yVals,s=0.001,c='b',marker='o')
plt.scatter(NxVals,NyVals,s=0.1,c='r',marker='o')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('RA vs DEC for NGC 6397 Cluster')

#Calcualte Radial Velocity Average
radV = []

for i in range(1,len(NGCvals)):
    if np.isnan(NGCvals[i][rV]):
        continue
    if NGCvals[i][rVE] <= 6:
        radV.append(NGCvals[i])

r = 0.0
re = 0.0

for i in range(0,len(radV)):
    r += radV[i][rV]
    re += radV[i][rVE]
r /= len(radV)
re /= len(radV)
print("Radial Velocity: ", r," +/- ",re)

#Tangential Velocity
Vt = []
VtE = []
mu = []
parallax = []

for i in range(1,len(NGCvals)):
    a = NGCvals[i][pmA]
    d = NGCvals[i][pmD]
    p = NGCvals[i][par]
    ea = NGCvals[i][pmAE]
    ed = NGCvals[i][pmDE]
    ep = NGCvals[i][parE]

    u = sqrt((a/1000)**2 + (d/1000)**2)
    Vtan = 4.74*u/(p/1000)
    VtanE = 4.74*sqrt((ea/1000)**2 + (ed/1000)**2)/(ep/1000)
    if abs(Vtan) > 1500:
        continue
    Vt.append(Vtan)
    VtE.append(VtanE)
    mu.append(u)
    parallax.append(p/1000)

Vtan = np.sum(Vt)/len(Vt)
VtanE = np.sum(VtE)/len(VtE)
muVal = np.sum(mu)/len(mu)
plax = np.sum(parallax)/len(parallax)

#Find the average velocity of stars in the cluster
VtVals = []
for i in range(0,len(Vt)):
    if abs(Vt[i] - Vtan) < 200:
        VtVals.append(Vt[i] - Vtan)

print("Star Velocity, 1STD: ",np.std(VtVals))

f13 = plt.figure(13)
plt.hist(VtVals, bins=100)
plt.xlabel('Velocity')
plt.ylabel('Amount')
plt.title('NGC 6397 Cluster: Star Velocity - Cluster Tangential Velocity')

#Assuming (quite inaccurately) that V_R = V_T = 4.74*mu*d, we can get a rough estimate for d
mu = sqrt( (NGCpmA/1000)**2 + (NGCpmD/1000)**2)
print("Tangential Velocity: ",Vtan," +/-",VtanE)
d = r/(4.74*muVal)
print("Distance using Radial: ", d)
d = 1/plax
print("Distance using Parallax: ", d)
d = Vtan/(4.74*muVal)
print("Distance Using Tangential: ", d)
print("")

################################################################
################################################################
################################################################
#       PLOTTING OF OTHER CLUSTER PROPER MOTION AND HR PLOT
################################################################
################################################################
################################################################

fileOther = 'E:\\Documents\\School\\3Y03\\1612831725320O-result.csv' #Other
Other = np.genfromtxt(fileOther, dtype=('float'), delimiter=',') #break up the file into each row

xVals = []
yVals = []

#points that have been established within the Other cluster
Othervals = []
rad = 2

for i in range(1,len(Other)):
    x = Other[i][pmA]
    y = Other[i][pmD]

    if sqrt( (x+2.7)**2 + (y+2.2)**2 ) <= rad: #selects stars within the cluster, using the 'centre' and approximate radius
        Othervals.append(Other[i])

    xVals.append(x)
    yVals.append(y)

#~~~CLUSTER PLOT FOR Other~~~: Centre Around (-2.7, -2.2)
f5 = plt.figure(5)
plt.subplot(121)
plt.scatter(xVals,yVals,s=0.01,c='b',marker='o')
circle1 = plt.Circle((-2.7, -2.2), rad, color='r', fill=False)
fig = plt.gcf()
ax = fig.gca()
ax.add_patch(circle1)
plt.axis([-20,20,-20,20])
plt.xlabel('Proper Motion - Right Ascension')
plt.ylabel('Proper Motion - Declination')
plt.title('Proper Motion of Other Cluster')

xVals = []
yVals = []

MxVals = []
MyVals = []

for i in range(1,len(Other)):
    x = Other[i][bprp]
    y = Other[i][phot]
    xVals.append(x)
    yVals.append(y)

for i in range(1,len(Othervals)):
    x = Othervals[i][bprp]
    y = Othervals[i][phot]
    MxVals.append(x)
    MyVals.append(y)

#~~~HR DIAGRAM FOR Other~~~
plt.subplot(122)
plt.scatter(xVals,yVals,s=0.05,c='b',marker='o',label='All Data') #plot all the data collected
plt.scatter(MxVals,MyVals,s=0.5,c='r',marker='o',label='Located Cluster') #overplot of the actual cluster
plt.axis([0,3.5,21,10])
plt.legend(loc='upper right', markerscale=10.0)
plt.legend.framealpha = 0.8
plt.xlabel('Colour in B-R Filters')
plt.ylabel('Brightness in G Magnitude')
plt.title('Colour-Magnitude Diagram for Other Cluster')

#######################################################################################################

#Was used to plot the RA vs DEC of the data to help find the correct stars

xVals = []
yVals = []

MxVals = []
MyVals = []
for i in range(1,len(Other)):
    x = Other[i][ra]
    y = Other[i][dec]
    xVals.append(x)
    yVals.append(y)

for i in range(1,len(Othervals)):
    x = Othervals[i][ra]
    y = Othervals[i][dec]
    MxVals.append(x)
    MyVals.append(y)

#Other cluster plot
f6 = plt.figure(6)
plt.scatter(xVals,yVals,s=0.01,c='b',marker='o')
plt.scatter(MxVals,MyVals,s=0.1,c='r',marker='o')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('RA vs DEC for Other Cluster')

plt.show()




################################################################
################################################################
################################################################
#       Plots RA pm vs DEC pm vs Tangential for NGC 6397
#
#       FOR SOME REASON THIS MUST BE IN A DIFFERENT FILE
#                      TO RUN PROPERLY
#
################################################################
################################################################
################################################################
'''
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import csv
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D

#define variables for easier indexing
src = 0 #source id
ra = 1 #right ascension
raE = 2 #right ascension error
dec = 3 #declination
decE = 4 #declination error
par = 5 #parallax
parE = 6 #parallax error
pmA = 7 #proper motion in right ascension
pmAE = 8 #error in proper motion, RA
pmD = 9 #proper motion in declination
pmDE = 10 #error in peroper motion, DEC
phot = 11 #brightness in G magnitude
bprp = 12 #colour in BR filters
rV = 13 #radial velocity
rVE = 14 #radial velocity error



fileNGC = 'E:\\Documents\\School\\3Y03\\1611340488926O-result.csv' #NGC 6397
NGC = np.genfromtxt(fileNGC, dtype=('float'), delimiter=',') #break up the file into each row

xVals = []
yVals = []
NGCvals = []
Nxvals = []
Nyvals = []

rad = 0.5 #radius of the cluster on plot, started with 2 and determined experimentally that 1.4 was a good fit for the data
avgPA = 0.0


for i in range(1,len(NGC)):
    x = NGC[i][pmA]
    y = NGC[i][pmD]

    if sqrt( (x-3.24)**2 + (y+17.58)**2 ) <= rad: #selects stars within the cluster, using the 'centre' and approximate radius
        if sqrt( (NGC[i][ra]-265.2)**2 + (NGC[i][dec]+53.67)**2 ) <= 0.5: #using 0.5 degree radius as suggested by (https://ui.adsabs.harvard.edu/abs/1997A%26AS..122....1K/abstract)
            NGCvals.append(NGC[i])

    xVals.append(x)
    yVals.append(y)

Vt = []
decVals = []
raVals = []


for i in range(1,len(NGCvals)):
    a = NGCvals[i][pmA]
    d = NGCvals[i][pmD]
    p = NGCvals[i][par]
    ea = NGCvals[i][pmAE]
    ed = NGCvals[i][pmDE]
    ep = NGCvals[i][parE]

    u = sqrt((a/1000)**2 + (d/1000)**2)
    Vtan = 4.74*u/(p/1000)
    if abs(Vtan) > 1500:
        continue
    Vt.append(Vtan)
    decVals.append(NGCvals[i][pmD])
    raVals.append(NGCvals[i][pmA])

Vtan = np.sum(Vt)/len(Vt)
        

f2 = plt.figure(2)
ax = Axes3D(f2)
ax.scatter(raVals,decVals,Vt)
ax.set_xlabel('RA Proper Motion')
ax.set_ylabel('DEC Proper Motion')
ax.set_zlabel('Tangential Velocity')
plt.title('RA Proper Motion vs DEC Proper Motion vs Tangential Velocity - NGC 6397')

plt.show()
'''

################################################################
################################################################
################################################################
#         Plots RA pm vs DEC pm vs Tangential for M35
#
#       FOR SOME REASON THIS MUST BE IN A DIFFERENT FILE
#                      TO RUN PROPERLY
#
################################################################
################################################################
################################################################

'''
import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import csv
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D

#define variables for easier indexing
src = 0 #source id
ra = 1 #right ascension
raE = 2 #right ascension error
dec = 3 #declination
decE = 4 #declination error
par = 5 #parallax
parE = 6 #parallax error
pmA = 7 #proper motion in right ascension
pmAE = 8 #error in proper motion, RA
pmD = 9 #proper motion in declination
pmDE = 10 #error in peroper motion, DEC
phot = 11 #brightness in G magnitude
bprp = 12 #colour in BR filters
rV = 13 #radial velocity
rVE = 14 #radial velocity error



fileM35 = 'E:\\Documents\\School\\3Y03\\1611340333922O-result.csv' #NGC
M35 = np.genfromtxt(fileM35, dtype=('float'), delimiter=',') #break up the file into each row

xVals = []
yVals = []
M35vals = []
MxVals = []
MyVals = []

rad = 0.5 #radius of the cluster on plot, started with 2 and determined experimentally that 1.4 was a good fit for the data
avgPA = 0.0


for i in range(1,len(M35)):
    x = M35[i][pmA]
    y = M35[i][pmD]

    if sqrt( (x-2.31)**2 + (y+2.99)**2 ) <= rad: #selects stars within the cluster, using the 'centre' and approximate radius
        if sqrt( (M35[i][ra]-92.25)**2 + (M35[i][dec]-24.33)**2 ) <= 0.5: #restricting to stars within a 0.5 degree radius following research paper (https://ui.adsabs.harvard.edu/abs/2009ApJ...695..679M/abstract)
            M35vals.append(M35[i])

    xVals.append(x)
    yVals.append(y)

Vt = []
decVals = []
raVals = []


for i in range(1,len(M35vals)):
    a = M35vals[i][pmA]
    d = M35vals[i][pmD]
    p = M35vals[i][par]
    ea = M35vals[i][pmAE]
    ed = M35vals[i][pmDE]
    ep = M35vals[i][parE]

    u = sqrt((a/1000)**2 + (d/1000)**2)
    Vtan = 4.74*u/(p/1000)
    if abs(Vtan) > 150:
        continue
    Vt.append(Vtan)
    decVals.append(M35vals[i][pmD])
    raVals.append(M35vals[i][pmA])

Vtan = np.sum(Vt)/len(Vt)

#Calcualte Radial Velocity Average
radV = []
radVT = []
radVals = []

for i in range(1,len(M35vals)):
    if np.isnan(M35vals[i][rV]):
        continue
    if M35vals[i][rVE] <= 6:
        radV.append(M35vals[i])
        radVT.append(Vt[i])
        

f2 = plt.figure(2)
ax = Axes3D(f2)
ax.scatter(raVals,decVals,Vt)
ax.set_xlabel('RA Proper Motion')
ax.set_ylabel('DEC Proper Motion')
ax.set_zlabel('Tangential Velocity')
plt.title('RA Proper Motion vs DEC Proper Motion vs Tangential Velocity - M35')

plt.show()
'''
