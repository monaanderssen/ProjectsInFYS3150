from matplotlib import pyplot as plt
import numpy as np
import sys
import pylab
from matplotlib.ticker import FormatStrFormatter

input = sys.argv[1]
file = open(input, "r")

#input2 = sys.argv[2]
#file2 = open(input2, "r")

numberOfIterations  = int(file.readline())
dim = int(file.readline()) 
numberOfPlanets = int(file.readline())



print numberOfPlanets

X = np.zeros(((numberOfIterations)/1000, numberOfPlanets))
Y = np.zeros(((numberOfIterations)/1000, numberOfPlanets))
Z = np.zeros(((numberOfIterations)/1000, numberOfPlanets))
"""
X = np.zeros((numberOfIterations+1, numberOfPlanets))
Y = np.zeros((numberOfIterations+1, numberOfPlanets))
Z = np.zeros((numberOfIterations+1, numberOfPlanets))
"""
t = 0
p = 0

for line in file: 
    line_ = line.split()
    x_i = float(line_[0])                
    y_i = float(line_[1])
    z_i = float(line_[2])
    
    X[t,p] = x_i
    Y[t,p] = y_i
    Z[t,p] = z_i 

    p += 1
    
    if(p == numberOfPlanets):
        t += 1
        p = 0

#for i in range(0,numberOfPlanets):
#	X[:,i] -= X[:,0]
#	Y[:,i] -= Y[:,0] 

#print X[:,0], Y[:,0]

"""
numberOfIterations_2  = int(file2.readline())
dim_2 = int(file2.readline()) 
numberOfPlanets_2 = int(file2.readline())

X2 = np.zeros((500, numberOfPlanets-1))
Y2 = np.zeros((500, numberOfPlanets-1))
Z2 = np.zeros((500, numberOfPlanets-1))
"""
"""
t2 = 0
p2 = 0


for line in file2: 
    line_ = line.split()
    x_i = float(line_[0])                
    y_i = float(line_[1])
    z_i = float(line_[2])
   
    X2[t2,p2] = x_i
    Y2[t2,p2] = y_i
    Z2[t2,p2] = z_i 

    p2 += 1
    
    if(p2 == numberOfPlanets-1):
        t2 += 1
        p2 = 0

"""

#print X2
fig, ax = plt.subplots()


plt.plot(X[:,0], Y[:,0], '-', color='black')
#plt.plot(X[:,6], Y[:,6], 'grey')
#plt.plot(X[:,4], Y[:,4], 'brown')
#plt.plot(X[:,1], Y[:,1], 'g')
#plt.plot(X2[:,1], Y2[:,1], 'r')
#plt.plot(X[:,3], Y[:,3], 'firebrick')
#plt.plot(X[:,2], Y[:,2], 'b')
#plt.plot(X[:,5], Y[:,5])
#plt.plot(X[:,7], Y[:,7], 'coral')
#plt.plot(X[:,8], Y[:,8], 'dodgerblue')
#plt.plot(X[:,9], Y[:,9], 'chocolate')

#plt.plot(X[:,0],t)

plt.legend(['Orbit of Earth with Jupiter', 'Orbit of Earth without Jupiter'])
#sun = plt.Circle((X[0,0],Y[0,0]), 0.005, color='black')
#orbit = plt.Circle((0,0), 1, color='k', fill=False)
plt.title('Position of Earth: n=%s' % numberOfIterations, fontsize = 22)
plt.title("Sun around center of mass", fontsize = 22)
plt.xlabel('x [AU]', fontsize = 22)
#plt.axis([-0.00005-7.2277e-1, -0.00001-(7.2277e-1), 0.000002+(6.71438e-1), 0.000010+(6.71438e-1)])
#plt.axis('equal')
ax.ticklabel_format(useOffset=False)

plt.ylabel('y [AU]', fontsize = 22)# change to v_tilde(x)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.legend(['Sun'])
#plt.legend(['Sun'])
#plt.legend(['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto'], fontsize = 16) 
#ax.add_artist(sun)
#ax.add_artist(orbit)
plt.show()


