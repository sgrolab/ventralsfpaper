# Mark Aronson
# Sgro Lab
# Model code for "Ventral Stress Fibers Induce Plasma Membrane Deformation in Human Dermal Fibroblasts"


#%% membrane shape models

from matplotlib import pyplot as plt
import numpy as np
from scipy.misc import derivative

# inputs an x value and outputs the y value modeling the height of the membrane 
# based on the radius of the fiber and how far below the resting height of the 
# membrane the centerpoint is 
def f(x):
    
    # centerpoint above resting part of membrane - regime 1
    if cy >= 0:
        if abs(x) > r:
            return 0
        elif -1*np.sqrt(r**2-x**2) + cy > 0:
            return 0
        else:
            return -1*np.sqrt(r**2-x**2) + cy
    
    # centerpoint below resting level of membrane - regime 2 
    else: 
        if x < -(r+rx):
            return 0
        elif x < -r:
            return np.sqrt(cy**2-cy**2/(rx**2)*(x+(r+rx))**2)+cy
        elif x > r + rx:
            return 0
        elif x > r:
            return np.sqrt(cy**2-cy**2/(rx**2)*(x-(r+rx))**2)+cy
        else:
            return -1*np.sqrt(r**2-x**2) + cy

# calculates amount of curvature around a given x value based on the 
# surrounding shape of the membrane as defeined by f(x) 
def curvature(x):
    return abs(derivative(f, x, n=2, dx=1e-3))*(1 + derivative(f, x, n=1, dx=1e-3)**2)**-1.5

# set length window for membrane bending [#nm]
xmin = -2000
xmax = 2000


n = 1000        # number of steps for considering fiber curvature 
W = 2000        # nm
L = 45000       # length of fiber [nm]
r = 500         # radius of fiber [nm]
cy_range = np.linspace(-r,500,100)      # array of heights of fiber above/below resting point of membrane
K = np.zeros(len(cy_range))
Dx = W/n        # nm
Kb = 15         # standard bending energy of a membrane [kBT]

# define arrays to hold values of bending energy 
Gbend = np.zeros(len(cy_range))
Gbend_smooth = np.zeros_like(Gbend) # for smoothed curvature plots 

# iterate through heights of centerpoint of fiber above/below membrane resting height 
for j in range(0,len(cy_range)):

    cy = cy_range[j]
    rx = 100            # define radius of curvature for hyperbola of fibers when centerpoint is below membrane 
    
    # calculate membrane height profile using defined function f(x) 
    x = np.linspace(xmin,xmax,n)
    y = np.zeros(len(x))
    for i in range(0,len(x)):
        y[i] = f(x[i])
    
    # calculate curvature value at each point along membrane profile 
    curve_y = [curvature(t) for t in x]
    curve_y_smooth = np.zeros_like(curve_y)
    
    # smooth curvature values with a window size of 10 
    window = 10
    for i in range(0,len(curve_y)):
        if i < window:
            curve_y_smooth[i] = np.mean([curve_y[i:i+2*window]])
        elif i > len(curve_y)-window:
            curve_y_smooth[i] = np.mean([curve_y[i-2*window:i]])
        else:
            curve_y_smooth[i] = np.mean([curve_y[i-window:i+window]])
    
    # create arrays to store bending energy along whole membrane 
    Ebend = np.zeros(len(curve_y))
    Ebend_smooth = np.zeros_like(Ebend)
    
    # calculate bending energy at each x value given amount of curvature and membrane properties 
    for i in range(0,len(Ebend)):
        Ebend[i] = Kb/2*L*Dx*curve_y[i]**2
        Ebend_smooth[i] = Kb/2*L*Dx*curve_y_smooth[i]**2
    
    # total membrane bending energy is sum of bending energy at each point 
    Gbend[j] = np.sum(Ebend)
    Gbend_smooth[j] = np.sum(Ebend_smooth)
    
    # store for sum of curvature 
    K[j] = np.sum(curve_y)


# plot values 
plt.plot(cy_range,Gbend_smooth,color='blue',linewidth=3)
plt.xlabel('Position of Cylinder Center (nm)',fontsize=25)
plt.ylabel('Free Energy of Bending (kBT)',fontsize=25)




#%% Contraction energy

import numpy as np
import math
import matplotlib.pyplot as plt


Rfiber = np.linspace(50,250,6)              # radius of fiber in nm
contractMax = 15000                         # maximum distance of contraction in nm
contractionDist = np.linspace(0,contractMax)

fig, ax1 = plt.subplots()

# iterate through range of fiber radii 
for i in range(0,len(Rfiber)):
    nATP = np.zeros_like(contractionDist)
    
    # iterate through range of contraction distances
    for j in range(0,len(contractionDist)):
        
        # calculate cross-sectional area of the stress fiber
        CSAfiber = math.pi * Rfiber[i]**2   # [nm^2]
        
        # use literature value for single polymer radius
        Rpoly = 8                           # [nm] source: PBoC
        
        # calculate cross sectional area of single polymer
        CSApoly = math.pi * Rpoly**2        # [nm^2]
        
        # calculate number of polymers in given fiber radius
        Npoly = CSAfiber / CSApoly
        
        # get step size of myosin motor
        stepSize = 5                        # [nm] Hundt et al. 2016
        
        # calculate number of steps for fiber contraction
        nSteps = Npoly * (contractionDist[j] / stepSize)
        
        ATPperStep = 1  # Riet et al. 2000
        
        # calculate amount of ATP needed for given fiber contraction 
        nATP[j] = nSteps * ATPperStep
    
    # plot relationship between fiber contraction length and ATP needed for given fiber radii    
    ax1.plot(contractionDist,nATP,color=[145/255,35/255,57/255],linewidth=5,alpha=.5+i*.1)

# plotting settings 
plt.legend(['$r_{fiber}$ = 50nm','$r_{fiber}$ = 100nm','$r_{fiber}$ = 150nm','$r_{fiber}$ = 200nm','$r_{fiber}$ = 250nm'],frameon=0,fontsize=20)
ax1.set_xlim([0,contractMax])
ax1.set_ylim([0,np.max(nATP)])
for axiss in ['x','y']:
    ax1.tick_params(axis=axiss,labelsize=20,width=3,length=10)
for axis in ['left','bottom','right','top']:
    ax1.spines[axis].set_linewidth(2)
for axis in ['top']:
    ax1.spines[axis].set_linewidth(0)
ax1.set_xlabel('Contraction Distance (nm)',fontsize=30)
ax1.set_ylabel('ATP required',fontsize=30)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Free Energy [$k_BT$]',fontsize=30)  # we already handled the x-label with ax1
ATPtokBT = 20
kBT = nATP / ATPtokBT
ax2.plot(contractionDist, kBT,linewidth=0)
ax2.set_ylim([0,np.max(kBT)])

for axiss in ['x','y']:
    ax2.tick_params(axis=axiss,labelsize=20,width=3,length=10)
for axis in ['left','bottom','right','top']:
    ax2.spines[axis].set_linewidth(2)


fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()