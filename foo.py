import numpy as np
from math import *
from pylab import *
from matplotlib import pyplot as plt
from matplotlib import animation

# Constants
isqrt = 2**(-0.5)
omega = np.sqrt(2-np.sqrt(2))   #Angular velocity
L=4                             #Length of the system

n = 1                         #Normal mode number  
if n==1:
    z = [isqrt,1,isqrt]             #mode 1
elif n==2:
    z = [1,0,-1]                   #mode 2
elif n==3:
    z = [isqrt,-1,isqrt]           #mode 3

ex = [1,2,3]                    #x-coordinates of scatter points

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, L), ylim=(-1.1, 1.1))

line, = ax.plot([], [], lw=2, color='b')
scat, = ax.plot([],[], linestyle='', marker='o', color='b')
# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    scat.set_data([], [])
    return [scat,line,]

# animation function.  This is called sequentially
def animate(t):
    xinterval = np.arange(0,10,0.05)
    wave = np.cos(0.1*omega*t)*np.sin(n*xinterval*np.pi/L)
    line.set_data(xinterval, wave)
    dots = z*real(np.exp(0+(omega*0.1*t)*1j))

    scat.set_data(ex, dots)
    return [scat,line,]

# call the animator. 
anim = animation.FuncAnimation(fig, animate,init_func=init, frames=range(200), interval=20,     blit=True)
plt.grid(True)
plt.show()