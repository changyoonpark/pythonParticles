import numpy as np
from matplotlib import pyplot as plt

# Differential equation
# diff = y'= y/x (or say x+y)
def diff(x,y):
    return x/y # try also x+y

x = np.linspace(-10,10,50)
y = np.linspace(-10,10,50)

# use x,y
for j in x:
    for k in y:
        slope = diff(j,k)
        domain = np.linspace(j-0.07,j+0.07,2)
        def fun(x1,y1):
            z = slope*(domain-x1)+y1
            return z
        plt.plot(domain,fun(j,k),solid_capstyle='projecting',solid_joinstyle='bevel')

plt.title("Slope field y'")
plt.grid(True)
plt.show()
    
print("End of the program")