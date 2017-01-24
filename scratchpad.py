import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


fig = plt.figure(figsize=(7, 7))
ax = fig.add_axes([0, 0, 20, 20], frameon=False)
ax.set_xlim(0, 20), ax.set_xticks([])
ax.set_ylim(0, 20)
n_drops = 2
rain_drops = np.zeros(n_drops, dtype=(float, 2))
rain_drops = np.random.uniform(0, 1, (n_drops, 2))
scat = ax.scatter(rain_drops[:, 0], rain_drops[:, 1])
print(type(rain_drops[:,0]))
print(rain_drops[:,1])
def update(frame_number):
    dat = []
    for _ in range(0,1):
      rain_drops = np.float16(np.random.uniform(0, 1, 2))
      dat.append(rain_drops)
    print(frame_number)
    # print("--")
    # scat.set_offsets(dat)


animation = FuncAnimation(fig, update, interval=500)
plt.show()