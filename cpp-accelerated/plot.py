import numpy as np, matplotlib.pyplot as plt
import struct
from pathlib import Path

fig, ax = plt.subplots()
X, Y = np.meshgrid(np.arange(20), np.arange(20))
iter : int = 0

while iter < dtTot.size:
    ax.clear()
    ax.set_aspect('equal')
    ax.quiver(X, Y, ux, uy)
    ax.set_title("Distribution at t: {:.3f} [s].".format(dtTot[iter]))
    plt.pause(0.01)
