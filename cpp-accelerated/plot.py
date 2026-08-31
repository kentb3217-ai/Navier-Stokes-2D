import numpy as np, matplotlib.pyplot as plt
import struct
from pathlib import Path

file_path = Path(__file__).resolve().parent / "data.bin"

with open(file_path, "rb") as file:
    count = struct.unpack("<Q", file.read(8))[0]

    ux = np.fromfile(file, dtype="<f8", count=count)
    uy = np.fromfile(file, dtype="<f8", count=count)
    dtTot = np.fromfile(file, dtype="<f8", count=count)

fig, ax = plt.subplots()
X, Y = np.meshgrid(np.arange(ux.size), np.arange(uy.size))
iter : int = 0

while iter < dtTot.size:
    ax.clear()
    ax.set_aspect('equal')
    ax.quiver(X, Y, ux[iter], uy[iter])
    ax.set_title("Distribution at t: {:.3f} [s].".format(dtTot[iter]))
    plt.pause(0.01)
