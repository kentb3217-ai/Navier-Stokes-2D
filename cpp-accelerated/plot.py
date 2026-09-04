import numpy as np, matplotlib.pyplot as plt
import struct

def read_exact(file, size):
    data = file.read(size)
    return data

def read_uint32(file):
    return struct.unpack("=I", read_exact(file, 4))[0]

def read_vec(file):
    count = read_uint32(file)
    return list(struct.unpack(f"={count}d", read_exact(file, count * 8)))

def read_mat(file):
    count_row = read_uint32(file)
    return [read_vec(file) for _ in range(count_row)]

def read_vec_mat(file):
    count_mat = read_uint32(file)
    return [read_mat(file) for _ in range(count_mat)]

path = r"C:\Users\burne\OneDrive\Desktop\Coding\Projects\Navier-Stokes-2D\cpp-accelerated\data.bin"

with open(path, "rb") as file:
    ux = np.array(read_vec_mat(file))
    uy = np.array(read_vec_mat(file))
    dtTot = np.array(read_vec(file))

fig, ax = plt.subplots()
Y, X = np.meshgrid(np.arange(ux.shape[1]), np.arange(ux.shape[1]), indexing="ij")
iterations : int = 0

while iterations < dtTot.size:
    ax.clear()
    ax.set_aspect('equal')
    ax.quiver(X, Y, ux[iterations], uy[iterations], scale=1)
    ax.set_title("Distribution at t: {:.3f} [s].".format(dtTot[iterations]))
    plt.pause(0.01)

    iterations += 1
