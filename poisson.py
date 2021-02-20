#%%
import numpy as np
import matplotlib.pyplot as plt
from animateFun import Fun2D

Lx = 2
Ly = 1
nx = 61
ny = 61
dx = Lx / (nx-1)
dy = Ly / (ny-1)
nt = 100
dt = 0.01

x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)

p = np.zeros((nx, ny))
pd = np.zeros((ny, nx))
b  = np.zeros((ny, nx))
p_store = p.copy()

b[int(ny / 4), int(nx / 4)]  = 100
b[int(3 * ny / 4), int(3 * nx / 4)] = -100

eta = 1e-4
error_s = []
error = 1

while error > eta:
    pd = p.copy()
    p[1:-1,1:-1] = (((pd[1:-1, 2:] + pd[1:-1, :-2]) * dy**2 +
                    (pd[2:, 1:-1] + pd[:-2, 1:-1]) * dx**2 -
                    b[1:-1, 1:-1] * dx**2 * dy**2) / 
                    (2 * (dx**2 + dy**2)))

    p[0, :] = 0
    p[ny-1, :] = 0
    p[:, 0] = 0
    p[:, nx-1] = 0

    p_store = np.dstack((p_store, p))

    error = (np.sum(np.abs(p[:]) - np.abs(pd[:]))) / np.sum(np.abs(pd[:]))

    error_s.append(error)

# %%
Fun2D(p_store, "poisson", x, y, np.linspace(0, len(error_s), 100).astype(int), dt, int(len(error_s)/100))
