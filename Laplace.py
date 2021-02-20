#%%
import numpy as np
import matplotlib.pyplot as plt
from animateFun import Fun2D

Lx = 2
Ly = 2
nx = 61
ny = 61
dx = Lx / (nx-1)
dy = Ly / (ny-1)
nt = 150
dt = 0.01

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)
X, Y = np.meshgrid(x, y)
p = np.zeros((nx, ny))
p[:,0] = 0
p[:,-1] = y
p[0, :] = p[1, :]
p[-1, :] = p[-2, :]

p_store = p.copy()
eta = 1e-4
error_s = []
error = 1

while error > eta:
    pn = p.copy()
    p[1:ny-1, 1:nx-1] = (dx ** 2 * (pn[2:ny, 1:nx-1] + pn[0:ny-2, 1:nx-1]) + \
                        dy ** 2 * (pn[1:ny-1, 2:nx] + pn[1:ny-1, 0:nx-2])) / (2 * (dx**2 + dy**2))
    p[:,0] = 0
    p[:,-1] = y
    p[0, :] = p[1, :]
    p[-1, :] = p[-2, :]

    p_store = np.dstack((p_store, p))

    error = (np.sum(np.abs(p[:]) - np.abs(pn[:]))) / np.sum(np.abs(pn[:]))

    error_s.append(error)

#%%
Fun2D(p_store, "Laplace", x, y, np.linspace(0, len(error_s), 100).astype(int), dt, int(len(error_s)/100))

# %%
