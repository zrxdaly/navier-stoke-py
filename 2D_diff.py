#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import time, sys
from animateFun import Fun2D

L = 2
nx = 161
ny = 161
dx = L / (nx-1)
dy = L / (ny-1)
nu = .05
nt = 100
sigma = .25
dt = sigma * dx * dy / nu

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)
X, Y = np.meshgrid(x, y) 

# %% 1D linear convection case
# %%timeit
u = np.ones((ny, nx))
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
u_store = u.copy()

for n in range(nt + 1):
    un = u.copy()
    row, col = u.shape 
    u[1:row-1,1:col-1] = un[1:row-1,1:col-1] + nu * dt / dx ** 2 * (un[1:row-1,2:col] - 2 * un[1:row-1,1:col-1] + un[1:row-1,0:col-2])  \
                        + nu * dt / dy ** 2 * (un[2:row,1:col-1] - 2 * un[1:row-1,1:col-1] + un[0:row-2,1:col-1])
    # with - sgin of second term the solution will be totally different
    u[[0,-1],:] = 1
    u[:,[0,-1]] = 1

    u_store = np.dstack((u_store, u))

# %% animate function
Fun2D(u_store, "2D_diff", x, y, nt, dt)
