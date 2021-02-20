#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import time, sys
from animateFun import Fun2D

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)

L = 2
nx = 161
ny = 161
dx = L / (nx-1)
dy = L / (ny-1)
nt = 100
sigma = .2
dt = sigma * dx
c = 1

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
    u[1:row,1:col] = un[1:row,1:col] - c * dt / dx * (u[1:row, 1:col] - u[1:row, 0:col-1]) - c * dt / dy * (u[1:row, 1:col] - u[0:row-1, 1:col])
    u[[0,-1],:] = 1
    u[:,[0,-1]] = 1

    u_store = np.dstack((u_store, u))

#%% using the ploting function to animate
Fun2D(u_store, "2D_con", x, y, nt, dt)
