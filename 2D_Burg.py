#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import time, sys
from matplotlib.animation import FuncAnimation
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
nt = 150
sigma = .001
nu = .05
dt = sigma * dx * dy / nu

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)
X, Y = np.meshgrid(x, y) 

# %% 1D linear convection case
# %%timeit
u = np.ones((ny, nx))
v = np.ones((ny, nx))
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
v[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
u_store = u.copy()
v_store = v.copy()

for n in range(nt + 1):
    un = u.copy()
    vn = v.copy()
    row, col = u.shape 
    u[1:row-1,1:col-1] = un[1:row-1,1:col-1] - dt / dx * un[1:row-1,1:col-1] * (un[1:row-1,1:col-1] - un[1:row-1,0:col-2]) \
                            - dt / dy * vn[1:row-1,1:col-1]  * (un[1:row-1,1:col-1] - un[0:row-2,1:col-1]) \
                            + nu * dt * ((un[1:row-1,2:col] - 2 * un[1:row-1,1:col-1] + un[1:row-1,0:col-2]) / dx ** 2 \
                            + (un[2:row,1:col-1] - 2 * un[1:row-1,1:col-1] + un[0:row-2,1:col-1]) / dy ** 2)
    v[1:row-1,1:col-1] = vn[1:row-1,1:col-1] - dt / dx * un[1:row-1,1:col-1] * (vn[1:row-1,1:col-1] - vn[1:row-1,0:col-2]) \
                            - dt / dy * vn[1:row-1,1:col-1]  * (vn[1:row-1,1:col-1] - vn[0:row-2,1:col-1]) \
                            + nu * dt * ((vn[1:row-1,2:col] - 2 * vn[1:row-1,1:col-1] + vn[1:row-1,0:col-2]) / dx ** 2 \
                            + (vn[2:row,1:col-1] - 2 * vn[1:row-1,1:col-1] + vn[0:row-2,1:col-1]) / dy ** 2)     
                                                    
    u[[0,-1],:] = 1
    u[:,[0,-1]] = 1
    v[[0,-1],:] = 1
    v[:,[0,-1]] = 1

    u_store = np.dstack((u_store, u))
    v_store = np.dstack((v_store, u))

#%%
Fun2D(u_store, "2D_Burg_U", x, y, nt, dt)
Fun2D(v_store, "2D_Burg_V", x, y, nt, dt)
