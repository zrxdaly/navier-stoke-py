#%%
import numpy as np
import matplotlib.pyplot as plt
import time, sys
from matplotlib.animation import FuncAnimation

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)

L = 2
nx = 41
dx = L/(nx-1)
nt = 25
dt = 0.025
c = 1

# %% 1D diffusion
u = np.ones(nx)
u[int(.5 / dx):int(1 / dx + 1)] = 2
u_store = u.copy()

x = np.linspace(0, 2, nx)
nu = 0.3
sigma = .2
dt = sigma * dx **2 / nu

for n in range(nt):
    un = u.copy()
    # for i in range(1, nx-1): ## you can try commenting this line and...
    #     u[i] = un[i] + nu * dt / dx / dx * (un[i+1] - 2 * un[i] + un[i-1])

    u[1:nx-1] = un[1:nx-1] + nu * dt / (dx ** 2) * (un[2:nx] - 2 * un[1:nx-1] + un[0:nx-2])
    u_store = np.vstack((u_store, u))

fig = plt.figure()
ax = plt.axes(xlim = (0, 2), ylim = (0.8, 2.2))
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return(line,)

def animate(i):
    # x = np.linspace(0, 2, 1000)
    y = u_store[i,:]
    line.set_data(x, y)
    return(line,)

anim = FuncAnimation(fig, animate, init_func=init,
                               frames=26, interval=1, blit=True)
anim.save('1D_Diff.mp4', dpi=300, fps=4, extra_args=['-vcodec', 'libx264'])
plt.show()


# %%
