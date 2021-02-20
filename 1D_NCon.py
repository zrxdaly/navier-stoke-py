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
# %% 1D non-linear convective  -- inviscid burger equation
u = np.ones(nx)
u[int(.5 / dx):int(1 / dx + 1)] = 2

u_store = u.copy()

x = np.linspace(0, 2, nx)

for n in range(nt):
    un = u.copy()
    # for i in range(1, nx): 
    #     u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i-1])
    u[1:nx] = un[1:nx] - un[1:nx] * dt /dx * (un[1:nx] - un[:nx-1])
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
anim.save('1D_NCon.mp4', dpi=300, fps=4, extra_args=['-vcodec', 'libx264'])
plt.show()


