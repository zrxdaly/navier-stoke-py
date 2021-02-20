#%%
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import time, sys

from sympy import init_printing
init_printing(use_latex=True)
from sympy.utilities.lambdify import lambdify

nx = 101
dx = 2 * np.pi /(nx - 1)
nt = 100
nu0 = .07
dt = dx * nu0

# %% initial condition
x, nu, t = sp.symbols('x nu t')
phi = (sp.exp(-(x - 4.0 * t)**2 / (4.0 * nu * (t + 1))) +
       sp.exp(-(x - 4.0 * t - 2.0 * sp.pi)**2 / (4.0 * nu * (t + 1))))

phiprime = phi.diff(x)

u = -2 * nu * (phiprime / phi) + 4

ufunc = lambdify((t, x, nu), u)
# %% declarations of U
x = np.linspace(0, 2 * np.pi, nx)
un = np.empty(nx)
t = 0.0

# what happened here is strange
u = np.asarray([ufunc(t, x0, nu0) for x0 in x])

u_store = u.copy()
u_ana_store = u.copy()
# %% 
for n in range(nt):
    un = u.copy()
    # for i in range(1, nx-1):
    #     u[i] = un[i] - un[i] * dt / dx *(un[i] - un[i-1]) + nu0 * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])
    u[1:nx-1] = un[1:nx-1] - un[1:nx-1] * dt / dx *(un[1:nx-1] - un[0:nx-2]) + nu0 * dt / dx**2 * (un[2:nx] - 2 * un[1:nx-1] + un[0:nx-2])
    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu0 * dt / dx**2 * (un[1] - 2 * un[0] + un[-2])
    # the left boundary condition is calculated using
    u[-1] = u[0]

    u_analytical = np.asarray([ufunc(n * dt, xi, nu0) for xi in x])

    u_store = np.vstack((u_store, u))
    u_ana_store = np.vstack((u_ana_store, u_analytical))

# %%
fig = plt.figure()
ax = plt.axes(xlim = (0, 2 * np.pi), ylim = (0, 10))
line, = ax.plot([], [], marker='o', lw=2, label='Computational')
line2, = ax.plot([], [], label='Analytical')
plt.legend()

def init():
    line.set_data([], [])
    line2.set_data([], [])
    return(line, line2,)

def animate(i):
    y = u_store[i,:]
    y2 = u_ana_store[i,:]
    line.set_data(x, y)
    line2.set_data(x, y2)
    return(line,line2,)

anim = FuncAnimation(fig, animate, init_func=init, frames=26, interval=1, blit=True)
anim.save('1D_Burg.mp4', dpi=300, fps=4, extra_args=['-vcodec', 'libx264'])
plt.show()


# %%
