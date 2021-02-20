#%%
import numpy as np
import matplotlib.pyplot as plt
from animateFun import quiver2D

Lx = 2
Ly = 2
nx = 101
ny = 101
dx = Lx / (nx-1)
dy = Ly / (ny-1)
nt = 2000

x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)

rho = 1
nu = .01
dt = .001

u = np.zeros((ny, nx))
un = u.copy()
u_store = u.copy()
v = np.zeros((ny, nx))
vn = v.copy()
v_store = v.copy()
p = np.zeros((ny, nx))
pn = np.zeros((ny, nx))
p_store = p.copy()

b  = np.zeros((ny, nx))
eta = 1e-4

nit = 50
#%% solver
def build_b(b, rho, dt, u, v, dx, dy):
    b[1:-1, 1:-1] = rho * (((u[1:-1,2:] -u[1:-1,0:-2]) / (2 * dx) + (v[2:,1:-1] - v[0:-2,1:-1]) / (2 * dy)) / dt - \
                    ((u[1:-1,2:] - u[1:-1,0:-2]) / (2 * dx)) ** 2 - \
                    2 * (v[1:-1,2:] - v[1:-1,0:-2]) / (2 * dx) * (u[2:,1:-1] - u[0:-2,1:-1]) / (2 * dy) - \
                    ((v[2:,1:-1] - v[0:-2,1:-1]) / (2 * dy)) ** 2)
    return(b)


def poisson_pre(p, dx, dy, b):
    error = 1
    pn = p.copy()
    while error > eta:
        pn = p.copy()
        p[1:-1,1:-1] = (((pn[1:-1, 2:] + pn[1:-1, :-2]) * dy**2 +
                        (pn[2:, 1:-1] + pn[:-2, 1:-1]) * dx**2 -
                        b[1:-1, 1:-1] * dx**2 * dy**2) / 
                        (2 * (dx**2 + dy**2)))

        p[0,:] = p[1,:] 
        p[:,0] = p[:,1]
        p[:,-1] = p[:,-2]
        p[-1,:] = 0
        error = (np.sum(np.abs(p[:]) - np.abs(pn[:]))) / np.sum(np.abs(pn[:]))
    return(p)

def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, p_store, u_store, v_store):
    un = np.empty_like(u)
    vn = np.empty_like(v)
    b = np.zeros((ny, nx))
    p_store = p.copy()
    u_store = u.copy()
    v_store = v.copy()
    for ii in range(nt):
        un = u.copy()
        vn = v.copy()
        b = build_b(b, rho, dt, u, v, dx, dy)
        p = poisson_pre(p, dx, dy, b)
        u[1:-1,1:-1] = un[1:-1,1:-1] - dt / dx * un[1:-1,1:-1] * (un[1:-1,1:-1] - un[1:-1,0:-2]) - \
                        dt / dy * vn[1:-1,1:-1] * (un[1:-1,1:-1] - un[0:-2,1:-1]) - \
                        dt / rho / dx / 2 * (p[1:-1,2:] - p[1:-1,0:-2]) + \
                        dt * nu * (un[1:-1,2:] - 2 * un[1:-1,1:-1] + un[1:-1,0:-2]) / (dx ** 2) + \
                        dt * nu * (un[2:,1:-1] - 2 * un[1:-1,1:-1] + un[0:-2,1:-1]) / (dy ** 2)
        v[1:-1,1:-1] = vn[1:-1,1:-1] - dt / dx * un[1:-1,1:-1] * (vn[1:-1,1:-1] - vn[1:-1,0:-2]) - \
                        dt / dy * vn[1:-1,1:-1] * (vn[1:-1,1:-1] - vn[0:-2,1:-1]) - \
                        dt / rho / dy / 2 * (p[2:,1:-1] - p[0:-2,1:-1]) + \
                        dt * nu * (vn[1:-1,2:] - 2 * vn[1:-1,1:-1] + vn[1:-1,0:-2]) / (dx ** 2) + \
                        dt * nu * (vn[2:,1:-1] - 2 * vn[1:-1,1:-1] + vn[0:-2,1:-1]) / (dy ** 2)
        u[-1,:] = 10
        u[0,:] = 0
        u[:,[0,-1]] = 0
        v[[0,-1],:] = 0
        v[:,[0,-1]] = 0
        p_store = np.dstack((p_store, p))
        u_store = np.dstack((u_store, u))
        v_store = np.dstack((v_store, v))

    return(p_store, u_store, v_store)

p_store, u_store, v_store = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, p_store, u_store, v_store)

#%%
quiver2D(p_store, u_store, v_store, "cavity", x, y, 100, dt, int(len(p_store[...,0])/100))
