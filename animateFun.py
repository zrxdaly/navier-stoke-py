import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.animation import FuncAnimation

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='20')
plt.rc('ytick', labelsize='20')
plt.rc('text', usetex=True)

def end2D(u, case_txt, x, y, nt, dt):
    X, Y = np.meshgrid(x, y) 
    fig = plt.figure(figsize=(8, 6.4), dpi=150)
    ax = plt.axes(xlabel='x', ylabel='y')
    con = plt.contourf(X, Y, u, cmap=cm.viridis)
    plt.colorbar()
    plt.tight_layout()

def Fun2D(u_store, case_txt, x, y, nt, dt, fps):
    X, Y = np.meshgrid(x, y) 
    fig = plt.figure(figsize=(8, 6.4), dpi=300)
    ax = plt.axes(xlabel='x', ylabel='y')
    con = plt.contourf(X, Y, u_store[:,:,0], cmap=cm.viridis)
    plt.colorbar()
    plt.tight_layout()

    def animate(i):
    # global con
        # for c in con.collections:
        #     c.remove() 
        con = plt.contourf(X, Y, u_store[:,:,i], cmap=cm.viridis)
        plt.title(r'$%s$ t = %.5f sec' % (case_txt, i * dt))
        return(con)

    anim = FuncAnimation(fig, animate, frames=nt, interval=1)
    anim.save('%s.mp4'%case_txt, dpi=300, fps=10, extra_args=['-vcodec', 'libx264'])
    plt.show()


def quiver2D(p_store, u_store, v_store, case_txt, x, y, nt, dt, fps):
    X, Y = np.meshgrid(x, y) 
    fig, ax = plt.subplots(figsize=(8, 6.4), dpi=300)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    con = ax.contourf(X, Y, p_store[:,:,0], cmap=cm.viridis)
    fig.colorbar(con)
    qax = ax.quiver(x[::5], y[::5], u_store[::5,::5,0], v_store[::5,::5,0])
    fig.tight_layout()

    def animate(i):
        con = ax.contourf(X, Y, p_store[:,:,i], cmap=cm.viridis)
        qax = ax.quiver(x[::5], y[::5], u_store[::5,::5,i], v_store[::5,::5,i])
        plt.title(r'$%s$ t = %.5f sec' % (case_txt, i * dt))


    anim = FuncAnimation(fig, animate, frames=np.linspace(0, len(p_store[...,0]), nt).astype(int), interval=1)
    anim.save('%s.mp4'%case_txt, dpi=300, fps=10, extra_args=['-vcodec', 'libx264'])
    plt.show()

