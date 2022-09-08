import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

def show_results(p):
    T = p.res

    Tmin = np.min(T)
    Tmax = np.max(T)





    X, Y = np.meshgrid(p.y, p.x) # attention à l'inversion

    fig = plt.figure()

    axtime = plt.axes([0.1, 0.01, 0.8, 0.03])
    time_slider = Slider(
        ax=axtime,
        label='Time',
        valstep = np.arange(0,T.shape[0]),
        valmin=0,
        valmax=T.shape[0]-1,
        valinit=0,
    )



    ax = plt.axes(projection='3d')

    def update(val):
        ax.cla()
        ax.set_xlabel('y [m]')
        ax.set_ylabel('x [m]')
        ax.set_zlabel(f'T @ t={int(val)*p.dtLog:.1f} [s]')
        ax.set_zlim3d(Tmin, Tmax)
        ax.plot_surface(X, Y, T[int(val)], rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        fig.canvas.draw_idle()

    time_slider.on_changed(update)

    update(0)
    plt.show()