import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

def show_results(p):
    T = p.res_T
    F = p.res_Y  # je change de variable ici simplement pour concerver Y pour la dimension d'espace.

    Tmin = np.min(T)
    Tmax = np.max(T)





    X, Y = np.meshgrid(p.y, p.x) # attention à l'inversion

    fig = plt.figure()

    axtime = plt.axes([0.1, 0.01, 0.8, 0.03])
    time_slider = Slider(
        ax=axtime,
        label='i=',
        valstep = np.arange(0,T.shape[0]),
        valmin=0,
        valmax=T.shape[0]-1,
        valinit=0,
    )



    axT = fig.add_subplot(1, 2, 1, projection='3d')
    axY = fig.add_subplot(1, 2, 2, projection='3d')

    def update(val):
        # Pour la température

    
        axT.cla()
        axT.set_xlabel('y [m]')
        axT.set_ylabel('x [m]')
        axT.set_zlabel('T [°C]')
        axT.set_zlim3d(Tmin, Tmax)
        axT.plot_surface(X, Y, T[int(val)], rstride=1, cstride=1, cmap='viridis', edgecolor='none')

        # Pour la fraction liquide
        axY.cla()
        axY.set_xlabel('y [m]')
        axY.set_ylabel('x [m]')
        axY.set_zlabel('Y [-]')
        axY.set_zlim3d(0, 1)
        axY.plot_surface(X, Y, F[int(val)], rstride=1, cstride=1, cmap='viridis', edgecolor='none')

        fig.canvas.draw_idle()

    time_slider.on_changed(update)

    update(0)
    plt.show()