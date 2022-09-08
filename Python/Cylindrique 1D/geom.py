import numpy as np

def geometry(p):
    # J'ajoute ici dans l'objet p les données qui vont m'être utiles

    # ------------------------------------
    # Distances...
    # ------------------------------------

    p.dx = dx = p.Lx/p.M    # notez ici le chainage des affectations
    p.dy = dy = p.Ly/p.N

    # j'en profite pour créer un vecteur contenant les x(m) et y(n) qui seront (peut-être) utiles par la suite
    p.x = dx*(np.arange(0,p.M)+0.5)
    p.y = dy * (np.arange(0, p.N) + 0.5)

    # ------------------------------------
    # Surfaces...
    # ------------------------------------
    p.Sx = dy*1
    p.Sy = dx * 1

    # ------------------------------------
    # Volume(s)...
    # ------------------------------------
    p.V = dx*dy*1

    return p


def solve():
    pass

def show_results():
    pass