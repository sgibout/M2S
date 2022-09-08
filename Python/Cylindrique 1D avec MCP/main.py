# =============================================================================
# Exemple 1D cylindrique - Python
# SG.10.2021
# =============================================================================

# Attention: par simplicité (et comme je sais ce que je fais), je laisse plusieurs variables dans la portée "globale",
# ce qui veut dire de 1) elles seront accessible partout (ça c'est bien) et 2) on peut les modifiers partout (ça c'est moins bien)
# Vous remarquerez que, comme python est un langage dynamique, je n'ai pas à proprement parler de "déclaration" de variables

from data import DataDict
from geom import geometry
from model import solve
from display import show_results


# on crée le DataDict qui portera toute les données du problème
param = DataDict()
# et on charge les données
param.load("params.txt")

if (param.M % 2) != 0:
    print('M doit être impair --> STOP')
    exit()





# Pour vérifier ce qui se passe, vous pouvez ajouter des  "print(param)" de temps en temps

# Calcul des distances, surfaces et volumes (qu'on ajoute dans le param)
# ici comme param est un objet (donc passé par référence, j'aurais pu me passer du param = ...
# mais cela permet de garder en tête que geometry MODIFIE param.
# Attention, dans ce cas, ne pas oublier le return dans geometry
param = geometry(param)

# résolution
param = solve(param)


# affichage
show_results(param)