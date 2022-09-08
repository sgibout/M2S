import numpy as np
import time as time_pkg # ici je renomme le package pour pouvoir utiliser ma variable time


# C'est ici qu'on va implémenter le modèle

def solve(p):
    # Déclaration des tableaux
    T = np.zeros((p.M, p.N))  # Le champs de température (cohérent) à l'état i
    EVOL = np.zeros((p.M, p.N))  # L'évolution de température (pour le champs) au cours de l'itération i -> i+1

    # Tableau pour sauvegarde
    # on a une durée D avec une sauvegarde tout les dtlog...
    # on va reserver un tableau avec int(D/dtlog)+1 valeurs
    res = np.zeros((max(1,int(p.D/p.dtLog)),p.M, p.N))
    # et on se garde un compteur pour savoir où sauvegarder
    iSave = 0

    # Condition initiale
    T[:, :] = p.T0

    time = 0  # le temps physique
    nextlog = 0  # l'instant de la prochaine sauvegarde du champs

    # j'en profiction pour vous donner la syntaxe d'une "expression lambda"
    # qui est une sorte de fonction réduite à sa plus simple expression

    K = lambda T1, T2 : p.K0 + p.alphaK * (0.5 * (T1 + T2) - p.TREF)
    C = lambda T : p.C0 + p.alphaC * (T - p.TREF)

    start_time = time_pkg.time() # pour calculer la durée d'exécution
    while time < p.D:  # Boucle temporelle

        # ------------------------------------------
        # gestion de la sauvegarde
        # ------------------------------------------
        if time >= nextlog:
            res[iSave, :,:] = T[:,:]
            nextlog += p.dtLog
            iSave += 1

        # ------------------------------------------
        # Boucles spatiales
        # ------------------------------------------
        for m in range(0, p.M):
            for n in range(0, p.N):

                # Flux gauche
                if m == 0:
                    FG = 0
                else:
                    FG = K(T[m - 1, n], T[m, n]) * p.Sx * (T[m - 1, n] - T[m, n]) / p.dx

                # Flux droit
                if m == p.M - 1:
                    tmp = 1 / (1 / p.H1 + 0.5 * p.dx / K(T[m, n], T[m, n]))
                    FD = tmp * p.Sx * (p.TINF - T[m, n])
                else:
                    FD = K(T[m + 1, n], T[m, n]) * p.Sx * (T[m + 1, n] - T[m, n]) / p.dx

                # Flux haut
                if n == p.N - 1:
                    tmp = 1 / (1 / p.H2 + 0.5 * p.dx / K(T[m, n], T[m, n]))
                    FH = tmp * p.Sy * (p.TINF - T[m, n])
                    if m>p.M/2-1:
                        # on ajoute le flux radiatif
                        FH += p.Sy*p.PHI
                else:
                    FH = K(T[m, n+1], T[m, n]) * p.Sy * (T[m, n+1] - T[m, n]) / p.dy

                # Flux bas
                if n == 0:
                    FB = 2* K(T[m,n],T[m,n])* p.Sy * (p.TB - T[m, n])/p.dy
                else:
                    FB = K(T[m, n-1], T[m, n]) * p.Sy * (T[m, n-1] - T[m, n]) / p.dy


                #bilan
                EVOL[m,n] = p.dt*(FG+FD+FH+FB)/p.rho/C(T[m,n])/p.V
        # ------------------------------------------
        # Préparation de l'itération suivante
        # ------------------------------------------
        T[:, :] = T[:, :] + EVOL[:, :]

        time += p.dt


    p.elapsed = time_pkg.time()-start_time
    p.res = res



    # on oublie par de retourner le dictionnaire paramètre
    return p
