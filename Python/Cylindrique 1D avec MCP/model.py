import time as time_pkg  # ici je renomme le package pour pouvoir utiliser ma variable time

import numpy as np


# C'est ici qu'on va implémenter le modèle

def solve(p):
    # Déclaration des tableaux
    T = np.zeros((p.M, p.N))  # Le champs de température (cohérent) à l'état i
    Y = np.zeros((p.M, p.N))  # Le champs de faction liquide (cohérent) à l'état i
    E = np.zeros((p.M, p.N))  # Le champs d'enthalpie massique (spécifique) (cohérent) à l'état i

    EVOL = np.zeros((p.M, p.N))  # L'évolution d'ENERGIE (pour le champs) au cours de l'itération i -> i+1

    # Tableaux pour sauvegarde (T et Y)
    # on a une durée D avec une sauvegarde tout les dtlog...
    # on va reserver un tableau avec int(D/dtlog)+1 valeurs
    res_T = np.zeros((max(1, int(p.D / p.dtLog)), p.M, p.N))
    res_Y = np.zeros((max(1, int(p.D / p.dtLog)), p.M, p.N))
    # et on se garde un compteur pour savoir où sauvegarder
    iSave = 0

    # Condition initiale
    T[:, :] = p.T0
    if p.T0 < p.TF:  # SOLIDE
        Y[:, :] = 0.0
        E[:, :] = p.CS * (p.T0 - p.TF)
    elif p.T0 > p.TF:  # Liquide
        Y[:, :] = 1.0
        E[:, :] = p.LF + p.CL * (p.T0 - p.TF)
    else:
        print("Initialisation à T0=TF impossible")
        exit()

    time = 0  # le temps physique
    nextlog = 0  # l'instant de la prochaine sauvegarde du champs

    # j'en profiction pour vous donner la syntaxe d'une "expression lambda"
    # qui est une sorte de fonction réduite à sa plus simple expression

    K = lambda Y1, Y2: p.KS + 0.5 * (Y1 + Y2) * (p.KL - p.KS)
    C = lambda Y: p.CS + Y * (p.CL - p.CS)

    def calc_TY(e):
        # retourne le couple T,Y correspondant à l'énergie spécifique e (via l'équation d'état, cf. cours)
        if e < 0:  # Solide
            return p.TF + e / p.CS, 0
        if e > p.LF:  # Liquide
            return p.TF + (e - p.LF) / p.CL, 1
        else:  # en cours de transformation
            return p.TF, e / p.LF

    start_time = time_pkg.time()  # pour calculer la durée d'exécution
    while time < p.D:  # Boucle temporelle

        # ------------------------------------------
        # gestion de la sauvegarde
        # ------------------------------------------
        if time >= nextlog:
            res_T[iSave, :, :] = T[:, :]
            res_Y[iSave, :, :] = Y[:, :]
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
                    FG = K(Y[m - 1, n], Y[m, n]) * p.Sx * (T[m - 1, n] - T[m, n]) / p.dx

                # Flux droit
                if m == p.M - 1:
                    tmp = 1 / (1 / p.H1 + 0.5 * p.dx / K(Y[m, n], Y[m, n]))
                    FD = tmp * p.Sx * (p.TINF - T[m, n])
                else:
                    FD = K(Y[m + 1, n], Y[m, n]) * p.Sx * (T[m + 1, n] - T[m, n]) / p.dx

                # Flux haut
                if n == p.N - 1:
                    tmp = 1 / (1 / p.H2 + 0.5 * p.dx / K(Y[m, n], Y[m, n]))
                    FH = tmp * p.Sy * (p.TINF - T[m, n])
                    if m > p.M / 2 - 1:
                        # on ajoute le flux radiatif
                        FH += p.Sy * p.PHI
                else:
                    FH = K(Y[m, n + 1], Y[m, n]) * p.Sy * (T[m, n + 1] - T[m, n]) / p.dy

                # Flux bas
                if n == 0:
                    FB = 2 * K(Y[m, n], Y[m, n]) * p.Sy * (p.TB - T[m, n]) / p.dy
                else:
                    FB = K(Y[m, n - 1], Y[m, n]) * p.Sy * (T[m, n - 1] - T[m, n]) / p.dy

                # bilan
                EVOL[m, n] = p.dt * (FG + FD + FH + FB) / p.rho / p.V  # attention, pas de C au dénominateur
        # ------------------------------------------
        # Préparation de l'itération suivante
        # ------------------------------------------

        E[:, :] = E[:, :] + EVOL[:, :]

        # ON calcule les nouveaux champs T et Y
        for m in range(0, p.M):
            for n in range(0, p.N):
                T[m, n], Y[m, n] = calc_TY(E[m, n])

        time += p.dt

    p.elapsed = time_pkg.time() - start_time
    p.res_T = res_T
    p.res_Y = res_Y

    # on oublie par de retourner le dictionnaire paramètre
    return p
