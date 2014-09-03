PROGRAM Cart2dExp

	! ***********************************************************
	! ** Exemple de code - Cours M2S - ENSGTI 2A - Energétique **
	! ** Stéphane Gibout - 2014                                **
	! **                                                       **
	! ** Cartésien 2D Instationnaire - Schéma Explicit         **
	! ** avec des paramètres non constants                     **
	! **                                                       **
	! ***********************************************************

	IMPLICIT NONE
	
	
	! ====================================================================================
	! Déclaration des variables
	! ====================================================================================

	! Par commodité, les variables contenant les paramètres du modèle et les
	! résultats sont déclarées globales. Cela permet de ne pas avoir trop de 
	! paramètres à passer aux fonctions/procédures.

	! Attention car Fortran n'est pas sensible à la casse !!
	
	! --- CONSTANTES ---------------------------------------------------------------------
	
	DOUBLE PRECISION, PARAMETER :: PI = 4*ATAN(1.D0)	! Rapport entre le périmètre d'un cercle et son diamètre (PI !) 
	INTEGER, PARAMETER :: FICH=100						! Identifiant pour le fichier
	
	! --- Données du modèles -------------------------------------------------------------
	INTEGER :: Mmax, Nmax				! Nombre de pas d'espaces 1..Mmax, 1..Nmax
	
	DOUBLE PRECISION :: D, LX,LY		! Durée [s] et dimensions [m]
	DOUBLE PRECISION :: rho				! Paramètres thermophysiques [SI]
	DOUBLE PRECISION :: TREF			! Température de référence pour k et C
	DOUBLE PRECISION :: k0, k1
	DOUBLE PRECISION :: C0, C1
	DOUBLE PRECISION :: H1, H2			! Coefficients échanges globaux [SI]
	
	DOUBLE PRECISION :: T0				! Température initiale [°C]
	DOUBLE PRECISION :: TINF,P,DTINF	! Température ambiante [°C]
	
	DOUBLE PRECISION :: dt				! Pas de temps de résolution [s]
	DOUBLE PRECISION :: dtLog			! Pas de temps pour la sauvegarde [s]
	
	INTEGER :: NQ						! Nombre de termes sources
	DOUBLE PRECISION :: Qnom			! Puissance volumique nominale
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tabQ		! Tableau de terme sources
	
	
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T	! Le champs de température
	! Déclaré global pour simplifier l'écriture de la routine d'EXPORT
	
	! ====================================================================================
	! Début du programme principal
	! ====================================================================================
	

	! --- Chargement des paramètres depuis le fichier sphere.data ------------------------ 
	CALL load()		
	PRINT*,'Chargement des données.................OK'

	! --- Construction du maillage -------------------------------------------------------
	! Ici je n'ai pas besoin d'une routine spéciale pour le maillage qui est très simple
	! Comme ces valeurs ne sont utilisées que pour le calcul, je préfère les calculer
	! directement dans la routine de résolution.
	
	! --- Calcul -------------------------------------------------------------------------
	PRINT*,'Calcul..............................DEBUT'	
	CALL run()
	PRINT*,'Calcul................................FIN'	


	! ====================================================================================
	! Fin du programme principal
	! ====================================================================================

	PRINT*,'Bye Bye'	


CONTAINS	


	! ====================================================================================
	! Début des procédures et fonctions
	! ====================================================================================


! Chargement des données depuis depuis le fichier sphere.data
SUBROUTINE load()
	IMPLICIT NONE
	
	INTEGER :: k,m,n
	DOUBLE PRECISION :: a
	
	! ------------------------------------------------------------------------------------
	! Ouverture du fichier en lecture
	! ------------------------------------------------------------------------------------

	! Je ne controle pas le succés donc le programme s'arrêtera en cas de problème
	OPEN(FICH,FILE="./cart2d.data", ACTION="READ");

	! ------------------------------------------------------------------------------------
	! Lecture des données
	! ------------------------------------------------------------------------------------
	
	READ(FICH,*)Mmax
	READ(FICH,*)Nmax
	
	READ(FICH,*)LX
	READ(FICH,*)LY
	READ(FICH,*)D
	
	READ(FICH,*)rho
	READ(FICH,*)TREF
	READ(FICH,*)k0
	READ(FICH,*)k1
	READ(FICH,*)C0
	READ(FICH,*)C1
	
	READ(FICH,*)H1
	READ(FICH,*)H2

	READ(FICH,*)T0
	READ(FICH,*)DTINF
	READ(FICH,*)P
	
	READ(FICH,*)dt
	READ(FICH,*)dtLog

	READ(FICH,*)NQ
	READ(FICH,*)Qnom
	
	ALLOCATE(tabQ(1:Mmax,1:Nmax))		! ici je n'ai pas d'autre choix que de déclarer ici..
	DO k=1,NQ
		READ(FICH,*)m,n,a
		tabQ(m,n)=a
	ENDDO

	
	! ------------------------------------------------------------------------------------
	! Fermeture du fichier
	! ------------------------------------------------------------------------------------
	
	CLOSE(FICH)	
	
END SUBROUTINE load


! Routine principale de résolution
SUBROUTINE run()
	IMPLICIT NONE

	! --- Stockage des grandeurs géométriques --------------------------------------------
	! Ici nous n'avons pas besoin de tableaux ici puisque tout est constant !
	DOUBLE PRECISION :: dx, dy, Sx, Sy, V						! paramètres géométriques

	! --- Stockage des résultats ---------------------------------------------------------
	! Je ne mémorise que le champs de température courant (en global cf. ci-dessus) et
	! la correction vers le pas de temps suivant... donc pas d'indice de temps !
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: EVOL
	
	
	INTEGER :: m,n
	INTEGER :: isave						! un indice pour gérer la sauvegarde	
	DOUBLE PRECISION :: k					! la conductivité (locale)
	DOUBLE PRECISION :: TINF				! la température ambiante
	DOUBLE PRECISION :: FG, FD, FH, FB		! les flux
	DOUBLE PRECISION :: time		! le temps physique
	DOUBLE PRECISION :: nextLog		! l'instant de la prochaine sauvegarde du champs

	! ------------------------------------------------------------------------------------
	! Allocation du tableau
	! ------------------------------------------------------------------------------------

	ALLOCATE(T(1:Mmax,1:Nmax),EVOL(1:Mmax,1:Nmax))
	
	! ------------------------------------------------------------------------------------
	! Calcul des grandeurs du maillage
	! ------------------------------------------------------------------------------------

	dx = LX/Mmax
	dy = LY/Nmax
	
	Sx = dy*1.D0			! le 1.D0 est pour l'épaisseur (arbitraire)
	Sy = dy*1.D0
	
	V = dx*dy*1.D0

	! ------------------------------------------------------------------------------------
	! Condition initiale
	! ------------------------------------------------------------------------------------

	T(:,:) = T0

	! ------------------------------------------------------------------------------------
	! Boucle de résolution
	! ------------------------------------------------------------------------------------
	
	! Je n'ai pas d'indice de temps (pas besoin !)
	
	time = 0.D0					! je commence par le commencement
	nextLog = 0.D0				! je veux sauvegarder à l'instant 0.D0
	isave = 0					! la première sauvegarde est la #0
	
	! Boucle sur le temps
	DO WHILE (time<=D)
	
		TINF = T0+DTINF*sin(2*PI*time/P)			! on mémorise pour 
		
		! Doit-on sauvegarder ?
		IF (time>=nextLog) THEN
			! on doit sauvegarder
			CALL export(isave)
			isave = isave+1
			nextLog = nextLog + dtLog
		ENDIF
		
		
		! Boucle sur l'espace
		DO m=1,Mmax				
			DO n=1,Nmax				
		
		
				! ==== Calcul des flux ===========================================================
		
				! GAUCHE FG
				IF (m==1) THEN
					FG = 0.D0			! Symétrie
				ELSE
					k = calcK(T(m-1,n),T(m,n))
					FG = k*Sx*(T(m-1,n)-T(m,n))/dx
				ENDIF
		
				! DROITE FD
				IF (m==Mmax) THEN
					FD = H1*Sx*(TINF-T(m,n))
				ELSE
					k = calcK(T(m+1,n),T(m,n))
					FD = k*Sx*(T(m+1,n)-T(m,n))/dx
				ENDIF

				! BAS FB
				IF (n==1) THEN
					FB = 0.D0			! Symétrie
				ELSE
					k = calcK(T(m,n-1),T(m,n))
					FB = k*Sy*(T(m,n-1)-T(m,n))/dy
				ENDIF
		
				! HAUT FH
				IF (n==Nmax) THEN
					FH = H2*Sy*(TINF-T(m,n))
				ELSE
					k = calcK(T(m,n+1),T(m,n))
					FH = k*Sy*(T(m,n+1)-T(m,n))/dy
				ENDIF

				! ==== Bilans/Evolution ==========================================================
				! Obligatoire parce que j'ai encore besoin de l'ancien champs de température
				EVOL(m,n) = dt*(FG+FD+FH+FB+tabQ(m,n))/(rho*V)/(C0+C1*(T(m,n)-TREF))
		
			ENDDO
		ENDDO
		
		! ==== Mise à jour des température ===============================================
		! Maintenant je n'ai plus besoin des anciennes températures donc je peux les
		! écraser pour préparer l'itération suivante.
		
		T(:,:) = T(:,:) + EVOL(:,:)
		
		! et on passe au temps suivant...
		time = time + dt
		
	END DO
END SUBROUTINE run

! Ecriture des résultats pour exploitation (dans sphere.result)
SUBROUTINE export(i)
	IMPLICIT NONE
	
	INTEGER, INTENT(IN) :: i
	INTEGER :: m,n
	
	CHARACTER(LEN=5) :: txt			! pour l'indice du fichier
	
	! ------------------------------------------------------------------------------------
	! Ouverture du fichier en lecture
	! ------------------------------------------------------------------------------------

	! On enregistre un fichier par pas de temps...
	! L'extension .z est imposée par l'outil GLE...

	! Je ne controle pas le succés donc le programme s'arrêtera en cas de problème
	WRITE(txt,"(i5)"),i					! Astuce (??) Fortran pour convertir un nombre en chaîne...
	OPEN(FICH,FILE="field."//TRIM(adjustl(txt))//".z", ACTION="WRITE");   ! // pour la concaténation


	! ------------------------------------------------------------------------------------
	! Ecriture des données
	! ------------------------------------------------------------------------------------
	
	WRITE(FICH,*)"! nx ",Mmax," ny ",Nmax, " xmin 0 xmax ",LX," ymin 0 ymax ",LY 
	DO n=1,Nmax
		WRITE(FICH,*)(T(m,n),m=1,Mmax)
	ENDDO
	

	! ------------------------------------------------------------------------------------
	! Fermeture du fichier
	! ------------------------------------------------------------------------------------
	
	CLOSE(FICH)	
	
END SUBROUTINE export

FUNCTION calcK(T1,T2) 
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN) :: T1,T2
	DOUBLE PRECISION :: calcK
	
	calcK = k0+k1*(0.5*(T1+T2)-TREF)

END FUNCTION calcK

FUNCTION calcC(T1,T2) 
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN) :: T1,T2
	DOUBLE PRECISION :: calcC
	
	calcC = C0+C1*(0.5*(T1+T2)-TREF)

END FUNCTION calcC


END PROGRAM Cart2dExp



	
	
	
	