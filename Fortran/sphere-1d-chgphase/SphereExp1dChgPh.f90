PROGRAM SphereExp1dChgPh

	! ***********************************************************
	! ** Exemple de code - Cours M2S - ENSGTI 2A - Energétique **
	! ** Stéphane Gibout - 2014                                **
	! **                                                       **
	! ** Sphère 1D Instationnaire - Schéma Explicit            **
	! ** Avec changement de phase !                            **
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
	INTEGER :: Mmax					! Nombre de pas d'espace 1..Mmax
	INTEGER :: Imax					! Nombre de pas d'espace 0..Imax	
	
	DOUBLE PRECISION :: D, R		! Durée [s] et rayon [m]
	DOUBLE PRECISION :: rho, k		! Paramètres thermophysiques [SI]
	DOUBLE PRECISION :: CL, CS		! Paramètres thermophysiques [SI]
	DOUBLE PRECISION :: TF, LF		! Paramètres thermophysiques [SI]
	DOUBLE PRECISION :: HG			! Coefficient échange global [SI]
	
	DOUBLE PRECISION :: T0, TINF	! Température initiale et ambiante [°C]
	
	
	! --- Stockage des grandeurs géométriques --------------------------------------------
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DR, S, V		! distances, surfaces et volumes
	DOUBLE PRECISION :: dt										! le pas de temps (constant)

	! --- Stockage des résultats ---------------------------------------------------------
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: T		! température
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Y		! fraction liquide
	
	! ====================================================================================
	! Début du programme principal
	! ====================================================================================
	
	

	! --- Chargement des paramètres depuis le fichier sphere.data ------------------------ 
	CALL load()		
	PRINT*,'Chargement des données.................OK'

	! --- Construction du maillage -------------------------------------------------------
 	CALL mesh()
	PRINT*,'Construction du maillage...............OK'

	! --- Calcul -------------------------------------------------------------------------
	CALL run()
	PRINT*,'Calcul.................................OK'	

	! --- Export des résultats (dans sphere.result) --------------------------------------
	CALL export()
	PRINT*,'Export.................................OK'	
	
	
PRINT*,S(Mmax),sum(V)	

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
	
	! ------------------------------------------------------------------------------------
	! Ouverture du fichier en lecture
	! ------------------------------------------------------------------------------------

	! Je ne controle pas le succés donc le programme s'arrêtera en cas de problème
	OPEN(FICH,FILE="./sphere.data", ACTION="READ");

	! ------------------------------------------------------------------------------------
	! Lecture des données
	! ------------------------------------------------------------------------------------
	
	READ(FICH,*)Mmax
	READ(FICH,*)Imax
	
	READ(FICH,*)R
	READ(FICH,*)D
	
	READ(FICH,*)rho
	READ(FICH,*)k
	READ(FICH,*)CS
	READ(FICH,*)CL
	READ(FICH,*)TF
	READ(FICH,*)LF
	
	
	READ(FICH,*)HG
	READ(FICH,*)T0
	READ(FICH,*)TINF

	! ------------------------------------------------------------------------------------
	! Fermeture du fichier
	! ------------------------------------------------------------------------------------
	
	CLOSE(FICH)	
	
END SUBROUTINE load

! Calcul des informations relatives à la géométrie
SUBROUTINE mesh()
	IMPLICIT NONE

	INTEGER :: m
	DOUBLE PRECISION :: deltaR
	
	! ------------------------------------------------------------------------------------
	! Allocation des tableaux
	! ------------------------------------------------------------------------------------
	! Pas de contrôle de succés... je vous laisse le faire
	
	! Distance entre deux noeuds (m=1..Mmax-1) ou entre le noeuds et la frontière (m=Mmax)
	! Ici on n'a pas besoin de la distance entre le centre de la sphère et le 1er noeud.
	! Je fais le choix de considérer directement le maillage non-uniforme juste pour vous
	! donner la technique utilisable...
	ALLOCATE(DR(1:Mmax))
	
	! Surface d'échange entre le noeud m et le noeuds m+1 (ou extérieur si m=M)
	! Ici encore, on n'a pas besoin de la surface en r=0 (puisqu'elle est nulle)
	ALLOCATE(S(1:Mmax))

	! Volume du noeud
	ALLOCATE(V(1:Mmax))

	
	! ------------------------------------------------------------------------------------
	! Calculs 
	! ------------------------------------------------------------------------------------
	
	! Discrétisation en temps
	dt = D/Imax

	! Discrétisation en espace
	
	deltaR = R/Mmax
	
	! Les distances m -> m+1
	DR(1:Mmax-1)=deltaR
	DR(Mmax)=deltaR/2

	
	! Les surfaces
	DO m=1,Mmax
		S(m) = 4*PI*m**2*deltaR**2  	! correspond au demi-pas supérieur
	ENDDO

	! Les volumes
	DO m=1,Mmax
		V(m) = 4*PI*((3*m-3)*m+1 )*deltaR**3/3
	ENDDO

END SUBROUTINE mesh

! Routine principale de résolution
SUBROUTINE run()
	IMPLICIT NONE
	
	INTEGER :: i,m
	DOUBLE PRECISION :: FG, FD		! les flux

	! On n'a pas besoin de stocker les enthalpies donc un tableau
	! mono-dimensionnel est suffisant
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: h

	! ------------------------------------------------------------------------------------
	! Allocation du tableau
	! ------------------------------------------------------------------------------------

	ALLOCATE(h(1:Mmax))

	ALLOCATE(T(1:Mmax,0:Imax))
	ALLOCATE(Y(1:Mmax,0:Imax))
		
	! ------------------------------------------------------------------------------------
	! Condition initiale
	! ------------------------------------------------------------------------------------

	T(:,0) = T0
	h(:) = calcH(T0)		! on travaille en massique
	Y(:,0) = calcY(calcH(T0))

	! ------------------------------------------------------------------------------------
	! Boucle de résolution
	! ------------------------------------------------------------------------------------

	DO i=0,Imax-1				! i représente le pas de temps courant i.e. connu
		DO m=1,Mmax				! et m le noeud courant
		
		
		! ==== Calcul des flux ===========================================================
		! Notez ici que je sépare systématiquement le traitement des différents flux
		
		! GAUCHE FG
		IF (m==1) THEN
			FG = 0.D0			! Symétrie
		ELSE
			! Vous remarquerez ici que le choix fait pour représenter les indices de 
			! surfaces et de distances conduit à des expressions simples et vérifiables.
			FG = k*S(m-1)*(T(m-1,i)-T(m,i))/DR(m-1)
		ENDIF
		
		! DROITE FD
		IF (m==Mmax) THEN
			FD = HG*S(m)*(TINF-T(m,i))
		ELSE
			FD = k*S(m)*(T(m+1,i)-T(m,i))/DR(m)
		ENDIF

		! ==== Bilans ====================================================================
		! Vous pouvez vérifier au passage que j'écris toujours à l'intérieur de mes 
		! tableaux (compilez avec un -fbounds-check pour vérifier)
		
		! C'est ici qu'on doit découpler le bilan en enthalpie et la mise à jours des
		! des températures et fractions liquides
		
		h(m) = h(m) + dt*(FG+FD)/(rho*V(m))
		
		T(m,i+1) = calcT(h(m))
		Y(m,i+1) = calcY(h(m))
		
		
		ENDDO
	ENDDO
END SUBROUTINE run



! Routine de calcul de h en fonction de T
FUNCTION calcH(T) 
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: T
	DOUBLE PRECISION :: calcH
	
	IF (T<TF) THEN
		! SOLIDE
		calcH = CS*(T-TF)
	ELSEIF (T>TF) THEN
		! LIQUIDE
		calcH = LF+CL*(T-TF)
	ELSE ! T=TF
		! FUSION
		WRITE(*,*)"ERREUR T=TF"
		STOP
	ENDIF
END FUNCTION

! Routine de calcul de T en fonction de h
FUNCTION calcT(h) 
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: h
	DOUBLE PRECISION :: calcT
	
	IF (h<0) THEN
		! SOLIDE
		calcT = TF + h/CS
	ELSEIF (h>LF) THEN
		! LIQUIDE
		calcT = TF+(h-LF)/CL
	ELSE
		! EQUILIBRE
		calcT = TF
	ENDIF
END FUNCTION

! Routine de calcul de Y en fonction de h
FUNCTION calcY(h) 
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: h
	DOUBLE PRECISION :: calcY
	
	IF (h<0) THEN
		! SOLIDE
		calcY = 0
	ELSEIF (h>LF) THEN
		! LIQUIDE
		calcY = 1
	ELSE
		! EQUILIBRE
		calcY = h/LF
	ENDIF
END FUNCTION


! Ecriture des résultats pour exploitation (dans sphere.result)
SUBROUTINE export()
	IMPLICIT NONE
	
	INTEGER :: i,m
	
	! ------------------------------------------------------------------------------------
	! Ouverture du fichier en lecture
	! ------------------------------------------------------------------------------------

	! Je ne controle pas le succés donc le programme s'arrêtera en cas de problème
	OPEN(FICH,FILE="./sphere.result", ACTION="WRITE");

	! ------------------------------------------------------------------------------------
	! Ecriture des données
	! ------------------------------------------------------------------------------------
	
	! Pour faciliter le traitement avec GLE (ou excel), je vais cracher le champs de
	! température pas te temps par pas de temps...
	! Avec en première colonne le temps physique
	
	DO i=0,Imax
		WRITE(FICH,*)i*dt,(T(m,i), Y(m,i),m=1,mMax)
	ENDDO
	

	! ------------------------------------------------------------------------------------
	! Fermeture du fichier
	! ------------------------------------------------------------------------------------
	
	CLOSE(FICH)	
	
END SUBROUTINE export


END PROGRAM SphereExp1dChgPh



	
	
	
	
