// ***********************************************************
// ** Exemple de code - Cours M2S - ENSGTI 2A - Energétique **
// ** Stéphane Gibout - 2014                                **
// **                                                       **
// ** Domaine cartésien 2D instationnaire  (implicite)      **
// **                                                       **
// ***********************************************************



function [ x , y , T ] = run(M,N,I,LX,LY,D,rho,C,k,T0,H1,H2,TINF,PHI)

    // -------------------------------------------------------
    // Valeurs calculées
    // -------------------------------------------------------

    // Pas de temps
    dt = D/I

    // Pas d'espace
    dx = LX/M
    dy = LY/N

    // Coefficient <<alpha>>
    bx = k*dt/(rho*C*dx*dx);
    by = k*dt/(rho*C*dy*dy);

    gx = dt/(rho*C*dx);
    gy = dt/(rho*C*dy);

    // -------------------------------------------------------
    // Fonction utilitaire pour le calcul des indices
    // -------------------------------------------------------

    function [k] = k(m,n) 
       k =  m+(n-1)*M
    endfunction

    // -------------------------------------------------------
    // Construction de la matrice [A] et du vecteur {C}
    // -------------------------------------------------------

    A = zeros(M*N,M*N);
    C = zeros(M*N,1);

    // Centre 
    for m=2:M-1
        for n=2:N-1
            // terme m,n-1
            A(k(m,n),k(m,n-1)) = -by;
            // terme m-1,n
            A(k(m,n),k(m-1,n)) = -bx;
            // terme m,n
            A(k(m,n),k(m,n)) = (1+2*bx+2*by);
            // terme m+1,n
            A(k(m,n),k(m+1,n)) = -bx;
            // terme m,n+1
            A(k(m,n),k(m,n+1)) = -by;
        end
    end

    // Côté gauche
    m=1;
    for n=2:N-1
        // terme m,n-1
        A(k(m,n),k(m,n-1)) = -by;
        // terme m,n
        A(k(m,n),k(m,n)) = (1+bx+2*by);
        // terme m+1,n
        A(k(m,n),k(m+1,n)) = -bx;
        // terme m,n+1
        A(k(m,n),k(m,n+1)) = -by;
    end

    // Côté droit
    m=M;
    for n=2:N-1
        // terme m,n-1
        A(k(m,n),k(m,n-1)) = -by;
        // terme m-1,n
        A(k(m,n),k(m-1,n)) = -bx;
        // terme m,n
        A(k(m,n),k(m,n)) = (1+bx+2*by+gx*H1);
        // terme m,n+1
        A(k(m,n),k(m,n+1)) = -by;

        C(k(m,n)) = gx*H1*TINF;
    end

    // Côté bas
    n=1
    for m=2:M-1
        // terme m-1,n
        A(k(m,n),k(m-1,n)) = -bx;
        // terme m,n
        A(k(m,n),k(m,n)) = (1+2*bx+by);
        // terme m+1,n
        A(k(m,n),k(m+1,n)) = -bx;
        // terme m,n+1
        A(k(m,n),k(m,n+1)) = -by;
    end

    // Côté haut
    n=N;
    for m=2:M-1
        // terme m,n-1
        A(k(m,n),k(m,n-1)) = -by;
        // terme m-1,n
        A(k(m,n),k(m-1,n)) = -bx;
        // terme m,n
        A(k(m,n),k(m,n)) = (1+2*bx+by+gy*H2);
        // terme m+1,n
        A(k(m,n),k(m+1,n)) = -bx;

        C(k(m,n)) = gy*(H2*TINF+PHI);
    end

    // Angle bas-gauche
    m=1;
    n=1;
    // terme m,n
    A(k(m,n),k(m,n)) = (1+bx+by);
    // terme m+1,n
    A(k(m,n),k(m+1,n)) = -bx;
    // terme m,n+1
    A(k(m,n),k(m,n+1)) = -by;

    // Angle bas-droite
    m=M;
    n=1;
    // terme m-1,n
    A(k(m,n),k(m-1,n)) = -bx;
    // terme m,n
    A(k(m,n),k(m,n)) = (1+bx+by+gx*H1);
    // terme m,n+1
    A(k(m,n),k(m,n+1)) = -by;

    C(k(m,n)) = gx*H1*TINF;

    // Angle haut-gauche
    m=1;
    n=N;
    // terme m,n-1
    A(k(m,n),k(m,n-1)) = -by;
    // terme m,n
    A(k(m,n),k(m,n)) = (1+bx+by+gy*H2);
    // terme m+1,n
    A(k(m,n),k(m+1,n)) = -bx;

    C(k(m,n)) = gy*(H2*TINF+PHI);


    // Angle haut-droit
    m=M;
    n=N;
    // terme m,n-1
    A(k(m,n),k(m,n-1)) = -by;
    // terme m-1,n
    A(k(m,n),k(m-1,n)) = -bx;
    // terme m,n
    A(k(m,n),k(m,n)) = (1+bx+by+gx*H1+gy*H2);

    C(k(m,n)) = gx*H1*TINF + gy*(H2*TINF+PHI);


    // -------------------------------------------------------
    // Déclaration du vecteur solution et CI
    // -------------------------------------------------------

    T = zeros(M*N,I+1);
    T(:,1) = T0;


    // -------------------------------------------------------
    // Résolution
    // -------------------------------------------------------

    // calcul de l'inverse de A

    iA = inv(A);

    // Résolution

    for i=1:I
        T(:,i+1) = iA*(C+T(:,i)); 
    end


    // -------------------------------------------------------
    // On réorganise la matrice (et on construit les x et y)
    // -------------------------------------------------------

    // On retourne une <<hyper-matrice>> à trois indices.
    // T(:,:,i) permet de retourner le champ à l'itération i

    T = matrix(T,[M N I+1]);

    x = ((1:M) -0.5)*dx;
    y = ((1:N) -0.5)*dy;



endfunction

M = 20;         // Nombre de noeuds en espace / x
N = 20;         // Nombre de noeuds en espace / y
LX = 0.1;       // Dimensions [m] /x
LY = 0.1;       // Dimensions [m] /y
I = 3600;       // Nombre de pas de temps
D = 3600;       // Durée totale [s]

rho = 1000;     // Masse volumique [kg/m3]
C = 1000;       // Capacité calorifique [J/(kg.K)]
k = 1;          // Conductivité thermique [W/(m.K)]

T0 = 20;        // Température initiale [°C]
H1 = 100;        // Coefficient d'échange côté droit [SI]
H2 = 100;       // Coefficient d'échange côté haut [SI]
TINF = 75;      // Température extérieure [°C]
PHI = 5000;     // Densité de flux surface sup [W/m2]


[ X,Y,T ] = run(M,N,I,LX,LY,D,rho,C,k,T0,H1,H2,TINF,PHI)


clf();
xset("colormap",jetcolormap(64))
for i=1:10:I+1
    drawlater
    clf
    plot3d1(X,Y,T(:,:,i),theta=-120,alpha=%pi/10)
    drawnow
end






