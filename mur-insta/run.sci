// ***********************************************************
// ** Exemple de code - Cours M2S - ENSGTI 2A - Energétique **
// ** Stéphane Gibout - 2014                                **
// **                                                       **
// ** Mur Instationnaire - Schéma de Crank-Nicholson        **
// **                                                       **
// ***********************************************************


function [ T ] = run(M,I,D,L,rho,C,k,T0,TG,TD)
    
    // -------------------------------------------------------
    // Valeurs calculées
    // -------------------------------------------------------
        
    // Pas de temps
    dt = D/I
    
    // Pas de temps
    dx = L/M
    
    // Coefficient <<alpha>>
    a = (dt*k)/(rho*C*dx*dx);
    
    // -------------------------------------------------------
    // Déclaration des vecteurs et matrices
    // avec remplissage . CI
    // -------------------------------------------------------

    // Matrice A : je profite ici de la structure particulière
    A = diag([(1+1.5*a) , (1+a)*linspace(1,1,M-2) , (1+1.5*a) ]);
    A = A + diag(-0.5*a*linspace(1,1,M-1),1) +  diag(-0.5*a*linspace(1,1,M-1),-1);

    // Matrice B : idem
    B = diag([(1-1.5*a) , (1-a)*linspace(1,1,M-2) , (1-1.5*a) ]);
    B = B + diag(0.5*a*linspace(1,1,M-1),1) +  diag(0.5*a*linspace(1,1,M-1),-1);

    // Vecteur C : on met des 0 partout puis on modifie les deux seuls valeurs non nulles
    C = zeros(M,1); 
    C([1 M])= 2*a*[TG TD]';
    
    // On peut maintenant calculer les matrices <<primes>>
    AA = inv(A)*B;
    CC = inv(A)*C;
    
    // Chaque colonne correspond au champs de température pour un pas
    // de temps donné. Puisque la numérotation des indices débute à 1
    // avec Scilab, on doit décaler les indices de colonne...
    T = zeros(M,I+1);
    T(:,1) = T0;
    
    
    // -------------------------------------------------------
    // Résolution
    // -------------------------------------------------------
    
    for i=1:I
        T(:,i+1) = AA*T(:,i) + CC; 
    end
    
    
    // -------------------------------------------------------
    // C'est tout !
    // -------------------------------------------------------

    
    
endfunction

M = 100;         // Nombre de noeuds en espace
I = 7200;        // Nombre de pas de temps
L = 0.1;        // Epaisseur du mur [m]
D = 5000;        // Durée totale [s]

rho = 1000;     // Masse volumique [kg/m3]
C = 1000;       // Capacité calorifique [J/(kg.K)]
k = 1;          // Conductivité thermique [W/(m.K)]

T0 = 20;        // Température initiale [°C]
TG = 30;        // Température imposée en x=0 [°C]
TD = 00;        // Température imposée en x=L [°C]

[ T ] = run(M,I,D,L,rho,C,k,T0,TG,TD)


clf();
plot2d(linspace(0,100*L,M),T(:,linspace(1,I+1,10)))
xlabel("x [cm]");
