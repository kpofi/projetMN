// fonction membrane
N0 = 80;
Nr = 40;
CFL = 0.50;
c = 1;
l01 = 2.40483;
l11 = 3.83171;
fin = 1;

function w = initialise(N0, Nr, l01, l11)
    w = zeros(Nr,N0 + 1);
    // dr represente detha
    dr = 1/(Nr - 1);
    // d0 represente dtheta
    d0 = 2*%pi/N0;
    for i=1:Nr-1
        for j=1:N0 + 1
            w(i,j) = besselj(0,l01*(i-1)*dr) + besselj(1,l11*(i-1)*dr)*cos((j-1)*d0);
        end
    end
    
endfunction

function membrane(N0, Nr, CFL, c, l01, l11, fin)
    // dr represente detha
    dr = 1/(Nr - 1);
    // d0 represente dtheta
    d0 = 2*%pi/N0;
    // dt represente dtho
    dt = CFL*dr*d0;
    // matrice de la position transversale w0 a l'instant 0 et 1
    w0 = initialise(N0,Nr,l01,l11);
    w1 = w0;
    w2 = zeros(Nr,N0 + 1);
    wex = w2;
    
    // coordonnées 
    x = zeros(Nr,(N0+1));
    y = zeros(Nr,(N0+1));
    zmax = abs(w0(1,1)*(1+50/100));
    zmin = -zmax;
    ri = 0; tj = 0;
    
    // variable pour les images
    p = 1000;
    
    // calcul des coodonnées
    for i=1:Nr
        for j=1:N0+1
            x(i,j) = ri*cos(tj);
            y(i,j) = ri*sin(tj);
            tj = tj + d0;
        end
        ri = ri + dr;
    end
    
    NF = fix(fin/dt);
    err = zeros(NF+1,1);
    errglob = zeros(NF+1,1);
    
    for n=0:NF
        drawlater;
        // calcul de w quand r=0
        w2(1,1) = 2*w1(1,1) - w0(1,1);
        for j=1:N0
            w2(1,1) = w2(1,1) + ( 4*(dt^2)/(N0*(dr^2)) )*(w1(2,j) - w1(1,1));
        end
         
        // r=0 donc invariance par rapport a theta
        for j=2:N0+1
            w2(1,j) = w2(1,1);
        end
        
        // calcul en utilisant le schema
        for i=2:Nr-1
            for j=2:N0
                w2(i,j) = 2*w1(i,j) - w0(i,j);
                w2(i,j) = w2(i,j) + ((dt/dr)^2)*(w1(i+1,j) - 2*w1(i,j) + w1(i-1,j));
                w2(i,j) = w2(i,j) + ((dt/((i-1)*dr*d0))^2)*(w1(i,j+1) - 2*w1(i,j) + w1(i,j-1));
                w2(i,j) = w2(i,j) + ((dt/dr)^2)*(1/(2*(i-1)))*(w1(i+1,j) - w1(i-1,j));
            end
        end
        
        // j = N0 + 1 pareil pour j = 1 pour i et t fixe car w est periodique 
        for i=2:Nr-1
            w2(i,N0+1) = 2*w1(i,N0+1) - w0(i,N0+1);
            w2(i,N0+1) = w2(i,N0+1) + ((dt/dr)^2)*(w1(i+1,N0+1) - 2*w1(i,N0+1) + w1(i-1,N0+1));
            w2(i,N0+1) = w2(i,N0+1) + ((dt/((i-1)*dr*d0))^2)*(w1(i,2) - 2*w1(i,N0+1) + w1(i,N0));
            w2(i,N0+1) = w2(i,N0+1) + ((dt/dr)^2)*(1/(2*(i-1)))*(w1(i+1,N0+1) - w1(i-1,N0+1));
            
            w2(i,1) = w2(i,N0+1);
        end
        
        for i=1:Nr
            for j=1:N0+1
                wex(i,j) = besselj(0,l01*(i-1)*dr)*cos(l01*n*dt) + besselj(1,l11*(i-1)*dr)*cos((j-1)*d0)*cos(l11*n*dt);
                if ( errglob(n+1) < abs( wex(i,j) - w0(i,j) ) ) then
                    errglob(n+1) = abs(wex(i,j) - w0(i,j));
                end
            end
        end
        
        err(n+1) = abs( wex(1,1) - w0(1,1) );
        
        w0 = w1;
        w1 = w2;
        
 
        clf;
        a=gca();
        a.data_bounds=[-1,-1,zmin; 1,1,zmax];
        surf(x,y,w0);
        
        p = p + 1;
        if ((modulo(n,10) == 0) & (n < 2000)) then
            winnum=winsid();
            nom_image1 ='image11_'+string(p)+'.gif';
            xs2gif(winnum($),nom_image1);
        end;
        drawnow;
        
    end
    clf;
    temps = (0:1:NF)';
    temps = dt*temps;
    plot2d(temps,[err,errglob],style=[2,5]);
    
    // erreur globale en rouge
    winnum=winsid();
    nom_image1 ='erreur1.gif';
    xs2gif(winnum($),nom_image1);
    
endfunction

