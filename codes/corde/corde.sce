   //fonction corde
   driver("Rec");
   SR = 44100;                     // taux d'échantillonage
   B = 0.001;                      // paramètre d'inharmonicité
   f = 110;                        // fréquence fondamentale
   TF = 4;                         // durée de la simulation
   x0 = 0.1;                       // position où la corde est pincée
   c0 = 1;                         // hauteur du pincement
   rp = [0.3 0.7];                 // positions des micros
   loss = [100, 10; 1000, 8];      // coouples de fréquence/temps décroissance
   
   function e = epsi(w, g, kp)
       e = (-(g^2) + ((g^4) + 4*((kp*w)^2) )^(1/2) )/2*(kp^2);
   endfunction;
   
   function hm = hmin(g, kp, k)
       hm = ( ( (g*k)^2 + ( (g*k)^4 + 16*(kp*k)^2 )^(1/2) )/2 )^(1/2);
   endfunction;
   
   function u = etat_initial (x0, c0, N)
       u = ones(N-1,1);
       a0 = c0/x0;
       a1 = c0/(1-x0);
       for i=1:N-1
           if ((i/N) < x0) then
               u(i) = a0*(i/N);
           else
               u(i) = a1*(1-(i/N));
           end
       end
   endfunction
   
   function k = indice(p, N)
       i = 0;
       while (i/N) <= p
           i = i + 1;
       end
       k = i - 1;
   endfunction
   
   function p = variation(U,p11,p1,N)
       a = 0; b = 0;
       if (p11 == 0) then
           a = N*U(1);
       else
           if (p11 == N - 1) then
               a = -N*U(p11);
               b = -a;
           else
               a = N*(U(p11 + 1) - U(p11));
               b = U(p11) - a*(p11/N);
           end
       end
       p = a*p1 + b;
   endfunction
   
   function corde(SR, B, f, TF, x0, c0, rp, loss)
       w1 = 2*%pi*loss(1);
       w2 = 2*%pi*loss(2);
       // gamma
       g = 2*f;
       // kappa
       kp = (g*B^(1/2))/%pi;
       // calcul de epsilon 1 et 2
       e1 = epsi(w1,g,kp);
       e2 = epsi(w2,g,kp);
       // calcul de sigma 0 et 1
       s0 = ((6*log(10))/(e2-e1))*((e2/loss(3))-(e1/loss(4)));
       s1 = ((6*log(10))/(e2-e1))*((-1/loss(3))+(1/loss(4)));
       // pas en temps
       k = 1/SR;
       // calcul de h
       N = fix(1/hmin(g,kp,k));
       h = 1/N;
       // matrice de Toeplitz
       n = N - 1;
       v = 0*[1:n];
       v(1) = -2; v(2) = 1;
       D = (1/(h^2))*toeplitz(v);
       D2 = D^2;
       // calcul de A B et C
       A = (1 + s0*k)*eye(n,n) - s1*k*D;
       B0 = -2*eye(n,n) -((g*k)^2)*D + ((kp*k)^2)*D2;
       C = (1 - s0*k)*eye(n,n) + s1*k*D;
       // calcul de U0 et U1
       U0 = etat_initial(x0,c0,N);
       U1 = U0;
       U2 = 0*ones(n,1);
       x = (1/N):(1/N):(1-(1/N));
       p = 0;
       // indice des points qui encadrent p1
       p11 = indice(rp(1),N);
       // indice des points qui encadrent p2
       p21 = indice(rp(2),N);
       // nombre total d'iterations
       NF = fix(TF/k);
       // vecteur out
       out = zeros(2,NF);
       // enregistrement de la variation
       i = 1;
       out(1,i) = variation(U0,p11,rp(1),N);
       out(2,i) = variation(U0,p21,rp(2),N);
       
       
       for t=2:NF
           drawlater;
           U2 = A\(-(B0*U1 + C*U0));
           U0 = U1;
           U1 = U2;
           
           //enregistrement de la variation
           i = i + 1;
           out(1,i) = variation(U0,p11,rp(1),N);
           out(2,i) = variation(U0,p21,rp(2),N);
           clf;
           plot(x,U0');
           a=gca();
           a.data_bounds=[0, -c0*(1 + 20/100); 1, c0*(1 + 20/100)];
           drawnow;
           
           
           // on fait 1000 + t pour pouvoir faire l'animation apres
           p = 1000 + t;
           //on enregistre 1/10 image
           if ((modulo(t,10) == 0) & (t < 2000) ) then
               nom_image='image_'+string(p)+'.gif';
               winnum=winsid();
               xs2gif(winnum($),nom_image);
           end
       end
       
       // on affiche l'evolution des points p1 et p2
       clf;
       temps = (0:1:NF-1)';
       temps = k*temps;
       plot2d(temps,[out(1,1:NF)', out(2,1:NF)'],style=[2,5]);
       // on enregistre l'image out1 en bleu et out2 en rouge
       nom_image='aa_son.gif';
       winnum=winsid();
       xs2gif(winnum($),nom_image);
       
       // transformées de Fourier
       clf;
       s1 = fftshift( out(1,1 : size(out, 2)) );
       s2 = fftshift( out(2,1 : size(out, 2)) );
       plot2d(temps,[s1', s2'],style=[2,5]);
       // on enregistre l'image s1 en bleu et s2 en rouge
       nom_image='aa_fourier.gif';
       winnum=winsid();
       xs2gif(winnum($),nom_image);
       
       // on joue et on enregistre le son 1
       playsnd(out(1,1:size(out,2)),SR);
       savewave('son1.wav', out(1,1 : size(out, 2)), SR);
       
       // on joue et on enregistre le son 2
       playsnd(out(2,1:size(out,2)),SR);
       savewave('son2.wav', out(2,1 : size(out, 2)), SR);
       
   endfunction;
   
