a = fibra.radio/7;
r = linspace(0, a, 400);
phi = linspace(0, 2*pi, 400);
[r,phi] = meshgrid(r, phi);
x = r.*cos(phi); y=r.*sin(phi);

r2 = linspace(a, 1.5*a, 400);
phi2 = linspace(0, 2*pi, 400);
[r2,phi2] = meshgrid(r2, phi2);
x2 = r2.*cos(phi2); y2=r2.*sin(phi2);

% funciones respecto a parámetros U y W 

beta = 2*pi./lambda;
k = beta;
U = a.*sqrt(k.^2*n1^2 - beta.^2);
W = a.*sqrt(-beta.^2 + k.^2.*n2^2);
psi_clad11 = abs( ( (besselj(1,U(1)) ./ besselk(1,W(1)) )).* (besselk(1,W(1)*r2./(a)))).*(cos(1.*phi2)) /(pi*a^2) ;
psi01 = besselj(0,(U(1)*r./(a))).*(cos(0.*phi)) /(pi*a^2) ;
psi11 = besselj(1,(U(1)*r./(a))).*(cos(1.*phi)) /(pi*a^2) ;

hold on;mesh(x,y,psi11.^2) ; mesh(x2,y2,psi_clad11.^2);

%%

clear r r2 phi phi2

psi11 = @(r,phi) besselj(1,(U(1)*r./(a))).*(cos(1.*phi));
psi_clad11 = @(r,phi) abs( ( (besselj(1,U(1)) ./ besselk(1,W(1)) )).* (besselk(1,W(1)*r./(a)))).*(cos(1.*phi));


Gamma_norm_intensity = ( integral2( @(r,phi) abs( psi11(r,phi).*r) , 0, a, 0, 2*pi) ...
                        + integral2( @(r,phi) abs( psi_clad11(r,phi).*r ) , a, inf, 0, 2*pi) ) /(pi*a^2) ;
                    
psi01 = @(r,phi) besselj(0,(U(1)*r./(a))).*(cos(0.*phi));
psi_clad01 = @(r,phi) abs( ( (besselj(0,U(1)) ./ besselk(0,W(1)) )).* (besselk(0,W(1)*r./(a)))).*(cos(0.*phi));

Gamma_norm_intensity = ( integral2( @(r,phi) abs( psi01(r,phi).*r) , 0, a, 0, 2*pi) ...
                        + integral2( @(r,phi) abs( psi_clad01(r,phi).*r ) , a, inf, 0, 2*pi) ) /(pi*a^2) 
                    
                    
%% funciones respecto a frecuencia normalizada 

a = fibra.radio;
v=(pi*a./lambda)*sqrt(fibra.n1^2-fibra.n2^2);
r = linspace(0, a, 400);
phi = linspace(0, 2*pi, 400);
[r,phi] = meshgrid(r, phi);
x = r.*cos(phi); y=r.*sin(phi);

r2 = linspace(a, 2*a, 400);
phi2 = linspace(0, 2*pi, 400);
[r2,phi2] = meshgrid(r2, phi2);
x2 = r2.*cos(phi2); y2=r2.*sin(phi2);

core11 = besselj(1,(v(1)*r./(a))).*(cos(1.*phi)) /(pi*a^2) ;
clad11 =  ( ( (besselj(1,v(1)) ./ besselk(1,v(1)) )).* (besselk(1,v(1)*r2./(a)))).*(cos(1.*phi2)) /(pi*a^2) ;


hold on;mesh(x,y,core11.^2) ; mesh(x2,y2,clad11.^2);
 
%%


clear r r2 phi phi2

core11 = @(r,phi) besselj(1,(v(1)*r./(a))).*(cos(1.*phi)) /(pi*a^2) ;
clad11 = @(r,phi) ( ( (besselj(1,v(1)) ./ besselk(1,v(1)) )).* (besselk(1,v(1)*r./(a)))).*(cos(1.*phi)) /(pi*a^2) ;


Gamma_norm_intensity = ( integral2( @(r,phi) abs( core11(r,phi).*r) , 0, a, 0, 2*pi) ...
                        + integral2( @(r,phi) abs( clad11(r,phi).*r ) , a, inf, 0, 2*pi) )  ;
                    
core01 = @(r,phi) besselj(1,(v(1)*r./(a))).*(cos(0.*phi)) /(pi*a^2) ;
clad01 = @(r,phi) ( ( (besselj(0,v(1)) ./ besselk(0,v(1)) )).* (besselk(0,v(1)*r./(a)))).*(cos(0.*phi)) /(pi*a^2) ;

Gamma_norm_intensity = ( integral2( @(r,phi) abs( core01(r,phi).*r) , 0, a, 0, 2*pi) ...
                        + integral2( @(r,phi) abs( clad01(r,phi).*r ) , a, inf, 0, 2*pi) ) 
