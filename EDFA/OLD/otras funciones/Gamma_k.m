function Gamma = Gamma_k(a,l,n_core,n_clad,lambda)
% rho : concentracion totales de iones de erbio
% rho = 3e24 ; %concentracion de iones de erbio
beta = 2*pi/lambda;
%k = 2*pi/lambda_0;
mu_0 = 4*pi*1e-7;   % Permeabilidad magnética del vacío
e_0 =  8.854*1e-12; % Permititividad eléctrica del vacío
k_0 = 2*pi*lambda * sqrt(mu_0*e_0);
%k = k_0;
k = beta;
U = a.*sqrt(k^2*n_core^2 - beta^2);
W = a.*sqrt(-beta^2 + k^2*n_clad^2);

% Psi(r,phi) = g(r) * h(phi)
% donde : g(r) = psi_r(r) + psi2_r(r)

psi_r = @(r)  besselj(l-1,(U.*r./a)).^2;                     % dentro del nucleo, valido para [0,a]
psi2_r = @(r) abs( ( (besselj(l-1,U) ./ besselk(l,W) ) .^2 ) ...
    .* (besselk(l,W*r/a)).^2 );                          % fuera del nucleo, valido para [a,inf]
radio = linspace(-a,a,5000);
plot(radio,psi_r(radio))

if(mod(l,2) == 0 )              % Modos pares
    f = @(phi) (cos(l*phi).^2);
elseif(mod(l,2)==1)             % Modos impares
    f = @(phi) (sin(l*phi).^2);
end

den = integral( @(r) (psi_r(r).*r) , 0 , a ) + integral( @(r) (psi2_r(r).*r) , a , inf )...
    * integral( @(phi) (f(phi)),0,2*pi );

% I(r,phi) = k * Psi(r,phi) , con k = 1/den

I_r = @(r) psi_r(r)/den;         
I_phi = @(phi) f(phi);

%Gamma = (1/rho) * integral(@(r) I_r(r).*r,0,a) * integral(@(phi) I_phi(phi),0,2*pi);
%Gamma = integral( @(r) I_r(r).*r,0,a ) * integral( @(phi) I_phi(phi),0,2*pi );


I_r2 = @(r) psi2_r(r);

Gamma = ( integral( @(r) I_r(r).*r,0,a ) + integral( @(r) I_r2(r).*r,a,inf ) ) * integral( @(phi) I_phi(phi),0,2*pi );
%IR1 = integral( @(r) I_r(r).*r,0,a );
%IR2 = integral( @(r) I_r2(r).*r,a,inf ) 
%IPHI = integral( @(phi) I_phi(phi),0,2*pi );
%fprintf(' Dentro del nucleo: %.4f \n Fuera del nucleo: %.4fe-14 \n' ,IR1 , IR2*1e14)

%fprintf(' Γ = %.4f \n' , Gamma)

% Gamma_k(25e-6,1,1.45,1.475,980e-9)