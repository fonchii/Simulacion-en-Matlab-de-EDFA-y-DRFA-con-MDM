clear all;clc
% Calculo y gráfico de frecuencias de corte para cada movo

% Parametros
a=25e-6 ; n1=1.46 ; n2=1.455 ; lambda=980e-9;

% Calculo de modos presentes en esa fibra
%tic;
%Modos = modos(a,n1,n2,lambda);
%t_end1=toc;

%freqs=modos.freqs_sort;
%isempty(find(modos.m2==freqs(5)));

%for i=1:length(Modos.LPsort)
%   Modos.LPPsort{i} = strcat('LP',Modos.LPsort{1});
%end

tic;
Modosv2 = modosv2(a,n1,n2,lambda);
t_end2=toc;

%fprintf('Tiempo 1: %.2f\nTiempo 2: %.2f\n',t_end1,t_end2)

%% Calculos superposicion
close all;
r = linspace(0, a, 400);
phi = linspace(0, 2*pi, 400);
[r,phi] = meshgrid(r, phi);
x = r.*cos(phi); y=r.*sin(phi);
for i = 1:5
    %v = Modosv2.v;
    v = Modosv2.freqs_sort(i);
    l = floor(str2double(Modosv2.LPsort(i))/10);
    m = mod(str2double(Modosv2.LPsort(i)),10);
    if l == 0
        z = besselj(l,v*r/a).*cos(l.*phi) / (pi*a^2);
        figure(i)
        mesh(x,y,z.^2); view(2); title(strcat('Modo LP',Modosv2.LPsort(i)));
    else 
        z_a = besselj(l,v*r/a).*cos(l.*phi) / (pi*a^2) ;
        z_b = besselj(l,v*r/a).*sin(l.*phi) / (pi*a^2) ;
        figure(i)
        subplot 211; mesh(x,y,z_a.^2); view(2); title(strcat(strcat('Modo LP',Modosv2.LPsort(i)),'a'));
        subplot 212; mesh(x,y,z_b.^2); view(2); title(strcat(strcat('Modo LP',Modosv2.LPsort(i)),'b'));
    end
end

%%
clear r phi ; %close all

v_01 = Modosv2.freqs_sort(1);
v_11 = Modosv2.freqs_sort(2);
v_21 = Modosv2.freqs_sort(3);
v_02 = Modosv2.freqs_sort(4);


M_01 = @(r,phi) ( besselj(0,v_01*r/a) ) .* ( cos(0*phi)) ;


M_11_a = @(r,phi) ( besselj(1,v_11*r/a) ) .* ( cos(1*phi) ) ;
M_11_b = @(r,phi) ( besselj(1,v_11*r/a) ) .* ( sin(1*phi) );

psi_21_a = @(r) ( besselj(2,v_21*r/a) ).^2 ;
f_21_a = @(phi) ( cos(1*phi) ).^2 ;
psi_21_b = @(r) ( besselj(2,v_21*r/a) ).^2 ;
f_21_b = @(phi) ( sin(1*phi) ).^2 ;

psi_02 = @(r) ( besselj(1,v_02*r/a) ).^2 ;
f_02 = @(phi) ( cos(0*phi) ).^2 ;


num = integral2( @(r,phi) abs (M_01(r,phi) ).* abs((M_11_b(r,phi) + M_11_a(r,phi) )) , 0, a, 0, 2*pi).^2;

den = integral2( @(r,phi) abs(M_01(r,phi)).^2 , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs(M_11_b(r,phi) + M_11_a(r,phi) ) ).^2 , 0, a, 0, 2*pi);


Gamma = (num/den)




