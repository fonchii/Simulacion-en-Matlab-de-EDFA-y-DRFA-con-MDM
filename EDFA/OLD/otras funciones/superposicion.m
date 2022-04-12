%function Superposicion = superposicion(a,l,l_s,n_core,n_clad,lambda_p,lambda_s)
    % superposicion(25e-7,1,1,1.475,1.45,980e-9,1550e-9);

close all; clc
a = 10e-7; l= 0; l_s=2; n_core=1.475; n_clad=1.45; lambda_p=980e-9; lambda_s=1550e-9;


% Señal
beta_s = 2*pi/lambda_s;
k_s = beta_s;
U_s = a.*sqrt(k_s^2*n_core^2 - beta_s^2);
W_s = a.*sqrt(-beta_s^2 + k_s^2*n_clad^2);

    % Funcion
psi_r_s = @(r)  besselj(l_s,(U_s.*r./a)).^2;                     % dentro del nucleo, valido para [0,a]
psi2_r_s = @(r) abs( ( (besselj(l_s,U_s) ./ besselk(l,W_s) ) .^2 ) ...
    .* (besselk(l,W_s*r/a)).^2 );                          % fuera del nucleo, valido para [a,inf]


% Bombeo
beta_p = 2*pi/lambda_p;
k_p = beta_p;
U_p = a.*sqrt(k_p^2*n_core^2 - beta_p^2);
W_p = a.*sqrt(-beta_p^2 + k_p^2*n_clad^2);

    % Funcion
psi_r_p = @(r)  besselj(l,(U_p.*r./a)).^2;                     % dentro del nucleo, valido para [0,a]
psi2_r_p = @(r) ( ( (besselj(l,U_p) ./ besselk(l,W_p) ) .^2 ) ...
    .* (besselk(l,W_p*r/a)).^2 );                          % fuera del nucleo, valido para [a,inf]

if(mod(l,2) == 0)    % Modos pares
    %f_p = @(phi) (cos(l*phi).^2);                            % Definición de la tesis
    f_p = @(phi) ( (cos((l+1)*phi).^2+cos((l-1)*phi).^2 ));
elseif(mod(l,2)==1)	% Modos impares
    %f_p = @(phi) (sin(l*phi).^2);                           % Definición de la tesis
    f_p = @(phi) ( (sin((l+1)*phi).^2+sin((l-1)*phi).^2 ));
end

if(mod(l,2) == 0 )	% Modos pares
    %f_s = @(phi) (cos(l_s*phi).^2);                          % Definición de la tesis
    f_s = @(phi) ( (cos((l_s+1)*phi).^2+cos((l_s-1)*phi).^2 ));
elseif(mod(l,2)==1)	% Modos impares
    %f_s = @(phi) (sin(l_s*phi).^2);                           % Definición de la tesis
    f_s = @(phi) ( (sin((l_s+1)*phi).^2+sin((l_s-1)*phi).^2 ));
end

B = integral( @(r) (psi_r_s(r).*psi_r_s(r)) , -a , a ) * integral( @(phi) (f_s(phi).*f_s(phi)),0,2*pi );
Superposicion = integral( @(r) (psi_r_p(r).*psi_r_s(r)) , -a , a )* integral( @(phi) (f_p(phi).*f_s(phi)),0,2*pi )/B;
fprintf('Modo de señal: %i \nModo de Bombeo: %i\nSuperposición: %.3f\n\n' ,l_s,l, Superposicion)

% Gráficos
radio = linspace(0,a,50000) ; radio2 = linspace(1.001*a,1.3*a,50000); radio3 = linspace(0,3*a,50000);
legend1 = strcat(strcat('Pump @',int2str(lambda_p*1e9)),'nm'); legend2 = strcat(strcat('Signal @',int2str(lambda_s*1e9)),'nm');
        
% Intensidad vs radio
f1 = figure; 
f1.WindowState = 'maximized';
subplot 221
plot(radio3,psi_r_p(radio3),'DisplayName',legend1,'LineWidth',1) ; hold on ; plot(radio3,psi_r_s(radio3),'DisplayName',legend2,'LineWidth',1); title('Dentro del nucleo')
xticks([-a -a/2 0 a/2 a 2*a 3*a]) ; xticklabels({'-a','-a/2','0','a/2','a','2a','3a'}) ;
%set(fig1,'Xticks',[-a -a/2 0 a/2 a] , 'Xticklabels',({'-a','-a/2','0','a/2','a'})) ;
xlim([min(radio3),max(radio3)]); legend();
xline(a,'--',{'Nucleo - Revestimiento'});

subplot 222
plot(radio2,psi2_r_p(radio2),'DisplayName',legend1 ,'LineWidth',1) ; hold on ; plot(radio2,psi2_r_s(radio2),'DisplayName',legend2 ,'LineWidth',1); title('Fuera del nucleo')
xlim([min(radio2),max(radio2)]); legend();
xticks([a 1.1*a 1.2*a 1.3*a]) ; xticklabels({'a','1.1 a','1.2 a', '1.3 a'}) ;

subplot (2,2,[3 4])
plot([radio radio2], [psi_r_p(radio) psi2_r_p(radio2)],'DisplayName',legend1 ,'LineWidth',1) ; hold on
plot([radio radio2], [psi_r_s(radio) psi2_r_s(radio2)],'DisplayName',legend2 ,'LineWidth',1) ; legend();
xticks([0 a/4 a/2 3*a/4 a 1.1*a 1.2*a 1.3*a]) ; xticklabels({'0','a/4','a/2','3a/4','a','1.1 a','1.2 a', '1.3 a'}) ;
xline(a,'--',{'Nucleo - Revestimiento'});


    % Grafico en colores
x = linspace(-a,a,200);
y = linspace(-pi/2,pi/2,200);
y_degree = linspace(-180,180,length(y));
[X,Y] = meshgrid(x,y);
psi  = @(r,phi) ( psi_r_s(r).*f_s(phi) );
PSI = psi(X,Y);

% f2=figure();
% f2.WindowState = 'maximized';
% subplot 121 ; pcolor(x,y,PSI)
% subplot 122 ; polarPcolor(y,y_degree,PSI)



