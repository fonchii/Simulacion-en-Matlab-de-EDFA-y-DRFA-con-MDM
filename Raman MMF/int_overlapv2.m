function fij = int_overlapv2(fibra,modoi,lambdai,modoj,lambdaj)


if nargin<1
    fibra.radio=25e-6; fibra.n1=1.46; fibra.n2=1.455;
    modoi='14'; lambdai = 1550e-9 ; modoj = '21' ; lambdaj = 1550e-9;
end

%% COMPROBACION ENTRADA DE MODOS
k1 = extract(modoi,1); k2 = extract(modoj,1);
if k1{1} == 'L'
    if length(modoi) == 4
        modoii = strcat( extract(modoi,3) , extract(modoi,4) );
        modoi = modoii{1};
    elseif length(modoi) == 5
        modoii = strcat( strcat( extract(modoi,3) , extract(modoi,4) , extract(modoi,5)  ) );
        modoi = modoii{1};
    end
end
if k2{1} == 'L'
    if length(modoj) == 4
        modojj = strcat( extract(modoj,3) , extract(modoj,4) );
        modoj = modojj{1};
    elseif length(modoj) == 5
        modojj = strcat( strcat( extract(modoj,3) , extract(modoj,4) , extract(modoj,5) ) );
        modoj = modojj{1};
    end
end

%% Comprobar si es posible transmitir el modo pedido
parami = modesparam(modoi,lambdai,fibra) ;
paramj = modesparam(modoj,lambdaj,fibra) ;
if isstruct(parami) == 0 || isstruct(paramj) == 0
    error('No es posible transmitir el modo en la fibra')
end


a = fibra.radio;
l_i = str2num(modoi(1)) ; l_j = str2num(modoj(1)) ; %l_i = extract(modoi,1); l_j = extract(modoj,1); 

Corei = @(r) besselj(l_i,parami.u*r/a).^2 ;
Cladi = @(r) ( (besselj(l_i,(parami.u)) / besselk(l_i,parami.w) )^2 ) .* besselk(l_i,parami.w*r/a).^2 ;
%Cladi = @(r) besselk(l_i,parami.w*r/a).^2 ;

Corej = @(r) besselj(l_j,paramj.u*r/a).^2 ;
Cladj = @(r) ( (besselj(l_j,paramj.u)/besselk(l_j,paramj.w))^2 ) .* besselk(l_j,paramj.w*r/a).^2 ;
%Cladj = @(r) besselk(l_j,paramj.w*r/a).^2 ;


% Definir polarizacion de los modos dependiendo del caso

% Funciones de Modo i
if modoi(end) == modoi(2)   %if (extract(modoi,strlength(modoi))) == extract(modoi,2)
    % Se toman ambas polarizaciones
    if l_i == 0
        fi = @(phi) ( cos(l_i*phi) )  ;

    else
        fi = @(phi) ( cos(l_i*phi) + sin(l_i*phi) )  ;
    end

else
    % Se toma solo una polarizacion
    pola = modoi(end); %pola = extract(modoi,strlength(modoi));

    if pola=="a"
        fi = @(phi) ( cos(l_i*phi) ) ;

    elseif pola=="b"
        fi = @(phi)  ( sin(l_i*phi) ) ;
    end
end

% Funciones de Modo j

if modoj(end) == modoj(2)%if (extract(modoj,strlength(modoj))) == extract(modoj,2)
    % Se toman ambas polarizaciones
    if l_j == 0
        fj = @(phi) ( cos(l_j*phi) )  ;

    else
        fj = @(phi) ( cos(l_j*phi) + sin(l_j*phi) )  ;
    end

else
    % Se toma solo una polarizacion
    pola = modoj(end) ; %pola = extract(modoi,strlength(modoi));

    if pola=="a"
        fj = @(phi) ( cos(l_j*phi) ) ;

    elseif pola=="b"
        fj = @(phi) ( sin(l_j*phi) ) ;
    end
end

den = integral2( @(r,phi) ( Corei(r).*r .*abs(fi(phi)) ) , 0,a, 0,2*pi ) * ...
    integral2( @(r,phi) ( Corej(r).*r .*abs(fj(phi)) ) , 0,a, 0,2*pi  )  +...
    integral2( @(r,phi) ( Cladi(r).*r .*abs(fi(phi)) ) , a,5*a, 0,2*pi )*...
    integral2( @(r,phi) ( Cladj(r).*r .*abs(fi(phi)) ) , a,5*a, 0,2*pi ); 

num = ( integral2( @(r,phi) ( Corei(r) .*abs(fi(phi)) )  .* ( Corej(r).*abs(fj(phi)) ) .*r , 0,a, 0,2*pi  ) )*pi*a^2  ;

fij = num / den;

% % Graficar
% r = linspace(0, a, 400); phi = linspace(0, 2*pi, 400); [r,phi] = meshgrid(r, phi);
% x = r.*cos(phi); y=r.*sin(phi); zi = Corei(r).*fi(phi) ; zj = Corej(r).*fj(phi) ; 
% subplot(1,2,1) ; mesh(x,y,abs(zi)); view(2); subplot(1,2,2) ; mesh(x,y,abs(zj)); view(2);
