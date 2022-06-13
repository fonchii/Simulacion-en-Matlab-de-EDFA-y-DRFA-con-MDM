function [Gamma,Beta_0] = norm_intensity2(fibra,modos,lambda_s)
% NOTA:
% las funciones de bessel se evaluan en la frecuencia normalizada de
% operacion, por lo qe para un lambda dado, el modo 01 y 02 son iguales.


if nargin<1
    fibra.radio=25e-6; fibra.n1=1.46; fibra.n2=1.455;
    modos="14"; lambda_s = 1550e-9 ; fibra.Nalpha = inf; fibra.M = 100; %fibra.Nalpha = 1.95
end

%% Comprobar si es posible transmitir el modo pedido
param = modesparam(modos,lambda_s,fibra) ;
if isstruct(param) == 0
    error('No es posible transmitir el modo en la fibra')
end

Beta_0=param.beta*1e6; %[1/m]

a = fibra.radio;
l_s = str2num(extract(modos,1)); m_s = str2num(extract(modos,2));


% Definir intensidad del modo en Core y Clad dependiendo del caso

if (extract(modos,strlength(modos))) == extract(modos,2)
    % Se toman ambas polarizaciones
    if l_s == 0

        Core = @(r,phi) ( besselj(l_s,param.u*r/a) ./besselj(l_s,param.u) ) .* ( cos(l_s*phi) )  ;
        Clad = @(r,phi) ( besselk(l_s, (param.w/a) *r) ./ besselk(l_s, param.w) ) .* ( cos(l_s*phi) ) ;

    else
        Core = @(r,phi) ( besselj(l_s,param.u*r/a) ./besselj(l_s,param.u) ) .* ( cos(l_s*phi) + sin(l_s*phi) )  ;
        Clad = @(r,phi) ( besselk(l_s, (param.w/a) *r) ./ besselk(l_s, param.w) ) .* ( cos(l_s*phi) + sin(l_s*phi) ) ;

    end

else
    % Se toma solo una polarizacion
    pola = extract(modos,strlength(modos));

    if pola=="a"

        Core = @(r,phi) ( besselj(l_s,param.u*r/a) ./besselj(l_s,param.u) ) .* ( cos(l_s*phi) ) ;
        Clad = @(r,phi) ( besselk(l_s, (param.w/a) *r) ./ besselk(l_s, param.w) ).*  ( cos(l_s*phi) ) ;
    elseif pola=="b"

        Core = @(r,phi) ( besselj(l_s,param.u*r/a) ./besselj(l_s,param.u) ) .* ( sin(l_s*phi) ) ;
        Clad = @(r,phi) ( besselk(l_s, (param.w/a) *r) ./ besselk(l_s, param.w) ).*  ( sin(l_s*phi) ) ;

    end

end

den = ( integral2( @(r,phi) abs( Core(r,phi).*r ) , 0, a, 0, 2*pi) + ...
    integral2( @(r,phi) abs( Clad(r,phi).*r ) , a, inf, 0, 2*pi) ) ; %/ (pi*a^2) ;


num = integral2 ( @(r,phi) abs( Core(r,phi).*r ) , 0, a, 0, 2*pi )  ;

Gamma = num / den;




