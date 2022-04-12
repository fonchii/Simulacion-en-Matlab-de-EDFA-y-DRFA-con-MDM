function fij = int_overlap(fibra,modoi,lambdai,modoj,lambdaj)


if nargin<1
    fibra.radio=25e-6; fibra.n1=1.46; fibra.n2=1.455;
    modoi="14"; lambdai = 1550e-9 ; modoj = "21" ; lambdaj = 1550e-9;
end

%% Comprobar si es posible transmitir el modo pedido
parami = modesparam(modoi,lambdai,fibra) ;
paramj = modesparam(modoj,lambdaj,fibra) ;
if isstruct(parami) == 0 || isstruct(paramj) == 0
    error('No es posible transmitir el modo en la fibra')
end


a = fibra.radio;
l_i = str2num(extract(modoi,1)); l_j = str2num(extract(modoj,1)); 


% Definir intensidad del modo en Core y Clad dependiendo del caso
% Funciones de Modo I
if (extract(modoi,strlength(modoi))) == extract(modoi,2)
    % Se toman ambas polarizaciones
    if l_i == 0
        Corei = @(r,phi) ( besselj(l_i,parami.u*r/a) ./besselj(l_i,parami.u) ) .* ( cos(l_i*phi) )  ;
        Cladi = @(r,phi) ( besselk(l_i, (parami.w/a) *r) ./ besselk(l_i, parami.w) ) .* ( cos(l_i*phi) ) ;

    else
        Corei = @(r,phi) ( besselj(l_i,parami.u*r/a) ./besselj(l_i,parami.u) ) .* ( cos(l_i*phi) + sin(l_i*phi) )  ;
        Cladi = @(r,phi) ( besselk(l_i, (parami.w/a) *r) ./ besselk(l_i, parami.w) ) .* ( cos(l_i*phi) + sin(l_i*phi) ) ;

    end

else
    % Se toma solo una polarizacion
    pola = extract(modoi,strlength(modoi));

    if pola=="a"

        Corei = @(r,phi) ( besselj(l_i,parami.u*r/a) ./besselj(l_i,parami.u) ) .* ( cos(l_i*phi) ) ;
        Cladi = @(r,phi) ( besselk(l_i, (parami.w/a) *r) ./ besselk(l_i, parami.w) ).*  ( cos(l_i*phi) ) ;
    elseif pola=="b"

        Corei = @(r,phi) ( besselj(l_i,parami.u*r/a) ./besselj(l_i,parami.u) ) .* ( sin(l_i*phi) ) ;
        Cladi = @(r,phi) ( besselk(l_i, (parami.w/a) *r) ./ besselk(l_i, parami.w) ).*  ( sin(l_i*phi) ) ;

    end

end

% Funciones de Modo J

if (extract(modoj,strlength(modoj))) == extract(modoj,2)
    % Se toman ambas polarizaciones
    if l_j == 0
        Corej = @(r,phi) ( besselj(l_j,paramj.u*r/a) ./besselj(l_j,paramj.u) ) .* ( cos(l_j*phi) )  ;
        Cladj = @(r,phi) ( besselk(l_j, (paramj.w/a) *r) ./ besselk(l_j, paramj.w) ) .* ( cos(l_j*phi) ) ;

    else
        Corej = @(r,phi) ( besselj(l_j,paramj.u*r/a) ./besselj(l_j,paramj.u) ) .* ( cos(l_j*phi) + sin(l_j*phi) )  ;
        Cladj = @(r,phi) ( besselk(l_j, (paramj.w/a) *r) ./ besselk(l_j, paramj.w) ) .* ( cos(l_j*phi) + sin(l_j*phi) ) ;

    end

else
    % Se toma solo una polarizacion
    pola = extract(modoi,strlength(modoi));

    if pola=="a"

        Corej = @(r,phi) ( besselj(l_j,paramj.u*r/a) ./besselj(l_j,paramj.u) ) .* ( cos(l_j*phi) ) ;
        Cladj = @(r,phi) ( besselk(l_j, (paramj.w/a) *r) ./ besselk(l_j, paramj.w) ).*  ( cos(l_j*phi) ) ;
    elseif pola=="b"

        Corej = @(r,phi) ( besselj(l_j,paramj.u*r/a) ./besselj(l_j,paramj.u) ) .* ( sin(l_j*phi) ) ;
        Cladj = @(r,phi) ( besselk(l_j, (paramj.w/a) *r) ./ besselk(l_j, paramj.w) ).*  ( sin(l_j*phi) ) ;

    end

end

Ii = @(r,phi) ( Corei(r,phi) ) ; % Ii = @(r,phi) ( Cladi(r,phi) + Corei(r,phi) ) ; % Tomar en cuenta core y clad
Ij = @(r,phi) ( Corej(r,phi) ) ; % Ij = @(r,phi) ( Cladj(r,phi) + Corej(r,phi) ) ; % Tomar en cuenta core y clad

den = ( integral2( @(r,phi) abs( Ii(r,phi).*r ) , 0,a, 0,2*pi ) * ...
    integral2( @(r,phi) abs( Ij(r,phi).*r ) , 0,a, 0,2*pi ) ) 


num = integral2 ( @(r,phi) abs( (Ii(r,phi).*r).*(Ij(r,phi).*r) ) , 0,a, 0,2*pi )  

fij = num / den;