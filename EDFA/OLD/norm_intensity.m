function Gamma_norm_intensity = norm_intensity(fibra,modos,lambda_s)
% NOTA:
% las funciones de bessel se evaluan en la frecuencia normalizada de
% operacion, por lo qe para un lambda dado, el modo 01 y 02 son iguales.

%% Comprobar si es posible transmitir el modo pedido
%vs_sort = [0 , 2.40482554431517,3.83170596389892,5.13562229949996,5.52007810983591,6.38017274704006,7.01558666864552,7.58834243446895,8.41724413965471,8.65372791248466,8.77148381537071];
%LP = ["01","11","21","02","31","12","41","22","03","51"];

if nargin<1
    fibra.radio=25e-6; fibra.n1=1.46; fibra.n2=1.455; 
    modos="01"; lambda_s = 1550e-9 ;
end

% ref: https://www.researchgate.net/publication/319971584_Complete_measurement_of_spatiotemporally_complex_multi-spatial-mode_ultrashort_pulses_from_multimode_optical_fibers_using_delay-scanned_wavelength-multiplexed_holography
vs_sort = [0, 2.048, 3.8317, 3.8317, 5.1356, 5.5201, 6.38, 7.016, 7.016, 7.588, 8.417, 8.654, 8.771, 9.761];
LP = ["01","11","21","02","31","12","41","22","03","51", "32", "13", "61", "42"];
if (extract(modos,strlength(modos))) == extract(modos,2)
    M = modos;
else
    M = split(modos,"_");
    M = M(1);
end
for i=1:length(LP)
    if M == LP(i)
        v_corte = vs_sort(i);
    end
end


a = fibra.radio;
if isfield(fibra,'AN')
    vs = (2*pi*a/lambda_s)*fibra.AN;
else
    n1 = fibra.n1; n2 = fibra.n2;
    vs = (2*pi*a/lambda_s)*sqrt(n1^2-n2^2);
end

if vs<v_corte
    err = strcat( strcat("No es posible transmitir este modo con estos parámetros de fibra, frecuencia normalizada de operacion: ", num2str(vs)) , strcat(" para el modo LP"),modos);
    error(err)
end

l_s = str2num(extract(modos,1)); m_s = str2num(extract(modos,2));

if (extract(modos,strlength(modos))) == extract(modos,2) 
    % Se toman ambas polarizaciones
        if l_s == 0
            Core = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( cos(l_s*phi) )  ;
            Clad = @(r,phi) ( (besselj(l_s, vs) ./ besselk(l_s, vs) )).* ( besselk(l_s, vs*r/a ) ) .* ( cos(l_s*phi) ) ;
            
            Gamma_norm_intensity = ( integral2( @(r,phi) abs( Core(r,phi).*r ) , 0, a, 0, 2*pi) + ...
                integral2( @(r,phi) abs( Clad(r,phi).*r ) , a, inf, 0, 2*pi) ) /(pi*a^2) ;
            
        else 
            Core_a = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( cos(l_s*phi) ) ;
            Core_b = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( sin(l_s*phi) ) ;
            Clad_a = @(r,phi) ( (besselj(l_s, vs) ./ besselk(l_s, vs) )).* ( besselk(l_s, vs*r/a ) ) .* ( cos(l_s*phi) ) ;
            Clad_b = @(r,phi) ( (besselj(l_s, vs) ./ besselk(l_s, vs) )).* ( besselk(l_s, vs*r/a ) ) .* ( sin(l_s*phi) ) ;
            
            Gamma_norm_intensity = ( integral2( @(r,phi) abs( Core_b(r,phi).*r + Core_a(r,phi).*r ) , 0, a, 0, 2*pi) ...
                        + integral2( @(r,phi) abs( Clad_a(r,phi).*r + Clad_b(r,phi).*r ) , a, inf, 0, 2*pi) ) /(pi*a^2) ;

        end
        
else 
    % Se toma solo una polarizacion
    
        pola = extract(modos,strlength(modos));
        
        if pola=="a"
            Core = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( cos(l_s*phi) ) ;
            Clad = @(r,phi) ( (besselj(l_s, vs) ./ besselk(l_s, vs) )).* ( besselk(l_s, vs*r/a ) ) .* ( cos(l_s*phi) ) ;
        elseif pola=="b"
            Core = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( sin(l_s*phi) ) ;
            Clad = @(r,phi) ( (besselj(l_s, vs) ./ besselk(l_s, vs) )).* ( besselk(l_s, vs*r/a ) ) .* ( sin(l_s*phi) ) ;
        end
        Gamma_norm_intensity = ( integral2( @(r,phi) abs( Core(r,phi).*r ) , 0, a, 0, 2*pi) +...
            integral2( @(r,phi) abs( Clad(r,phi).*r ) , a, inf, 0, 2*pi) ) /(pi*a^2) ;

end

%fprintf('Intensidad normalizada para el modo %s : %f\n',modos,Gamma_norm_intensity)




