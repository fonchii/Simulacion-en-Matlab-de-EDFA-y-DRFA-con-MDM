function Coupling = coupling_coef(fibra,modos,modop,lambda_s,lambda_p)
% NOTA:
% las funciones de bessel se evaluan en la frecuencia normalizada de
% operacion, por lo qe para un lambda dado, el modo 01 y 02 son iguales.

if nargin<1
    fibra.radio=25e-6; fibra.n1=1.46; fibra.n2=1.455; modos="01"; modop= "01";
    lambda_s = 1550e-9 ; lambda_p = 980e-9;
end
a = fibra.radio;
if isfield(fibra,'AN')
    vs = (pi*a/lambda_s)*fibra.AN;
    vp = (pi*a/lambda_p)*fibra.AN;
else
    n1 = fibra.n1; n2 = fibra.n2; 
    vs = (pi*a/lambda_s)*sqrt(n1^2-n2^2);
    vp = (pi*a/lambda_p)*sqrt(n1^2-n2^2);
end


l_s = str2num(extract(modos,1)); m_s = str2num(extract(modos,2));
l_p = str2num(extract(modop,1)); m_p = str2num(extract(modop,2));

if (extract(modos,strlength(modos)))==extract(modos,2)
    if (extract(modop,strlength(modop)))==extract(modop,2)
        % se estan tomando ambas polarizaciones (a y b) para señal y bombeo
        if l_s == 0
            Ms = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( cos(l_s*phi) ) ;
        else 
            Ms_a = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( cos(l_s*phi) ) ;
            Ms_b = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( sin(l_s*phi) ) ;
        end

        if l_p == 0
            Mp = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( cos(l_p*phi) ) ;
        else 
            Mp_a = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( cos(l_p*phi) ) ;
            Mp_b = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( sin(l_p*phi) ) ;
        end

        if l_s==0 
            if l_p==0
                num = integral2( @(r,phi) abs( Ms(r,phi) ).* abs((Mp(r,phi))) , 0, a, 0, 2*pi).^2;
                den = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs(Mp(r,phi)) ).^2 , 0, a, 0, 2*pi);
            else
                num = integral2( @(r,phi) abs( Ms(r,phi) ).* abs( Mp_b(r,phi) + Mp_a(r,phi) ) , 0, a, 0, 2*pi).^2 ;
                den = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs( Mp_b(r,phi).^2 + Mp_a(r,phi).^2 ) ) , 0, a, 0, 2*pi);
            end
        else 
            if l_p==0
                num = integral2( @(r,phi) abs( Ms_b(r,phi) + Ms_a(r,phi) ).* abs( Mp(r,phi) ) , 0, a, 0, 2*pi).^2 ;
                den = integral2( @(r,phi) abs( Ms_b(r,phi).^2 + Ms_a(r,phi).^2 ) , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs( Mp(r,phi) ) ).^2 , 0, a, 0, 2*pi);
            else
                num =  integral2( @(r,phi) abs( Ms_b(r,phi) + Ms_a(r,phi) ).* abs( Mp_b(r,phi) + Mp_a(r,phi) ) , 0, a, 0, 2*pi).^2 ;
                den = integral2( @(r,phi) abs( Ms_b(r,phi).^2 + Ms_a(r,phi).^2 ) , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs( Mp_b(r,phi).^2 + Mp_a(r,phi).^2 ) ) , 0, a, 0, 2*pi);
            end
        end
        
    else
        % bombeo esta tomando solo una polarización, señal ambas
        % señal puede ser 0j o ij; bombeo solo ij_a o ij_b. con i,j!=0 
        if l_s == 0
            Ms = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( cos(l_s*phi) ) ;
        else 
            Ms_a = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( cos(l_s*phi) ) ;
            Ms_b = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( sin(l_s*phi) ) ;
        end

        pola = extract(modop,strlength(modop));
        
        if pola=="a"
            Mp = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( cos(l_p*phi) ) ;
        elseif pola=="b"
            Mp = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( sin(l_p*phi) ) ;
        end

        if l_s==0 
            num = integral2( @(r,phi) abs( Ms(r,phi) ).* abs((Mp(r,phi))) , 0, a, 0, 2*pi).^2;
            den = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs(Mp(r,phi)) ).^2 , 0, a, 0, 2*pi);

        else 
            num = integral2( @(r,phi) abs( Ms_b(r,phi) + Ms_a(r,phi) ).* abs( Mp(r,phi) ) , 0, a, 0, 2*pi).^2 ;
            den = integral2( @(r,phi) abs( Ms_b(r,phi).^2 + Ms_a(r,phi).^2 ) , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs( Mp(r,phi) ) ).^2 , 0, a, 0, 2*pi);
        end
    end
    
else
    % señal esta tomando solo una polarización, 
    % señal puede ser solo ij_a o ij_b con i,j!=0
    pola_s = extract(modos,strlength(modos));
    
    if (extract(modop,strlength(modop)))==extract(modop,2)
        % bombeo esta tomando ambas polarizaciones
        if pola_s == "a"
            Ms = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( cos(l_s*phi) ) ;
        elseif pola_s == "b" 
            Ms = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( sin(l_s*phi) ) ;
        end

        if l_p == 0
            Mp = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( cos(l_p*phi) ) ;
        else 
            Mp_a = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( cos(l_p*phi) ) ;
            Mp_b = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( sin(l_p*phi) ) ;
        end

        if l_p==0
                num = integral2( @(r,phi) abs( Ms(r,phi) ).* abs( Mp(r,phi) ) , 0, a, 0, 2*pi).^2 ;
                den = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs( Mp(r,phi) ) ).^2 , 0, a, 0, 2*pi);
        else
                num =  integral2( @(r,phi) abs( Ms(r,phi)).* abs( Mp_b(r,phi) + Mp_a(r,phi) ) , 0, a, 0, 2*pi).^2 ;
                den = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs( Mp_b(r,phi).^2 + Mp_a(r,phi).^2 ) ) , 0, a, 0, 2*pi);
        end
    else
        % bombeo y señal estan tomando ambas polarizaciones
        % ambas pueden ser solo de la forma ij_a o ij_b con i,j!=0
        pola = extract(modop,strlength(modop));
        pola_s = extract(modos,strlength(modos));
        
        if pola=="a"
            Mp = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( cos(l_p*phi) ) ;
        elseif pola=="b"
            Mp = @(r,phi) ( besselj(l_p,vp*r/a) ) .* ( sin(l_p*phi) ) ;
        end
        
        if pola_s == "a"
            Ms = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( cos(l_s*phi) ) ;
        elseif pola_s == "b" 
            Ms = @(r,phi) ( besselj(l_s,vs*r/a) ) .* ( sin(l_s*phi) ) ;
        end
        
        num = integral2( @(r,phi) abs( Ms(r,phi) ).* abs((Mp(r,phi))) , 0, a, 0, 2*pi).^2;
        den = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , 0, a, 0, 2*pi) * integral2( @(r,phi) ( abs(Mp(r,phi)) ).^2 , 0, a, 0, 2*pi);
    end
end

Coupling = num/den;