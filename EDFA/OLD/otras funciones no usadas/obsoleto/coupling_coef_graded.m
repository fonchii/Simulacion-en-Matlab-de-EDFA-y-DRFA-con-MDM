function Coupling = coupling_coef_graded(fibra,modos,modop,lambda_s,lambda_p)
% NOTA:
% las funciones de bessel se evaluan en la frecuencia normalizada de
% operacion, por lo qe para un lambda dado, el modo 01 y 02 son iguales.

if nargin<1
    fibra.radio=25e-6; fibra.n1=1.46; fibra.n2=1.455; modos="01"; modop= "01";
    lambda_s = 1550e-9 ; lambda_p = 980e-9 ; fibra.Nalpha = 1.95 ; fibra.M = 100;
end
a = fibra.radio;
aa = linspace(0,a,fibra.M);
alpha = fibra.Nalpha ; % α característico del indice de refraccion
delta = ( fibra.n1^2 - fibra.n2^2)/2*fibra.n1^2;

l_s = str2num(extract(modos,1)); m_s = str2num(extract(modos,2));
l_p = str2num(extract(modop,1)); m_p = str2num(extract(modop,2));

for m = 1:fibra.M-1
    if m==1
        n1 = fibra.n1; n2 = n1*sqrt( 1 - 2*delta*( aa(m+1)/a )^alpha );
        vs = (2*pi*a/lambda_s)*sqrt(n1^2-n2^2);
        vp = (2*pi*a/lambda_p)*sqrt(n1^2-n2^2);
        
    else
        n1 = fibra.n1*sqrt( 1 - 2*delta*( aa(m)/a )^alpha );
        n2 = n1*sqrt( 1 - 2*delta*( aa(m+1)/a )^alpha );
        vs = (2*pi*a/lambda_s)*sqrt(n1^2-n2^2);
        vp = (2*pi*a/lambda_p)*sqrt(n1^2-n2^2);
    end
    
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
                    num(m) = integral2( @(r,phi) abs( Ms(r,phi) ).* abs((Mp(r,phi))) , aa(m), aa(m+1), 0, 2*pi).^2;
                    den(m) = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , aa(m), aa(m+1), 0, 2*pi) * integral2( @(r,phi) ( abs(Mp(r,phi)) ).^2 , aa(m), aa(m+1), 0, 2*pi);
                else
                    num(m) = integral2( @(r,phi) abs( Ms(r,phi) ).* abs( Mp_b(r,phi) + Mp_a(r,phi) ) , aa(m), aa(m+1), 0, 2*pi).^2 ;
                    den(m) = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , aa(m), aa(m+1), 0, 2*pi) * integral2( @(r,phi) ( abs( Mp_b(r,phi).^2 + Mp_a(r,phi).^2 ) ) , aa(m), aa(m+1), 0, 2*pi);
                end
            else
                if l_p==0
                    num(m) = integral2( @(r,phi) abs( Ms_b(r,phi) + Ms_a(r,phi) ).* abs( Mp(r,phi) ) , aa(m), aa(m+1), 0, 2*pi).^2 ;
                    den(m) = integral2( @(r,phi) abs( Ms_b(r,phi).^2 + Ms_a(r,phi).^2 ) , aa(m), aa(m+1), 0, 2*pi) * integral2( @(r,phi) ( abs( Mp(r,phi) ) ).^2 , aa(m), aa(m+1), 0, 2*pi);
                else
                    num(m) =  integral2( @(r,phi) abs( Ms_b(r,phi) + Ms_a(r,phi) ).* abs( Mp_b(r,phi) + Mp_a(r,phi) ) , aa(m), aa(m+1), 0, 2*pi).^2 ;
                    den(m) = integral2( @(r,phi) abs( Ms_b(r,phi).^2 + Ms_a(r,phi).^2 ) , aa(m), aa(m+1), 0, 2*pi) * integral2( @(r,phi) ( abs( Mp_b(r,phi).^2 + Mp_a(r,phi).^2 ) ) , aa(m), aa(m+1), 0, 2*pi);
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
                num(m) = integral2( @(r,phi) abs( Ms(r,phi) ).* abs((Mp(r,phi))) , aa(m), aa(m+1), 0, 2*pi).^2;
                den(m) = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , aa(m), aa(m+1), 0, 2*pi) * integral2( @(r,phi) ( abs(Mp(r,phi)) ).^2 , aa(m), aa(m+1), 0, 2*pi);
                
            else
                num(m) = integral2( @(r,phi) abs( Ms_b(r,phi) ).* abs( Mp(r,phi) ) , aa(m), aa(m+1), 0, 2*pi).^2 + integral2( @(r,phi) abs( Ms_a(r,phi) ).* abs( Mp(r,phi) ) , aa(m), aa(m+1), 0, 2*pi).^2 ;
                den(m) = integral2( @(r,phi) abs( Ms_b(r,phi).^2 + Ms_a(r,phi).^2 ) , aa(m), aa(m+1), 0, 2*pi) * integral2( @(r,phi) ( abs( Mp(r,phi) ) ).^2 , aa(m), aa(m+1), 0, 2*pi);
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
                num(m) = integral2( @(r,phi) abs( Ms(r,phi) ).* abs( Mp(r,phi) ) , aa(m), aa(m+1), 0, 2*pi).^2 ;
                den(m) = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , aa(m), aa(m+1), 0, 2*pi) * integral2( @(r,phi) ( abs( Mp(r,phi) ) ).^2 , aa(m), aa(m+1), 0, 2*pi);
            else
                num(m) =  integral2( @(r,phi) abs( Ms(r,phi)).* abs( Mp_b(r,phi) + Mp_a(r,phi) ) , aa(m), aa(m+1), 0, 2*pi).^2 ;
                den(m) = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , aa(m), aa(m+1), 0, 2*pi) * integral2( @(r,phi) ( abs( Mp_b(r,phi).^2 + Mp_a(r,phi).^2 ) ) , aa(m), aa(m+1), 0, 2*pi);
            end
        else
            % bombeo y señal estan tomando ambas polarizaciones
            % ambas pueden(m) ser solo de la forma ij_a o ij_b con i,j!=0
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
            
            num(m) = integral2( @(r,phi) abs( Ms(r,phi) ).* abs((Mp(r,phi))) , aa(m), aa(m+1), 0, 2*pi).^2;
            den(m) = integral2( @(r,phi) abs( Ms(r,phi) ).^2 , aa(m), aa(m+1), 0, 2*pi) * integral2( @(r,phi) ( abs(Mp(r,phi)) ).^2 , aa(m), aa(m+1), 0, 2*pi);
        end
    end
end

Coupling = num./den;