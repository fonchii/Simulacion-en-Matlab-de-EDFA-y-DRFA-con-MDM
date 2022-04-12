function graf_modos = graficar_modos(fibra,modos,lambdas)
% Modos se grafican respecto a frecuencia de operación
a = fibra.radio ;
r = linspace(0, a, 400);
phi = linspace(0, 2*pi, 400);
[r,phi] = meshgrid(r, phi);
x = r.*cos(phi); y=r.*sin(phi);


Nmodos = 0;
parfor i = 1:length(modos)
    if (extract(modos(i),strlength(modos(i))))==extract(modos(i),2) % se toman ambas polarizaciones
        if (str2num(extract(modos(i),1)) == 0) 
            Nmodos = Nmodos + 1;
        else
            %Nmodos = Nmodos+2; % Los modos con dos polarizaciones suman 2
            Nmodos = Nmodos + 1; % Se grafican en un mismo gráfico 
        end
    else % Se toma una de las dos polarizaciones
        Nmodos = Nmodos + 1 ; 
    end
end
gvert = 1;
if Nmodos>1
    gvert = 2;
end
if mod(Nmodos,2)==0
    Nmodos = Nmodos/2;
else
    Nmodos = round((Nmodos+1)/2);
end
grafico=1;


for i=1:length(modos)
    %v=(pi*fibra.radio/lambdas(1))*sqrt(fibra.n1^2-fibra.n2^2);
    if isfield(fibra,'AN')
        v = (pi*fibra.radio/lambdas(1))*fibra.AN;
    else
        n1 = fibra.n1; n2 = fibra.n2;
        v = (pi*fibra.radio/lambdas(1))*sqrt(n1^2-n2^2);
    end
    
    if (extract(modos(i),strlength(modos(i))))==extract(modos(i),2) % se toman ambas polarizaciones
        l = str2num(extract(modos(i),1));        
        if l == 0
            z = besselj(l,v*r/a).*cos(l.*phi) / (pi*a^2);
            subplot(gvert,(Nmodos),grafico); grafico = grafico + 1;
            mesh(x,y,z.^2); view(2); title(strcat('Modo LP',modos(i)));
        else
            z_a = besselj(l,v*r/a).*cos(l.*phi) / (pi*a^2) ;
            z_b = besselj(l,v*r/a).*sin(l.*phi) / (pi*a^2) ;
            subplot(gvert,(Nmodos),grafico); grafico = grafico + 1;
            mesh(x,y,z_a.^2); view(2); title(strcat('Modo LP',modos(i))); hold on
            %subplot(2,(Nmodos),grafico); grafico = grafico + 1;
            mesh(x,y,z_b.^2); view(2); %title(strcat(strcat('Modo LP',modos(i)),'b'));
            
%             z_a = besselj(l,v*r/a).*cos(l.*phi) / (pi*a^2) ;
%             z_b = besselj(l,v*r/a).*sin(l.*phi) / (pi*a^2) ;
%             subplot(2,(Nmodos),grafico); grafico = grafico + 1;
%             mesh(x,y,(z_a.^2 + z_b.^2)); view(2); title(strcat('Modo LP',modos(i)));
        end
    else % se toma solo una polarización, lp!=0
        L = split(modos(i),"_");
        l = str2num(extract(L(1),1));
        pola = L(2);
        if pola=="a"
            z = ( besselj(l,v*r/a) ) .* ( cos(l*phi) ) ;
        elseif pola=="b"
            z = ( besselj(l,v*r/a) ) .* ( sin(l*phi) ) ;
        end
        subplot(gvert,(Nmodos),grafico); grafico = grafico + 1;
        mesh(x,y,z.^2); view(2); title(strcat('Modo LP',modos(i)));
    end
end
