function graf_modos = graficar_modos(fibra,modos)
% Respecto a frecuencia de corte
a = fibra.radio ; n1 = fibra.n1 ; n2 = fibra.n2 ; lambda=980e-9;

r = linspace(0, a, 400);
phi = linspace(0, 2*pi, 400);
[r,phi] = meshgrid(r, phi);
x = r.*cos(phi); y=r.*sin(phi);

vs_sort = [2.40482554431517,3.83170596389892,5.13562229949996,5.52007810983591,6.38017274704006,7.01558666864552,7.58834243446895,8.41724413965471,8.65372791248466,8.77148381537071];
LP = ["01","11","21","02","31","12","41","22","03","51"];

Nmodos = 0;
for i = 1:length(modos)
    if (extract(modos(i),strlength(modos(i))))==extract(modos(i),2) % se toman ambas polarizaciones
        if (str2num(extract(modos(i),1)) == 0) 
            Nmodos = Nmodos + 1;
        else
            Nmodos = Nmodos+2; % Los modos con dos polarizaciones suman 2
        end
    else
        Nmodos = Nmodos + 1 ; % Todos suman 1
    end
end
        
Nmodos = round((Nmodos+1)/2);
grafico=1;            

for i=1:length(modos)
    if (extract(modos(i),strlength(modos(i))))==extract(modos(i),2) % se toman ambas polarizaciones
        l = str2num(extract(modos(i),1));
        for j = 1:length(LP)
            if modos(i) == LP(j)
                v = vs_sort(j);
            end
        end
        if l == 0
            z = besselj(l,v*r/a).*cos(l.*phi) / (pi*a^2);
            subplot(2,(Nmodos),grafico); grafico = grafico + 1;
            mesh(x,y,z.^2); view(2); title(strcat('Modo LP',modos(i)));
        else
            z_a = besselj(l,v*r/a).*cos(l.*phi) / (pi*a^2) ;
            z_b = besselj(l,v*r/a).*sin(l.*phi) / (pi*a^2) ;
            %subplot(2,(round(length(modos)+1))/2,i);
            subplot(2,(Nmodos),grafico); grafico = grafico + 1;
            mesh(x,y,z_a.^2); view(2); title(strcat(strcat('Modo LP',modos(i)),'a'));
            %subplot(2,(round(length(modos)+1))/2,i);
            subplot(2,(Nmodos),grafico); grafico = grafico + 1;
            mesh(x,y,z_b.^2); view(2); title(strcat(strcat('Modo LP',modos(i)),'b'));
        end
    else % se toma solo una polarización, lp!=0
        L = split(modos(i),"_");
        l = str2num(extract(L(1),1));
        for j = 1:length(LP)
            if L(1) == LP(j)
                v = vs_sort(j)
            end
        end
        pola = extract(modos,strlength(modos));
        if pola=="a"
            z = ( besselj(l,v*r/a) ) .* ( cos(l*phi) ) ;
        elseif pola=="b"
            z = ( besselj(l,v*r/a) ) .* ( sin(l*phi) ) ;
        end
        %subplot(2,(round(length(modos)+1))/2,i);
        subplot(2,(Nmodos),grafico); grafico = grafico + 1;
        mesh(x,y,z.^2); view(2); title(strcat('Modo LP',modos(i)));
    end
end
        