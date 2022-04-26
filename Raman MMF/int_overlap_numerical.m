function fij = int_overlap_numerical(fibra,modoi,lambdai,modoj,lambdaj)


if nargin<1
    fibra.radio=25e-6; fibra.n1=1.46; fibra.n2=1.450;
    modoi='22_a' ; lambdai = 1550e-9 ; modoj = '02' ; lambdaj = 1550e-9;
end
k1 = extract(modoi,1); k2 = extract(modoj,1);
if k1{1} == 'L'
    if length(modoi) == 4
        modoii = strcat( extract(modoi,3) , extract(modoi,4) );
        modoi = modoii{1};
    elseif length(modoi) == 6
        modoii = strcat( strcat( extract(modoi,3) , extract(modoi,4) ) , strcat(extract(modoi,5) , extract(modoi,6)) );
        modoi = modoii{1};
    end
end
if k2{1} == 'L'
    if length(modoj) == 4
        modojj = strcat( extract(modoj,3) , extract(modoj,4) );
        modoj = modojj{1};
    elseif length(modoj) == 6
        modojj = strcat( strcat( extract(modoj,3) , extract(modoj,4) ) , strcat(extract(modoj,5) , extract(modoj,6)) );
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
%l_i = str2num(extract(modoi,1)); l_j = str2num(extract(modoj,1)); 
l_i = str2num(modoi(1)) ; l_j = str2num(modoj(1));
gridSize = 500;             % tamaño de grilla para calculo numerico
maxPlotRadius = a*2;        % limite de grilla para calculo numerico 

[modeCosi, modeSini, profilei, pixelSizei] = plot_LP_mode(l_i, parami.u, parami.w, a, maxPlotRadius, gridSize) ;
[modeCosj, modeSinj, profilej, pixelSizej] = plot_LP_mode(l_j, paramj.u, paramj.w, a, maxPlotRadius, gridSize) ;

 %Ii = ( abs(modeCosi) + abs(modeSini) ) ; Ij = ( abs(modeCosj) + abs(modeSinj) ) ;


% Definir intensidad del modo en Core y Clad dependiendo del caso
% Funciones de Modo I
%if (extract(modoi,strlength(modoi))) == extract(modoi,2)
if modoi( strlength(modoi) ) == modoi(2)
    % Se toman ambas polarizaciones
    if l_i == 0
        Ii = ( abs(modeCosi) ) ;

    else
        Ii = ( abs(modeCosi) + abs(modeSini) ) ;
        
    end

else
    % Se toma solo una polarizacion
    %pola = extract(modoi,strlength(modoi));
    pola = modoi(strlength(modoi)) ;

    if pola=="a"
        Ii = ( abs(modeCosi) ) ;


    elseif pola=="b"
        Ii = ( abs(modeSini) ) ;

    end

end

% Funciones de Modo J

%if (extract(modoj,strlength(modoj))) == extract(modoj,2)
if modoj(strlength(modoj)) == modoj(2)
    % Se toman ambas polarizaciones
    if l_j == 0
        Ij = ( abs(modeCosj) ) ;

        
    else
        Ij = ( abs(modeCosj) + abs(modeSinj) ) ;

    end

else
    % Se toma solo una polarizacion
    %pola = extract(modoj,strlength(modoj));
    pola = modoj(strlength(modoj));

    if pola=="a"
        Ij = ( abs(modeCosj) ) ;
        
    elseif pola=="b"
        Ij = ( abs(modeSinj) ) ;
        
    end

end

num = sum( sum(Ii.*Ij) ) ; den = a*sum( sum(Ii) ).* sum( sum(Ij) ) ;
fij = num/den;

% % Graficar
% centers = [gridSize/2,gridSize/2] ; radii = gridSize/4; 
% close all; subplot(2,1,1);imagesc(Ii);hold on; viscircles(centers,125,'Color','black') ;
% subplot(2,1,2);imagesc(Ij);hold on; viscircles(centers,radii,'Color','black')


