%% Indicaciones:
% Las cantidad de canales por modo deben ser iguales

close all; clear all; clc

    % Modos
signal.modos = ["01" "11"] ;
signal.lambda.LP_01 = [1530e-9];                            P0_signal.LP_01 = [-15];
signal.lambda.LP_11_a = [1550e-9];                          P0_signal.LP_11_a = [-15];
signal.lambda.LP_11_b = [1550e-9];                          P0_signal.LP_11_b = [-15];
signal.lambda.LP_11 = [1530e-9 ];                           P0_signal.LP_11 = [-15];              
signal.lambda.LP_02 = [1540e-9 1560e-9 1570e-9];            P0_signal.LP_01 = [-15 -15 -15];


%pump.modos = ["01","11_a","11_b"] ;
pump.modos = ["01"] ;
pump.lambda.LP_01 = [980e-9];                           %P0_pump.LP_01 = [250e-3]  ;  %% Comentar para iterar por Potencias
pump.lambda.LP_11_a = [980e-9];                         %P0_pump.LP_11_a = [250e-3 5e-3] ;  
pump.lambda.LP_11_b = [982e-9];                         %P0_pump.LP_11_b = [250e-3 5e-3]  ;
pump.lambda.LP_11 = [980e-9];                           %P0_pump.LP_11 = [250e-3]  ;
pump.lambda.LP_02 = [980e-9];                           %P0_pump.LP_02 = [250e-3]  ;
pump.lambda.LP_21 = [980e-9];                           %P0_pump.LP_21 = [250e-3]  ;
    %POTENCIAS
for i=1:length(signal.modos)        % Potencia de señal a mW
    for j=1:length( P0_signal.(strcat("LP_",signal.modos(i))) )
        P0_signal.( strcat("LP_",signal.modos(i) ) )(j) = 1e-3*10^(P0_signal.( strcat("LP_",signal.modos(i) ) )(j)/10);
    end
end ;clear i j;
signal.P0 = P0_signal; 
Pase = -50;                                                  %dBm

% pump.P0 = P0_pump;% COMENTAR PARA ITERAR POR POTENCIAS

fibra.largo = 30;  % COMENTAR PARA ITERAR POR LARGOS

    % Datos de la fibra
fibra.nucleos = 1;  % Numero de nucleos
fibra.radio = 8e-6 ; fibra.N = 1e24; 
%fibra.AN = 0.1 ;   % Se prioriza el uso de AN en caso de entregar n1,n2,AN
fibra.n1 = 1.47; fibra.n2 = 1.46;
%fibra.n1 = 1.47 ; fibra.n2 = 1.44; % PAPER
fibra.M = 10; fibra.Nalpha = inf;
fibra.ASEFlag = 1; % EVITA CALCULO DE ESPECTRO ASE

% Iteraciones por potencia
bombeos = [];
for i=1:5
    P0_pump = (i)*8e-3  ;                                          %mW
    pump.P0.(strcat("LP_",pump.modos(1))) = P0_pump;
    tic; EDFA.(strcat("bombeo",int2str(50*i))) = EDFA_MM_radialv2(fibra,signal,pump,Pase); t = toc;
    clc;fprintf('Iteración %i terminada en %.2f segundos',i,t)
    bombeos = [bombeos P0_pump];
end

% Iteraciones por largo de fibra
% largos = [];
% for i=1:10
%     fibra.largo = i*10  ;
%     tic; EDFA.(strcat("Largo",int2str(10*i))) = EDFA_MMv2(fibra,signal,pump,Pase); t = toc;
%     clc;fprintf('Iteración %i terminada en %.2f segundos',i,t)
%     largos = [largos i*10];
% end

%%

% Ganancias
close all; ylab = 'Potencia [dBm]';
figure()
% Ganancia vs Potencia bombeo

for s = 1:length(signal.modos)      
    G = [];
    for i = 1:length(bombeos)
        G = [G EDFA.(strcat("bombeo",int2str(50*i))).Nucleo1.salida.ganancias.(strcat("LP_",signal.modos(s)))(1,1)];
    end
    graf.ganancias.(strcat("LP_",signal.modos(s))) = G;
    ejex = bombeos.*1e3;
    leyenda = strcat( "Modo "  , strcat( "LP",signal.modos(s) ) ); titulo = strcat( "Ganancia vs Potencia de Bombeo" );
    plot(ejex , graf.ganancias.(strcat("LP_",signal.modos(s))) , '-o', 'DisplayName', leyenda) ; hold on ; title(titulo) ; legend(); grid on ;
    xlabel('Potencia de bombeo [mW]') ; ylabel(ylab); %ylim([-10 30]);
end

% Ganancia vs Largo del amplificador, 

% for s = 1:length(signal.modos)      
%     G = [];
%     for i = 1:length(largos)
%         G = [G EDFA.(strcat("Largo",int2str(10*i))).Nucleo1.salida.ganancias(1,1,s)];
%     end
%     graf.ganancias.(strcat("LP_",signal.modos(s))) = G;
%     ejex = largos;
%     leyenda = strcat( "Modo "  , strcat( "LP",signal.modos(s) ) ); titulo = strcat( "Ganancia vs Largo del Amplificador" );
%     plot(ejex , graf.ganancias.(strcat("LP_",signal.modos(s))) , '-o', 'DisplayName', leyenda) ; hold on ; title(titulo) ; legend(); grid on ;
%     xlabel('Largo del Amplificador [m]') ; ylabel(ylab);
% end






% Ganancia vs Lambda, usar varios canales por modo

% for i = 1:10     
%     for s = 1:length(signal.modos)
%         graf.ganancias.(strcat("LP_",signal.modos(s))).(strcat("bombeo",int2str(50*i))) = EDFA.(strcat("bombeo",int2str(50*i))).Nucleo1.salida.ganancias(:,:,s);
%         ejex = signal.lambda.LP_01.*1e9;
%         leyenda = strcat( strcat(" bombeo " , int2str(50*i))  , "mW"); titulo = strcat( strcat("LP", signal.modos(s) ) );
%         ax = subplot (1,length(signal.modos),s);
%         plot(ejex , graf.ganancias.(strcat("LP_",signal.modos(s))).(strcat("bombeo",int2str(50*i))) , '-o', 'DisplayName', leyenda) ; hold on ; title(titulo) ; legend(); grid on ; 
%         xlabel('Potencia de bombeo [mW]') ; ylabel(ylab);
%     end
%     sgtitle('Ganancias')
% end


