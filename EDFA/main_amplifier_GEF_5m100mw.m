% Simulador EDFA MM
close all; 
clear all; clc

%% Parámetros de entrada

    % Señal : Modos y Canales
NCh = 31;

Signal.modos = ["01" "11" ] ;

Pin = -20;
Signal.lambda.LP_01     = linspace(1530e-9,1560e-9,NCh);        P0_signal.LP_01     = Pin*ones(1,length(Signal.lambda.LP_01));
Signal.lambda.LP_11     = linspace(1530e-9,1560e-9,NCh);        P0_signal.LP_11     = Pin*ones(1,length(Signal.lambda.LP_11));


    % Bombeo : Modos y Canales
Pump.modos = ["01" ]   ;

Pump.lambda.LP_01   = 980e-9;                                 P0_pump.LP_01   = 100e-3  ;  


    % POTENCIAS

for i=1:length(Signal.modos)        % Potencia: W -> mW
    for j=1:length( P0_signal.(strcat("LP_",Signal.modos(i))) )
        P0_signal.( strcat("LP_",Signal.modos(i) ) )(j) = 1e-3*10^(P0_signal.( strcat("LP_",Signal.modos(i) ) )(j)/10);
    end
end ;clear i j;
Signal.P0 = P0_signal; 
Pump.P0 = P0_pump;
ASE = -200;                                                  %dBm  -50
Signal.NumberOfChannels=NCh;

    % Datos de la fibra

Fibra.nucleos = 1;                                           % Numero de nucleos
Fibra.largo = 5     ; Fibra.radio = 5e-6   ; Fibra.N = 7e24; % fibra.N = 3e24; 

Fibra.dvk=300e9;
Signal.NumberOfChannels=NCh;

Fibra.CrossSectionParameters = "OptiSystem"; % "OptiSystem" , "VPI"
Fibra.n1 = 1.45 ;   
Fibra.n2 = 1.4354 ;
%Fibra.dvk= P.OpticalBW; % diferencia : max_lambda - min_lambda 

Fibra.WaitBar = 1; Fibra.Avance = 1;    % Despliegue de info
Fibra.ASEFlag = 1;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)


%%
tic;
% ---------- Filtrado de equalizacion de ganancias ---------- %

% v3 : mas rapido
%EDFA = EDFA_MM_GEF_Calcv3(Fibra,Signal,Pump,ASE); 

% v4 : lento, puede lograr mejores resultados - recomendado usar v3 primero
% EDFA = EDFA_MM_GEF_Calcv4(Fibra,Signal,Pump,ASE); 

Fibra.GEF = load('GEF_Filters/EDFA_5m_100mw_OptiSystem') ;
EDFA = EDFA_MM_GEF_Use(Fibra,Signal,Pump,ASE)  ;
t_end = toc; fprintf('Tiempo de cómputo: %.2f segundos\n', t_end);

% --------- SaveFilter ---------- %
%Filtro = EDFA.Nucleo1.GEF.best_weight_Function;

    %% Graficos

close all
graf.Nc = length(fieldnames(EDFA)); 
graf.z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';


% % --------- Just Gain --------- %
% 
%load("EDFA_Gain_OptiSystem_SinGEF.mat")
%load("EDFA_Gain_VPI_SinGEF.mat")

for s = 1:length(Signal.modos)
    figure(1)
    graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.Nucleo1.salida.ganancias.(strcat("LP_",Signal.modos(s)));
    leyenda = strcat((strcat("LP",Signal.modos(s))));
    
    ejex = Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9;
    plot(ejex,graf.ganancias.(strcat("LP_",Signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
    xlabel('Longitud de onda [nm]'); ylabel(ylab) %xlabel('λ [nm]') 
    
    legend('Location','southoutside','Box','off','NumColumns',2)
    ylim([14.5 19.5])
    hold on
    set(gca,"ColorOrderIndex",s)
    plot(ejex,EDFA.Nucleo1.salida.ganancias_sinGEF.(strcat("LP_",Signal.modos(s))) , '--*' , 'DisplayName',strcat(leyenda,' sin GEF' )   )
end ; clear s ejex leyenda;
dim = [.132 .63 .4 .3]; str = strcat("Ripple: ",num2str(round(EDFA.Nucleo1.GEF.Ripple,2)) , ' dB en Modo LP01');
annotation('textbox',dim,'String',str,'FitBoxToText','on');

% % % ----- Save ----- %
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'GananciaEqualizada_5m-100mw'  

% figure(2)
% plot(Signal.lambda.(strcat("LP_",Signal.modos(1))).*1e9 , EDFA.Nucleo1.GEF.best_weight_Function ) 
% title('Filtro utilizado') ; xlabel('Longitud de onda [nm]') ; ylabel('Magnitud')
% 
% % % %----- Save ----- %
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'FiltroEqualizador_5m-100mw' 


% plot(EDFA.Nucleo1.Pase.LP_01(2,:)) ; ylim([-60 -20])
% plot(EDFA.Nucleo1.salida.OSNR.LP_01(2,:)) ; ylim([25 50])

%% Graficos 3D de Modos 

% % ----- Señal  ----- %
% figure()
% graficar_modos(Fibra,Signal.modos,Signal.lambda.LP_01)
% sgtitle('Modos de Señal') 

% % ----- Bombeo ----- %
% figure()
% graficar_modos(Fibra,Pump.modos,Pump.lambda.LP_01)
% sgtitle('Modos de Bombeo') ; 

%% Graficar potencia vs frecuencia
% figure()
% xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
% graf.Nc = length(fieldnames(EDFA)); 
% for n = 1:graf.Nc
%     for s = 1:length(Signal.modos)
%         graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
%         leyenda = strcat("LP",Signal.modos(s));
%         ejex = fliplr((3*10^8./(Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9)));
%         plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",Signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%         %ax.ColorOrderIndex = s;
%         xlabel('Frecuencia [THz]') ; ylabel(ylab)
% 
%     end ; clear s ejex leyenda n;
% end


%% Graficar potencia vs wavelenght
% n=n+1; figure()
% xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
% graf.Nc = length(fieldnames(EDFA)); 
% for n = 1:graf.Nc
%    for s = 1:length(Signal.modos)
%        graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
%        leyenda = strcat((strcat("LP",Signal.modos(s))));
%        ejex = fliplr((3*10^8./(Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9)));
%        plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",Signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%        %ax.ColorOrderIndex = s;
%        xlabel('Frecuencia [THz]') ; ylabel(ylab)
%        set(gca,'xdir','reverse')
%    end ; clear s ejex leyenda;
% end

