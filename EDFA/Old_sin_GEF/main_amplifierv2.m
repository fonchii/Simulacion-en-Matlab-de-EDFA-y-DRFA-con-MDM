% Simulador EDFA MM
close all; 
clear all; clc

%% Parámetros de entrada

    % Señal : Modos y Canales
NCh = 31;

Signal.modos = ["01" "11" ] ;

Signal.lambda.LP_01     = linspace(1530e-9,1560e-9,NCh);        P0_signal.LP_01     = -15*ones(1,length(Signal.lambda.LP_01));
Signal.lambda.LP_11     = linspace(1530e-9,1560e-9,NCh);        P0_signal.LP_11     = -15*ones(1,length(Signal.lambda.LP_11));


    % Bombeo : Modos y Canales
Pump.modos = ["01" ]   ;

Pump.lambda.LP_01   = 980e-9;                                 P0_pump.LP_01   = 300e-3  ;  


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
%EDFA = EDFA_MM(fibra,signal,pump,ASE);         % Sin efecto acomplamiento de Potencia intermodal

%EDFA = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);      % Con efecto acomplamiento de Potencia intermodal

%EDFAVPI= EDFA_MMvpi2(Fibra,Signal,Pump,ASE);

EDFA = EDFA_MM_GEF_Calcv4(Fibra,Signal,Pump,ASE); % Filtrado de equalizacion de ganancias
%Original_EDFA = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);   
t_end = toc; fprintf('Tiempo de cómputo: %.2f segundos\n', t_end);


    %% Graficos

close all
graf.Nc = length(fieldnames(EDFA)); 
graf.z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';


% for n = 1:graf.Nc
%     figure(n)
%     
%     % Señal
%     LineMode = ["-" , "--","-*"];
%     if Fibra.ASEFlag == 0
%         for s = 1:length(Signal.modos) % Grafico
%             graf.signal.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
%             subplot 221
%             for f = 1:6:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
%                 plot(graf.z,graf.signal.(strcat("LP_",Signal.modos(s)))(f,:) , LineMode(s) ,'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
%             end
%             legend('Location','southoutside','Box','off','NumColumns',4)
%         end; clear s f ;
%     else % No calcula espectro ASE
%         for s = 1:length(Signal.modos) % Grafico
%             graf.signal.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
%             subplot 121
%             for f = 1:6:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
%                 plot(graf.z,graf.signal.(strcat("LP_",Signal.modos(s)))(f,:) , LineMode(s) , 'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
%             end
%             legend('Location','southoutside','Box','off','NumColumns',4) ; 
%         end; clear s f ;
%     end
%     
% %     % Bombeo
% %     for s = 1:length(Pump.modos) % Grafico
% %         graf.pump.(strcat("LP_",Pump.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).pump.Potencia_dBm(:,:,s);
% %         figure(n)
% %         subplot 222
% %         for f = 1:length(Pump.lambda.LP_01)
% %             plot(graf.z,graf.pump.(strcat("LP_",Pump.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Pump.modos(s))," @"),strcat(int2str(Pump.lambda.(strcat("LP_",Pump.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Pump}') ; legend(); grid on
% %         end
% %     end
% legend()
%     % Ganancias
%     if Fibra.ASEFlag == 0
%         for s = 1:length(Signal.modos)
%             graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
%             leyenda = strcat((strcat("LP",Signal.modos(s))));
%             ax = subplot (2,2,2);
%             
%             ejex = Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9;
%             plot(ejex,graf.ganancias.(strcat("LP_",Signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%             %ax.ColorOrderIndex = s;
%             xlabel('λ [nm]') ; ylabel(ylab)
% 
%             legend('Location','southoutside','Box','off','NumColumns',5)
%         end ; clear s ejex leyenda;
% 
%     else % No calcula espectro ASE
%         for s = 1:length(Signal.modos)
%             graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
%             leyenda = strcat((strcat("LP",Signal.modos(s))));
%             ax = subplot (1,2,2);
%             
%             ejex = Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9;
%             plot(ejex,graf.ganancias.(strcat("LP_",Signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%             %ax.ColorOrderIndex = s;
%             xlabel('λ [nm]') ; ylabel(ylab)
%             ylim([10 25])
%             
%             legend('Location','southoutside','Box','off','NumColumns',5)
%         end ; clear s ejex leyenda;
%     end
%     
%     % Espectro de ganancia
%     if Fibra.ASEFlag == 0
%         for ase = 1:length(Signal.modos)
%             graf.freq_ase = EDFA.(strcat('Nucleo',int2str(n))).ASE_Spectrum.lambdas*1e9;
%             graf.ASE_Spectrum.(strcat("LP_",Signal.modos(ase))) = EDFA.(strcat("Nucleo",int2str(n))).ASE_Spectrum.mag.(strcat("LP_",Signal.modos(ase)));
%             subplot(2,2,[3,4])
%             plot( graf.freq_ase, graf.ASE_Spectrum.(strcat("LP_",Signal.modos(ase))) , 'DisplayName', strcat( "LP",Signal.modos(ase) ) ) ; hold on;
%             xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Potencia a la salida del EDFA') ; legend()
%         end ; clear ase xlab ylab t_end ;
%     end
%     
% end


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
    %ax.ColorOrderIndex = s;
    xlabel('Longitud de onda [nm]'); ylabel(ylab) %xlabel('λ [nm]') 
    
    legend('Location','southoutside','Box','off','NumColumns',2)
    %ylim([15 20])
    hold on
    set(gca,"ColorOrderIndex",s)
    plot(ejex,EDFA.Nucleo1.salida.ganancias_sinGEF.(strcat("LP_",Signal.modos(s))) , '--*' , 'DisplayName',strcat(leyenda,' sin GEF' )   )
end ; clear s ejex leyenda;
dim = [.132 .63 .4 .3]; str = strcat("Ripple: ",num2str(round(EDFA.Nucleo1.GEF.Ripple,2)) , ' dB en Modo LP01');
annotation('textbox',dim,'String',str,'FitBoxToText','on');
 
figure(2)
plot(Signal.lambda.(strcat("LP_",Signal.modos(1))).*1e9 , EDFA.Nucleo1.GEF.best_weight_Function ) 
title('Filtro utilizado') ; xlabel('Longitud de onda [nm]') ; ylabel('Magnitud')

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

