% Simulador EDFA MM
close all; 
clear all; clc

%% Parámetros de entrada

    % Señal : Modos y Canales
NCh = 41;%36;
Signal.modos = ["01"] ;

%Signal.lambda.LP_01     = linspace(1540e-9,1575e-9,NCh);              P0_signal.LP_01     = -15*ones(1,length(Signal.lambda.LP_01));
Signal.lambda.LP_01     = linspace(1530e-9,1570e-9,NCh);              P0_signal.LP_01     = -15*ones(1,length(Signal.lambda.LP_01));

    % Bombeo : Modos y Canales
Pump.modos = ["01" ]   ;

Pump.lambda.LP_01   = 980e-9;                                 P0_pump.LP_01   = [250e-3]  ;  

    % POTENCIAS

for i=1:length(Signal.modos)        % Potencia: dBm -> mW
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
Fibra.largo = 1     ; Fibra.radio = 5e-6   ; Fibra.N = 7e24; % fibra.N = 3e24; 

Fibra.dvk=300e9;
Signal.NumberOfChannels=30;


Fibra.n1 = 1.45 ;   
Fibra.n2 = 1.4354 ;
%Fibra.dvk= P.OpticalBW; % diferencia : max_lambda - min_lambda 

% Despliegue de Info
Fibra.WaitBar = 1; Fibra.Avance = 1;    

% Calculo ASE?
Fibra.ASEFlag = 1;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)


%% 


tic;
% --------- Iterar por Largos --------- %
% for largos=[3 5 7 10 15 20]%largos=1:2:20
%         for pump = [150,200,250,300,400,500,700,1000,1500]
%             %P0_pump.LP_01   = [150e-3]  ; Pump.P0 = P0_pump;
%             P0_pump.LP_01   = pump*1e-3  ; Pump.P0 = P0_pump;
%             Fibra.largo = largos;
%             fprintf('Iniciando Fibra %.0f de %.0f\n', largos,30);
%         
%             %Largos.(strcat("EDFA_",num2str(largos),'m')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);       % Con efecto acomplamiento de Potencia intermodal
%             %Largos.(strcat("EDFA_",num2str(largos),'m')) = EDFA_MMvpi2(Fibra,Signal,Pump,ASE);       
%             Largos_v2.(strcat("EDFA_",num2str(largos),'m')).(strcat('Pump',num2str(pump),'mw')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);
%         end
%     
% end
% largos_time=toc; fprintf('Tiempo Total de cómputo: %.2f segundos\n', largos_time )



% --------- Iterar por Potencias --------- %

%for largos=3:2:20
% cont = 1;
% for largos=[3 5 7 10 15 20]
%     for potencias=150:50:1000
%         P0_pump.LP_01 = potencias.*1e-3; Pump.P0 = P0_pump;
%         fprintf('Iniciando Fibra %.0f de %.0f\n', cont,(length(150:50:1000)*length([3 5 7 10 15 20]))); cont = cont+1;
%         Fibra.largo = largos; % 3,5,7,10
%         Potencias.(strcat('L',num2str(largos),'m')).(strcat("EDFA_",num2str(potencias),'mw')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);       % Con efecto acomplamiento de Potencia intermodal
%     
%     end
% end
% largos_time=toc; fprintf('Tiempo Total de cómputo: %.2f segundos\n', largos_time )

% --------- Ganancias Espectral --------- %
% Cambiar arriba a estos datos:
%Signal.lambda.LP_01     = linspace(1530e-9,1570e-9,41);              P0_signal.LP_01     = -15*ones(1,length(Signal.lambda.LP_01)); Signal.NumberOfChannels=41;
for largos=3:2:20
    cont = 1;
    for largos=[3 5 7 10 15 20]
        for potencias=150:50:1000
            P0_pump.LP_01 = potencias.*1e-3; Pump.P0 = P0_pump;
            fprintf('Iniciando Fibra %.0f de %.0f\n', cont,(length(150:50:1000)*length([3 5 7 10 15 20]))); cont = cont+1;
            Fibra.largo = largos; % 3,5,7,10
            GainSpectrum.(strcat('L',num2str(largos),'m')).(strcat("EDFA_",num2str(potencias),'mw')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);       % Con efecto acomplamiento de Potencia intermodal
        end
    end
end
spectrumG_time=toc; fprintf('Tiempo Total de cómputo: %.2f segundos\n', spectrumG_time )




% -------- OLD ---------
% % for potencias=150:50:1000       
% %     P0_pump.LP_01 = potencias.*1e-3; Pump.P0 = P0_pump;
% %     fprintf('Iniciando Fibra %.0f de %.0f\n', potencias,1000);
% %     Fibra.largo = 20; % 3,5,7,10
% %     Potencias.(strcat("EDFA_",num2str(potencias),'mw')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);       % Con efecto acomplamiento de Potencia intermodal
% % 
% % end
% % largos_time=toc; fprintf('Tiempo Total de cómputo: %.2f segundos\n', largos_time )



    %% Graficos

% close all
% graf.Nc = length(fieldnames(EDFA)); 
% graf.z = EDFA.(strcat('Nucleo',int2str(1))).z;
% xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
% for n = 1:graf.Nc
%     figure(n)
%     
%     % Señal
%     if Fibra.ASEFlag == 0
%         for s = 1:length(Signal.modos) % Grafico
%             graf.signal.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
%             subplot 221
%             for f = 1:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
%                 plot(graf.z,graf.signal.(strcat("LP_",Signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
%             end
%         end; clear s f ;
%     else % No calcula espectro ASE
%         for s = 1:length(Signal.modos) % Grafico
%             graf.signal.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
%             subplot 121
%             for f = 1:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
%                 plot(graf.z,graf.signal.(strcat("LP_",Signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
%             end
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
% 
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
%         end ; clear s ejex leyenda;
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
%             
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
% 
% %% Graficos 3D de Modos 
% 
% % Señal
% figure()
% graficar_modos(Fibra,Signal.modos,Signal.lambda.LP_01)
% sgtitle('Modos de Señal') 
% 
% % Bombeo
% figure()
% graficar_modos(Fibra,Pump.modos,Pump.lambda.LP_01)
% sgtitle('Modos de Bombeo') ; 
% 
% %% Graficar potencia vs frecuencia
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

