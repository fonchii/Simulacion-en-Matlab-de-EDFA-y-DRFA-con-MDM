% Simulador EDFA MM
close all; 
clear all; clc

%% Parámetros de entrada

    % Señal : Modos y Canales

%Signal.modos = ["01" "11_a" "11_b" "02" "21_a" "21_b" "32"] ;
Signal.modos = ["01" "11_a" "11_b"] ;

Signal.lambda.LP_01     = linspace(1524.6e-9,1570.1e-9,30);         P0_signal.LP_01     = -15*ones(1,length(Signal.lambda.LP_01));
Signal.lambda.LP_11_a   = linspace(1524.6e-9,1570.1e-9,30);         P0_signal.LP_11_a   = -15*ones(1,length(Signal.lambda.LP_11_a));
Signal.lambda.LP_11_b   = linspace(1524.6e-9,1570.1e-9,30);         P0_signal.LP_11_b   = -15*ones(1,length(Signal.lambda.LP_11_b));
Signal.lambda.LP_11     = [1550e-9];                                P0_signal.LP_11     = -15*ones(1,length(Signal.lambda.LP_11));
Signal.lambda.LP_02     = [1550e-9];                                P0_signal.LP_02     = -15*ones(1,length(Signal.lambda.LP_02));
Signal.lambda.LP_21_a   = [1550e-9];                                P0_signal.LP_21_a   = -15*ones(1,length(Signal.lambda.LP_21_a));
Signal.lambda.LP_21_b   = [1550e-9];                                P0_signal.LP_21_b   = -15*ones(1,length(Signal.lambda.LP_21_b));
Signal.lambda.LP_23     = linspace(1524.6e-9,1570.1e-9,20);         P0_signal.LP_23     = -15*ones(1,length(Signal.lambda.LP_23));
Signal.lambda.LP_12_a   = [1550e-9];                                P0_signal.LP_12_a   = -15*ones(1,length(Signal.lambda.LP_12_a));
Signal.lambda.LP_12_b   = [1550e-9];                                P0_signal.LP_12_b   = -15*ones(1,length(Signal.lambda.LP_12_b));
Signal.lambda.LP_31_a   = [1550e-9];                                P0_signal.LP_31_a   = -15*ones(1,length(Signal.lambda.LP_31_a));
Signal.lambda.LP_31_b   = [1550e-9];                                P0_signal.LP_31_b   = -15*ones(1,length(Signal.lambda.LP_31_b));
Signal.lambda.LP_32   = [1550e-9];                                P0_signal.LP_32   = -15*ones(1,length(Signal.lambda.LP_31_b));

    % Bombeo : Modos y Canales
Pump.modos = ["01" ]   ;

Pump.lambda.LP_01   = [980e-9];                                 P0_pump.LP_01   = [250e-3]  ;  
Pump.lambda.LP_11_a = [980e-9];                                 P0_pump.LP_11_a = [150e-3] ;  
Pump.lambda.LP_11_b = [980e-9];                                 P0_pump.LP_11_b = [150e-3]  ;
Pump.lambda.LP_11   = [980e-9];                                 P0_pump.LP_11   = [150e-3]  ;
Pump.lambda.LP_14   = [980e-9];                                 P0_pump.LP_14   = [150e-3]  ;
Pump.lambda.LP_02   = [980e-9];                                 P0_pump.LP_02   = [40e-3]  ;
Pump.lambda.LP_21   = [980e-9];                                 P0_pump.LP_21   = [1000e-3]  ;
Pump.lambda.LP_21_a = [980e-9];                                 P0_pump.LP_21_a = [110e-3]  ;
Pump.lambda.LP_21_b = [980e-9];                                 P0_pump.LP_21_b = [115e-3]  ;
Pump.lambda.LP_12_a = [980e-9];                                 P0_pump.LP_12_a = [590e-3]  ;
Pump.lambda.LP_12_b = [980e-9];                                 P0_pump.LP_12_b = [590e-3]  ;
Pump.lambda.LP_31_a = [980e-9];                                 P0_pump.LP_31_a = [180e-3]  ;
Pump.lambda.LP_31_b = [980e-9];                                 P0_pump.LP_31_b = [190e-3]  ;
Pump.lambda.LP_12_a = [980e-9];                                 P0_pump.LP_12_a = [250e-3 250e-3]  ;
Pump.lambda.LP_22_a = [980e-9];                                 P0_pump.LP_22_a = [110e-3]  ;
Pump.lambda.LP_22_b = [980e-9];                                 P0_pump.LP_22_b = [115e-3]  ;
Pump.lambda.LP_42_a = [980e-9];                                 P0_pump.LP_42_a = [250e-3 250e-3]  ;

    % POTENCIAS

for i=1:length(Signal.modos)        % Potencia: W -> mW
    for j=1:length( P0_signal.(strcat("LP_",Signal.modos(i))) )
        P0_signal.( strcat("LP_",Signal.modos(i) ) )(j) = 1e-3*10^(P0_signal.( strcat("LP_",Signal.modos(i) ) )(j)/10);
    end
end ;clear i j;
Signal.P0 = P0_signal; 
Pump.P0 = P0_pump;
ASE = -41;                                                  %dBm  -50

%h=6.626*10^(-34);
%ASE= 10*log10(P.Np*h*P.Fc*P.OpticalBW*1e3); % [dBm] % thesis 2016-Phd Thesis Spatially Integrated Erbium-Doped 

    % Datos de la fibra

Fibra.nucleos = 1;                                           % Numero de nucleos
Fibra.largo = 3     ; Fibra.radio = 5e-6   ; Fibra.N = 7e24; % fibra.N = 3e24; 
%Fibra.n1 = 1.47    ; Fibra.n2 = 1.44       ; % PAPER
%Fibra.AN = 0.2     ;  % Se prioriza el uso de AN en caso de entregar n1,n2,AN
%Fibra.n1 = 1.44959 ; Fibra.n2 = 1.43502    ; %VPI
Fibra.dvk=300e9;
Signal.NumberOfChannels=30;

%Fibra.AN = 0.185;   % = 0.2 ; 
Fibra.n1 = 1.45 ;   %Fibra.n2 = sqrt((Fibra.n1^2-Fibra.AN^2));
Fibra.n2 = 1.4354 ;
%Fibra.dvk= P.OpticalBW; % diferencia : max_lambda - min_lambda 

Fibra.WaitBar = 1; Fibra.Avance = 1;    % Despliegue de info
Fibra.ASEFlag = 0;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)


%%
tic;
%EDFA = EDFA_MM(fibra,signal,pump,ASE);         % Sin efecto acomplamiento de Potencia intermodal
EDFA = EDFA_MMvPCCv2(Fibra,Signal,Pump,ASE);      % Con efecto acomplamiento de Potencia intermodal
t_end = toc; fprintf('Tiempo de cómputo: %.2f segundos\n', t_end);


    %% Graficos

close all
graf.Nc = length(fieldnames(EDFA)); 
graf.z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
for n = 1:graf.Nc
    figure(n)
    
    % Señal
    if Fibra.ASEFlag == 0
        for s = 1:length(Signal.modos) % Grafico
            graf.signal.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
            subplot 221
            for f = 1:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
                plot(graf.z,graf.signal.(strcat("LP_",Signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
            end
        end; clear s f ;
    else % No calcula espectro ASE
        for s = 1:length(Signal.modos) % Grafico
            graf.signal.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
            subplot 121
            for f = 1:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
                plot(graf.z,graf.signal.(strcat("LP_",Signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
            end
        end; clear s f ;
    end
    
%     % Bombeo
%     for s = 1:length(Pump.modos) % Grafico
%         graf.pump.(strcat("LP_",Pump.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).pump.Potencia_dBm(:,:,s);
%         figure(n)
%         subplot 222
%         for f = 1:length(Pump.lambda.LP_01)
%             plot(graf.z,graf.pump.(strcat("LP_",Pump.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Pump.modos(s))," @"),strcat(int2str(Pump.lambda.(strcat("LP_",Pump.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Pump}') ; legend(); grid on
%         end
%     end

    % Ganancias
    if Fibra.ASEFlag == 0
        for s = 1:length(Signal.modos)
            graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
            leyenda = strcat((strcat("LP",Signal.modos(s))));
            ax = subplot (2,2,2);
            
            ejex = Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9;
            plot(ejex,graf.ganancias.(strcat("LP_",Signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
            %ax.ColorOrderIndex = s;
            xlabel('λ [nm]') ; ylabel(ylab)
            
        end ; clear s ejex leyenda;
    else % No calcula espectro ASE
        for s = 1:length(Signal.modos)
            graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
            leyenda = strcat((strcat("LP",Signal.modos(s))));
            ax = subplot (1,2,2);
            
            ejex = Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9;
            plot(ejex,graf.ganancias.(strcat("LP_",Signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
            %ax.ColorOrderIndex = s;
            xlabel('λ [nm]') ; ylabel(ylab)
            
        end ; clear s ejex leyenda;
    end
    
    % Espectro de ganancia
    if Fibra.ASEFlag == 0
        for ase = 1:length(Signal.modos)
            graf.freq_ase = EDFA.(strcat('Nucleo',int2str(n))).ASE_Spectrum.lambdas*1e9;
            graf.ASE_Spectrum.(strcat("LP_",Signal.modos(ase))) = EDFA.(strcat("Nucleo",int2str(n))).ASE_Spectrum.mag.(strcat("LP_",Signal.modos(ase)));
            subplot(2,2,[3,4])
            plot( graf.freq_ase, graf.ASE_Spectrum.(strcat("LP_",Signal.modos(ase))) , 'DisplayName', strcat( "LP",Signal.modos(ase) ) ) ; hold on;
            xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Potencia a la salida del EDFA') ; legend()
        end ; clear ase xlab ylab t_end ;
    end
    
end

%% Graficos 3D de Modos 

% Señal
figure()
graficar_modos(Fibra,Signal.modos,Signal.lambda.LP_01)
sgtitle('Modos de Señal') 

% Bombeo
figure()
graficar_modos(Fibra,Pump.modos,Pump.lambda.LP_01)
sgtitle('Modos de Bombeo') ; 

%% Graficar potencia vs frecuencia
figure()
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
graf.Nc = length(fieldnames(EDFA)); 
for n = 1:graf.Nc
    for s = 1:length(Signal.modos)
        graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
        leyenda = strcat("LP",Signal.modos(s));
        ejex = fliplr((3*10^8./(Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9)));
        plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",Signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
        %ax.ColorOrderIndex = s;
        xlabel('Frecuencia [THz]') ; ylabel(ylab)

    end ; clear s ejex leyenda n;
end


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

