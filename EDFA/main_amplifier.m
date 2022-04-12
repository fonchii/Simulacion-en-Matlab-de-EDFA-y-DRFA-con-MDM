% Simulador EDFA MM
close all; clear all; clc

% Parámetros de entrada

    % Modos y canales de señal y bombeo

%signal.modos = ["01" "11_a" "11_b" "02" "21_a" "21_b" "12_a" "12_b" "31_a" "31_b"] ;
signal.modos = ["01" "11_a" "11_b" "02" "21_a" "21_b"] ;

%signal.lambda.LP_01     = [1555e-9 1560e-9 1570e-9];                P0_signal.LP_01     = [-15 -15 -15];
%signal.lambda.LP_11_a   = [1550e-9 1560e-9 1570e-9];                P0_signal.LP_11_a   = [-15 -15 -15];
%signal.lambda.LP_01     = [1520e-9:2.5e-9:1570e-9];                 P0_signal.LP_01     = -15*ones(1,length(signal.lambda.LP_01));
signal.lambda.LP_01     = [1520e-9:1e-9:1570e-9];                                  P0_signal.LP_01     = -15*ones(1,length(signal.lambda.LP_01));
signal.lambda.LP_11_a   = 1550e-9;                                  P0_signal.LP_11_a   = -15*ones(1,length(signal.lambda.LP_11_a));
signal.lambda.LP_11_b   = [1550e-9];                                P0_signal.LP_11_b   = -15*ones(1,length(signal.lambda.LP_11_b));
signal.lambda.LP_11     = [1550e-9];                                P0_signal.LP_11     = -15*ones(1,length(signal.lambda.LP_11));
signal.lambda.LP_02     = [1550e-9];                                P0_signal.LP_02     = -15*ones(1,length(signal.lambda.LP_02));
signal.lambda.LP_21_a   = [1550e-9];                                P0_signal.LP_21_a   = -15*ones(1,length(signal.lambda.LP_21_a));
signal.lambda.LP_21_b   = [1550e-9];                                P0_signal.LP_21_b   = -15*ones(1,length(signal.lambda.LP_21_b));
signal.lambda.LP_23     = [1520e-9:2.5e-9:1570e-9];                 P0_signal.LP_23     = -15*ones(1,length(signal.lambda.LP_23));
signal.lambda.LP_12_a   = [1550e-9];                                P0_signal.LP_12_a   = -15*ones(1,length(signal.lambda.LP_12_a));
signal.lambda.LP_12_b   = [1550e-9];                                P0_signal.LP_12_b   = -15*ones(1,length(signal.lambda.LP_12_b));
signal.lambda.LP_31_a   = [1550e-9];                                P0_signal.LP_31_a   = -15*ones(1,length(signal.lambda.LP_31_a));
signal.lambda.LP_31_b   = [1550e-9];                                P0_signal.LP_31_b   = -15*ones(1,length(signal.lambda.LP_31_b));

%signal.lambda.LP_01 = [1530e-9 1535e-9 1540e-9 1545e-9 1550e-9 1555e-9 1560e-9 1565e-9 1570e-9]; P0_signal.LP_01 = [-10 -10 -10 -10 -10 -10 -10 -10 -10];
%signal.lambda.LP_11_a = [1530e-9 1535e-9 1540e-9 1545e-9 1550e-9 1555e-9 1560e-9 1565e-9 1570e-9]; P0_signal.LP_11_a = [-10 -10 -10 -10 -10 -10 -10 -10 -10];

%pump.modos = ["11_a" "11_b" "02" "21_a" "21_b" "12_a" "12_b" "31_a" "31_b"]   ;
pump.modos = ["01" "14"]   ;

pump.lambda.LP_01   = [980e-9];                         P0_pump.LP_01   = [150e-3]  ;  
pump.lambda.LP_11_a = [980e-9];                         P0_pump.LP_11_a = [150e-3] ;  
pump.lambda.LP_11_b = [980e-9];                         P0_pump.LP_11_b = [150e-3]  ;
pump.lambda.LP_11   = [980e-9];                         P0_pump.LP_11   = [150e-3]  ;
pump.lambda.LP_14   = [980e-9];                         P0_pump.LP_14   = [150e-3]  ;
pump.lambda.LP_02   = [980e-9];                         P0_pump.LP_02   = [40e-3]  ;
pump.lambda.LP_21   = [980e-9];                         P0_pump.LP_21   = [1000e-3]  ;
pump.lambda.LP_21_a = [980e-9];                         P0_pump.LP_21_a = [110e-3]  ;
pump.lambda.LP_21_b = [980e-9];                         P0_pump.LP_21_b = [115e-3]  ;
pump.lambda.LP_12_a = [980e-9];                         P0_pump.LP_12_a = [590e-3]  ;
pump.lambda.LP_12_b = [980e-9];                         P0_pump.LP_12_b = [590e-3]  ;
pump.lambda.LP_31_a = [980e-9];                         P0_pump.LP_31_a = [180e-3]  ;
pump.lambda.LP_31_b = [980e-9];                         P0_pump.LP_31_b = [190e-3]  ;
pump.lambda.LP_12_a = [980e-9];                         P0_pump.LP_12_a = [250e-3 250e-3]  ;
pump.lambda.LP_42_a = [980e-9];                         P0_pump.LP_42_a = [250e-3 250e-3]  ;

    % POTENCIAS
for i=1:length(signal.modos)        % Potencia de señal a mW
    for j=1:length( P0_signal.(strcat("LP_",signal.modos(i))) )
        P0_signal.( strcat("LP_",signal.modos(i) ) )(j) = 1e-3*10^(P0_signal.( strcat("LP_",signal.modos(i) ) )(j)/10);
    end
end ;clear i j;
signal.P0 = P0_signal; 
pump.P0 = P0_pump;
ASE = -50;                                                  %dBm

%h=6.626*10^(-34);
%ASE= 10*log10(P.Np*h*P.Fc*P.OpticalBW*1e3); % [dBm] % thesis 2016-Phd Thesis Spatially Integrated Erbium-Doped 

    % Datos de la fibra
fibra.nucleos = 1;                                           % Numero de nucleos
fibra.largo = 3     ; fibra.radio = 15e-6   ; fibra.N = 7e24; % fibra.N = 3e24; 
%fibra.n1 = 1.47    ; fibra.n2 = 1.44       ; % PAPER
%fibra.AN = 0.2     ;  % Se prioriza el uso de AN en caso de entregar n1,n2,AN
%fibra.n1 = 1.44959 ; fibra.n2 = 1.43502    ; %VPI
% PAPER
fibra.AN = 0.185;   %fibra.AN = 0.2 ; 
fibra.n1 = 1.45 ;   fibra.n2 = sqrt((fibra.n1^2-fibra.AN^2));
%fibra.dvk= P.OpticalBW; % diferencia : max_lambda - min_lambda 

fibra.ASEFlag = 1; % EVITA CALCULO DE ESPECTRO ASE


%%
tic;
EDFA = EDFA_MM(fibra,signal,pump,ASE); 
t_end = toc; fprintf('Tiempo de cómputo: %.2f segundos\n', t_end);


    %% Graficos

close all
graf.Nc = length(fieldnames(EDFA)); 
graf.z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
for n = 1:graf.Nc
    figure(n)
    
    % Señal
    for s = 1:length(signal.modos) % Grafico
        graf.signal.(strcat("LP_",signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",signal.modos(s)));
        subplot 221
        for f = 1:length(signal.lambda.(strcat("LP_",signal.modos(s))))
            plot(graf.z,graf.signal.(strcat("LP_",signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",signal.modos(s))," @"),strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
        end
    end; clear s f ;
    
%     % Bombeo
%     for s = 1:length(pump.modos) % Grafico
%         graf.pump.(strcat("LP_",pump.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).pump.Potencia_dBm(:,:,s);
%         figure(n)
%         subplot 222
%         for f = 1:length(pump.lambda.LP_01)
%             plot(z,graf.pump.(strcat("LP_",pump.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",pump.modos(s))," @"),strcat(int2str(pump.lambda.(strcat("LP_",pump.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Pump}') ; legend(); grid on
%         end
%     end

    % Ganancias
    for s = 1:length(signal.modos)
        graf.ganancias.(strcat("LP_",signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",signal.modos(s)));
        leyenda = strcat((strcat("LP",signal.modos(s))));
        ax = subplot (2,2,2);
        
        ejex = signal.lambda.(strcat("LP_",signal.modos(s))).*1e9;
        plot(ejex,graf.ganancias.(strcat("LP_",signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
        %ax.ColorOrderIndex = s;
        xlabel('λ [nm]') ; ylabel(ylab)
        
    end ; clear s ejex leyenda;
    
    % Espectro de ganancia
    for ase = 1:length(signal.modos)
        graf.freq_ase = EDFA.(strcat('Nucleo',int2str(n))).ASE_Spectrum.lambdas*1e9;
        graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase))) = EDFA.(strcat("Nucleo",int2str(n))).ASE_Spectrum.mag.(strcat("LP_",signal.modos(ase)));
        subplot(2,2,[3,4])
        plot( graf.freq_ase, graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase)))(:,2) , 'DisplayName', strcat( "LP",signal.modos(ase) ) ) ; hold on;
        xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Potencia a la salida del EDFA') ; legend()
    end ; clear ase xlab ylab t_end ;
    
end

%% Graficos Modos

% Señal
n=n+1 ; figure(n)
graficar_modos(fibra,signal.modos,signal.lambda.LP_01)
sgtitle('Modos de Señal') 

% Bombeo
n=n+1 ; figure(n)
graficar_modos(fibra,pump.modos,pump.lambda.LP_01)
sgtitle('Modos de Bombeo') ; %clear n ;

%% Graficar potencia vs frecuencia
n=n+1; figure(n)
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
graf.Nc = length(fieldnames(EDFA)); 
for n = 1:graf.Nc
    for s = 1:length(signal.modos)
        graf.ganancias.(strcat("LP_",signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",signal.modos(s)));
        leyenda = strcat("LP",signal.modos(s));
        ejex = fliplr((3*10^8./(signal.lambda.(strcat("LP_",signal.modos(s))).*1e9)));
        plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
        %ax.ColorOrderIndex = s;
        xlabel('Frecuencia [THz]') ; ylabel(ylab)

    end ; clear s ejex leyenda;
end


%% Graficar potencia vs wavelenght
%n=n+1; figure(n)
%xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
%graf.Nc = length(fieldnames(EDFA)); 
%for n = 1:graf.Nc
%    for s = 1:length(signal.modos)
%        graf.ganancias.(strcat("LP_",signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",signal.modos(s)));
%        leyenda = strcat((strcat("LP",signal.modos(s))));
%        ejex = fliplr((3*10^8./(signal.lambda.(strcat("LP_",signal.modos(s))).*1e9)));
%        plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%        %ax.ColorOrderIndex = s;
%        xlabel('Frecuencia [THz]') ; ylabel(ylab)
%        set(gca,'xdir','reverse')
%    end ; clear s ejex leyenda;
%end

