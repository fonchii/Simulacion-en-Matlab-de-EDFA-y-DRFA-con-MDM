%% Indicaciones:
% Las cantidad de canales por modo deben ser iguales

close all; clear all; clc

    % Modos
signal.modos = ["01" "11_a" "11_b"] ;
signal.lambda.LP_01 = [1540e-9 1560e-9 1570e-9];
signal.lambda.LP_11_a = [1540e-9 1560e-9 1570e-9];
signal.lambda.LP_11_b = [1540e-9 1560e-9 1570e-9];
signal.lambda.LP_11 = [1540e-9 1560e-9 1570e-9];
signal.lambda.LP_02 = [1540e-9 1560e-9 1570e-9];


%pump.modos = ["01","11_a","11_b"] ;
pump.modos = ["11" "21"] ;
pump.lambda.LP_01 = [980e-9 990e-9];
pump.lambda.LP_11_a = [980e-9 990e-9];
pump.lambda.LP_11_b = [982e-9 960e-9];
pump.lambda.LP_11 = [980e-9 990e-9];
pump.lambda.LP_02 = [980e-9 990e-9];
pump.lambda.LP_21 = [980e-9 990e-9];

    %POTENCIAS
P0_signal = - 20 ;                                           %dBm
P0_signal = 1e-3*10.^(P0_signal/10);                         %mW
%P0_pump = 100e-3  ;                                         %mW
P0_pump = 250e-3  ;                                          %mW
signal.P0 = P0_signal; 
pump.P0 = P0_pump;
Pase = -50;                                                  %dBm

    % Datos de la fibra
fibra.nucleos = 1;                                           % Numero de nucleos
fibra.largo = 10;  fibra.AN = 0.2 ; fibra.radio = 10e-6 ; fibra.N = 3e24; 
%fibra.n1 = 1.46; fibra.n2 = 1.4462;
fibra.n1 = 1.47; fibra.n2 = 1.44;

tic;
EDFA = EDFA_MMv2(fibra,signal,pump,Pase);
t_end = toc; fprintf('Tiempo de iteración: %.2f segundos\n', t_end);


%% Graficos

close all
Nch = length(fieldnames(EDFA)); 
z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
for n = 1:Nch
    figure(n)
    
    % Señal
    for s = 1:length(signal.modos) % Grafico
        graf.signal.(strcat("LP_",signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm(:,:,s);
        subplot 221
        for f = 1:length(signal.lambda.LP_01)
            plot(z,graf.signal.(strcat("LP_",signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",signal.modos(s))," @"),strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
        end
    end
    
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
        graf.ganancias.(strcat("LP_",signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias(:,:,s);
        ejex = signal.lambda.LP_01.*1e9;
        ax = subplot (2,2,2);
        plot(ejex,graf.ganancias.(strcat("LP_",signal.modos(s))) , 'DisplayName',strcat((strcat("LP",signal.modos(s))))) ; hold on ; title('Ganancias') ; legend(); grid on ; 
        ax.ColorOrderIndex = s; scatter(ejex,graf.ganancias.(strcat("LP_",signal.modos(s))),'filled', 'DisplayName','') ; xlabel('λ [nm]') ; ylabel(ylab)
    end
    
    % Espectro de ganancia
    for ase = 1:length(signal.modos)
        freq_ase = EDFA.(strcat('Nucleo',int2str(n))).ASE_Spectrum.lambdas*1e9;
        graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase))) = EDFA.(strcat("Nucleo",int2str(n))).ASE_Spectrum.mag(:,:,ase);
        subplot(2,2,[3,4])
        plot( freq_ase, graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase)))(:,2) , 'DisplayName', strcat( "LP",signal.modos(ase) ) ) ; hold on;
        xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Espectro de ganancia') ; legend()
    end
    
end

%% Graficos Modos
n=1;
n=n+1 ; figure(n)
% Señal
graficar_modosv2(fibra,signal.modos,signal.lambda.LP_01)
sgtitle('Modos de Señal') 

n=n+1 ; figure(n)
% Bombeo
graficar_modosv2(fibra,pump.modos,pump.lambda.LP_01)
sgtitle('Modos de Bombeo') 
    

