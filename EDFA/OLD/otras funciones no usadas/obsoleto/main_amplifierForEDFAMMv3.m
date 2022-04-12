%% Indicaciones:
% Las cantidad de canales por modo deben ser iguales
% Adaptación con coeficiente de acoplamiento del VPI

close all; clear all; clc

    % Modos y canales de señal y bombeo
signal.modos = ["01" "11_a"] ;
%signal.lambda.LP_01 = [1550e-9 1560e-9 1570e-9];
%signal.lambda.LP_11_a = [1550e-9 1560e-9 1570e-9];
signal.lambda.LP_01 = [1520e-9 1525e-9 1530e-9 1535e-9 1540e-9 1545e-9 1550e-9 1555e-9 1560e-9 1565e-9 1570e-9 1575e-9 1580e-9];
signal.lambda.LP_11_a = [1520e-9 1525e-9 1530e-9 1535e-9 1540e-9 1545e-9 1550e-9 1555e-9 1560e-9 1565e-9 1570e-9 1575e-9 1580e-9];

signal.Gamma.LP_01 = 0.9565;
signal.Gamma.LP_11_a = 0.871;


pump.modos = ["01"] ;
pump.lambda.LP_01 = [980e-9];
pump.lambda.LP_11_a = [980e-9 990e-9];
pump.lambda.LP_11_b = [982e-9 960e-9];
pump.lambda.LP_11 = [980e-9 990e-9];
pump.lambda.LP_02 = [980e-9 990e-9];
pump.lambda.LP_21 = [980e-9 990e-9];

pump.Gamma.LP_01 = 0.987209;

    %POTENCIAS
P0_signal = - 15;                                           %dBm
P0_signal = 1e-3*10.^(P0_signal/10);                         %mW
P0_pump = 250e-3  ;                                          %mW
signal.P0 = P0_signal; 
pump.P0 = P0_pump;
Pase = -50;                                                  %dBm

    % Datos de la fibra
fibra.nucleos = 1;                                           % Numero de nucleos
fibra.largo = 2;  fibra.AN = 0.2 ; fibra.radio = 5e-6 ; fibra.N = 7e24; % fibra.N = 3e24; 
%fibra.n1 = 1.46; fibra.n2 = 1.4462;
fibra.n1 = 1.47 ; fibra.n2 = 1.44; % PAPER
fibra.acoplamiento = 1e-4 * fibra.largo;
%fibra.n1 = 1.44959; fibra.n2 = 1.43502; %VPI

tic;
EDFA = EDFA_MMv2(fibra,signal,pump,Pase);
t_end = toc; fprintf('Tiempo de iteración: %.2f segundos\n', t_end);


    %% Graficos

close all
Nc = length(fieldnames(EDFA)); 
z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
% for n = 1:Nc
%     figure(n)
%     
%     % Señal
%     for s = 1:length(signal.modos) % Grafico
%         graf.signal.(strcat("LP_",signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm(:,:,s);
%         subplot 221
%         for f = 1:length(signal.lambda.LP_01)
%             plot(z,graf.signal.(strcat("LP_",signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",signal.modos(s))," @"),strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
%         end
%     end
%     
% %     % Bombeo
% %     for s = 1:length(pump.modos) % Grafico
% %         graf.pump.(strcat("LP_",pump.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).pump.Potencia_dBm(:,:,s);
% %         figure(n)
% %         subplot 222
% %         for f = 1:length(pump.lambda.LP_01)
% %             plot(z,graf.pump.(strcat("LP_",pump.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",pump.modos(s))," @"),strcat(int2str(pump.lambda.(strcat("LP_",pump.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Pump}') ; legend(); grid on
% %         end
% %     end
% 
%     % Ganancias
%     for s = 1:length(signal.modos)
%         graf.ganancias.(strcat("LP_",signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias(:,:,s);
%         ejex = signal.lambda.LP_01.*1e9;
%         ax = subplot (2,2,2);
%         plot(ejex,graf.ganancias.(strcat("LP_",signal.modos(s))) , 'DisplayName',strcat((strcat("LP",signal.modos(s))))) ; hold on ; title('Ganancias') ; legend(); grid on ; 
%         ax.ColorOrderIndex = s; scatter(ejex,graf.ganancias.(strcat("LP_",signal.modos(s))),'filled', 'DisplayName','') ; xlabel('λ [nm]') ; ylabel(ylab)
%     end
%     
%     % Espectro de ganancia
%     for ase = 1:length(signal.modos)
%         freq_ase = EDFA.(strcat('Nucleo',int2str(n))).ASE_Spectrum.lambdas*1e9;
%         graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase))) = EDFA.(strcat("Nucleo",int2str(n))).ASE_Spectrum.mag(:,:,ase);
%         subplot(2,2,[3,4])
%         plot( freq_ase, graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase)))(:,2) , 'DisplayName', strcat( "LP",signal.modos(ase) ) ) ; hold on;
%         xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Potencia a la salida del EDFA') ; legend()
%     end
%     
% end

%% Graficos Modos

% % Señal
% n=n+1 ; figure(n)
% graficar_modosv2(fibra,signal.modos,signal.lambda.LP_01)
% sgtitle('Modos de Señal') 
% 
% % Bombeo
% n=n+1 ; figure(n)
% graficar_modosv2(fibra,pump.modos,pump.lambda.LP_01)
% sgtitle('Modos de Bombeo') 


%%
for s = 1:length(signal.modos)
    graf.ganancias.(strcat("LP_",signal.modos(s))) = EDFA1.(strcat("Nucleo",int2str(1))).salida.ganancias(:,:,s);
    ejex = signal.lambda.LP_01.*1e9;
    ax = subplot (1,1,1);
    leyenda = strcat( strcat( strcat("LP", signal.modos(s)) ), " Modelo CoefAcop" );
    plot(ejex , graf.ganancias.(strcat("LP_",signal.modos(s))) ,'-o', 'DisplayName' , leyenda) ; hold on ; title('Ganancias') ; legend(); grid on ;
    %ax.ColorOrderIndex = s;  
    xlabel('λ [nm]') ; ylabel(ylab)
end

for s = 1:length(signal.modos)
    graf.ganancias.(strcat("LP_",signal.modos(s))) = EDFA2.(strcat("Nucleo",int2str(1))).salida.ganancias(:,:,s);
    ejex = signal.lambda.LP_01.*1e9;
    ax = subplot (1,1,1);
    leyenda = strcat((strcat("LP",signal.modos(s)))) ;
    plot(ejex , graf.ganancias.(strcat("LP_",signal.modos(s))) ,'-o', 'DisplayName' , leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
    %ax.ColorOrderIndex = s; 
    xlabel('λ [nm]') ; ylabel(ylab)
end
