% close all; 
clear all; clc ; close all

% Parámetros de entrada

    % Modos y canales de señal y bombeo
signal.modos = ["11_a" "21_a"] ;

Signal.NumberOfChannels=41;
Wavelength_gridS=linspace(1530,1570,Signal.NumberOfChannels).*1e-9; % Banda C
c=299.792458e6; % [m/s]

Pin=0; %[dBm]

signal.lambda.LP_11_a     = Wavelength_gridS;                                   P0_signal.LP_11_a     = Pin*ones(1,length(signal.lambda.LP_11_a));
signal.lambda.LP_21_a   = Wavelength_gridS;                                     P0_signal.LP_21_a   = Pin*ones(1,length(signal.lambda.LP_21_a));

pump.modos = "01" ;
Wavelength_gridP=980e-9;
Ppump= 300e-3; %[W]

pump.lambda.LP_01   = Wavelength_gridP;                         P0_pump.LP_01   = Ppump  ;  

ModoS=strcat("LP_",signal.modos(:));
ModoP=strcat("LP_",pump.modos(:));

    % POTENCIAS

for i=1:length(signal.modos)        % Potencia de señal a W
    for j=1:length(P0_signal.(ModoS(i)))
        P0_signal.(ModoS(i))(j) = 1e-3*10^(P0_signal.(ModoS(i))(j)/10);
    end
end ;clear i j;
signal.P0 = P0_signal; 
pump.P0 = P0_pump;
h=6.62607015*10^(-34);
P.Np=2; 
P.Fc=c/Wavelength_gridS(ceil(length(Wavelength_gridS)/2)); P.Fb = 50e9; 
ASE= -200;

    % Datos de la fibra
fibra.nucleos = 1;
fibra.largo = 5; fibra.radio = 5.5e-6 ; fibra.N = 7e24; 
fibra.n1 = 1.45 ;   fibra.IndexContrast=0.01;
fibra.AN=fibra.n1*sqrt(2*fibra.IndexContrast);
fibra.n2 =sqrt((fibra.n1^2-fibra.AN^2));
fibra.dvk=P.Fb;
fibra.PumpMode = "reverse";

fibra.WaitBar = 1; fibra.Avance = 1;    % Despliegue de info
fibra.ASEFlag = 1;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)

load("EDFA.mat")
load("EDFA_RP.mat")


    %% Graficos
% SAVE:
set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
%print -dpdf 'NAME'


close all 
Nc = length(fieldnames(EDFA)); 
z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
smodos = ["11a","21a"];

% % % % % % Señal

% for s = 1:length(signal.modos) % Grafico
%     figure(s)
% 
%     graf.signal = EDFA.("Nucleo1").signal.Potencia_dBm;
%     grafReverse.signal = EDFA_RP.("Nucleo1").signal.Potencia_dBm;
%     var = 2;
%     for f = 21:2:31
%         plot(z , graf.signal.(strcat("LP_",signal.modos(s)))(f,:) , 'DisplayName', strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm') ) ; 
%         hold on ; 
%     end
%     set(gca,"ColorOrderIndex",1,'FontSize',8)
%     for f = 21:2:31
%         plot(z , grafReverse.signal.(strcat("LP_",signal.modos(s)))(f,:) , '--' , 'DisplayName', strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm') )
%     end
%     xlabel(xlab,'FontSize',14) ; ylabel(ylab,'FontSize',14); title(strcat('Distribución Axial de la Señal LP',smodos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 6,'FontSize', 9);
%     annotation('textbox', [0.1254, 0.148, 0, 0], 'string', 'ForwardPump')
%     annotation('textbox', [0.1254, 0.128, 0, 0], 'string', 'BackwardPump')
% end



% % % % % % Bombeo
% for s = 1:length(pump.modos) % Grafico
%     graf.pump = EDFA.("Nucleo1").pump.Potencia_dBm;
%     grafReverse.pump = EDFA_RP.Nucleo1.pump.Potencia_dBm;
%     figure(3)
%     for f = 1:length(pump.lambda.(strcat("LP_",pump.modos(s))))
%         plot(z,graf.pump.(strcat("LP_",pump.modos(s)))(f,:) , 'DisplayName', strcat("Forward ", int2str(pump.lambda.(strcat("LP_",pump.modos(s)))(f)*1e9) ,' nm') ) ; 
%         hold on ; 
%     end
%     set(gca,'FontSize',8)
%     for f = 1:length(pump.lambda.(strcat("LP_",pump.modos(s))))
%         plot(z,grafReverse.pump.(strcat("LP_",pump.modos(s)))(f,:) , 'DisplayName', strcat("Backward " , int2str(pump.lambda.(strcat("LP_",pump.modos(s)))(f)*1e9) ,' nm') )
%     end
%     xlabel(xlab,'FontSize', 14) ; ylabel(ylab,'FontSize', 14); title('Distribución Axial del Bombeo','FontSize', 14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize', 9);
% end


% % % % % % Ganancias

% for s = 1:length(signal.modos)
%     graf.ganancias = EDFA.("Nucleo1").salida.ganancias;
%     ejex = signal.lambda.LP_11_a.*1e9;
%     plot(ejex,graf.ganancias.(strcat("LP_",signal.modos(s))) , '-o' , 'DisplayName',strcat("Forward Pump LP" , smodos(s)) )  ; hold on ;
% end
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% for s = 1:length(signal.modos)
%     
%     grafReverse.ganancias = EDFA_RP.("Nucleo1").salida.ganancias;
%     ejex = signal.lambda.LP_11_a.*1e9;
%     plot(ejex,grafReverse.ganancias.(strcat("LP_",signal.modos(s))) , '-*' , 'DisplayName', strcat("Backward Pump LP" , smodos(s)) )
%      
% end
% title('Distribución Espectral de Ganancias','FontSize',14) ; xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel("Magnitud [dB]",'FontSize',14) ;
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize', 9);
 

% % % % % % OSNR

% % MODO LP11a
% s = 1; graf.OSNR.LP_11_a = EDFA.("Nucleo1").OSNR.LP_11_a(:,end);
% ejex = signal.lambda.LP_11_a.*1e9;
% plot(ejex,graf.OSNR.(ModoS(s)) , '-o' , 'DisplayName',"Forward Pump" ) ; hold on ; 
% 
% set(gca,'FontSize',8)
% 
% grafReverse.OSNR.LP_11_a = EDFA_RP.("Nucleo1").OSNR.LP_11_a(:,end);
% ejex = signal.lambda.LP_11_a.*1e9;
% plot(ejex,grafReverse.OSNR.(ModoS(s)) , '-*' , 'DisplayName', "Backward Pump"  )
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize', 9)
% title(strcat('OSNR en Modo LP',smodos(s)),'FontSize', 14) ; %grid on ; 
% xlabel('Longitud de Onda [nm]','FontSize', 14) ; ylabel("Magnitud [dB]",'FontSize', 14) ; 

% % MODO LP11a
% s = 2; graf.OSNR.LP_21_a = EDFA.("Nucleo1").OSNR.LP_21_a(:,end);
% ejex = signal.lambda.LP_21_a.*1e9;
% plot(ejex,graf.OSNR.(ModoS(s)) , '-o' , 'DisplayName','Forward Pump' ); hold on ; 
% 
% set(gca,'FontSize',8)
% 
% grafReverse.OSNR.LP_21_a = EDFA_RP.("Nucleo1").OSNR.LP_21_a(:,end);
% ejex = signal.lambda.LP_21_a.*1e9;
% plot(ejex,grafReverse.OSNR.(ModoS(s)) , '-*' , 'DisplayName', "Backward Pump" )
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize', 9)
% title('OSNR en Modo LP21a','FontSize', 14) ; %grid on ; 
% xlabel('Longitud de Onda [nm]','FontSize', 14) ; ylabel("Magnitud [dB]",'FontSize', 14) ; 



% % % % % % Figuras de Ruido

% for s = 1:length(signal.modos)
%     graf.NF = EDFA.("Nucleo1").NF;
%     ejex = signal.lambda.LP_11_a.*1e9;
%     plot(ejex,graf.NF.(strcat("LP_",signal.modos(s))) , '-o' , 'DisplayName',strcat("Forward Pump LP",smodos(s)) ) ; hold on ; 
% end
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% for s = 1:length(signal.modos)
%     grafReverse.NF = EDFA_RP.("Nucleo1").NF;
%     ejex = signal.lambda.LP_21_a.*1e9;
%     plot(ejex,grafReverse.NF.(strcat("LP_",signal.modos(s))) , '-*' , 'DisplayName', strcat( "Backward Pump " ,"LP",smodos(s)) )
% end
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize', 9)
% title('Figuras de Ruido','FontSize', 14) ; %grid on ; 
% xlabel('Longitud de Onda [nm]','FontSize', 14) ; ylabel("Magnitud [dB]",'FontSize', 14) ; 




% % % % % ASE
% smodos = ["11a","21a"];
% for s = 1:length(signal.modos) % Grafico
%     figure(s)
%     graf.signal = EDFA.("Nucleo1").Pase;
%     grafReverse.signal = EDFA_RP.("Nucleo1").Pase;
%     var = 2;
%     for f = 21:2:31
%         plot(z , graf.signal.(strcat("LP_",signal.modos(s)))(f,:) , 'DisplayName', strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm') ) ; 
%         hold on ; 
%     end
%     set(gca,"ColorOrderIndex",1,'FontSize',8)
%     for f = 21:2:31
%         plot(z , grafReverse.signal.(strcat("LP_",signal.modos(s)))(f,:) , '--' , 'DisplayName', strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm') )
%     end
%     ylim([-50 -30])
%     xlabel(xlab,'FontSize',14) ; ylabel(ylab,'FontSize',14); title(strcat('Distribución Axial de la potencia ASE para el modo LP',smodos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 6,'FontSize', 9);
%     annotation('textbox', [0.1254, 0.148, 0, 0], 'string', 'ForwardPump')
%     annotation('textbox', [0.1254, 0.128, 0, 0], 'string', 'BackwardPump')
% end

