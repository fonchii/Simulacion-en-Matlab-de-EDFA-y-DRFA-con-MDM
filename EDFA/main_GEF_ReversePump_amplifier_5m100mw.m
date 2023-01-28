% close all; 
clear all; clc ; close all

% Parámetros de entrada

c=299.792458e6; % [m/s]
h=6.62607015*10^(-34);

    % Modos y canales de señal y bombeo

signal.modos = ["01" "11_a"] ;

Signal.NumberOfChannels=31;
Wavelength_gridS=linspace(1530,1560,Signal.NumberOfChannels).*1e-9;



Pin=-20; %[dBm]

signal.lambda.LP_01     = Wavelength_gridS;                                   P0_signal.LP_01     = Pin*ones(1,length(signal.lambda.LP_01));
signal.lambda.LP_11_a   = Wavelength_gridS;                                     P0_signal.LP_11_a   = Pin*ones(1,length(signal.lambda.LP_11_a));

pump.modos = "01" ;
Wavelength_gridP=980e-9;
Ppump= 100e-3; %[W]

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
ASE= -200;

    % Datos de la fibra
fibra.nucleos = 1;
fibra.largo = 5; fibra.radio = 5e-6 ; fibra.N = 7e24; 

fibra.n1 = 1.45 ;   fibra.IndexContrast=0.01;
fibra.AN=fibra.n1*sqrt(2*fibra.IndexContrast);
fibra.n2 =sqrt((fibra.n1^2-fibra.AN^2));

P.Np=2; P.Fc=c/Wavelength_gridS(ceil(length(Wavelength_gridS)/2)); P.Fb = 50e9; 
fibra.dvk=P.Fb;
fibra.PumpMode = "reverse";
fibra.GEF = load('GEF_Filters/EDFA_5m_100mw_OptiSystem') ; % GEF_Filters/EDFA_10m_150mw_OptiSystem

fibra.SigmaSpectrum = "VPI";
fibra.WaitBar = 1; fibra.Avance = 1;    % Despliegue de info
fibra.ASEFlag = 0;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)

%%
tic;
EDFA_RP = EDFA_ReversePump_MMvPCCv3_GEF(fibra,signal,pump,ASE);
EDFA = EDFA_MMvPCCv3_GEF(fibra,signal,pump,ASE);
t_end = toc; fprintf('Tiempo de cómputo: %.2f segundos\n', t_end);


    %% Graficos
% SAVE:
set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'NAME'


close all 
Nc = length(fieldnames(EDFA)); 
z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';


% Señal
smodos = ["11a","21a"];
for s = 1:length(signal.modos) % Grafico

    figure(s)

    graf.signal = EDFA.("Nucleo1").signal.Potencia_dBm;
    grafReverse.signal = EDFA_RP.("Nucleo1").signal.Potencia_dBm;
    var = 2;
    for f = 21:2:31
        plot(z , graf.signal.(strcat("LP_",signal.modos(s)))(f,:) , 'DisplayName', strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm') ) ; 
        hold on ; 
    end
    set(gca,"ColorOrderIndex",1,'FontSize',8)
    for f = 21:2:31
        plot(z , grafReverse.signal.(strcat("LP_",signal.modos(s)))(f,:) , '--' , 'DisplayName', strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm') )
    end
    xlabel(xlab,'FontSize',14) ; ylabel(ylab,'FontSize',14); title(strcat('Distribución Axial de la Señal LP',smodos(s)),'FontSize',14) ; 
    legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7,'FontSize', 9);
    annotation('textbox', [0.1254, 0.148, 0, 0], 'string', 'ForwardPump')
    annotation('textbox', [0.1254, 0.128, 0, 0], 'string', 'BackwardPump')
end



% Bombeo
% for s = 1:length(pump.modos) % Grafico
%     graf.pump = EDFA.("Nucleo1").pump.Potencia_dBm;
%     grafReverse.pump = EDFA_RP.(strcat("Nucleo",int2str(n))).pump.Potencia_dBm;
%     figure(3)
%     for f = 1:length(pump.lambda.(strcat("LP_",pump.modos(s))))
%         plot(z,graf.pump.(strcat("LP_",pump.modos(s)))(f,:) , 'DisplayName', strcat("Forward ", int2str(pump.lambda.(strcat("LP_",pump.modos(s)))(f)*1e9) ,' nm') ) ; 
%         hold on ; 
%     end
%     set(gca,'FontSize',8)
%     for f = 1:length(pump.lambda.(strcat("LP_",pump.modos(s))))
%         plot(z,grafReverse.pump.(strcat("LP_",pump.modos(s)))(f,:) , 'DisplayName', strcat("Backward " , int2str(pump.lambda.(strcat("LP_",pump.modos(s)))(f)*1e9) ,' nm') )
%     end
%     xlabel(xlab,'FontSize', 14) ; ylabel(ylab,'FontSize', 14); title('Distribución Axial del Bombeo LP12a','FontSize', 14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize', 9);
% end


% % Ganancias
figure(4)
for s = 1:length(signal.modos)
    graf.ganancias = EDFA.("Nucleo1").salida.ganancias;
    ejex = signal.lambda.(strcat('LP_',signal.modos(1))).*1e9;
    plot(ejex,graf.ganancias.(strcat("LP_",signal.modos(s))) , '-o' , 'DisplayName',strcat("LP",signal.modos(s)) ) ; hold on ;
end
set(gca,'ColorOrderIndex',1,'FontSize',8)
for s = 1:length(signal.modos)
    
    grafReverse.ganancias = EDFA_RP.("Nucleo1").salida.ganancias;
    ejex = signal.lambda.(strcat('LP_',signal.modos(1))).*1e9;
    plot(ejex,grafReverse.ganancias.(strcat("LP_",signal.modos(s))) , '-*' , 'DisplayName', strcat("Backward Pump " , "LP",signal.modos(s)) )
     
end
title('Ganancias','FontSize',14) ; xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel("Magnitud [dB]",'FontSize',14) ;
legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize', 9);


% % Figuras de Ruido
%figure(n+2+length(signal.modos)-1)
figure(5)

for s = 1:length(signal.modos)
    graf.NF = EDFA.("Nucleo1").NF;
    ejex = signal.lambda.(strcat('LP_',signal.modos(1))).*1e9;
    plot(ejex,graf.NF.(strcat("LP_",signal.modos(s))) , '-o' , 'DisplayName',strcat("LP",smodos(s)) ) ; hold on ; 
end
set(gca,'ColorOrderIndex',1,'FontSize',8)
for s = 1:length(signal.modos)
    grafReverse.NF = EDFA_RP.("Nucleo1").NF;
    ejex = signal.lambda.(strcat('LP_',signal.modos(1))).*1e9;
    plot(ejex,grafReverse.NF.(strcat("LP_",signal.modos(s))) , '-*' , 'DisplayName', strcat( "Backward Pump " ,"LP",smodos(s)) )
end
legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize', 9)
title('Figuras de Ruido','FontSize', 14) ; %grid on ; 
xlabel('Longitud de Onda [nm]','FontSize', 14) ; ylabel("Magnitud [dB]",'FontSize', 14) ; 






%%    
%     % Espectro de ganancia
%     for ase = 1:length(Signal.modos)
%         freq_ase = EDFA.(strcat('Nucleo',int2str(n))).ASE_Spectrum.lambdas*1e9;
%         graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase))) = EDFA.(strcat("Nucleo",int2str(n))).ASE_Spectrum.mag(:,:,ase);
%         subplot(2,2,[3,4])
%         plot( freq_ase, graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase)))(:,2) , 'DisplayName', strcat( "LP",signal.modos(ase) ) ) ; hold on;
%         xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Potencia a la salida del EDFA') ; legend()
%     end
    

% %% Graficos Modos
% 
% % Señal
% n=n+1 ; figure(n)
% graficar_modos(fibra,signal.modos,signal.lambda.LP_01)
% sgtitle('Modos de Señal') 
% 
% % Bombeo
% n=n+1 ; figure(n)
% graficar_modos(fibra,pump.modos,pump.lambda.LP_01)
% sgtitle('Modos de Bombeo') ; %clear n ;

%% Graficar potencia vs frecuencia
% n=n+1; figure(n)
% figure(2)
% xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
% graf.Nc = length(fieldnames(EDFA)); 
% for n = 1:graf.Nc
%     for s = 1:length(signal.modos)
%         graf.ganancias.(ModoS(s)) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(ModoS(s));
%         leyenda = strcat((ModoS(s)));
%         ejex = fliplr((3*10^8./(signal.lambda.(ModoS(s)).*1e9)));
%         plot(ejex/1000,fliplr(graf.ganancias.(ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%         %ax.ColorOrderIndex = s;
%         xlabel('Frecuencia [THz]') ; ylabel(ylab)
% 
%     end ; clear s ejex leyenda;
% end


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% archivos de salida
% aux_out1=[graf.freq_ase' graf.ASE_Spectrum.LP_01(:,1) graf.ASE_Spectrum.LP_01(:,2) graf.ASE_Spectrum.LP_11_a(:,1) graf.ASE_Spectrum.LP_11_a(:,2)]; 
% aux_out1=[Frequency_gridS' Wavelength_gridS' EDFA.Nucleo1.salida.signal.potencia_dBm.LP_01 EDFA.Nucleo1.salida.signal.potencia_dBm.LP_11_a];
% aux_out2= [graf.ganancias.LP_01' graf.ganancias.LP_11_a'];
% 
% save 'C:\Users\HP -\Documents\UTFSM\Investigación Doctorado\Simulaciones Matlab\registro de resultados\Datos importados del simulador FMF\noise_simul.dat' aux_out1 -ascii
% save 'C:\Users\HP -\Documents\UTFSM\Investigación Doctorado\Simulaciones Matlab\registro de resultados\Datos importados del simulador FMF\signal_simul.dat' aux_out1 -ascii
% save 'C:\Users\HP -\Documents\UTFSM\Investigación Doctorado\Simulaciones Matlab\registro de resultados\Datos importados del simulador FMF\gain_simul.dat' aux_out2 -ascii



%% GUARDAR datos

% S01 = EDFA.Nucleo1.signal.Potencia_dBm.(ModoS(1));
% S11_a = EDFA.Nucleo1.signal.Potencia_dBm.(ModoS(2));
% 
% frecuencias = Frequency_gridS';
% lambdas = Wavelength_gridS';
% G1 = graf.ganancias.LP_01';
% G2 = graf.ganancias.LP_11_a';
% NF01  = EDFA.Nucleo1.NF.LP_01; 
% NF11a = EDFA.Nucleo1.NF.LP_11_a;
% 
% DAT = [lambdas G1 lambdas G2 frecuencias NF01 NF11a];

