% close all; 
clear all; clc

% Parámetros de entrada

    % Modos y canales de señal y bombeo
signal.NumberOfChannels=50;
signal.modos = ["01" "11_a"] ;
Frequency_gridS=linspace(191.19421875e12,193.64421875e12,signal.NumberOfChannels);
c=299.792458e6; % [m/s]
Wavelength_gridS=c./Frequency_gridS;

Pin=-10; %[dBm]
%Pin=0; %[dBm]

signal.lambda.LP_01     = Wavelength_gridS;                                  P0_signal.LP_01     = Pin*ones(1,length(signal.lambda.LP_01));
signal.lambda.LP_11_a   = Wavelength_gridS;                                  P0_signal.LP_11_a   = Pin*ones(1,length(signal.lambda.LP_11_a));


pump.modos = "12_a" ;%"01" ;
Wavelength_gridP=980e-9;
Ppump= 250e-3; %[W] 
%Ppump= 1000e-3; %[W]

pump.lambda.LP_12_a   = Wavelength_gridP;                         P0_pump.LP_12_a   = Ppump  ;  

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
P.Fc=c/Wavelength_gridS(ceil(length(Wavelength_gridS)/2)); P.Fb = 300e9; 
% Ruido ASE Entrada
ASE= -58;%-200;

    % Datos de la fibra
fibra.nucleos = 1;                                           % Numero de nucleos
%fibra.largo = 3; fibra.radio = 5e-6 ; fibra.N = 7e24; % fibra.N = 3e24; 
fibra.largo = 5; fibra.radio = 5.5e-6 ; fibra.N = 7e24; 
%fibra.AN = 0.2;   %fibra.AN = 0.2 ; %fibra.AN=sqrt(fibra.n1^2-fibra.n2^2);
fibra.n1 = 1.45 ;   fibra.IndexContrast=0.01;
fibra.AN=fibra.n1*sqrt(2*fibra.IndexContrast);
fibra.n2 =sqrt((fibra.n1^2-fibra.AN^2));
%fibra.M = 10; fibra.Nalpha = inf;
%fibra.ASEFlag = 1; % EVITA CALCULO DE ESPECTRO ASE
fibra.dvk=P.Fb;
fibra.ASEFlag = 0;
fibra.Avance = 1; fibra.WaitBar = 1;

%%
tic;
EDFA = EDFA_MMvPCCv3_Compare(fibra,signal,pump,ASE); %21.45s
% EDFA = EDFA_MM(fibra,signal,pump,ASE); %197.82 segundos
t_end = toc; fprintf('Tiempo de cómputo: %.2f segundos\n', t_end);


    %% Graficos

close all

%% Graficar potencia vs frecuencia
% n=n+1; figure(n)
figure(1)
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
graf.Nc = length(fieldnames(EDFA)); 
for n = 1:graf.Nc
    for s = 1:length(signal.modos)
        graf.ganancias.(ModoS(s)) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(ModoS(s));
        leyenda = strcat("LP",signal.modos(s));%strcat((ModoS(s)));
        
        %THz
        %ejex = fliplr((3*10^8./(signal.lambda.(ModoS(s)).*1e9))); 
        %plot(ejex/1000,fliplr(graf.ganancias.(ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
        % xlabel('Frecuencia [THz]') ; ylabel(ylab)

        % nm
        ejex = signal.lambda.(ModoS(s)).*1e9; %nm
        plot(ejex,fliplr(graf.ganancias.(ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
        %ax.ColorOrderIndex = s;
        xlabel('Longitud de onda [nm]') ; ylabel(ylab)

    end ; clear s ejex leyenda;
end


%% % % NF Y OSNR

% SE AJUSTO Ppump DE 1W A 250mw !! y Pin signal a -10dbm
%close all
figure(2)
% 
% Pan = EDFA.Nucleo1.Pan;
Pap = EDFA.Nucleo1.Pap;
%Pase = EDFA.Nucleo1.Pase;

OSNR_01  = EDFA.Nucleo1.OSNR.LP_01;
OSNR_11a = EDFA.Nucleo1.OSNR.LP_11_a;
NF_01  = OSNR_01(:,3) - OSNR_01(:,end);
NF_11a  = OSNR_11a(:,3) - OSNR_11a(:,end);

ejex = signal.lambda.(ModoS(1)).*1e9;

% NoiseFigure
plot(ejex , NF_01 , "DisplayName" , "Modo LP01" ) , hold on
plot(ejex , NF_11a , "DisplayName" , "Modo LP11a" )
set(gca,'ColorOrderIndex',1,'FontSize',13)
legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize',13)
xlabel('Longitud de Onda [nm]','FontSize',16) ; ylabel('NF [dB]','FontSize',16); title('Figura de Ruido','FontSize',18)
set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'Pump12a_NF_Ajustada'

figure(3)
%close all
% Optical Signal Noise Ratio
for i=1:8:length(OSNR_01)
    temp_osnr = OSNR_01(i,:);
    plot(EDFA.Nucleo1.z(3:end) , temp_osnr(3:end) , "DisplayName" , strcat( int2str(ejex(i)) , " nm") ) ; hold on
end
set(gca,'ColorOrderIndex',1,'FontSize',13)
legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 5,'FontSize',13)
xlabel('Longitud de Onda [nm]','FontSize',16) ; ylabel('OSNR [dB]','FontSize',16); title('Distribución axial de la OSNR','FontSize',18)
set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'Pump12a_OSNR_Ajustada'




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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Datos VPI
% p_pump_vpi=vpiedfamm.axialdist.Pump;
% p_sig_LP01_vpi=vpiedfamm.axialdist.Sig.LP_01;
% p_sig_LP11a_vpi=vpiedfamm.axialdist.Sig.LP_11a;
% p_ase_LP01_vpi=vpiedfamm.axialdist.ASE.LP_01;
% p_ase_LP11a_vpi=vpiedfamm.axialdist.ASE.LP_11a;
% 
% Z=linspace(0,3,101);
% colorlist=[1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0 0 0; 0.85 0.3250 0.980; 0.4940 0.184 0.5560; 0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840];
% 
% % Figuras
% % Pump
% figure(1)
% plot(Z,EDFA.Nucleo1.pump.Potencia_dBm.LP_01,'DisplayName','Matlab'); hold on
% plot(Z,p_pump_vpi,'DisplayName','VPI')
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend()
% title('Power axial distibution Pump (LP01,980 [nm], 250 [mW] ')
% 
% %Signal primeras 10 frecuencias LP01
% figure(2)
% for aux_i=1:length(colorlist)-1
%     plot(Z,EDFA.Nucleo1.signal.Potencia_dBm.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(Frequency_gridS(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i,:));hold on
%     plot(Z,p_sig_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Signal (LP01,-10[dBm])')
% 
% %Signal últimas 10 frecuencias LP01
% figure(3)
% for aux_i=1+10:length(colorlist)-1+10
%     plot(Z,EDFA.Nucleo1.signal.Potencia_dBm.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(Frequency_gridS(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-10,:));hold on
%     plot(Z,p_sig_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-10,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Signal (LP01,-10[dBm])')
% 
% % ASE forward LP01
% figure(4)
% for aux_i=1:length(colorlist)
%     plot(Z,EDFA.Nucleo1.PaseF.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i,:));hold on
%     plot(Z,p_ase_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP01,-200 [dBm])')
% 
% figure(5)
% for aux_i=1+11:length(colorlist)+11
%     plot(Z,EDFA.Nucleo1.PaseF.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-11,:));hold on
%     plot(Z,p_ase_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-11,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP01,-200 [dBm])')
% 
% figure(6)
% for aux_i=1+11*2:length(colorlist)+11*2
%     plot(Z,EDFA.Nucleo1.PaseF.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-11*2,:));hold on
%     plot(Z,p_ase_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-11*2,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP01,-200 [dBm])')
% 
% % LP_11a
% %Signal LP01a
% figure(7)
% for aux_i=1:length(colorlist)-1
%     plot(Z,EDFA.Nucleo1.signal.Potencia_dBm.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(Frequency_gridS(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i,:));hold on
%     plot(Z,p_sig_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Signal (LP11a,-10[dBm])')
% 
% %Signal últimas 10 frecuencias LP01
% figure(8)
% for aux_i=1+10:length(colorlist)-1+10
%     plot(Z,EDFA.Nucleo1.signal.Potencia_dBm.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(Frequency_gridS(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-10,:));hold on
%     plot(Z,p_sig_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-10,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Signal (LP11a,-10[dBm])')
% 
% % ASE forward LP01
% figure(9)
% for aux_i=1:length(colorlist)
%     plot(Z,EDFA.Nucleo1.PaseF.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i,:));hold on
%     plot(Z,p_ase_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP11a,-200 [dBm])')
% 
% figure(10)
% for aux_i=1+11:length(colorlist)+11
%     plot(Z,EDFA.Nucleo1.PaseF.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-11,:));hold on
%     plot(Z,p_ase_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-11,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP11a,-200 [dBm])')
% 
% figure(11)
% for aux_i=1+11*2:length(colorlist)+11*2
%     plot(Z,EDFA.Nucleo1.PaseF.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-11*2,:));hold on
%     plot(Z,p_ase_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-11*2,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP11a,-200 [dBm])')
% 
