clear all; close all; clc

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])


% SAVE:
% print -dpdf 'NAME'

load("RamanVPIE3.mat")


%load("Raman_Matlab_EJ1_Gain3.mat")
%load("Raman_Matlab_EJ1_Gain3_390mw.mat")
%load("Raman_Matlab_EJ1_Gain3_ASE300mw.mat")
%load("Raman_Matlab_EJ1_Gain3backward.mat")
%load("Raman_Matlab_EJ1_Gain3backward_400mw.mat")

%load("Raman_MatlabE3.mat")
load("Raman_MatlabE3_ASE.mat")

c = 299792458;
wavelenghts_matlab = ( c./linspace(c/(1.66551366e-6) , c/(1.50046275e-6) , 100)).*1e9;

GainDiff = abs( RamanVPIE3.GainOnOff - Raman.Sig.GainOnOFF.LP01 );


% % Ganancias

% plot(wavelenghts_matlab , Raman.Sig.GainOnOFF.LP01,"DisplayName","Matlab") ; hold on
% plot(RamanVPIE3.wavelength.*1e9 , RamanVPIE3.GainOnOff,"DisplayName","VPI") 
% 
% set(gca,'FontSize',8)
% title("Distribución espectral de Ganancias",'FontSize',14) ; xlabel("Longitud de onda [nm]",'FontSize',14) ; ylabel("Ganancia [dBm]",'FontSize',14)
% 
%     % Diferencia de Ganancias
% yyaxis right
% ylabel('Diferencia de Ganancias','FontSize',14)
% plot(wavelenghts_matlab , GainDiff ,  "-o" , DisplayName = "Diferencia de Ganancias" , Color='#0072BD')
% ylim([0, 5]) ; %grid minor
% legend('Location', 'southoutside','Orientation','horizontal','Box','off' , "NumColumns",2,"FontSize",9)

%Analisis de forma:
%plot((Raman.Sig.GainOnOFF.LP01./max(Raman.Sig.GainOnOFF.LP01)).*max(RamanVPIE3.GainOnOff),"DisplayName","Matlab_ajustado") ; legend()

%% Señales

% for i=0:6
%     plot(Raman.z , Raman.Sig.Power.LP01(1+15*i,:)  ,"DisplayName",strcat(num2str(round(wavelenghts_matlab(1+15*i))), 'nm') ); hold on
% end
% 
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% 
% for i=0:6
%     plot(RamanVPIE3.z , (1e-3.*10.^( RamanVPIE3.Signals(:,1+15*i)./10)) , '--' ,"DisplayName",strcat( num2str(round(RamanVPIE3.wavelength(1+15*i)*1e9) ), 'nm')  ) 
% end
% 
% title('Distribución Axial de la Potencia de Señal','FontSize',14) ; xlabel('Posición en fibra [km]','FontSize',14) ; ylabel('Potencia [mW]','FontSize',14)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off','NumColumns',7,'FontSize',9)
% annotation('textbox', [0.13, 0.148, 0, 0], 'string', 'Matlab')
% annotation('textbox', [0.13, 0.127, 0, 0], 'string', 'VPIphotonics')


% % % Pump
% 
% plot(Raman.z , Raman.Pump.forward.LP01 , "DisplayName", "Matlab  Forward") ; hold on
% plot(Raman.z , Raman.Pump.backward.LP01 , "DisplayName", "Matlab Backward")
% 
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% 
% plot(RamanVPIE3.z , (1e-3.*10.^( RamanVPIE3.PumpFwd./10)) , '--' , "DisplayName","VPI Forward" ) 
% plot(RamanVPIE3.z , (1e-3.*10.^( RamanVPIE3.PumpBwd./10)) , '--' , "DisplayName","VPI Backward" )
% 
% title('Distribución Axial de la Potencia de Bombeo','FontSize',14) ; xlabel('Posición en fibra [km]','FontSize',14) ; ylabel('Potencia [mW]','FontSize',14)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off' , "NumColumns",2,"FontSize",9)


% % % ASE

% for i=0:6
%     plot(Raman.z , Raman.ASE.LP01(1+15*i,:)*1000  ,"DisplayName",strcat(num2str(round(wavelenghts_matlab(1+15*i))), 'nm') ); hold on
% end
% % 
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% ase_wl = c./RamanVPIE3.ASE_freqs;
% for i=0:6
%     plot(RamanVPIE3.z , RamanVPIE3.ASE_Fwd_mW(:,1+15*i) , '--' ,"DisplayName",strcat( num2str(round(ase_wl(RamanVPIE3.ASE_freqs_idx(1+15*i))*1e9) ), 'nm')  ) 
% end
% 
% title('Distribución Axial de la Potencia ASE','FontSize',14) ; xlabel('Posición [km]','FontSize',14) ; ylabel('Potencia [mW]','FontSize',14)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off','NumColumns',7,'FontSize',9)
% annotation('textbox', [0.13, 0.148, 0, 0], 'string', 'Matlab')
% annotation('textbox', [0.13, 0.127, 0, 0], 'string', 'VPIphotonics')



% % % OSNR
% OSNR_Matlab = 10*log10( Raman.Sig.Power.LP01(:,:)./Raman.ASE.LP01(:,:) );
% NF_Matlab = OSNR_Matlab(:,12)-OSNR_Matlab(:,end);
% 
% %signals_VPI = 1e-3.*10.^( RamanVPIE3.Signals(:,:)./10);
% OSNR_VPI = 10*log10( 1e-3.*10.^( RamanVPIE3.Signals(:,:)./10)./(RamanVPIE3.ASE_Fwd_mW(:,:)./1000) );
% NF_VPI = OSNR_VPI(3,:)-OSNR_VPI(end,:);
%  
% for i=0:6
%     plot(Raman.z , OSNR_Matlab(1+15*i,:)  ,"DisplayName",strcat(num2str(round(wavelenghts_matlab(1+15*i))), 'nm') ); hold on
% end
% % 
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% ase_wl = c./RamanVPIE3.ASE_freqs;
% for i=0:6
%     plot(RamanVPIE3.z , OSNR_VPI(:,1+15*i) , '--' ,"DisplayName",strcat( num2str(round(ase_wl(RamanVPIE3.ASE_freqs_idx(1+15*i))*1e9) ), 'nm')  ) 
% end
% 
% title('Distribución Axial de la OSNR','FontSize',14) ; xlabel('Posición [km]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off','NumColumns',7,'FontSize',9)
% annotation('textbox', [0.13, 0.148, 0, 0], 'string', 'Matlab')
% annotation('textbox', [0.13, 0.127, 0, 0], 'string', 'VPIphotonics')


% % % NF
OSNR_Matlab = 10*log10( Raman.Sig.Power.LP01(:,:)./Raman.ASE.LP01(:,:) );
NF_Matlab = OSNR_Matlab(:,13)-OSNR_Matlab(:,end);

%signals_VPI = 1e-3.*10.^( RamanVPIE3.Signals(:,:)./10);
OSNR_VPI = 10*log10( 1e-3.*10.^( RamanVPIE3.Signals(:,:)./10)./(RamanVPIE3.ASE_Fwd_mW(:,:)./1000) );
NF_VPI = OSNR_VPI(3,:)-OSNR_VPI(end,:);
 
plot(wavelenghts_matlab , NF_Matlab  ,"DisplayName",strcat('Matlab') ); hold on

plot(RamanVPIE3.wavelength.*1e9 , NF_VPI , "DisplayName",strcat( 'VPIPhotonics')  ) 


title('Figura de Ruido','FontSize',14) ; xlabel('Longitud de onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
legend('Location', 'southoutside','Orientation','horizontal','Box','off','NumColumns',7,'FontSize',9)








%% Cargar resultado vpi.. desde columna ec
% VPISignal = resultadovpipropagacionaxial;
% 
% for i=0:19
%     VPISignal_01(:,i+1) = VPISignal(:,2+4*i)/10e13;
%     VPISignal_11a(:,i+1) = VPISignal(:,4+4*i)/10e13;
% end
% VPISignal_01([1,2],:) = VPISignal_01([1,2],:).*10e13;
% VPISignal_11a([1,2],:) = VPISignal_11a([1,2],:).*10e13;



% load("Raman_vpi_0Param.mat")
% % Frequency	Wavelength	 OSNR_Out 	 Gain 	 NoiseFigureStandard 
% RamanVPIE3_0P.freq = RamanVPI_E1_0Param(:,1);
% RamanVPIE3_0P.wavelength = RamanVPI_E1_0Param(:,2);
% RamanVPIE3_0P.Gain = RamanVPI_E1_0Param(:,4);
% RamanVPIE3_0P.GainOnOff = RamanVPI_E1_0Param(:,4)+20;





% load("Raman_vpi.mat")
% % Frequency ,	Wavelength,	 SignalPowerIn, 	SignalPowerOut,	 NoisePowerOut,     OSNR_Out, 	 Gain, 	 NoiseFigureStandard, GainOnOff
% 
% RamanVPIE3.freq = Ejemplo1RamanVPI(:,1);
% RamanVPIE3.wavelength = Ejemplo1RamanVPI(:,2);
% RamanVPIE3.signal.In = Ejemplo1RamanVPI(:,3);
% RamanVPIE3.signal.Out = Ejemplo1RamanVPI(:,4);
% RamanVPIE3.ASE = Ejemplo1RamanVPI(:,5);
% RamanVPIE3.OsnrOut= Ejemplo1RamanVPI(:,6);
% RamanVPIE3.Gain = Ejemplo1RamanVPI(:,7);
% RamanVPIE3.NF = Ejemplo1RamanVPI(:,8);
% RamanVPIE3.GainOnOff= Ejemplo1RamanVPI(:,9);


% Get ASE signals

% ASE_freqs vs freqs idx get
% freqs_idx = zeros(1,length(RamanVPIE3.freqs));
% for j = 1:length(RamanVPIE3.freqs)
%     prev_diff = 1000000000000;
%     for i = 1:length(RamanVPIE3.ASE_freqs)
%         diff = abs(RamanVPIE3.ASE_freqs(i) - RamanVPIE3.freqs(j));
%         if diff < prev_diff
%             i
%             prev_diff = diff;          
%             freqs_idx(j) = i;
%         end
%     end
% end
% 
% for i=3:102 
%     vpi_aseFwd(:,i-2) = 10.^( RamanVPIE3.ASE_Fwd(:,i) ./10 );
% end


