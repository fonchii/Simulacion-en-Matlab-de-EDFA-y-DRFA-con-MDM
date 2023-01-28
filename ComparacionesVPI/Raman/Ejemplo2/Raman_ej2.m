clear all; close all; clc

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])


% SAVE:
% print -dpdf 'NAME'

load("RamanVPIE2.mat")


%load("BestRaman_Matlab_EJ1_v0.mat")

c = 299792458;
wavelenghts_matlab = ( c./linspace(c/(1.66551366e-6) , c/(1.50046275e-6) , 100)).*1e9;

%GainDiff = abs( RamanVPIE1.GainOnOff - Raman.Sig.GainOnOFF.LP01 );

% % Ganancias

% plot(wavelenghts_matlab , Raman.Sig.GainOnOFF.LP01,"DisplayName","Matlab") ; hold on
%  plot((c./(RamanVPIE2.freqs)).*1e9,RamanVPIE2.GainOnOff,"DisplayName","VPI") 
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
%plot((Raman.Sig.GainOnOFF.LP01./max(Raman.Sig.GainOnOFF.LP01)).*max(RamanVPIE1.GainOnOff),"DisplayName","Matlab_ajustado") ; legend()

%% Señales

for i=0:6
    plot(Raman.z , Raman.Sig.Power.LP01(1+15*i,:)  ,"DisplayName",strcat(num2str(round(wavelenghts_matlab(1+15*i))), 'nm') ); hold on
end

set(gca,'ColorOrderIndex',1,'FontSize',8)

for i=0:6
    plot(RamanVPIE1.z , (1e-3.*10.^( RamanVPIE1.Signals(:,1+15*i)./10)) , '--' ,"DisplayName",strcat( num2str(round(RamanVPIE1.wavelength(1+15*i)*1e9) ), 'nm')  ) 
end

title('Distribución Axial de la Potencia de Señal','FontSize',14) ; xlabel('Posición en fibra [km]','FontSize',14) ; ylabel('Potencia [mW]','FontSize',14)
legend('Location', 'southoutside','Orientation','horizontal','Box','off','NumColumns',7,'FontSize',9)
annotation('textbox', [0.13, 0.148, 0, 0], 'string', 'Matlab')
annotation('textbox', [0.13, 0.127, 0, 0], 'string', 'VPIphotonics')


% % % Potencia
% 
% plot(RamanVPIE1.z , (1e-3.*10.^( RamanVPIE1.Pump./10)) , "DisplayName","VPI" ) ; hold on
% set(gca,'FontSize',8)
% plot(Raman.z , Raman.Pump.backward.LP01 , "DisplayName", "Matlab")
% title('Distribución Axial de la Potencia de Bombeo','FontSize',14) ; xlabel('Posición en fibra [km]','FontSize',14) ; ylabel('Potencia [mW]','FontSize',14)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off' , "NumColumns",2,"FontSize",9)


% % ASE


%plot(RamanVPIE1.z , (1e-3.*10.^( RamanVPIE1.ASE_Fwd./10)) ) ; hold on

% set(gca,'ColorOrderIndex',1,'FontSize',8)
% for i=0:6
%     plot(Raman.z , (1e-3.*10.^( Raman.Sig.Power.ASE.LP01(1+15*i,:)./10)) ) ; hold on
% end

%title('Distribución Axial de la Potencia ASE','FontSize',14) ; xlabel('Posición [km]','FontSize',14) ; ylabel('Potencia [mW]','FontSize',14)











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
% RamanVPIE1_0P.freq = RamanVPI_E1_0Param(:,1);
% RamanVPIE1_0P.wavelength = RamanVPI_E1_0Param(:,2);
% RamanVPIE1_0P.Gain = RamanVPI_E1_0Param(:,4);
% RamanVPIE1_0P.GainOnOff = RamanVPI_E1_0Param(:,4)+20;





% load("Raman_vpi.mat")
% % Frequency ,	Wavelength,	 SignalPowerIn, 	SignalPowerOut,	 NoisePowerOut,     OSNR_Out, 	 Gain, 	 NoiseFigureStandard, GainOnOff
% 
% RamanVPIE1.freq = Ejemplo1RamanVPI(:,1);
% RamanVPIE1.wavelength = Ejemplo1RamanVPI(:,2);
% RamanVPIE1.signal.In = Ejemplo1RamanVPI(:,3);
% RamanVPIE1.signal.Out = Ejemplo1RamanVPI(:,4);
% RamanVPIE1.ASE = Ejemplo1RamanVPI(:,5);
% RamanVPIE1.OsnrOut= Ejemplo1RamanVPI(:,6);
% RamanVPIE1.Gain = Ejemplo1RamanVPI(:,7);
% RamanVPIE1.NF = Ejemplo1RamanVPI(:,8);
% RamanVPIE1.GainOnOff= Ejemplo1RamanVPI(:,9);

