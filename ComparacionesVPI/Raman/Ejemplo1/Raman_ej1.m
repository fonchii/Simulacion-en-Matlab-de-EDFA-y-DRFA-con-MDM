clear all; close all; clc

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])


% SAVE:
% print -dpdf 'NAME'

load("RamanVPIE1.mat")



%load("Raman_Matlab_EJ1_Gain3.mat")
%load("Raman_Matlab_EJ1_Gain3_390mw.mat")
%load("Raman_Matlab_EJ1_Gain3_ASE300mw.mat")
%load("Raman_Matlab_EJ1_Gain3backward.mat")
%load("Raman_Matlab_EJ1_Gain3backward_400mw.mat")


load("BestRaman_Matlab_EJ1_v0.mat")
%load("BestRaman_Matlab_EJ1_ASE.mat")

c = 299792458;
wavelenghts_matlab = ( c./linspace(c/(1.66551366e-6) , c/(1.50046275e-6) , 100)).*1e9;

GainDiff = abs( RamanVPIE1.GainOnOff - Raman.Sig.GainOnOFF.LP01 );

% % Ganancias

plot(wavelenghts_matlab , Raman.Sig.GainOnOFF.LP01,"DisplayName","Matlab") ; hold on
plot(RamanVPIE1.wavelength.*1e9 , RamanVPIE1.GainOnOff,"DisplayName","VPI") 

set(gca,'FontSize',8)
title("Distribución espectral de Ganancias",'FontSize',14) ; xlabel("Longitud de onda [nm]",'FontSize',14) ; ylabel("Ganancia [dBm]",'FontSize',14)

    % Diferencia de Ganancias
yyaxis right
ylabel('Diferencia de Ganancias','FontSize',14)
plot(wavelenghts_matlab , GainDiff ,  "-o" , DisplayName = "Diferencia de Ganancias" , Color='#0072BD')
ylim([0, 5]) ; %grid minor
legend('Location', 'southoutside','Orientation','horizontal','Box','off' , "NumColumns",2,"FontSize",9)

%Analisis de forma:
%plot((Raman.Sig.GainOnOFF.LP01./max(Raman.Sig.GainOnOFF.LP01)).*max(RamanVPIE1.GainOnOff),"DisplayName","Matlab_ajustado") ; legend()

%% Señales

% for i=0:6
%     plot(Raman.z , Raman.Sig.Power.LP01(1+15*i,:)  ,"DisplayName",strcat(num2str(round(wavelenghts_matlab(1+15*i))), 'nm') ); hold on
% end
% 
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% 
% for i=0:6
%     plot(RamanVPIE1.z , (1e-3.*10.^( RamanVPIE1.Signals(:,1+15*i)./10)) , '--' ,"DisplayName",strcat( num2str(round(RamanVPIE1.wavelength(1+15*i)*1e9) ), 'nm')  ) 
% end
% 
% title('Distribución Axial de la Potencia de Señal','FontSize',14) ; xlabel('Posición en fibra [km]','FontSize',14) ; ylabel('Potencia [mW]','FontSize',14)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off','NumColumns',7,'FontSize',9)
% annotation('textbox', [0.13, 0.148, 0, 0], 'string', 'Matlab')
% annotation('textbox', [0.13, 0.127, 0, 0], 'string', 'VPIphotonics')


% % % Potencia
% 
% plot(RamanVPIE1.z , (1e-3.*10.^( RamanVPIE1.Pump./10)) , "DisplayName","VPI" ) ; hold on
% set(gca,'FontSize',8)
% plot(Raman.z , Raman.Pump.backward.LP01 , "DisplayName", "Matlab")
% title('Distribución Axial de la Potencia de Bombeo','FontSize',14) ; xlabel('Posición en fibra [km]','FontSize',14) ; ylabel('Potencia [mW]','FontSize',14)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off' , "NumColumns",2,"FontSize",9)


% % % ASE
% plot(Raman.z , ( Raman.ASE.LP01(1,:))*1000 ,"DisplayName",strcat(num2str(round(wavelenghts_matlab(1))), 'nm - Matlab') ) ; hold on
% 
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% 
% plot(RamanVPIE1.z , (1e-3.*10.^( RamanVPIE1.ASE_Fwd./10)) , '--',"DisplayName",strcat( num2str(round(RamanVPIE1.wavelength(1)*1e9)) , 'nm - VPIphotonics' ) ) ; hold on

% for i=0:6
%     %plot(Raman.z , (1e-3.*10.^( Raman.Sig.Power.ASE.LP01(1+15*i,:)./10)) ) ; hold on
%     plot(Raman.z , ( Raman.ASE.LP01(1+15*i,:)) ) ; hold on
% end

% title('Distribución Axial de la Potencia ASE','FontSize',14) ; xlabel('Posición [km]','FontSize',14) ; ylabel('Potencia [mW]','FontSize',14)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off','NumColumns',1,'FontSize',9)


% % % OSNR
% OSNR_Matlab = 10*log10( Raman.Sig.Power.LP01(1,:)./Raman.ASE.LP01(1,:) );
% 
% OSNR_VPI = 10*log10( 1e-3.*10.^( RamanVPIE1.Signals(:,1)./10)./((1e-3.*10.^( RamanVPIE1.ASE_Fwd./10))./1000) );
%  
% plot(Raman.z , OSNR_Matlab(1,:)  ,"DisplayName",strcat(num2str(round(wavelenghts_matlab(1))), 'nm - Matlab') ); hold on
% 
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% 
% plot(RamanVPIE1.z , OSNR_VPI(:,1) , '--' ,"DisplayName",strcat( num2str(round(RamanVPIE1.wavelength(1)*1e9) ), 'nm - VPIphotonics')  ) 
% 
% title('Distribución Axial de la OSNR','FontSize',14) ; xlabel('Posición [km]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off','NumColumns',1,'FontSize',9)






%% GANANCIAS

% figure(1)

% plot(DAT(:,1).*1e9,DAT(:,2), "DisplayName","LP01 - Simulador Matlab" , Color= '#0072BD') ; hold on
% plot(AmplifiedSignals(:,1),AmplifiedSignals(:,2),"DisplayName","LP01 - VPIPhotonics" , LineStyle="--", Color= '#0072BD' ) ; 
% plot(DAT(:,1).*1e9,DAT(:,4),  "DisplayName","LP11a - Simulador Matlab", "Color", '#D95319' )
% plot(AmplifiedSignals(:,1),AmplifiedSignals(:,4),"DisplayName","LP11a - VPIPhotonics", LineStyle="--" , Color='#D95319') 
% xlabel('Longitud de Onda [nm]') ; ylabel('Ganancia [dB]') ; title("Ganancia EDFA con Bombeo en LP12a") ; 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4)
% 
% DG01 =  abs(DAT(:,2) - AmplifiedSignals(:,2)) ;
% DG11 = abs(DAT(:,4) - AmplifiedSignals(:,4)) ;
%     % Diferencia de Ganancias
% yyaxis right
% ylabel('Diferencia de Ganancias')
% plot(DAT(:,1).*1e9 , DG01 ,  "-o" , DisplayName = "DG - LP01" , Color='#0072BD')
% plot(DAT(:,1).*1e9 , DG11 ,  "-o" , DisplayName = "DG - LP11a" , Color='#D95319')
% ylim([0, 2]) ; %grid minor




%% NOISE FIGURE
% figure(2)
% 
% plot(DAT(:,1).*1e9,DAT(:,6), "DisplayName","LP01 - Simulador" , Color= '#0072BD') ; 
% hold on
% plot(AmplifiedSignals(:,1),fliplr(ModeGainsandNFs(:,4)),"DisplayName","LP01 - VPIPhotonics" , LineStyle="--", Color= '#0072BD' ) ; 
% plot(DAT(:,1).*1e9,DAT(:,7),  "DisplayName","LP11a - Simulador", "Color", '#D95319' )
% plot(AmplifiedSignals(:,1),fliplr(ModeGainsandNFs(:,8)),"DisplayName","LP11a - VPIPhotonics", LineStyle="--" , Color='#D95319') 
% xlabel('Longitud de Onda [nm]') ; ylabel('Noise Figure [dB]') ; title("Comparación de NF obtenidas mediante Simulador y VPIphotonics") ; 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off')
% 
% 
% 
% VPISignal(end,1) = VPISignal(end,1)/1e10;
% zVPI = VPISignal(:,1);
% %lambda = DAT(1,1);  % ORDENADOS DE MAYOR A MENOR, IGUAL QE LOS DE VPISIGNAL
% 
% 
% colors = ["#0072BD" , "#D95319" , "#EDB120" , "#7E2F8E" , "#77AC30" , "#4DBEEE" , "#A2142F"];
% 
% for i = 0
%     lambda = DAT(1+i,1);
%     plot(z,S01(1+i,:), "DisplayName",strcat("LP01 - Simulador @", int2str(lambda*1e9) , "nm"), Color= colors(i+1)) ; 
%     hold on
%     plot(zVPI,VPISignal_01(:,1+i),"DisplayName",strcat("LP01 - VPIPhotonics",int2str(lambda*1e9) , "nm" ) , LineStyle="--", Color= colors(i+1) ) ;
% end
% %plot(DAT(:,1).*1e9,DAT(:,7),  "DisplayName","LP11a - Simulador", "Color", '#D95319' )
% %plot(AmplifiedSignals(:,1),fliplr(ModeGainsandNFs(:,8)),"DisplayName","LP11a - VPIPhotonics", LineStyle="--" , Color='#D95319') 
% xlabel('Posición en EDFA [m]') ; ylabel('Potencia [dBm]') ; title("Comparación de Evolución Potencias obtenidas mediante Simulador y VPIphotonics") ; 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off')


%% Propagación

% % %              Propagación SEÑAL LP01

% z = linspace(0,5,length(S01(1,:)));
% zvpi = linspace(0,5,length(vpiedfamm.axialdist.Sig.LP_01(:,1)));
% 
% 
% for i=0:6
%     plot(z , EDFA.Nucleo1.signal.Potencia_dBm.LP_01(1+8*i,:) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm")) ; hold on
% end
% set(gca,'ColorOrderIndex',1,'FontSize', 8)
% for i=0:6
%     plot(zvpi , vpiedfamm.axialdist.Sig.LP_01(:,1+8*i) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm") , LineStyle="--" )
% end
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7,'FontSize', 9)
% xlabel('Posición Axial [m]','FontSize', 14) ; ylabel('Potencia [dBm]','FontSize', 14); title('Distribución Axial de la Potencia de Señal LP01','FontSize', 14)


%               Propagacion SEÑAL LP11

% z = linspace(0,5,length(S01(1,:)));
% zvpi = linspace(0,5,length(vpiedfamm.axialdist.Sig.LP_01(:,1)));
% 
% 
% for i=0:6
%     plot(z , EDFA.Nucleo1.signal.Potencia_dBm.LP_11_a(1+8*i,:) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm")) ; hold on
% end
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% for i=0:6
%     plot(zvpi , vpiedfamm.axialdist.Sig.LP_11a(:,1+8*i) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm") , LineStyle="--" )
% end
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7,'FontSize', 9)
% xlabel('Posición Axial [m]','FontSize', 14) ; ylabel('Potencia [dBm]','FontSize', 14); title('Distribución Axial de la Potencia de Señal LP11a','FontSize', 14)

%               Propagación Bombeo

% z = linspace(0,5,length(S01(1,:)));
% zvpi = linspace(0,5,length(vpiedfamm.axialdist.Sig.LP_01(:,1)));
% 
% 
% plot(z , EDFA.Nucleo1.pump.Potencia_dBm.LP_12_a(1,:) , "DisplayName", "Matlab") ; hold on
% 
% plot(zvpi , vpiedfamm.axialdist.Pump(:,1) , "DisplayName", "VPIphotonics" )% ,LineStyle="--" )
% 
% set(gca,'FontSize', 8)
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7,'FontSize', 9)
% xlabel('Posición Axial [m]','FontSize', 14) ; ylabel('Potencia [dBm]','FontSize', 14); title('Distribución Axial de la Potencia de Bombeo LP12a','FontSize', 14)
% ylim([8 31])


%% %               Propagación ASE , OSNR , NF

% %         LP01
% % % ASE+ 
%  %%  Matlab
% for i=0:6
%     plot(z , Pap.LP_01(1+8*i,:) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm")) ; hold on
% end
% % % %%  VPI
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% for i=1:length(vpiasePos)
%     plot(zvpi , vpiedfamm.axialdist.ASE.LP_01(:,vpiasePos(i)) , "DisplayName", strcat(int2str(vpiedfamm.ASElam(1,vpiasePos(i))*1e9) , " nm") , LineStyle="--" )
% end
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7,'FontSize', 9)
% xlabel('Posición Axial [m]','FontSize', 14) ; ylabel('Potencia [dBm]','FontSize', 14); title('Distribución Axial de la Potencia de Ruido ASE+ en modo LP01','FontSize', 14)
% ylim([-90 -30])
% 
% maxAseDif = -100000;
% for i = 0:6
%     if maxAseDif < abs( Pap.LP_01(1+8*i,end) - vpiedfamm.axialdist.ASE.LP_01(end,vpiasePos(i+1)) )
%         maxAseDif = abs( Pap.LP_01(1+8*i,end) - vpiedfamm.axialdist.ASE.LP_01(end,vpiasePos(i+1)) ) ;
%         maxdiffLambda = (DAT(1+8*i,1)*1e9);
%     end
% end
% %dB
% for i = 0:6
%     if maxAseDif < abs( 1e-3*10.^(Pap.LP_01(1+3*i,end)) - 1e-3*10.^(vpiedfamm.axialdist.ASE.LP_01(end,vpiasePos(i+1))) )
%         maxAseDifdB = 10*log10( abs( 1e-3*10.^(Pap.LP_01(1+3*i,end)) - 1e-3*10.^(vpiedfamm.axialdist.ASE.LP_01(end,vpiasePos(i+1))) )/1e3) ;
%         maxdiffLambda = (DAT(1+3*i,1)*1e9);
%     end
% end
% annotation('textbox', [0.1254, 0.148, 0, 0], 'string', 'Matlab')
% annotation('textbox', [0.1254, 0.128, 0, 0], 'string', 'VPIphotonics')

        % LP11a
%%%% ASE+ Matlab
% for i=0:6
%     plot(z , Pap.LP_11_a(1+8*i,:) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm")) ; hold on
% end
% 
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% 
% %%%% ASE+ VPI
% for i=1:length(vpiasePos)
%     plot(zvpi , vpiedfamm.axialdist.ASE.LP_11a(:,vpiasePos(i)) , "DisplayName", strcat(int2str(vpiedfamm.ASElam(1,vpiasePos(i))*1e9) , " nm") , LineStyle="--" )
% end
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7,'FontSize', 9)
% xlabel('Posición Axial [m]','FontSize', 14) ; ylabel('Potencia [dBm]','FontSize', 14); title('Distribución Axial de la Potencia de Ruido ASE+ en modo LP11a','FontSize', 14)
% ylim([-90 -30])
% 
% maxAseDif = -100000;
% for i = 0:6
%     if maxAseDif < abs( Pap.LP_11_a(1+8*i,end) - vpiedfamm.axialdist.ASE.LP_11a(end,vpiasePos(i+1)) )
%         maxAseDif = abs( Pap.LP_11_a(1+8*i,end) - vpiedfamm.axialdist.ASE.LP_11a(end,vpiasePos(i+1)) ) ;
%         maxdiffLambda = (DAT(1+8*i,1)*1e9);
%     end
% end
% annotation('textbox', [0.1254, 0.148, 0, 0], 'string', 'Matlab')
% annotation('textbox', [0.1254, 0.128, 0, 0], 'string', 'VPIphotonics')


%               Propagación ASE + -


% %% LP01
% % ASE+ Matlab
% for i=0:6
%     plot(z , Pap.LP_01(1+8*i,:) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm")) ; hold on
% end
% 
% % ASE- Matlab
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% for i=0:6
%     plot(z , Pan.LP_01(1+8*i,:) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm") ) ; hold on
% end
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7,'FontSize',9)
% xlabel('Posición Axial [m]','FontSize',14) ; ylabel('Potencia [dBm]','FontSize',14); title('Distribución Axial de Ruido ASE Total en modo LP01','FontSize',14)
% ylim([-90 -30])


% %% LP11a

% for i=0:6
%     plot(z , Pap.LP_11_a(1+8*i,:) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm")) ; hold on
% end
% 
% % ASE- Matlab
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% for i=0:6
%     plot(z , Pan.LP_11_a(1+8*i,:) , "DisplayName", strcat(int2str(DAT(1+8*i,1)*1e9) , " nm") ) ; hold on
% end
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7,'FontSize',9)
% xlabel('Posición Axial [m]','FontSize',14) ; ylabel('Potencia [dBm]','FontSize',14); title('Distribución Axial de Ruido ASE Total en modo LP11a','FontSize',14)
% ylim([-90 -30])




% % % NF Y OSNR
 
% Pan = EDFA.Nucleo1.Pan;
% Pap = EDFA.Nucleo1.Pap;
%Pase = EDFA.Nucleo1.Pase;

% Osnr01Matlab = EDFA.Nucleo1.signal.Potencia_dBm.LP_01 - Pap.LP_01;
% NF01Matlab = Osnr01Matlab(:,2) - Osnr01Matlab(:,end);
% Osnr11aMatlab = EDFA.Nucleo1.signal.Potencia_dBm.LP_11_a - Pap.LP_11_a;
% NF11aMatlab = Osnr11aMatlab(:,2) - Osnr11aMatlab(:,end);
% 
% 
% for i=1:50
%     Osnr01VPI(:,i) = vpiedfamm.axialdist.Sig.LP_01(:,i) - vpiedfamm.axialdist.ASE.LP_01(:,i+3);
%     Osnr11aVPI(:,i) = vpiedfamm.axialdist.Sig.LP_11a(:,i) - vpiedfamm.axialdist.ASE.LP_11a(:,i+3);
% end
% NF01VPI = Osnr01VPI(3,:) - Osnr01VPI(end,:);
% NF11aVPI = Osnr11aVPI(3,:) - Osnr11aVPI(end,:);
% 
% plot(DAT(:,1).*1e9 , NF01Matlab , "DisplayName" , "Matlab LP01" ) , hold on
% plot(DAT(:,1).*1e9 , NF11aMatlab , "DisplayName" , "MatlabLP11a" )
% 
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% 
% plot(DAT(:,1).*1e9 , NF01VPI , "DisplayName" , "VPIphotonics LP01" , LineStyle="--"  )
% plot(DAT(:,1).*1e9 , NF11aVPI , "DisplayName" , "VPIphotonics LP11a" , LineStyle="--"  )
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2,'FontSize',9)
% xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('NF [dB]','FontSize',14); title('Figura de Ruido','FontSize',14)
% %ylim([8 15])
% 
% for i=1:length(NF01VPI)
%     DiffNF01(1,i) = (abs(NF01VPI(i) - NF01Matlab(i)));
%     DiffNF01(2,i) = DAT(i,1);
%     DiffNF11(1,i) = (abs(NF11aVPI(i) - NF11aMatlab(i)));
%     DiffNF11(2,i) = DAT(i,1);
% end
% 
% annotation('textbox', [0.28, 0.148, 0, 0], 'string', 'Matlab')
% annotation('textbox', [0.28, 0.128, 0, 0], 'string', 'VPIphotonics')
% ylim([10 15])


%% N1 y N2

% plot(z,EDFA.Nucleo1.N1 , "DisplayName" , "Iones en estado basal (N1)") ; hold on
% plot(z,EDFA.Nucleo1.N2 , "DisplayName" , "Iones en estado basal (N2)")
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 1,'FontSize',9)
% xlabel('Posición Axial [m]','FontSize',14) ; ylabel('Cantidad de Iones','FontSize',14); title('Distribución Axial de Densidades Poblacionenales de Iones','FontSize',14)

% %% Inversion Fraccional

% plot(z,(EDFA.Nucleo1.N2./EDFA.Nucleo1.Nt)*100 , "DisplayName" , "Matlab") ; hold on
% plot(VPI_InversionFraction(:,1),VPI_InversionFraction(:,2) , "DisplayName" , "VPIphotonics")
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 1,'FontSize',9)
% xlabel('Posición Axial [m]','FontSize',14) ; ylabel('Iones Excitados (N2) [%]','FontSize',14); title('Distribución Axial de Iones Excitados','FontSize',14)


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

