clear all; close all; clc

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% SAVE:
% print -dpdf 'NAME'

% DAT = [lambdas Gain01 lambdas Gain11a frecuencias NF01 NF11a];
load("MatlabSigGains.mat")
load("ModeGainsandNFs.mat")
load("AmplifiedSignals.mat");
load('MatlabEDFA_Pump01.mat');
load('MatlabS01.mat') ; load('MatlabS11_a.mat');

%% GANANCIAS

% figure(1)

% plot(DAT(:,1).*1e9,DAT(:,2), "DisplayName","LP01 - Simulador Matlab" , Color= '#0072BD') ; hold on
% plot(AmplifiedSignals(:,3),fliplr(ModeGainsandNFs(:,2)),"DisplayName","LP01 - VPIPhotonics" , LineStyle="--", Color= '#0072BD' ) ; 
% plot(DAT(:,1).*1e9,DAT(:,4),  "DisplayName","LP11a - Simulador Matlab", "Color", '#D95319' )
% plot(AmplifiedSignals(:,1),fliplr(ModeGainsandNFs(:,6)),"DisplayName","LP11a - VPIPhotonics", LineStyle="--" , Color='#D95319') 
% xlabel('Longitud de Onda [nm]') ; ylabel('Magnitud [dB]') ; title("Ganancia EDFA con Bombeo en LP01") ; 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off')
% 
% DG01 =  abs(DAT(:,2) - fliplr(ModeGainsandNFs(:,2))) ;
% DG11 = abs(DAT(:,4) - fliplr(ModeGainsandNFs(:,6))) ;
%     % Diferencia de Ganancias
% yyaxis right
% ylabel('Diferencia de Ganancias')
% plot(DAT(:,1).*1e9 , DG01 ,  "-o" , DisplayName = "DG - LP01" , Color='#0072BD')
% plot(DAT(:,1).*1e9 , DG11 ,  "-o" , DisplayName = "DG - LP11a" , Color='#D95319')
% ylim([0, 2]) ; grid minor
% 



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

load('VPI_EDFAmm_axialdist.mat') % Bombeo LP01


%               Propagación SEÑAL LP01
% figure(3)

% z = linspace(0,3,length(S01(1,:)));
% zvpi = linspace(0,3,length(vpiedfamm.axialdist.Sig.LP_01(:,1)));
% 
% 
% for i=0:6
%     plot(z , EDFA.Nucleo1.signal.Potencia_dBm.LP_01(1+3*i,:) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm")) ; hold on
%     % plot(z , S01(2,:)) 
%     % plot(z , S01(3,:)) 
% end
% set(gca,'ColorOrderIndex',1)
% for i=0:6
%     plot(zvpi , vpiedfamm.axialdist.Sig.LP_01(:,1+3*i) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm") , LineStyle="--" )
%     % plot(zvpi , vpiedfamm.axialdist.Sig.LP_01(:,2) , LineStyle="--")
%     % plot(zvpi , vpiedfamm.axialdist.Sig.LP_01(:,3) , LineStyle="--")
% end
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7)
% xlabel('Posición Axial [m]') ; ylabel('Potencia [dBm]'); title('Distribución Axial de la Potencia de Señal LP01')


%               Propagacion SEÑAL LP11

% z = linspace(0,3,length(S01(1,:)));
% zvpi = linspace(0,3,length(vpiedfamm.axialdist.Sig.LP_01(:,1)));
% 
% 
% for i=0:6
%     plot(z , EDFA.Nucleo1.signal.Potencia_dBm.LP_11_a(1+3*i,:) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm")) ; hold on
%     % plot(z , S01(2,:)) 
%     % plot(z , S01(3,:)) 
% end
% set(gca,'ColorOrderIndex',1)
% for i=0:6
%     plot(zvpi , vpiedfamm.axialdist.Sig.LP_11a(:,1+3*i) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm") , LineStyle="--" )
%     % plot(zvpi , vpiedfamm.axialdist.Sig.LP_01(:,2) , LineStyle="--")
%     % plot(zvpi , vpiedfamm.axialdist.Sig.LP_01(:,3) , LineStyle="--")
% end
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7)
% xlabel('Posición Axial [m]') ; ylabel('Potencia [dBm]'); title('Distribución Axial de la Potencia de Señal LP11a')


%               Propagación Bombeo

% z = linspace(0,3,length(S01(1,:)));
% zvpi = linspace(0,3,length(vpiedfamm.axialdist.Sig.LP_01(:,1)));
% 
% 
% plot(z , EDFA.Nucleo1.pump.Potencia_dBm.LP_01(1,:) , "DisplayName", "Matlab") ; hold on
% 
% plot(zvpi , vpiedfamm.axialdist.Pump(:,1) , "DisplayName", "VPIphotonics" )% ,LineStyle="--" )
% 
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7)
% xlabel('Posición Axial [m]') ; ylabel('Potencia [dBm]'); title('Distribución Axial de la Potencia de Bombeo LP01')


%% %               Propagación ASE , OSNR , NF
vpiasePos = [7 , 10 , 13, 17, 20, 22, 25];

Pan = EDFA.Nucleo1.Pan;
Pap = EDFA.Nucleo1.Pap;
z = linspace(0,3,length(S01(1,:)));
zvpi = linspace(0,3,length(vpiedfamm.axialdist.Sig.LP_01(:,1)));

% % % % ASE+ Matlab

% for i=0:6
%     plot(z , Pap.LP_01(1+3*i,:) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm")) ; hold on
% end
% % ASE+ VPI
% 
% set(gca,'ColorOrderIndex',1)
% for i=1:length(vpiasePos)
%     plot(zvpi , vpiedfamm.axialdist.ASE.LP_01(:,vpiasePos(i)) , "DisplayName", strcat(int2str(vpiedfamm.ASElam(1,vpiasePos(i))*1e9) , " nm") , LineStyle="--" )
% end
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7)
% xlabel('Posición Axial [m]') ; ylabel('Potencia [dBm]'); title('Distribución Axial de la Potencia de Ruido ASE+ en modo LP01')
% ylim([-100 -10])
% 
% maxAseDif = -100000;
% for i = 0:6
%     if maxAseDif < abs( Pap.LP_01(1+3*i,end) - vpiedfamm.axialdist.ASE.LP_01(end,vpiasePos(i+1)) )
%         maxAseDif = abs( Pap.LP_01(1+3*i,end) - vpiedfamm.axialdist.ASE.LP_01(end,vpiasePos(i+1)) ) ;
%         maxdiffLambda = (DAT(1+3*i,1)*1e9);
%     end
% end
%dB
% for i = 0:6
%     if maxAseDif < abs( 1e-3*10.^(Pap.LP_01(1+3*i,end)) - 1e-3*10.^(vpiedfamm.axialdist.ASE.LP_01(end,vpiasePos(i+1))) )
%         maxAseDifdB = 10*log10( abs( 1e-3*10.^(Pap.LP_01(1+3*i,end)) - 1e-3*10.^(vpiedfamm.axialdist.ASE.LP_01(end,vpiasePos(i+1))) )/1e3) ;
%         maxdiffLambda = (DAT(1+3*i,1)*1e9);
%     end
% end

    % LP11a

% for i=0:6
%     plot(z , Pap.LP_11_a(1+3*i,:) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm")) ; hold on
% end
% % ASE+ VPI
% 
% set(gca,'ColorOrderIndex',1)
% for i=1:length(vpiasePos)
%     plot(zvpi , vpiedfamm.axialdist.ASE.LP_11a(:,vpiasePos(i)) , "DisplayName", strcat(int2str(vpiedfamm.ASElam(1,vpiasePos(i))*1e9) , " nm") , LineStyle="--" )
% end
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7)
% xlabel('Posición Axial [m]') ; ylabel('Potencia [dBm]'); title('Distribución Axial de la Potencia de Ruido ASE+ en modo LP11a')
% ylim([-100 -10])
% 
% maxAseDif = -100000;
% for i = 0:6
%     if maxAseDif < abs( Pap.LP_11_a(1+3*i,end) - vpiedfamm.axialdist.ASE.LP_11a(end,vpiasePos(i+1)) )
%         maxAseDif = abs( Pap.LP_11_a(1+3*i,end) - vpiedfamm.axialdist.ASE.LP_11a(end,vpiasePos(i+1)) ) ;
%         maxdiffLambda = (DAT(1+3*i,1)*1e9);
%     end
% end

%               Propagación ASE + -

Pan = EDFA.Nucleo1.Pan;
Pap = EDFA.Nucleo1.Pap;
Pase = EDFA.Nucleo1.Pase;

z = linspace(0,3,length(S01(1,:)));
zvpi = linspace(0,3,length(vpiedfamm.axialdist.Sig.LP_01(:,1)));

% %%    LP01
% ASE+ Matlab
% for i=0:6
%     plot(z , Pap.LP_01(1+3*i,:) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm")) ; hold on
% end
% 
% % ASE- Matlab
% set(gca,'ColorOrderIndex',1,'FontSize',8)
% for i=0:6
%     plot(z , Pan.LP_01(1+3*i,:) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm") ) ; hold on
% end
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 7,'FontSize',9)
% xlabel('Posición Axial [m]','FontSize',14) ; ylabel('Potencia [dBm]','FontSize',14); title('Distribución Axial de Ruido ASE Total en modo LP01','FontSize',14)
% ylim([-100 -10])

% %%    LP11a

for i=0:6
    plot(z , Pap.LP_11_a(1+3*i,:) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm")) ; hold on
end

% ASE- Matlab
set(gca,'ColorOrderIndex',1,'FontSize',8)
for i=0:6
    plot(z , Pan.LP_11_a(1+3*i,:) , "DisplayName", strcat(int2str(DAT(1+3*i,1)*1e9) , " nm") ) ; hold on
end

legend('Location', 'southoutside','Orientation','horizontal','Box','on', "NumColumns" , 7,'FontSize',9)
xlabel('Posición Axial [m]','FontSize',14) ; ylabel('Potencia [dBm]','FontSize',14); title('Distribución Axial de Ruido ASE Total en modo LP11a','FontSize',14)
ylim([-100 -10])




% % % NF Y OSNR
% 
% Pan = EDFA.Nucleo1.Pan;
% Pap = EDFA.Nucleo1.Pap;
% Pase = EDFA.Nucleo1.Pase;
% 
% Osnr01Matlab = EDFA.Nucleo1.signal.Potencia_dBm.LP_01 - Pap.LP_01;
% NF01Matlab = Osnr01Matlab(:,2) - Osnr01Matlab(:,end);
% Osnr11aMatlab = EDFA.Nucleo1.signal.Potencia_dBm.LP_11_a - Pap.LP_11_a;
% NF11aMatlab = Osnr11aMatlab(:,2) - Osnr11aMatlab(:,end);
% 
% 
% for i=7:26
%     Osnr01VPI(:,i-6) = vpiedfamm.axialdist.Sig.LP_01(:,i-6) - vpiedfamm.axialdist.ASE.LP_01(:,i);
%     Osnr11aVPI(:,i-6) = vpiedfamm.axialdist.Sig.LP_11a(:,i-6) - vpiedfamm.axialdist.ASE.LP_11a(:,i);
% end
% NF01VPI = Osnr01VPI(2,:) - Osnr01VPI(end,:);
% NF11aVPI = Osnr11aVPI(2,:) - Osnr11aVPI(end,:);
% 
% plot(DAT(:,1).*1e9 , NF01Matlab , "DisplayName" , "Matlab LP01" ) , hold on
% plot(DAT(:,1).*1e9 , NF11aMatlab , "DisplayName" , "MatlabLP11a" )
% 
% set(gca,'ColorOrderIndex',1)
% 
% plot(vpiedfamm.SigLam.*1e9 , NF01VPI , "DisplayName" , "VPIphotonics LP01" , LineStyle="--"  )
% plot(vpiedfamm.SigLam.*1e9 , NF11aVPI , "DisplayName" , "VPIphotonics LP11a" , LineStyle="--"  )
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 2)
% xlabel('Longitud de Onda [nm]') ; ylabel('NF [dB]'); title('Figura de Ruido')
% ylim([8 15])
% 
% max(abs(NF01VPI' - NF01Matlab))
% max(abs(NF11aVPI' - NF11aMatlab))




%% N1 y N2

% plot(z,EDFA.Nucleo1.N1 , "DisplayName" , "Iones en estado basal (N1)") ; hold on
% plot(z,EDFA.Nucleo1.N2 , "DisplayName" , "Iones en estado basal (N2)")
% 
% legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 1)
% xlabel('Posición Axial [m]') ; ylabel('Cantidad de Iones'); title('Distribución Axial de Densidades Poblacionenales de Iones')

%% Cargar resultado vpi.. desde columna ec
% VPISignal = resultadovpipropagacionaxial;
% 
% for i=0:19
%     VPISignal_01(:,i+1) = VPISignal(:,2+4*i)/10e13;
%     VPISignal_11a(:,i+1) = VPISignal(:,4+4*i)/10e13;
% end
% VPISignal_01([1,2],:) = VPISignal_01([1,2],:).*10e13;
% VPISignal_11a([1,2],:) = VPISignal_11a([1,2],:).*10e13;
