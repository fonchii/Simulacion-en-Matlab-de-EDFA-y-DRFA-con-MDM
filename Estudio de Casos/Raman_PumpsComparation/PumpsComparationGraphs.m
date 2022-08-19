close all ; clc ; clear all


%% Datos Entrada
c = 299792458;
% FIBRA
In.Fibra.RamanMethod              = 'Forward&Backward';                   % 'Forward', 'Backward', 'Forward&Backward'
In.Fibra.AttenuationMethod        = 'Dynamic';                    % 'Dynamic' , 'Static'
In.Fibra.Length                   = 50;                          % fibre length (km)
In.Fibra.T                        = 25;                           % Temperatura Fibra (ambiente)
In.Fibra.PolarizationFactor       = 0.5;                          % C_R_max
In.Fibra.n1=1.46;  In.Fibra.n2=1.42; In.Fibra.radio=25e-6; In.Fibra.area=pi*(In.Fibra.radio)^2;


% BOMBEOS : 
%     % LP01

In.Pump.LP01.Wavelengths       = 1450;               %c/(200e12) ;                    % [nm]
In.Pump.LP01.Powers            = 200*1e-3;                                                % [mW]

% SEÑALES : 
Nch = 100;
    % LP01  
In.Signal.LP11a.Wavelengths          = linspace(1500,1600,Nch) ;
In.Signal.LP11a.Powers               = -30*ones( 1,length(In.Signal.LP11a.Wavelengths) );                 %[dBm]
In.ASE.LP11a                         = -200*ones( 1,length(In.Signal.LP11a.Wavelengths) );
    % LP21a
In.Signal.LP11a.Wavelengths          = linspace(1500,1600,Nch) ;
In.Signal.LP11a.Powers               = -30*ones( 1,length(In.Signal.LP11a.Wavelengths) );                 %[dBm]
In.ASE.LP11a                         = -200*ones( 1,length(In.Signal.LP11a.Wavelengths) );


%% Graficar
clc;
load('RamanBackward.mat')
load('RamanForward.mat')
load('RamanHibrido2.mat')

close all; 
z = RamanBackward.z;
for mp = 1:length(RamanBackward.ModoP)
    pump.(RamanBackward.ModoP{mp}).forward = RamanBackward.Pump.forward.(RamanBackward.ModoP{mp});
    pump.(RamanBackward.ModoP{mp}).backward = RamanBackward.Pump.backward.(RamanBackward.ModoP{mp});
end
signal = RamanBackward.Sig.Power;
gain = RamanBackward.Sig.GainOnOFF;
osnr = RamanBackward.OSNR;

fs = 1; fp =1;
ModoS = RamanForward.ModoS; ModoP = RamanForward.ModoP;

% SAVE:
set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
%print -dpdf 'NAME'


% for ms = 1:length(Raman.ModoS)% Potencia Señal
%     for lS=1:10:length(In.Signal.(Raman.ModoS{ms}).Wavelengths)
%         figure(1); 
%         strlambda = strcat( num2str( In.Signal.(Raman.ModoS{ms}).Wavelengths(lS)) , "nm");
%         subplot(1,length(Raman.ModoS),fs) ; 
%         plot(z,10.*log10(signal.(Raman.ModoS{ms})(lS,:)./1e-3) , 'DisplayName',strlambda) ,ylabel("Potencia [dBm]") ; hold on
%         %subplot(1,length(Raman.ModoS),fs) ; plot(z,10.*log10(signal.Off.(Raman.ModoS{ms})./1e-3) ) ,ylabel("Potencia [dBm]") ; hold on
%         %subplot(1,fs,fs) ; plot( z,(signal.(Raman.ModoS{ms}) ) ,ylabel("Potencia [mW]"); hold on
%         xlabel("Largo [km]") ; title(strcat("Signal Modo ", Raman.ModoS{ms})) ; legend('Location','southoutside','NumColumns',10) 
%     end
%     fs = fs+1;
% end
% set(gca,'ColorOrderIndex',1)
% fs = 1;
% for ms = 1:length(Raman.ModoS)% Potencia Señal Off
%     for lS=1:10:length(In.Signal.(Raman.ModoS{ms}).Wavelengths)
%         %figure(1); 
%         strlambda = strcat( num2str( In.Signal.(Raman.ModoS{ms}).Wavelengths(lS)) , "nm OFF");
%         %h=subplot(1,length(Raman.ModoS),fs) ; h.ColorOrderIndex=1;
%         plot(z,10.*log10(signal.Off.(Raman.ModoS{ms})(lS,:)./1e-3) , '--', 'DisplayName',strlambda) ,ylabel("Potencia [dBm]") ; hold on
%         %subplot(1,fs,fs) ; plot( z,(signal.Off.(Raman.ModoS{ms})) ) ,ylabel("Potencia [mW]"); hold on
%         xlabel("Largo [km]") ; title(strcat("Signal Modo ", Raman.ModoS{ms})) ; legend('Location','southoutside','NumColumns',10)
%     end
%     fs = fs+1;
% end
% 
% switch In.Fibra.RamanMethod
%     case 'Forward&Backward'
%         for mp = 1:length(Raman.ModoP)% Potencia Bombeo
%             figure(2); 
%             for i = 1:length(In.Pump.(Raman.ModoP{mp}).Wavelengths)
%                 DispName1 = strcat("Forward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
%                 DispName2 = strcat("Backward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
%                 %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).forward(i,:)./1e-3),"DisplayName","Forward") ,ylabel("Potencia [dBm]") , hold on
%                 %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).backward(i,:)./1e-3),"DisplayName","Backward") ,ylabel("Potencia [dBm]")
%                 subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).forward(i,:),"DisplayName",DispName1) ,ylabel("Potencia [W]") , hold on
%                 subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).backward(i,:),"DisplayName",DispName2) ,ylabel("Potencia [W]")
%                 xlabel("Largo [km]") ; title(strcat("Pump Modo ", Raman.ModoP{mp})) , legend("Location","best")
%             end
%             fp = fp+1;
%         end
%     case 'Backward'
%         for mp = 1:length(Raman.ModoP)% Potencia Bombeo
%             figure(2); 
%             for i = 1:length(In.Pump.(Raman.ModoP{mp}).Wavelengths)
%                 DispName2 = strcat("Backward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
%                 %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).backward(i,:)./1e-3),"DisplayName",DispName2) ,ylabel("Potencia [dBm]"), hold on
%                 subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).backward(i,:),"DisplayName",DispName2) ,ylabel("Potencia [W]"), hold on
%                 xlabel("Largo [km]") ; title(strcat("Pump Modo ", Raman.ModoP{mp})) , legend("Location","best")
%             end
%             fp = fp+1;
%         end
%     case 'Forward'
%         for mp = 1:length(Raman.ModoP)% Potencia Bombeo
%             figure(2); 
%             for i = 1:length(In.Pump.(Raman.ModoP{mp}).Wavelengths)
%                 DispName1 = strcat("Forward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
%                 subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).forward(i,:),"DisplayName",DispName1) ,ylabel("Potencia [W]") , hold on
%                 xlabel("Largo [km]") ; title(strcat("Pump Modo ", Raman.ModoP{mp})) , legend("Location","best")
%             end
%             fp = fp+1;
%         end
% end

% % Ganancias


% % LP11a

% ms=1;
% plot(In.Signal.(RamanForward.ModoS{ms}).Wavelengths , RamanForward.Sig.GainOnOFF.(ModoS{ms}), '-o','DisplayName','Bombeo Forward')  ; hold on
% plot(In.Signal.(RamanBackward.ModoS{ms}).Wavelengths , RamanBackward.Sig.GainOnOFF.(ModoS{ms}), '-o','DisplayName','Bombeo Backward')
% plot(In.Signal.(RamanHibrido.ModoS{ms}).Wavelengths , RamanHibrido.Sig.GainOnOFF.(ModoS{ms}), '-o','DisplayName','Bombeo Hibrido')
% set(gca,"ColorOrderIndex",1,'FontSize',8)
% xlabel("Longitud de Onda [nm]",'FontSize',14) ,ylabel("Ganancia [dB]",'FontSize',14); title("Ganancias On-Off modo LP11a",'FontSize',14) ; 
% legend("Location","southoutside",'FontSize', 9,'Orientation','horizontal','Box','off')

%set(gca,'ColorOrderIndex',1)


% % LP21a

% ms = 2;
% plot(In.Signal.(RamanForward.ModoS{1}).Wavelengths , RamanForward.Sig.GainOnOFF.(ModoS{ms}), '-o','DisplayName','Bombeo Forward')  ; hold on
% plot(In.Signal.(RamanBackward.ModoS{1}).Wavelengths , RamanBackward.Sig.GainOnOFF.(ModoS{ms}), '-o','DisplayName','Bombeo Backward')
% plot(In.Signal.(RamanHibrido.ModoS{1}).Wavelengths , RamanHibrido.Sig.GainOnOFF.(ModoS{ms}), '-o','DisplayName','Bombeo Hibrido')
% set(gca,"ColorOrderIndex",1,'FontSize',8)
% xlabel("Longitud de Onda [nm]",'FontSize',14) ,ylabel("Ganancia [dB]",'FontSize',14); title("Ganancias On-Off modo LP21a",'FontSize',14) ; 
% legend("Location","southoutside",'FontSize', 9,'Orientation','horizontal','Box','off')


% % OSNR

% % LP11a

% ms = 1;
% plot(In.Signal.(ModoS{ms}).Wavelengths , RamanForward.OSNR.(ModoS{ms}), '-o','DisplayName','Bombeo Forward')  ; hold on
% plot(In.Signal.(ModoS{ms}).Wavelengths , RamanBackward.OSNR.(ModoS{ms}), '-o','DisplayName','Bombeo Backward')  ; 
% plot(In.Signal.(ModoS{ms}).Wavelengths , RamanHibrido.OSNR.(ModoS{ms}), '-o','DisplayName','Bombeo Hibrido')  ; 
% set(gca,'FontSize',8)
% xlabel("Longitud de Onda [nm]",'FontSize',14) ,ylabel("Magnitud [dB]",'FontSize',14) ; title("OSNR en modo LP11a",'FontSize',14) ; 
% legend("Location","southoutside",'FontSize', 9,'Orientation','horizontal','Box','off')

% % LP21a

% ms = 2;
% plot(In.Signal.(ModoS{1}).Wavelengths , RamanForward.OSNR.(ModoS{ms}), '-o','DisplayName','Bombeo Forward')  ; hold on
% plot(In.Signal.(ModoS{1}).Wavelengths , RamanBackward.OSNR.(ModoS{ms}), '-o','DisplayName','Bombeo Backward')  ; 
% plot(In.Signal.(ModoS{1}).Wavelengths , RamanHibrido.OSNR.(ModoS{ms}), '-o','DisplayName','Bombeo Hibrido')  ; 
% set(gca,'FontSize',8)
% xlabel("Longitud de Onda [nm]",'FontSize',14) ,ylabel("Magnitud [dB]",'FontSize',14) ; title("OSNR en modo LP21a",'FontSize',14) ; 
% legend("Location","southoutside",'FontSize', 9,'Orientation','horizontal','Box','off')


% % NF

% % LP11a

% ms = 1;
% plot(In.Signal.(ModoS{ms}).Wavelengths , RamanForward.NF.(ModoS{ms}), '-o','DisplayName','Bombeo Forward')  ; hold on
% plot(In.Signal.(ModoS{ms}).Wavelengths , RamanBackward.NF.(ModoS{ms}), '-o','DisplayName','Bombeo Backward')  ; 
% plot(In.Signal.(ModoS{ms}).Wavelengths , RamanHibrido.NF.(ModoS{ms}), '-o','DisplayName','Bombeo Hibrido')  ; 
% set(gca,'FontSize',8)
% xlabel("Longitud de Onda [nm]",'FontSize',14) ,ylabel("Magnitud [dB]",'FontSize',14) ; title("Figura de Ruido en modo LP11a",'FontSize',14) ; 
% legend("Location","southoutside",'FontSize', 9,'Orientation','horizontal','Box','off')

% % LP21a

ms = 2;
plot(In.Signal.(ModoS{1}).Wavelengths , RamanForward.NF.(ModoS{ms}), '-o','DisplayName','Bombeo Forward')  ; hold on
plot(In.Signal.(ModoS{1}).Wavelengths , RamanBackward.NF.(ModoS{ms}), '-o','DisplayName','Bombeo Backward')  ; 
plot(In.Signal.(ModoS{1}).Wavelengths , RamanHibrido.NF.(ModoS{ms}), '-o','DisplayName','Bombeo Hibrido')  ; 
set(gca,'FontSize',8)
xlabel("Longitud de Onda [nm]",'FontSize',14) ,ylabel("Magnitud [dB]",'FontSize',14) ; title("Figura de Ruido en modo LP21a",'FontSize',14) ; 
legend("Location","southoutside",'FontSize', 9,'Orientation','horizontal','Box','off')






% % ASE
% fs = 1;
% for ms = 1:10:length(Raman.ModoS)
%     for lS=1:10:length(In.Signal.(Raman.ModoS{ms}).Wavelengths)
%         figure(5); 
%         strlambda = strcat( num2str( In.Signal.(Raman.ModoS{ms}).Wavelengths(lS)) , "nm");
%         subplot(1,length(Raman.ModoS),fs) ; plot(z,10.*log10(Raman.ASE.(Raman.ModoS{ms})(lS,:)./1e-3) , 'DisplayName',strlambda) ,ylabel("Potencia [dBm]") ; hold on
%         xlabel("Largo [km]") ; title(strcat("Signal Modo ", Raman.ModoS{ms})) ; legend('Location','southoutside','NumColumns',10)
%     end
%     fs = fs+1;
% end

clear fs fp ms mp i strlambda DispName1 DispName2;

