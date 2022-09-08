% close all; 
clear all; clc ; close all

% SPAN para fibra con Pin -35 dBm ; 10m ; 150mw Pump


% Parámetros de entrada
c = 299.792458e6; % [m/s]
h = 6.62607015*10^(-34);

    % Modos y canales de señal y bombeo



Signal.NumberOfChannels=31;
Wavelength_gridS=linspace(1530,1560,Signal.NumberOfChannels).*1e-9; % Banda C



Pin = -20 ; %-35; %[dBm]
Signal.modos = ["01" "11_a" ];%"21_a"] ;

Signal.lambda.LP_01     = Wavelength_gridS;                                  P0_signal.LP_01     = Pin*ones(1,length(Signal.lambda.LP_01));
Signal.lambda.LP_11_a   = Wavelength_gridS;                                  P0_signal.LP_11_a   = Pin*ones(1,length(Signal.lambda.LP_11_a));
Signal.lambda.LP_21_a   = Wavelength_gridS;                                  P0_signal.LP_21_a   = Pin*ones(1,length(Signal.lambda.LP_21_a));
Signal.lambda.LP_02     = Wavelength_gridS;                                  P0_signal.LP_02     = Pin*ones(1,length(Signal.lambda.LP_02));


%Pump.modos = "12_a" ;
Pump.modos = "01" ;
Wavelength_gridP = 980e-9;
Ppump= 100e-3; %150e-3; %[W]

Pump.lambda.LP_01   = Wavelength_gridP;                         P0_pump.LP_01   = Ppump  ;  
%Pump.lambda.LP_12_a   = Wavelength_gridP;                         P0_pump.LP_12_a   = Ppump  ;  

ModoS=strcat("LP_",Signal.modos(:));
ModoP=strcat("LP_",Pump.modos(:));


    % POTENCIAS

for i=1:length(Signal.modos)        % Potencia de señal a W
    for j=1:length(P0_signal.(ModoS(i)))
        P0_signal.(ModoS(i))(j) = 1e-3*10^(P0_signal.(ModoS(i))(j)/10);
    end
end ;clear i j;
Signal.P0 = P0_signal; 
Pump.P0 = P0_pump;

P.Np=2; 
P.Fc=c/Wavelength_gridS(ceil(length(Wavelength_gridS)/2)); P.Fb = 50e9; 
ASE= -200;

    % Datos de la fibra
Fibra.nucleos = 1;                                           % Numero de nucleos
Fibra.largo = 5; %10; 
Fibra.radio = 5e-6 ; Fibra.N = 7e24; 
Fibra.n1 = 1.45 ;   Fibra.IndexContrast=0.01;
Fibra.AN=Fibra.n1*sqrt(2*Fibra.IndexContrast);
Fibra.n2 =sqrt((Fibra.n1^2-Fibra.AN^2));
%Fibra.GEF = load('GEF_Filters/EDFA_10m_150mw_OptiSystem') ;


Fibra.dvk=P.Fb;

Fibra.WaitBar = 1; Fibra.Avance = 1;    % Despliegue de info
Fibra.ASEFlag = 1;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)


Nspans = 3;
LargoFibra = 90; %175;        % [km]



%% Graficos

clc
load("Datos_Span5m_100mw.mat")
Nspans = length(fieldnames(Span))-1;  
ModoS = Signal.modos ; 
Pump.modos = "01" ; ModoP = "01" ; 


GrafSpans.z_total = [Span.EDFA1.Nucleo1.z]; 
GrafSpans.Signal_total = Span.EDFA1.Nucleo1.signal.Potencia_dBm;
GrafSpans.ASE_total = Span.EDFA1.Nucleo1.Pap;

z = Span.EDFA1.Nucleo1.z;

for span = 1:Nspans
    GrafSpans.z_total = [GrafSpans.z_total (z+(0.1*LargoFibra+GrafSpans.z_total(end)))];
    for i=1:length(Signal.modos)
        GrafSpans.Signal_total.(strcat("LP_",ModoS(i))) = [GrafSpans.Signal_total.(strcat("LP_",ModoS(i))) Span.(strcat("EDFA",num2str(span+1))).Nucleo1.signal.Potencia_dBm.(strcat("LP_",ModoS(i)))];
        GrafSpans.ASE_total.(strcat("LP_",ModoS(i))) = [GrafSpans.ASE_total.(strcat("LP_",ModoS(i))) Span.(strcat("EDFA",num2str(span+1))).Nucleo1.Pap.(strcat("LP_",ModoS(i)))];
    end ; clear i;
end

    % % Graficos

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% % SAVE:
% % print -dpdf 'NAME'
% 
close all
 
graf.Nc = length(fieldnames(Span.EDFA1)); 
graf.z = Span.EDFA1.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';

Name_Modos = ["LP01" "LP11a" "LP21a"];

% %%% Power Out
% s=1;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.SignalPump_OUT.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.salida.signal.potencia_dBm.(strcat("LP_",ModoS(s)))(:,end);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.SignalPump_OUT.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Potencia de Salida en Modo LP',Name_Modos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Potencia [dBm]','FontSize',14)
%     ylim([-4 1]);
% 
%     % % % Save:
%     if span==Nspans+1
%     set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
%     if s==1 ; print -dpdf 'Cascada_PotenciaSalida01'
%     elseif s==2; print -dpdf 'Cascada_PotenciaSalida11a';
%     elseif s==3 ; print -dpdf 'Cascada_PotenciaSalida21a'; 
%     end;end
% end ; clear s ejex leyenda;




 % Ganancias
 % LP01
% s=1;
% figure(s)
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.ganancias.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.salida.ganancias.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.ganancias.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Ganancias Por Amplificador en Modo LP',Signal.modos(s)),'FontSize',14) ; legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Ganancia [dB]','FontSize',14)
%     ylim([15 20])%ylim([30.5 35.5])
% end ; clear s ejex leyenda;
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'Cascada_Ganancias_LP01'

 % LP11a
% s=2;
% figure
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.ganancias.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.salida.ganancias.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.ganancias.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title('Ganancias Por Amplificador en Modo LP 11a','FontSize',14) ; legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Ganancia [dB]','FontSize',14)
%     ylim([13.5 18.5])%ylim([29 34])
% end ; clear s ejex leyenda;
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'Cascada_Ganancias_LP11a'

 % LP21a
% s=3;
% figure
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.ganancias.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.salida.ganancias.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.ganancias.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Ganancias Por Amplificador en Modo LP 21a'),'FontSize',14) ; legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Ganancia [dB]','FontSize',14)
%     ylim([24.5 29.5])
% end ; clear s ejex leyenda;
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'Cascada_Ganancias_LP21a'


%  % NF
% LP01
% s=1;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.NF.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.NF.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.NF.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title('Figuras de Ruido Por Amplificador en Modo LP 01','FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
%     ylim([0 10])
% end ; clear s ejex leyenda;
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'Cascada_NF01'

% % LP11a
% s=2;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.NF.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.NF.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.NF.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title('Figuras de Ruido Por Amplificador en Modo LP11a','FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
%     ylim([0 10])
% end ; clear s ejex leyenda;
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'Cascada_NF11a'

% % LP21a
% s=3;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.NF.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.NF.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.NF.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title('Figuras de Ruido Por Amplificador en Modo LP 21a','FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
%     ylim([0 10])
% end ; clear s ejex leyenda;
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'Cascada_NF21a'


% % OSRN OUT

%LP01
% s=1;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.all_ASE.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.Pase.(strcat("LP_",ModoS(s)));
%     graf.all_Signal.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.signal.Potencia_dBm.(strcat("LP_",ModoS(s)));
%     graf.OSNR = graf.all_Signal.(strcat("LP_",ModoS(s))) - graf.all_ASE.(strcat("LP_",ModoS(s)));
%     graf.OSNR_OUT.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.OSNR.(strcat("LP_",ModoS(s)))(:,end);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.OSNR_OUT.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title('OSNR en Modo LP01','FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
%     %ylim([0 10])
% end ; clear s ejex leyenda;
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'Cascada_OSNR01'

% %LP11a
% s=2;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.OSNR_OUT.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.OSNR.(strcat("LP_",ModoS(s)))(:,end);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.OSNR_OUT.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title('OSNR en Modo LP11a','FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'Cascada_OSNR11a'

%LP21a
% s=3;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.OSNR_OUT.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.OSNR.(strcat("LP_",ModoS(s)))(:,end);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.OSNR_OUT.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title('OSNR en Modo LP21a','FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% print -dpdf 'Cascada_OSNR21a'


% % % RUIDO ASE

% s=1;
% figure(s)
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.ASE.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.Pase.(strcat("LP_",ModoS(s)))(26,:);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Span.(strcat("EDFA",num2str(span))).Nucleo1.z;
%     plot(ejex,graf.ASE.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Distribución axial del ruido ASE a 1555 nm en el Modo LP',Name_Modos(s)),'FontSize',14) ; legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Posición en la fibra [m]','FontSize',14) ; ylabel('Potencia [dBm]','FontSize',14)
%     ylim([-60 -20])%ylim([30.5 35.5])
% 
%     % % % Save:
%     if span==Nspans+1
%     set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
%     if s==1 ; print -dpdf 'Cascada_PotenciaASE01'
%     elseif s==2; print -dpdf 'Cascada_PotenciaASE11a';
%     elseif s==3 ; print -dpdf 'Cascada_PotenciaASE21a'; 
%     end;end
%         
% end ; clear s ejex leyenda;


% % % % Señal

% s=2;
% figure(1)
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.Signal.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.signal.Potencia_dBm.(strcat("LP_",ModoS(s)))(26,:);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Span.(strcat("EDFA",num2str(span))).Nucleo1.z;
%     plot(ejex,graf.Signal.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Distribución axial del canal a 1555 nm en el Modo LP',Name_Modos(s)),'FontSize',14) ; legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Posición en la fibra [m]','FontSize',14) ; ylabel('Potencia [dBm]','FontSize',14)
%     %ylim([-60 -20])%ylim([30.5 35.5])
% 
%     % % % Save:
%     if span==Nspans+1
%     set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
%     if s==1 ; print -dpdf 'Cascada_PotenciaSeñalenEDFA01'
%     elseif s==2; print -dpdf 'Cascada_PotenciaSeñalenEDFA11a';
%     elseif s==3 ; print -dpdf 'Cascada_PotenciaSeñalenEDFA21a'; 
%     end;end
% end ; clear s ejex leyenda;

