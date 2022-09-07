% close all; 
clear all; clc ; close all

% Parámetros de entrada

    % Modos y canales de señal y bombeo
Signal.NumberOfChannels=41;
Signal.modos = ["11_a" "21_a" "02"];

% Ejemplo Pump01 - channels 20
%Frequency_gridS=linspace(191.07234e12,196.7723e12,Signal.NumberOfChannels);
% Ejemplo Pump12 - channels 50

%Frequency_gridS=linspace(191.19421875e12,193.64421875e12,Signal.NumberOfChannels);
c=299.792458e6; % [m/s]
%Wavelength_gridS=c./Frequency_gridS;

Wavelength_gridS=linspace(1530,1570,Signal.NumberOfChannels).*1e-9; % Banda C

Pin=0; %[dBm]

Signal.lambda.LP_11_a   = Wavelength_gridS;                                  P0_signal.LP_11_a   = Pin*ones(1,length(Signal.lambda.LP_11_a));
Signal.lambda.LP_21_a   = Wavelength_gridS;                                  P0_signal.LP_21_a   = Pin*ones(1,length(Signal.lambda.LP_21_a));
Signal.lambda.LP_02     = Wavelength_gridS;                                  P0_signal.LP_02     = Pin*ones(1,length(Signal.lambda.LP_02));

 Pump.modos = "01" ;
Wavelength_gridP=980e-9;
Ppump= 1000e-3; %[W]

Pump.lambda.LP_01   = Wavelength_gridP;                         P0_pump.LP_01   = Ppump  ;  

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
h=6.62607015*10^(-34);
P.Np=2; 
P.Fc=c/Wavelength_gridS(ceil(length(Wavelength_gridS)/2)); P.Fb = 50e9; 
ASE= -200;

    % Datos de la fibra
Fibra.nucleos = 1;                                           % Numero de nucleos
Fibra.largo = 3; Fibra.radio = 5.5e-6 ; Fibra.N = 7e24; 
Fibra.n1 = 1.45 ;   Fibra.IndexContrast=0.01;
Fibra.AN=Fibra.n1*sqrt(2*Fibra.IndexContrast);
Fibra.n2 =sqrt((Fibra.n1^2-Fibra.AN^2));

Fibra.dvk=P.Fb;

Fibra.WaitBar = 1; Fibra.Avance = 1;    % Despliegue de info
Fibra.ASEFlag = 1;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)

%% Primer Amplificador
% tic;
%EDFA = EDFA_MMvpi2(Fibra,Signal,Pump,ASE);%EDFA_MMvPCCv3(fibra,signal,pump,ASE);
% EDFA = EDFA_MM(fibra,signal,pump,ASE); %197.82 segundos
% t_end = toc; fprintf('Tiempo de cómputo: %.2f segundos\n', t_end);


%% Iterar en varios Spans de fibra-Amplificación

Nspans = 3;
LargoFibra = 40;        % [km]

%alp = load('Dynamic_Attenuation.dat');
%Attenuation = @(f) interp1( (alp(:,1).*1e-9) , (alp(:,2)) ,f);

%alpha = Attenuation(Wavelength_gridS);               % [dB/km]
%Att = alpha*LargoFibra;

%Span.EDFA1 = EDFA;

%tic;
%Fibra.Nspans = Nspans+1;
%for span=1:Nspans
%     Fibra.span = span+1;
%     fprintf('Iniciando Span %.0f de %.0f\n', span+1,Nspans+1);
% 
%     %  %% ---- Actualizar señales de entrada ---- %% %
%     % Señal y ASE
%     for i=1:length(Signal.modos)        % Potencia de señal y ASE a W
%             P0_signal.(ModoS(i)) = 1e-3*10.^(( (Span.(strcat("EDFA",num2str(span))).Nucleo1.salida.signal.potencia_dBm.(ModoS(i)))' -Att)/10);
%             P0_ASE.(ModoS(i)) = 1e-3*10.^(( (Span.(strcat("EDFA",num2str(span))).Nucleo1.Pap.(ModoS(i))(:,end)) -Att' )/10);
%     end ;clear i;
%     Signal.P0 = P0_signal; 
%     EDFA = Span_EDFA_MMvpi2(Fibra,Signal,Pump,P0_ASE);
% 
%     %  %% ---- Guardar Datos de Iteracion  ---- %% %
%     Span.(strcat("EDFA",num2str(span+1))) = EDFA;
% end
% span_time=toc; fprintf('Tiempo Total de cómputo: %.2f segundos\n', t_end+span_time )


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc
load("NewSpan.mat")
Nspans = length(fieldnames(Span))-1; LargoFibra = 40;
Signal.modos = ["11a" "21a" "02"] ; ModoS = ["11_a" "21_a" "02"] ; 
Pump.modos = "01" ; ModoP = "01" ; 
Signal.NumberOfChannels=41 ; Wavelength_gridS=linspace(1530,1570,Signal.NumberOfChannels).*1e-9; 
Signal.lambda.LP_11_a = Wavelength_gridS;  Signal.lambda.LP_21_a   = Wavelength_gridS; Signal.lambda.LP_02   = Wavelength_gridS; 



GrafSpans.z_total = [Span.EDFA1.Nucleo1.z]; 
GrafSpans.Signal_total = Span.EDFA1.Nucleo1.signal.Potencia_dBm;
GrafSpans.ASE_total = Span.EDFA1.Nucleo1.Pap;

z = Span.EDFA1.Nucleo1.z;

% for span = 1:Nspans
%     GrafSpans.z_total = [GrafSpans.z_total (z+(0.1*LargoFibra+GrafSpans.z_total(end)))];
%     for i=1:length(Signal.modos)
%         GrafSpans.Signal_total.(strcat("LP_",ModoS(i))) = [GrafSpans.Signal_total.(strcat("LP_",ModoS(i))) Span.(strcat("EDFA",num2str(span+1))).Nucleo1.signal.Potencia_dBm.(strcat("LP_",ModoS(i)))];
%         GrafSpans.ASE_total.(strcat("LP_",ModoS(i))) = [GrafSpans.ASE_total.(strcat("LP_",ModoS(i))) Span.(strcat("EDFA",num2str(span+1))).Nucleo1.Pap.(strcat("LP_",ModoS(i)))];
%     end ; clear i;
% end

    % % Graficos

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% % SAVE:
% % print -dpdf 'NAME'
% 
close all
 
graf.Nc = length(fieldnames(Span.EDFA1)); 
graf.z = Span.EDFA1.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';


% %%% Power Out
% s=1;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.SignalPump_OUT.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.salida.signal.potencia_dBm.(strcat("LP_",ModoS(s)))(:,end);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.SignalPump_OUT.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Potencia de Salida en Modo LP',Signal.modos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;



 % Ganancias
 % LP11a
 s=1;
figure(s)
for span = 1:Nspans+1
    set(gca,'FontSize',8)
    graf.ganancias.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.salida.ganancias.(strcat("LP_",ModoS(s)));
    leyenda = strcat(" EDFA ",num2str(span) );
    ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
    plot(ejex,graf.ganancias.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
    title(strcat('Ganancias Por Amplificador en Modo LP',Signal.modos(s)),'FontSize',14) ; legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
    xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
end ; clear s ejex leyenda;

 % LP21a
% s=2;
% figure
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.ganancias.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.salida.ganancias.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.ganancias.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Ganancias Por Amplificador en Modo LP',Signal.modos(s)),'FontSize',14) ; legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;

 % LP02
% s=3;
% figure
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.ganancias.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.salida.ganancias.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.ganancias.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Ganancias Por Amplificador en Modo LP',Signal.modos(s)),'FontSize',14) ; legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;


%  % NF
% LP11a
% s=1;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.NF.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.NF.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.NF.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Figuras de Ruido Por Amplificador en Modo LP',Signal.modos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;

% % LP21a
% s=2;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.NF.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.NF.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.NF.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Figuras de Ruido Por Amplificador en Modo LP',Signal.modos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;

% % LP02
% s=3;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.NF.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.NF.(strcat("LP_",ModoS(s)));
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.NF.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('Figuras de Ruido Por Amplificador en Modo LP',Signal.modos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;


% % OSRN OUT

%LP11a
% s=1;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.OSNR_OUT.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.OSNR.(strcat("LP_",ModoS(s)))(:,end);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.OSNR_OUT.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('OSNR en Modo LP',Signal.modos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;


% %LP21a
% s=2;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.OSNR_OUT.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.OSNR.(strcat("LP_",ModoS(s)))(:,end);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.OSNR_OUT.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('OSNR en Modo LP',Signal.modos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;

%LP02
% s=3;
% for span = 1:Nspans+1
%     set(gca,'FontSize',8)
%     graf.OSNR_OUT.(strcat("LP_",ModoS(s))) = Span.(strcat("EDFA",num2str(span))).Nucleo1.OSNR.(strcat("LP_",ModoS(s)))(:,end);
%     leyenda = strcat(" EDFA ",num2str(span) );
%     ejex = Signal.lambda.(strcat("LP_",ModoS(s))).*1e9;
%     plot(ejex,graf.OSNR_OUT.(strcat("LP_",ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     title(strcat('OSNR en Modo LP',Signal.modos(s)),'FontSize',14) ; 
%     legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 4,'FontSize',9); 
%     xlabel('Longitud de Onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
% end ; clear s ejex leyenda;


% span_ver=4;
% for n = 1:graf.Nc
%     figure(n+1)
%     
%     % Propagación Señal
%     if Fibra.ASEFlag == 0
%         for s = 1:length(Signal.modos) % Grafico
%             graf.Signal.(strcat("LP_",Signal.modos(s))) = Span.(strcat("EDFA",num2str(span_ver))).(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
%             subplot 221
%             for f = 1:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
%                 plot(graf.z,graf.Signal.(strcat("LP_",Signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
%             end
%         end; clear s f ;
%     else % No calcula espectro ASE
%         for s = 1:length(Signal.modos) % Grafico
%             graf.Signal.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
%             subplot 121
%             for f = 1:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
%                 plot(graf.z,graf.Signal.(strcat("LP_",Signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
%             end
%         end; clear s f ;
%     end
% 
%     % Ganancias
%     if Fibra.ASEFlag == 0
%         for s = 1:length(Signal.modos)
%             graf.ganancias.(strcat("LP_",Signal.modos(s))) = Span.(strcat("EDFA",num2str(span_ver))).(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
%             leyenda = strcat((strcat("LP",Signal.modos(s))));
%             ax = subplot (2,2,2);
%             
%             ejex = Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9;
%             plot(ejex,graf.ganancias.(strcat("LP_",Signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%             %ax.ColorOrderIndex = s;
%             xlabel('λ [nm]') ; ylabel(ylab)
%             
%         end ; clear s ejex leyenda;
%     else % No calcula espectro ASE
%         for s = 1:length(Signal.modos)
%             graf.ganancias.(strcat("LP_",Signal.modos(s))) = Span.(strcat("EDFA",num2str(span_ver))).(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
%             leyenda = strcat((strcat("LP",Signal.modos(s))));
%             ax = subplot (1,2,2);
%             
%             ejex = Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9;
%             plot(ejex,graf.ganancias.(strcat("LP_",Signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%             %ax.ColorOrderIndex = s;
%             xlabel('λ [nm]') ; ylabel(ylab)
%             
%         end ; clear s ejex leyenda;
%     end
%     
%     % Espectro de ganancia
%     if Fibra.ASEFlag == 0
%         for ase = 1:length(Signal.modos)
%             graf.freq_ase = Span.(strcat("EDFA",num2str(span_ver))).(strcat('Nucleo',int2str(n))).ASE_Spectrum.lambdas*1e9;
%             graf.ASE_Spectrum.(strcat("LP_",Signal.modos(ase))) = Span.(strcat("EDFA",num2str(span_ver))).(strcat("Nucleo",int2str(n))).ASE_Spectrum.mag.(strcat("LP_",Signal.modos(ase)));
%             subplot(2,2,[3,4])
%             plot( graf.freq_ase, graf.ASE_Spectrum.(strcat("LP_",Signal.modos(ase))) , 'DisplayName', strcat( "LP",Signal.modos(ase) ) ) ; hold on;
%             xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Potencia a la salida del EDFA') ; legend()
%         end ; clear ase xlab ylab t_end ;
%     end
%     
% end

% %% Graficos 3D de Modos 
% 
% % Señal
% figure()
% graficar_modos(Fibra,Signal.modos,Signal.lambda.LP_01)
% sgtitle('Modos de Señal') 
% 
% % Bombeo
% figure()
% graficar_modos(Fibra,Pump.modos,Pump.lambda.LP_01)
% sgtitle('Modos de Bombeo') ; 
% 
% %% Graficar potencia vs frecuencia
% figure()
% xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
% graf.Nc = length(fieldnames(EDFA)); 
% for n = 1:graf.Nc
%     for s = 1:length(Signal.modos)
%         graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
%         leyenda = strcat("LP",Signal.modos(s));
%         ejex = fliplr((3*10^8./(Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9)));
%         plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",Signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%         %ax.ColorOrderIndex = s;
%         xlabel('Frecuencia [THz]') ; ylabel(ylab)
% 
%     end ; clear s ejex leyenda n;
% end
% 
% 
% %% Graficar potencia vs wavelenght
% % n=n+1; figure()
% % xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
% % graf.Nc = length(fieldnames(EDFA)); 
% % for n = 1:graf.Nc
% %    for s = 1:length(Signal.modos)
% %        graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
% %        leyenda = strcat((strcat("LP",Signal.modos(s))));
% %        ejex = fliplr((3*10^8./(Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9)));
% %        plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",Signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
% %        %ax.ColorOrderIndex = s;
% %        xlabel('Frecuencia [THz]') ; ylabel(ylab)
% %        set(gca,'xdir','reverse')
% %    end ; clear s ejex leyenda;
% % end
% 
