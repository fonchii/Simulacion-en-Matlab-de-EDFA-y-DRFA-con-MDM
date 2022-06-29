% close all; 
clear all; clc ; close all

% Parámetros de entrada

    % Modos y canales de señal y bombeo
signal.NumberOfChannels=20;
signal.modos = ["01" "11_a"] ;
% Ejemplo Pump01 - channels 20
%Frequency_gridS=linspace(191.07234e12,196.7723e12,signal.NumberOfChannels);
% Ejemplo Pump12 - channels 50
Frequency_gridS=linspace(191.19421875e12,193.64421875e12,signal.NumberOfChannels);
c=299.792458e6; % [m/s]
Wavelength_gridS=c./Frequency_gridS;

Pin=0; %[dBm]

signal.lambda.LP_01     = Wavelength_gridS;                                  P0_signal.LP_01     = Pin*ones(1,length(signal.lambda.LP_01));
signal.lambda.LP_11_a   = Wavelength_gridS;                                  P0_signal.LP_11_a   = Pin*ones(1,length(signal.lambda.LP_11_a));

pump.modos = "12_a" ;
Wavelength_gridP=980e-9;
Ppump= 1000e-3; %[W]

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
P.Fc=c/Wavelength_gridS(ceil(length(Wavelength_gridS)/2)); P.Fb = 50e9; 
ASE= -200;

    % Datos de la fibra
fibra.nucleos = 1;                                           % Numero de nucleos
fibra.largo = 5; fibra.radio = 5.5e-6 ; fibra.N = 7e24; 
fibra.n1 = 1.45 ;   fibra.IndexContrast=0.01;
fibra.AN=fibra.n1*sqrt(2*fibra.IndexContrast);
fibra.n2 =sqrt((fibra.n1^2-fibra.AN^2));

fibra.dvk=P.Fb;

fibra.WaitBar = 1; fibra.Avance = 1;    % Despliegue de info
fibra.ASEFlag = 1;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)

%%
tic;
EDFA = EDFA_MMvpi2(fibra,signal,pump,ASE);%EDFA_MMvPCCv3(fibra,signal,pump,ASE);
% EDFA = EDFA_MM(fibra,signal,pump,ASE); %197.82 segundos
t_end = toc; fprintf('Tiempo de cómputo: %.2f segundos\n', t_end);


SPAN.EDFA = EDFA;

%% Iterar en varios Spans de figra
Nspans = 3;
LargoFibra = 100;        % [km]
alpha = 0.2*ones(length(Wavelength_gridS),1);               % [Np/km]
Att = -alpha*LargoFibra;


%for i=1:Nspans
P0_signal.LP_01     = Pin*ones(1,length(signal.lambda.LP_01));

    %% Graficos

close all
graf.Nc = length(fieldnames(EDFA)); 
graf.z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
for n = 1:graf.Nc
    figure(n)
    
    % Señal
    if Fibra.ASEFlag == 0
        for s = 1:length(Signal.modos) % Grafico
            graf.signal.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
            subplot 221
            for f = 1:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
                plot(graf.z,graf.signal.(strcat("LP_",Signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
            end
        end; clear s f ;
    else % No calcula espectro ASE
        for s = 1:length(Signal.modos) % Grafico
            graf.signal.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(strcat("LP_",Signal.modos(s)));
            subplot 121
            for f = 1:length(Signal.lambda.(strcat("LP_",Signal.modos(s))))
                plot(graf.z,graf.signal.(strcat("LP_",Signal.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Signal.modos(s))," @"),strcat(int2str(Signal.lambda.(strcat("LP_",Signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
            end
        end; clear s f ;
    end
    
%     % Bombeo
%     for s = 1:length(Pump.modos) % Grafico
%         graf.pump.(strcat("LP_",Pump.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).pump.Potencia_dBm(:,:,s);
%         figure(n)
%         subplot 222
%         for f = 1:length(Pump.lambda.LP_01)
%             plot(graf.z,graf.pump.(strcat("LP_",Pump.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",Pump.modos(s))," @"),strcat(int2str(Pump.lambda.(strcat("LP_",Pump.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Pump}') ; legend(); grid on
%         end
%     end

    % Ganancias
    if Fibra.ASEFlag == 0
        for s = 1:length(Signal.modos)
            graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
            leyenda = strcat((strcat("LP",Signal.modos(s))));
            ax = subplot (2,2,2);
            
            ejex = Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9;
            plot(ejex,graf.ganancias.(strcat("LP_",Signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
            %ax.ColorOrderIndex = s;
            xlabel('λ [nm]') ; ylabel(ylab)
            
        end ; clear s ejex leyenda;
    else % No calcula espectro ASE
        for s = 1:length(Signal.modos)
            graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
            leyenda = strcat((strcat("LP",Signal.modos(s))));
            ax = subplot (1,2,2);
            
            ejex = Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9;
            plot(ejex,graf.ganancias.(strcat("LP_",Signal.modos(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
            %ax.ColorOrderIndex = s;
            xlabel('λ [nm]') ; ylabel(ylab)
            
        end ; clear s ejex leyenda;
    end
    
    % Espectro de ganancia
    if Fibra.ASEFlag == 0
        for ase = 1:length(Signal.modos)
            graf.freq_ase = EDFA.(strcat('Nucleo',int2str(n))).ASE_Spectrum.lambdas*1e9;
            graf.ASE_Spectrum.(strcat("LP_",Signal.modos(ase))) = EDFA.(strcat("Nucleo",int2str(n))).ASE_Spectrum.mag.(strcat("LP_",Signal.modos(ase)));
            subplot(2,2,[3,4])
            plot( graf.freq_ase, graf.ASE_Spectrum.(strcat("LP_",Signal.modos(ase))) , 'DisplayName', strcat( "LP",Signal.modos(ase) ) ) ; hold on;
            xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Potencia a la salida del EDFA') ; legend()
        end ; clear ase xlab ylab t_end ;
    end
    
end

%% Graficos 3D de Modos 

% Señal
figure()
graficar_modos(Fibra,Signal.modos,Signal.lambda.LP_01)
sgtitle('Modos de Señal') 

% Bombeo
figure()
graficar_modos(Fibra,Pump.modos,Pump.lambda.LP_01)
sgtitle('Modos de Bombeo') ; 

%% Graficar potencia vs frecuencia
figure()
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
graf.Nc = length(fieldnames(EDFA)); 
for n = 1:graf.Nc
    for s = 1:length(Signal.modos)
        graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
        leyenda = strcat("LP",Signal.modos(s));
        ejex = fliplr((3*10^8./(Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9)));
        plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",Signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
        %ax.ColorOrderIndex = s;
        xlabel('Frecuencia [THz]') ; ylabel(ylab)

    end ; clear s ejex leyenda n;
end


%% Graficar potencia vs wavelenght
% n=n+1; figure()
% xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
% graf.Nc = length(fieldnames(EDFA)); 
% for n = 1:graf.Nc
%    for s = 1:length(Signal.modos)
%        graf.ganancias.(strcat("LP_",Signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",Signal.modos(s)));
%        leyenda = strcat((strcat("LP",Signal.modos(s))));
%        ejex = fliplr((3*10^8./(Signal.lambda.(strcat("LP_",Signal.modos(s))).*1e9)));
%        plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",Signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%        %ax.ColorOrderIndex = s;
%        xlabel('Frecuencia [THz]') ; ylabel(ylab)
%        set(gca,'xdir','reverse')
%    end ; clear s ejex leyenda;
% end

