close all ; clc ; clear all


%% Datos Entrada

In.Fibra.Length                 = 50;                                                       % fibre length (km)
In.Fibra.T                      = 25;                                                       % Temperatura Fibra (ambiente)
In.Fibra.PolarizationFactor     = 0.5;                                                      % C_R_max
%In.Pump.LP01.Wavelengths        = [1400 1420 1440 1460] ;                                  % [nm]
%In.Pump.LP11.Wavelengths        = [1460 1470] ;                                            % [nm]
%In.Pump.LP11.Powers             = 100*1e-3*ones( 1,length(In.Pump.LP11.Wavelengths) );             % [mW]
%In.Pump.LP11.Alpha              = [0.25];                                                  % [dB/km]

In.Pump.LP01.Wavelengths        = [1450 ] ;                                             % [nm]
In.Pump.LP01.Powers             = 300*1e-3*ones( 1,length(In.Pump.LP01.Wavelengths) );      % [mW]
In.Pump.LP01.Alpha              = [0.28];   

In.Fibra.n1=1.46;  In.Fibra.n2=1.450; In.Fibra.radio=25e-6;

In.Fibra.RamanMethod            = 'Forward&Backward';            % 'Forward', 'Backward', 'Forward&Backward'
%P.Fibra.RamanNoiseMethod       = P.Fibre.RamanMethod;
%P.Fibra.RamanNoise             = 0;                    % 1: add noise every step, 0 add noise at the end of span
%P.Fibra.RamanRIN               = 0;        
%P.Fibra.RamanPumpRIN           = -130;                 % Raman pump RIN level (dBm/carrier)

In.Signal.LP01.Wavelengths      = [1500:10:1600];%[1500 1540 1550 1560 1600];
In.Signal.LP01.Powers           = -15*ones( 1,length(In.Signal.LP01.Wavelengths) );           %[dBm]
In.Signal.LP01.Alpha            = 0.25;                          % Fibre attenuation @Signal WL (dB/km)
In.Signal.LP11.Wavelengths      = [1560 1580 1600 ];
In.Signal.LP11.Powers           = -15*ones( 1,length(In.Signal.LP11.Wavelengths) );           %[dBm]
In.Signal.LP11.Alpha            = 0.2;                          % Fibre attenuation @Signal WL (dB/km)
%signal.modos.LP02              = [1540:10:1570];
%signal.modos.LP11              = [1540:10:1570];

%% Calculo de amplificación
tic;
Raman = RamanMM(In) ; tend = toc; fprintf("Tiempo de cómputo: %.2f",tend);

%% Graficar
close all; 
z = Raman.z;
for mp = length(Raman.ModoP)
    pump.(Raman.ModoP{mp}).forward = Raman.Pump.forward.(Raman.ModoP{mp});
    pump.(Raman.ModoP{mp}).backward = Raman.Pump.backward.(Raman.ModoP{mp});
end
signal = Raman.Sig.Power;
gain = Raman.Sig.Gain;

fs = 1; fp =1;


for ms = 1:length(Raman.ModoS)% Potencia Señal
    figure(1); 
    strlambda = strcat( num2str( In.Signal.(Raman.ModoS{ms}).Wavelengths(:)) , "nm");
    subplot(1,length(Raman.ModoS),fs) ; plot(z,10.*log10(signal.(Raman.ModoS{ms})./1e-3) ) ,ylabel("Potencia [dBm]") ; hold on
    %subplot(1,fs,fs) ; plot( z,(signal.(Raman.ModoS{ms}) ) ,ylabel("Potencia [mW]"); hold on
    xlabel("Largo [km]") ; title(strcat("Signal Modo ", Raman.ModoS{ms})) ; legend(strlambda,'Location','southwest')
    fs = fs+1;
end

for mp = 1:length(Raman.ModoP)% Potencia Bombeo
    figure(2); 
    for i = 1:length(In.Pump.(Raman.ModoP{mp}).Wavelengths)
        DispName1 = strcat("Forward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
        DispName2 = strcat("Backward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
        %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).forward(i,:)./1e-3),"DisplayName","Forward") ,ylabel("Potencia [dBm]") , hold on
        %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).backward(i,:)./1e-3),"DisplayName","Backward") ,ylabel("Potencia [dBm]")
        subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).forward(i,:),"DisplayName",DispName1) ,ylabel("Potencia [W]") , hold on
        subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).backward(i,:),"DisplayName",DispName2) ,ylabel("Potencia [W]")
        xlabel("Largo [km]") ; title(strcat("Pump Modo ", Raman.ModoP{mp})) , legend()
    end
    fp = fp+1;
end

for ms = 1:length(Raman.ModoS)% Ganancias
    figure(3)
    plot(In.Signal.(Raman.ModoS{ms}).Wavelengths , gain.(Raman.ModoS{ms}))  ; hold on
    %plot(In.Signal.(Raman.ModoS{ms}).Wavelengths(i) , gain.(Raman.ModoS{ms})(i) ) ,ylabel("Potencia [mW]"); hold on
end
    figure(3) ; xlabel("Longitud de Onda nm") ,ylabel("Ganancia [dB]") ; title("Ganancias") ; legend(Raman.ModoS{:})


clear fs fp ms mp i strlambda;


%% Otros Archivos
% eta = load('RADynamic_Rayleigh.dat') ; gr = load('RamanGainEfficiency_SMF28.dat');
% figure(1);plot(eta(:,1),eta(:,2)) ; figure(2) ; plot(gr(:,1),gr(:,2))