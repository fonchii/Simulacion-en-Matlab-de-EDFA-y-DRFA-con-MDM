close all ; clc ; clear all


%% Datos Entrada

In.Fibra.Length                 = 100;                                                       % fibre length (km)
In.Fibra.C_Rmax                 = 0.5;
%In.Pump.LP01.Wavelengths        = [1400 1420 1440 1460] ;                                   % [nm]
In.Pump.LP01.Wavelengths        = [1460] ;                                             % [nm]
In.Pump.LP01.Powers             = 400*1e-3*ones( 1,length(In.Pump.LP01.Wavelengths) );             % [mW]
In.Pump.LP01.Alpha              = [0.25];                                                   % [dB/km]

In.Fibra.RamanMethod            = 'Forward&Backward';            % 'Forward', 'Backward', 'Forward&Backward'
%P.Fibra.RamanNoiseMethod       = P.Fibre.RamanMethod;
%P.Fibra.RamanNoise             = 0;                    % 1: add noise every step, 0 add noise at the end of span
%P.Fibra.RamanRIN               = 0;        
%P.Fibra.RamanPumpRIN           = -130;                 % Raman pump RIN level (dBm/carrier)

In.Signal.LP01.Wavelengths      = [1500 1540 1550 1560 1600];
In.Signal.LP01.Powers           = -15*ones( 1,length(In.Signal.LP01.Wavelengths) );           %[dBm]
In.Signal.LP01.Alpha            = 0.2;                          % Fibre attenuation @Signal WL (dB/km)
In.Signal.LP11.Wavelengths      = [1560 1600];
In.Signal.LP11.Powers           = -15*ones( 1,length(In.Signal.LP11.Wavelengths) );           %[dBm]
In.Signal.LP11.Alpha            = 0.2;                          % Fibre attenuation @Signal WL (dB/km)
%signal.modos.LP02              = [1540:10:1570];
%signal.modos.LP11              = [1540:10:1570];

%% Calculo de amplificación
tic;
Raman = RamanMM(In) ; toc

%% Graficar
z = Raman.z;
for mp = length(Raman.ModoP)
    pump.(Raman.ModoP{mp}).forward = Raman.Pump.forward.(Raman.ModoP{mp});
    pump.(Raman.ModoP{mp}).backward = Raman.Pump.backward.(Raman.ModoP{mp});
end
signal = Raman.Sig;


fs = 1; fp =1;
for ms = 1:length(Raman.ModoS)
    figure(1); 
    subplot(1,length(Raman.ModoS),fs) ; plot(z,10.*log10(signal.(Raman.ModoS{ms})./1e-3)) ,ylabel("Potencia [dBm]")
    %subplot(1,fs,fs) ; plot(z,(signal.(Raman.ModoS{ms})) ,ylabel("Potencia [mW]")
    xlabel("Largo [km]") ; title(strcat("Signal Modo ", Raman.ModoS{ms}))
    fs = fs+1;
end

for mp = 1:length(Raman.ModoP)
    figure(2); 
    %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).forward./1e-3),"DisplayName","Forward") ,ylabel("Potencia [dBm]") , hold on
    %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).backward./1e-3),"DisplayName","Backward") ,ylabel("Potencia [dBm]")
    subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).forward,"DisplayName","Forward") ,ylabel("Potencia [W]") , hold on
    subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).backward,"DisplayName","Backward") ,ylabel("Potencia [W]")
    xlabel("Largo [km]") ; title(strcat("Pump Modo ", Raman.ModoP{mp})) , legend()
    fp = fp+1;
end
clear fs fp;
