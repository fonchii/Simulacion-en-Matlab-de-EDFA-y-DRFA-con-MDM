close all ; clc ; clear all


%% Datos Entrada

In.Fibra.Length                 = 100;                                                       % fibre length (km)
In.Fibra.C_Rmax                 = 0.5;
In.Pump.Wavelengths             = [1400 1420 1440 1460] ;                                   % [nm]
In.Pump.Powers                  = 200*1e-3*ones(1,length(In.Pump.Wavelengths));             % [mW]
In.Pump.Alpha                   = [0.25];                                                   % [dB/km]

In.Fibra.RamanMethod            = 'Forward&Backward';            % 'Forward', 'Backward', 'Forward&Backward','Ideal','NoRaman'
%P.Fibra.RamanNoiseMethod       = P.Fibre.RamanMethod;
%P.Fibra.RamanNoise             = 0;                    % 1: add noise every step, 0 add noise at the end of span
%P.Fibra.RamanRIN               = 0;        
%P.Fibra.RamanPumpRIN           = -130;                 % Raman pump RIN level (dBm/carrier)

In.Signal.Wavelenghts.LP01      = [1500 1540 1550 1560 1600];
In.Signal.Powers                = -15*ones(1,length(In.Signal.Wavelenghts.LP01));           %[dBm]
In.Signal.Alpha                 = 0.2;                  % Fibre attenuation @Signal WL (dB/km)
%signal.modos.LP02              = [1540:10:1570];
%signal.modos.LP11              = [1540:10:1570];

%% Calculo de amplificación
tic;
Raman = RamanSMFWDM(In);toc

%% Graficar
z = Raman.z; pumpLP01 = Raman.Pump;
for mp = 1 length()
pumpf = Raman.Pump.forward' ; pumpb = Raman.Pump.backward';
sig = Raman.Sig;
figure(1) ; subplot(1,2,1);%plot(z , 10*log10(pumpf/1e-3),'DisplayName', "Forward") , hold on;
                            plot(z , (pumpf/1e-3),'DisplayName', "Forward") , hold on; ylabel("Potencia [mW]") ; xlabel("Largo [km]");
            %plot(z , 10*log10(pumpb/1e-3),'DisplayName', "Backward") , title('Pump'); legend()
            plot(z , (pumpb/1e-3),'DisplayName', "Backward") , 
subplot(1,2,2); plot(z , 10*log10(sig/1e-3)) , title('Signal');ylabel("Potencia [dBm]") ; xlabel("Largo [km]");