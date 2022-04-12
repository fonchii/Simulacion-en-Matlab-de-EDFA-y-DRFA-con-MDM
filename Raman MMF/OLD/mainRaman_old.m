close all ; clc ; clear all


%% Datos Fibra
P.Att    = 0.2;             % Fibre attenuation (dB/km)
P.Length = 80;            % fibre length (km)
P.Fibre.C_Rmax = 0.5;
P.Fibre.RamanWavelengths = [1420 1450 1470 1500];   % [nm]
P.Fibre.RamanPowers      = [300 700 200 200]*1e-3;  % [mW]
P.Fibre.PumpAlpha        = [0.25];                  % [dB/km]
P.Fibre.RamanMethod = 'Forward';            % 'Forward', 'Backward', 'Forward&Backward','Ideal','NoRaman'
P.Fibre.RamanNoiseMethod = P.Fibre.RamanMethod;
P.Fibre.RamanNoise        = 0;        % 1: add noise every step, 0 add noise at the end of span
P.Fibre.RamanRIN          = 0;        %
P.Fibre.RamanPumpRIN      = -130;     % Raman pump RIN level (dBm/carrier)

signal.Ps0 = -15;
signal.modos.LP01 = [1540:10:1570];
signal.modos.LP02 = [1540:10:1570];
signal.modos.LP11 = [1540:10:1570];

%%

Raman = RamanMM(signal,P);

%%
z = Raman.z; pumpLP01 = Raman.Pump
pump = Raman.Pump.LP01(2,:) ; sig = Raman.Sig.LP01(2,:);
figure(1) ; subplot(1,2,1);plot(z , 10*log10(pump/1e-3)) , title('Pump');
subplot(1,2,2); plot(z , 10*log10(sig/1e-3)) , title('Signal');