% MAIN
    % Frecuencias
close all; clear all; clc
freq_central_signal = 1550;
numero_freqs_signal = 3; 
freq_central_pump = 980;
numero_freqs_pump = 1;

    % Modos
signal.modos = ["01","02"] ;
signal.lambda = wavelength(numero_freqs_signal , freq_central_signal);
pump.modos = ["01","11_a","11_b"] ;
pump.lambda = wavelength(numero_freqs_pump , freq_central_pump);

    %POTENCIAS
P0_signal = - 20 * ones(1,numero_freqs_signal);     %dBm
P0_signal = 1e-3*10.^(P0_signal/10);                %mW
P0_pump = 100e-3 * ones(1,numero_freqs_pump);       %mW
%P0_pump = [1e-3,2e-3,5e-3,10e-3,15e-3,20e-3,25e-3,50e-3,75e-3,100e-3];
signal.P0 = P0_signal; 
pump.P0 = P0_pump;
Pase = -50;                                         %dBm

    % Datos de la fibra
%fibra.nucleos = 2;      % Numero de nucleos
fibra.nucleos = 1;
fibra.largo = 10;  fibra.AN = 0.2 ; fibra.radio = 25e-6 ; fibra.N = 3e24; 
fibra.n1 = 1.46; fibra.n2 = 1.4462;

tic;
EDFA = EDFA_MM(fibra,signal,pump,Pase);
t_end = toc; fprintf('Tiempo de iteración: %.2f segundos\n', t_end);

%% Graficos de Potencias

close all ; xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
Nch = length(fieldnames(EDFA));
N_fig = 0;
tic;
for j = 1:1:Nch    % GRAFICAR TODOS LOS CANALES ESPACIALES
    Nsign = length(EDFA.(strcat('Nucleo',int2str(j))).signal.lambdas);
    Npump = length(EDFA.(strcat('Nucleo',int2str(j))).pump.lambdas);
    Pp_dbm = EDFA.(strcat('Nucleo',int2str(j))).pump.Potencia_dBm;
    Ps_dbm = EDFA.(strcat('Nucleo',int2str(j))).signal.Potencia_dBm;
    Pase_dbm = EDFA.(strcat('Nucleo',int2str(j))).Pase;
    OSNR = EDFA.(strcat('Nucleo',int2str(j))).OSNR;
    lambda_s = EDFA.(strcat('Nucleo',int2str(j))).signal.lambdas;
    lambda_p = EDFA.(strcat('Nucleo',int2str(j))).pump.lambdas;
    z = EDFA.(strcat('Nucleo',int2str(j))).z; N_fig = N_fig+1;
    figure(j)
    subplot 221
    for i = 1:1:Npump
        plot(z,Pp_dbm(i,:),'DisplayName',strcat(int2str(lambda_p(i)*1e9),' nm')); hold on; xlabel(xlab) ; ylabel(ylab) ;  title('P_{Pump}') ;  legend(); grid on ;
    end
    subplot 222
    for i = 1:2:Nsign
        plot(z,Ps_dbm(i,:),'DisplayName',strcat(int2str(lambda_s(i)*1e9) ,' nm')) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on
    end
    subplot 223
    for i = 1:2:Nsign
        plot(z,Pase_dbm(i,:),'DisplayName',strcat(int2str(lambda_s(i)*1e9) ,' nm')) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{ASE}') ; legend(); grid on;
    end
    subplot 224
    for i = 1:2:Nsign
        plot(z,OSNR(i,:),'DisplayName',strcat(int2str(lambda_s(i)*1e9) ,' nm')) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('OSNR') ; legend(); grid on
    end
end

figure (N_fig+1)
plot(EDFA.(strcat('Nucleo',int2str(1))).ASE_Spectrum.lambdas*1e9, EDFA.Nucleo1.ASE_Spectrum.mag(:,2)); xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Espectro de ganancia')
t_plot = toc; fprintf('Tiempo de ploteo: %.2f segundos\n', t_plot);
