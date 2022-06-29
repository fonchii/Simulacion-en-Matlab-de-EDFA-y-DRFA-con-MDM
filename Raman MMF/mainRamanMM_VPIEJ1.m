close all ; clc ; clear all


%% Datos Entrada

% FIBRA
In.Fibra.RamanMethod              = 'Forward';                   % 'Forward', 'Backward', 'Forward&Backward'
In.Fibra.AttenuationMethod        = 'Dynamic';                    % 'Dynamic' , 'Static'
In.Fibra.Length                   = 100;                          % fibre length (km)
In.Fibra.T                        = 25;                           % Temperatura Fibra (ambiente)
In.Fibra.PolarizationFactor       = 0.5;                          % C_R_max
In.Fibra.n1=1.45;  In.Fibra.n2=1.4354; In.Fibra.radio=5e-6;

% BOMBEOS : 
%     % LP01
% In.Pump.LP01.Wavelengths       = [1430 1450 1470] ;                                               % [nm]
% In.Pump.LP01.Powers            = [200 20 80]*1e-3;                                                % [mW]
% In.Pump.LP01.Alpha             = [0.29]; 
    % LP11
In.Pump.LP11.Wavelengths          = [1430 1445 1455 1470 1490] ;                                               % [nm]
In.Pump.LP11.Powers               = [50 30 30 50 90]*1e-3;                                                  % [mW]
In.Pump.LP11.Alpha                = 0.29;                                                           % [dB/km]
    % LP21
In.Pump.LP21.Wavelengths          = [1430 1445 1455 1470 1490] ;                                               % [nm]
In.Pump.LP21.Powers               = [50 30 30 50 90]*1e-3;                                                  % [mW]
In.Pump.LP21.Alpha                = 0.29;                                                           % [dB/km]

%     % LP02
% In.Pump.LP02.Wavelengths        = [1430 1450 1470 ] ;                                             % [nm]
% In.Pump.LP02.Powers             = 50.*1e-3*ones( 1,length(In.Pump.LP02.Wavelengths) );            % [mW]
% In.Pump.LP02.Alpha              = [0.29];                                                         % [dB/km]


% SEÑALES : 
Nch = 20;
    % LP01  
In.Signal.LP01.Wavelengths        = linspace(1525,1570,Nch);
In.Signal.LP01.Powers             = -15*ones( 1,length(In.Signal.LP01.Wavelengths) );                 %[dBm]
In.Signal.LP01.Alpha              = 0.2;                                                              % [dB/km]
%     % LP11
% In.Signal.LP11.Wavelengths      = [1530:1:1570];
% In.Signal.LP11.Powers           = -15*ones( 1,length(In.Signal.LP11.Wavelengths) );               %[dBm]
% In.Signal.LP11.Alpha            = 0.2;                          % Fibre attenuation @Signal WL (dB/km)

    % LP11_b
In.Signal.LP11_a.Wavelengths      = linspace(1525,1570,Nch);
In.Signal.LP11_a.Powers           = -15*ones( 1,length(In.Signal.LP11_a.Wavelengths) );             %[dBm]
In.Signal.LP11_a.Alpha            = 0.2;                                                            % [dB/km]

    % LP11_b
In.Signal.LP11_b.Wavelengths      = linspace(1525,1570,Nch);
In.Signal.LP11_b.Powers           = -15*ones( 1,length(In.Signal.LP11_b.Wavelengths) );             %[dBm]
In.Signal.LP11_b.Alpha            = 0.2;                                                            % [dB/km]

%     % LP21
% In.Signal.LP21.Wavelengths      = [1530:5:1565];
% In.Signal.LP21.Powers           = -15*ones( 1,length(In.Signal.LP21.Wavelengths) );               %[dBm]
% In.Signal.LP21.Alpha            = 0.25;                                                           % [dB/km]
%     % LP02
% In.Signal.LP02.Wavelengths      = [1530:5:1565];
% In.Signal.LP02.Powers           = -15*ones( 1,length(In.Signal.LP02.Wavelengths) );               %[dBm]
% In.Signal.LP02.Alpha            = 0.25;                                                           % [dB/km]


%% Calculo de amplificación
tic;
%Raman = RamanMM_comp(In) ; tend = toc; fprintf("Tiempo de cómputo: %.2f",tend);
Raman = RamanMMv3(In) ; tend = toc; fprintf("Tiempo de cómputo: %.2fs\n",tend);

%% Graficar

close all; 
z = Raman.z;
for mp = 1:length(Raman.ModoP)
    pump.(Raman.ModoP{mp}).forward = Raman.Pump.forward.(Raman.ModoP{mp});
    pump.(Raman.ModoP{mp}).backward = Raman.Pump.backward.(Raman.ModoP{mp});
end
signal = Raman.Sig.Power;
gain = Raman.Sig.GainOnOFF;
osnr = Raman.OSNR;

fs = 1; fp =1;


for ms = 1:length(Raman.ModoS)% Potencia Señal
    figure(1); 
    strlambda = strcat( num2str( In.Signal.(Raman.ModoS{ms}).Wavelengths(ms)) , "nm");
    subplot(1,length(Raman.ModoS),fs) ; plot(z,10.*log10(signal.(Raman.ModoS{ms})./1e-3) , 'DisplayName',strlambda) ,ylabel("Potencia [dBm]") ; hold on
    %subplot(1,length(Raman.ModoS),fs) ; plot(z,10.*log10(signal.Off.(Raman.ModoS{ms})./1e-3) ) ,ylabel("Potencia [dBm]") ; hold on
    %subplot(1,fs,fs) ; plot( z,(signal.(Raman.ModoS{ms}) ) ,ylabel("Potencia [mW]"); hold on
    xlabel("Largo [km]") ; title(strcat("Signal Modo ", Raman.ModoS{ms})) ; legend('Location','best')
    fs = fs+1;
end
fs = 1;
for ms = 1:length(Raman.ModoS)% Potencia Señal Off
    %figure(1); 
    strlambda = strcat( num2str( In.Signal.(Raman.ModoS{ms}).Wavelengths(ms)) , "nm OFF");
    h=subplot(1,length(Raman.ModoS),fs) ; h.ColorOrderIndex=1;
    plot(z,10.*log10(signal.Off.(Raman.ModoS{ms})./1e-3) , 'DisplayName',strlambda) ,ylabel("Potencia [dBm]") ; hold on
    %subplot(1,fs,fs) ; plot( z,(signal.Off.(Raman.ModoS{ms})) ) ,ylabel("Potencia [mW]"); hold on
    xlabel("Largo [km]") ; title(strcat("Signal Modo ", Raman.ModoS{ms})) ; legend('Location','best')
    fs = fs+1;
end

switch In.Fibra.RamanMethod
    case 'Forward&Backward'
        for mp = 1:length(Raman.ModoP)% Potencia Bombeo
            figure(2); 
            for i = 1:length(In.Pump.(Raman.ModoP{mp}).Wavelengths)
                DispName1 = strcat("Forward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
                DispName2 = strcat("Backward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
                %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).forward(i,:)./1e-3),"DisplayName","Forward") ,ylabel("Potencia [dBm]") , hold on
                %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).backward(i,:)./1e-3),"DisplayName","Backward") ,ylabel("Potencia [dBm]")
                subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).forward(i,:),"DisplayName",DispName1) ,ylabel("Potencia [W]") , hold on
                subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).backward(i,:),"DisplayName",DispName2) ,ylabel("Potencia [W]")
                xlabel("Largo [km]") ; title(strcat("Pump Modo ", Raman.ModoP{mp})) , legend("Location","best")
            end
            fp = fp+1;
        end
    case 'Backward'
        for mp = 1:length(Raman.ModoP)% Potencia Bombeo
            figure(2); 
            for i = 1:length(In.Pump.(Raman.ModoP{mp}).Wavelengths)
                DispName2 = strcat("Backward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
                %subplot(1,length(Raman.ModoP),fp) ; plot(z,10.*log10(pump.(Raman.ModoP{mp}).backward(i,:)./1e-3),"DisplayName",DispName2) ,ylabel("Potencia [dBm]"), hold on
                subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).backward(i,:),"DisplayName",DispName2) ,ylabel("Potencia [W]"), hold on
                xlabel("Largo [km]") ; title(strcat("Pump Modo ", Raman.ModoP{mp})) , legend("Location","best")
            end
            fp = fp+1;
        end
    case 'Forward'
        for mp = 1:length(Raman.ModoP)% Potencia Bombeo
            figure(2); 
            for i = 1:length(In.Pump.(Raman.ModoP{mp}).Wavelengths)
                DispName1 = strcat("Forward " ,  num2str( In.Pump.(Raman.ModoP{mp}).Wavelengths(i)) , "nm" );
                subplot(1,length(Raman.ModoP),fp) ; plot(z,pump.(Raman.ModoP{mp}).forward(i,:),"DisplayName",DispName1) ,ylabel("Potencia [W]") , hold on
                xlabel("Largo [km]") ; title(strcat("Pump Modo ", Raman.ModoP{mp})) , legend("Location","best")
            end
            fp = fp+1;
        end
end

for ms = 1:length(Raman.ModoS)% Ganancias
    figure(3)
    plot(In.Signal.(Raman.ModoS{ms}).Wavelengths , gain.(Raman.ModoS{ms}), '-o')  ; hold on
    %plot(In.Signal.(Raman.ModoS{ms}).Wavelengths(i) , gain.(Raman.ModoS{ms})(i) ) ,ylabel("Potencia [mW]"); hold on
end
    figure(3) ; xlabel("Longitud de Onda nm") ,ylabel("Ganancia [dB]") ; title("Ganancias On-Off") ; legend(Raman.ModoS{:},"Location","best")

for ms = 1:length(Raman.ModoS)% OSNR
    figure(4)
    plot(In.Signal.(Raman.ModoS{ms}).Wavelengths , osnr.(Raman.ModoS{ms}), '-o')  ; hold on
    %plot(In.Signal.(Raman.ModoS{ms}).Wavelengths(i) , osnr.(Raman.ModoS{ms})(i) ) ,ylabel("Potencia [mW]"); hold on
end
    figure(4) ; xlabel("Longitud de Onda nm") ,ylabel("Magnitud [dB]") ; title("OSNR") ; legend(Raman.ModoS{:},"Location","best")


clear fs fp ms mp i strlambda DispName1 DispName2;

% for mp = length(Raman.ModoP)
%    pumpv2.(Raman.ModoP{mp}).forward = Ramanv2.Pump.forward.(Raman.ModoP{mp});
%    pumpv2.(Raman.ModoP{mp}).backward = Ramanv2.Pump.backward.(Raman.ModoP{mp});
% end
% plot(z,pump.(Raman.ModoP{1}).backward(1,:)) ; hold on ; plot(z,pumpv2.(Raman.ModoP{1}).backward(1,:)) 

% %%% Atenuaciones

% alp = load('RADynamic_Attenuation.dat');
% plot(alp(:,1),alp(:,2))


%% Otros Archivos
% eta = load('RADynamic_Rayleigh.dat') ; gr = load('RamanGainEfficiency_SMF28.dat');
% figure(1);plot(eta(:,1),eta(:,2)) ; figure(2) ; plot(gr(:,1),gr(:,2))