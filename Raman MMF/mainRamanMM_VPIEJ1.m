close all ; clc ; clear all

%% Datos Entrada

% FIBRA
In.Fibra.RamanMethod              = 'Backward';                   % 'Forward', 'Backward', 'Forward&Backward'
In.Fibra.AttenuationMethod        = 'Static';                    % 'Dynamic' , 'Static'
In.Fibra.Length                   = 100;                          % fibre length (km)
In.Fibra.T                        = 25;                           % Temperatura Fibra (ambiente)
In.Fibra.PolarizationFactor       = 0.66;                          % C_R_max
In.Fibra.n1=1.47;  In.Fibra.n2=1.42; In.Fibra.radio=5.04626e-6; In.Fibra.area=pi*(In.Fibra.radio)^2;

% BOMBEOS : 
    % LP01
In.Pump.LP01.Wavelengths       = 1500 ;                                                    % [nm]
In.Pump.LP01.Powers            = 300 *1e-3;                                                % [mW]
In.Pump.LP01.Alpha             = 0.20; 



% SEÑALES : 
Nch = 100; c = 299792458;
    % LP01  
%In.Signal.LP01.Wavelengths        = linspace(1.66551366e3 , 1.50046275e3 , Nch);                  % [nm]
In.Signal.LP01.Wavelengths        = ( c./linspace(c/(1.66551366e-6) , c/(1.50046275e-6) , Nch)).*1e9;                  % [nm]
In.Signal.LP01.Powers             = -30*ones( 1,length(In.Signal.LP01.Wavelengths) );                 %[dBm]
In.Signal.LP01.Alpha              = 0.2;                                                              % [dB/km]



In.Monomodo = true;
In.Rayleight = false; In.ase_calc = false;

load("RamanVPIE1.mat")

%% Calculo de amplificación
bestGain = 1000000; worstGain = 0;
factor = [1,1]; %best , worst
tic;

Raman = RamanMMv3_VPI(In) ; 
% for fct=1:51
%     clc; fprintf('%.0f /51 \n',fct)
%     In.Fibra.PolarizationFactor       = fct / 50;
%     
%     Raman = RamanMMv3_VPI(In) ; 
%     
%     GainDiff = sum((RamanVPIE1.Gain - Raman.Sig.Gain.LP01).^2);
%     if GainDiff < bestGain
%         bestGain = GainDiff ; factor(1) = fct;
%         BestRaman = Raman;
%     end
%     if GainDiff > worstGain
%         worstGain = GainDiff ; factor(2) = fct;
%         WorstRaman = Raman;
%     end
% end
tend = toc; fprintf("Tiempo de cómputo: %.2fs\n",tend);

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


% for ms = 1:length(Raman.ModoS)% Potencia Señal
%     figure(1); 
%     strlambda = strcat( num2str( In.Signal.(Raman.ModoS{ms}).Wavelengths(ms)) , "nm");
%     subplot(1,length(Raman.ModoS),fs) ; plot(z,10.*log10(signal.(Raman.ModoS{ms})./1e-3) , 'DisplayName',strlambda) ,ylabel("Potencia [dBm]") ; hold on
%     %subplot(1,length(Raman.ModoS),fs) ; plot(z,10.*log10(signal.Off.(Raman.ModoS{ms})./1e-3) ) ,ylabel("Potencia [dBm]") ; hold on
%     %subplot(1,fs,fs) ; plot( z,(signal.(Raman.ModoS{ms}) ) ,ylabel("Potencia [mW]"); hold on
%     xlabel("Largo [km]") ; title(strcat("Signal Modo ", Raman.ModoS{ms})) ; legend('Location','best')
%     fs = fs+1;
% end
% fs = 1;
% for ms = 1:length(Raman.ModoS)% Potencia Señal Off
%     %figure(1); 
%     strlambda = strcat( num2str( In.Signal.(Raman.ModoS{ms}).Wavelengths(ms)) , "nm OFF");
%     h=subplot(1,length(Raman.ModoS),fs) ; h.ColorOrderIndex=1;
%     plot(z,10.*log10(signal.Off.(Raman.ModoS{ms})./1e-3) , 'DisplayName',strlambda) ,ylabel("Potencia [dBm]") ; hold on
%     %subplot(1,fs,fs) ; plot( z,(signal.Off.(Raman.ModoS{ms})) ) ,ylabel("Potencia [mW]"); hold on
%     xlabel("Largo [km]") ; title(strcat("Signal Modo ", Raman.ModoS{ms})) ; legend('Location','best')
%     fs = fs+1;
% end

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

for ms = 1:length(Raman.ModoS)% Ganancias
    figure(3)
    plot(In.Signal.(Raman.ModoS{ms}).Wavelengths , gain.(Raman.ModoS{ms}), '-o')  ; hold on
    %plot(In.Signal.(Raman.ModoS{ms}).Wavelengths(i) , gain.(Raman.ModoS{ms})(i) ) ,ylabel("Potencia [mW]"); hold on
end
    figure(3) ; xlabel("Longitud de Onda nm") ,ylabel("Ganancia [dB]") ; title("Ganancias On-Off") ; legend(Raman.ModoS{:},"Location","best")

% for ms = 1:length(Raman.ModoS)% OSNR
%     figure(4)
%     plot(In.Signal.(Raman.ModoS{ms}).Wavelengths , osnr.(Raman.ModoS{ms}), '-o')  ; hold on
%     %plot(In.Signal.(Raman.ModoS{ms}).Wavelengths(i) , osnr.(Raman.ModoS{ms})(i) ) ,ylabel("Potencia [mW]"); hold on
% end
%     figure(4) ; xlabel("Longitud de Onda nm") ,ylabel("Magnitud [dB]") ; title("OSNR") ; legend(Raman.ModoS{:},"Location","best")


clear fs fp ms mp i strlambda DispName1 DispName2;

% for mp = length(Raman.ModoP)
%    pumpv2.(Raman.ModoP{mp}).forward = Ramanv2.Pump.forward.(Raman.ModoP{mp});
%    pumpv2.(Raman.ModoP{mp}).backward = Ramanv2.Pump.backward.(Raman.ModoP{mp});
% end
% plot(z,pump.(Raman.ModoP{1}).backward(1,:)) ; hold on ; plot(z,pumpv2.(Raman.ModoP{1}).backward(1,:)) 

% %%% Atenuaciones

% alp = load('RADynamic_Attenuation.dat');
% plot(alp(:,1),alp(:,2))


% A                   =load('RamanGain3.dat');
% B                   = load('RamanGains/RamanGain_PSCF.mat'); B = B.RamanGainPSCF;
% C                   = load('RamanGains/RamanGain_DCF.mat'); C = C.RamanGainDCF;
% D                   = load('RamanGains/RamanGainSLA.mat'); D = D.RamanGainSLA;
% plot(A(:,2),"DisplayName","A") ; hold on ; plot(B(:,2),"DisplayName","B") ; plot(C(:,2)./5,"DisplayName","C") ; plot(D(:,2).*1000,"DisplayName","D")  ; legend()


%% Otros Archivos
% eta = load('RADynamic_Rayleigh.dat') ; gr = load('RamanGainEfficiency_SMF28.dat');
% figure(1);plot(eta(:,1),eta(:,2)) ; figure(2) ; plot(gr(:,1),gr(:,2))