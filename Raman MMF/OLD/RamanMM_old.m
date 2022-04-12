% evaluates the Raman gain for given signal structure assuming a pumpt at 1445 nm
function [Raman] = RamanMM(signal,P)
%% Load Raman responce for SMF 28
A = load('RamanResponse.txt');
g_R = smooth(A(:,2),20);

% Max_g_R = P.gR;                          % max Raman gain coeffiencit (G_R/Aeff)
C_R = P.Fibre.C_Rmax*(g_R/max(g_R));
FF_gR = A(:,1)*1e12;
% figure(10)
% plot(A(:,1),C_R)
% xlabel('Frequency difference (THz)')
% ylabel('C_{R} (1/W/km)')


%%
PumpWavelengths    = P.Fibre.RamanWavelengths;     % Pump Wavelengths (nm)
Ppf0                = P.Fibre.RamanPowers;          % Pump forward powers (W)
Ppb0                = P.Fibre.RamanPowers;          % Pump backward powers (W)
Ps00 = 10^( (signal.Ps0-30)/10 ) ;                  % Total Signal input power (W)
PumpAlpha           = P.Fibre.PumpAlpha;

%% functions
Gr_fun = @(f)interp1(FF_gR,C_R,f);

%% Calculate power profile for each channel and pump wavelength
    
    %[Ps00]    = PowerMeter(SignalIn);        % Total Signal input power (W)
    
    Z=(0:0.1:P.Length);                     % zed (km)
    zed = repmat(Z,length(PumpWavelengths),1);
    alpha_Ss   = (P.Att*log(10)/10);                          % fibre loss at the signal frequency (np/km)

                                   
    alpha_Pf   = (P.Fibre.PumpAlpha.*log(10)/10);               % fibre loss at the pump frequency Foward propagation (np/km)
    alpha_Pp_b = alpha_Pf*-1;                                   % fibre loss at the pump frequency Backward propagation np/km)
    
    
    Leff_f     = [1-exp(-alpha_Pf.*Z)]./alpha_Pf;             % Efective length in forward propagation (m)
    Leff_b     = [1-exp(-alpha_Pp_b.*Z)]./alpha_Pf;           % Efective length in Backward propagation (m)
    
    Ayy        = Ppb0.*exp(-alpha_Pf.*P.Length);                % backward pump power at z = 0 km
 

    %%
    Modos = fieldnames(signal.modos);

    for m = 1:length(Modos) % Recorrido por modo 
        modo = Modos{m};
        fprintf("modo: %s \n", modo)
        Nchan = length(signal.modos.(modo)) ;
        for ch = 1:Nchan
            fprintf("lambda: %f nm\n", signal.modos.(modo)(ch))
            F_s = 3e8/( signal.modos.(modo)(ch)*1e-9 );
            Delta_Omega = 3e8./(PumpWavelengths*1e-9)- F_s;
            gR      = Gr_fun(Delta_Omega); 
            
        
        evolution of pumps
        fig = m+1;
        switch P.Fibre.RamanMethod
            case 'Forward'
                % Raman Forward

                for g_index = 1:length(gR)
                    g_index
                    [z,pwr_evolution] = ode45(@(t,PP) evolution_f(t,PP,alpha_Ss,alpha_Pf,gR(g_index),gR(g_index),F_s,PumpWavelengths(g_index)), Z,[Ps00,Ppf0(g_index)] );
                    sig.(modo)(g_index,:) = pwr_evolution(:,1);
                    pump.(modo)(g_index,:) = pwr_evolution(:,2);
                end
%                     figure(fig), hold on, box on, grid on
%                     plot(z,10*log10(sig.(modo)/1e-3), 'DisplayName', int2str(signal.modos.(modo)(ch) ) )
%                     xlabel('distance (km)') ; ylabel('signal power (dBm)') ; title(strcat('Modo ', Modos{m}) )
%                     legend()
%                     
% 
%                     figure(fig+1),hold on, box on, grid on
%                     plot(z,10*log10(pump.(modo)/1e-3),'r')
%                     xlabel('distance (km)') ; ylabel('pump power (dBm)')
%                     title(strcat('Modo ', Modos{m}) )
                    
                
            case 'Backward'
                % Raman Backward
                for g_index = 1:length(gR)
                    [z,pwr_evolution]=ode45(@(t,PP) evolution_b(t,PP,alpha_Ss,alpha_Pf,gR(g_index),gR(g_index),F_s,PumpWavelengths(g_index)) , Z,[Ps00,Ayy(g_index)] );
                    sig.(modo)(g_index,:) = pwr_evolution(:,1);
                    pump.(modo)(g_index,:) = pwr_evolution(:,2);
                end
%                     figure(4), hold on, box on, grid on
%                     plot(z,10*log10(sig/1e-3))
%                     xlabel('distance (km)')
%                     ylabel('signal power (dBm)')
%                     
%                     figure(5), hold on, box on, grid on
%                     plot(z,10*log10(pump/1e-3),'r')
%                     xlabel('distance (km)')
%                     ylabel('pump power (dBm)')
        end % END SwitchCase
          
        end
    end

Raman.z = z';
Raman.Pump = pump;
Raman.Sig = sig ;

function [ Pout ] = evolution_f(t,PP,a,ap,Gr1,Gr2,Fs,Wp)
    % PP(1): Signal ; PP(2): Pump
        Fp = 3e8./(Wp*1e-9);
        Pout=zeros(length(PP),1);
       Pout(1)= -a*PP(1)+Gr1*PP(1).*PP(2);
       Pout(2)= -ap*PP(2)-Fp/Fs*Gr2*PP(1).*PP(2);
end

function [ Pout ] = evolution_b(t,PP,a,ap,Gr1,Gr2,Fs,Wp)
        Fp = 3e8./(Wp*1e-9);
        Pout=zeros(length(PP),1);
       Pout(1)= -a*PP(1)+Gr1*PP(1).*PP(2);
       Pout(2)= +ap*PP(2)-Fp/Fs*Gr2*PP(1).*PP(2);
end

end

function [Potencias] = RamanPot(X,alphaS,alphaP,lambdaS,lambdaP,Gr)
    S = X(1) ; Pf = X(2) ; Pb = X(3);
    Potencias = [-alphaS*S + Gr*(Pf+Pb)*S;
                 -alphaP*Pf -(lambdaS/lambdaP)*Gr*S*Pf;
                 alphaP*Pb + (lambdaS/lambdaP)*Gr*S*Pb]
end
