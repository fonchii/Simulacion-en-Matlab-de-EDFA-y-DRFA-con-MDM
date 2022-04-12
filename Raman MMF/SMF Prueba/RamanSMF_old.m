% evaluates the Raman gain for given signal structure assuming a pumpt at 1445 nm
function [Raman] = RamanSMF_old(signal,P) %with ODE
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
PpbL                = P.Fibre.RamanPowers;          % Pump backward powers (W)
Ps0 = 10^( (signal.Ps0-30)/10 ) ;                  % Total Signal input power (W)
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
    
    Ppb0        = PpbL*exp(-alpha_Pf.*P.Length);                % backward pump power at z = 0 km
 
    
    %%
    Modos = fieldnames(signal.modos);
    
%    for m = 1:length(Modos) % Recorrido por modo 
%        for ch = 1:Nchan
    modo = Modos{1};
    F_s = 3e8/( signal.modos.(modo)(1)*1e-9 );
    Delta_Omega = 3e8/(PumpWavelengths*1e-9)- F_s;
    gR      = Gr_fun(Delta_Omega);
    alphaS = alpha_Ss; alphaP = alpha_Pf ; lambdaP = PumpWavelengths ; lambdaS=signal.modos.(modo)(1)*1e-9;
    Z=(0:0.1:P.Length);
    [z,Pout] = ode45(@(L,X) ( RamanPot(L,X,alphaS,alphaP,lambdaS,lambdaP,gR) ) , Z,[Ps0,Ppb0,Ppf0]  ) ;
    %[z,Pout] = ode45(@(L,X) ( RamanPot(L,X,alphaS,alphaP,lambdaS,lambdaP,gR) ) , Z,[Ps0,Ppb0]  ) ;            
        
    Raman.z = z';
    Raman.Pump.forward = Pout(:,2);   %pump;
    Raman.Pump.backward = Pout(:,3);
    Raman.Sig = Pout(:,1);             %sig ;

end

function [Potencias] = RamanPot(L,X,alphaS,alphaP,lambdaS,lambdaP,Gr)
    S = X(1) ; Pf = X(2) ; Pb = X(3);
    Potencias = zeros(length(X),1);
    Potencias(1) = -alphaS*S + Gr.*(Pf+Pb)*S;
    %Potencias(1) = -alphaS*S + Gr.*(Pf)*S;
    Potencias(2) = -alphaP*Pf -(lambdaS/lambdaP)*Gr*S*Pf;
    Potencias(3) =  alphaP*Pb - (lambdaS/lambdaP)*Gr*S*Pb;

end

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




