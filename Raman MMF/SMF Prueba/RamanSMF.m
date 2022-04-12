% evaluates the Raman gain for given signal structure assuming a pumpt at 1445 nm
function [Raman] = RamanSMF(signal,P)
%% Load Raman responce for SMF 28
A = load('RamanResponse.txt');
g_R = smooth(A(:,2),20);

% Max_g_R = P.gR;                          % max Raman gain coeffiencit (G_R/Aeff)
C_R = P.Fibre.C_Rmax*(g_R/max(g_R));
FF_gR = A(:,1)*1e12;

% % functions
Gr_fun = @(f)interp1(FF_gR,C_R,f);

%% Calculate power profile for each channel and pump wavelength

alpha_Pf   = (P.Fibre.PumpAlpha.*log(10)/10) ;             % fibre loss at the pump frequency Foward propagation (np/km)
%alpha_Pf = P.Fibre.PumpAlpha/8.686 ;
alpha_Pp_b = alpha_Pf*-1;                                   % fibre loss at the pump frequency Backward propagation np/km)
alpha_Ss   = (P.Att*log(10)/10);                            % fibre loss at the signal frequency (np/km)
%alpha_Ss   = P.Att / 8.686;
PumpWavelengths    = P.Fibre.RamanWavelengths;     % Pump Wavelengths (nm)
Ppf0                = P.Fibre.RamanPowers;          % Pump forward powers (W)
PpbL                = P.Fibre.RamanPowers;          % Pump backward powers (W)
Ppb0        = PpbL*exp(-alpha_Pf.*P.Length);                % backward pump power at z = 0 km

Ps0 = 10^( (signal.Ps0-30)/10 ) ;                  % Total Signal input power (W)

Z=(0:0.1:P.Length);                     % zed (km)
%zed = repmat(Z,length(PumpWavelengths),1);
                          





%Leff_f     = [1-exp(-alpha_Pf.*Z)]./alpha_Pf;             % Efective length in forward propagation (m)
%Leff_b     = [1-exp(-alpha_Pp_b.*Z)]./alpha_Pf;           % Efective length in Backward propagation (m)
%Leff_f     = [1-exp(-alpha_P.*P.Length)]./alpha_Pf;


deltaZ = 0.1;
Z=(0:0.1:P.Length);
alphaS = alpha_Ss; alphaP = alpha_Pf ; lambdaP = PumpWavelengths ; 


Modos = fieldnames(signal.modos);
modo = Modos{1};
lambdaS=signal.modos.(modo)(1)*1e-9;
F_s = 3e8/( signal.modos.(modo)(1)*1e-9 );
Delta_Omega = 3e8/(PumpWavelengths*1e-9)- F_s;
gR      = Gr_fun(Delta_Omega);


% Calculo potencia forward/backward en ausencia de otros bombeo/señales
L = P.Length;
Leff_P     = [1-exp(-alphaP.*Z)]./alphaP;
Pb = zeros(length(Z),1);
Pf = zeros(length(Z),1);
Ps = zeros(length(Z),1);
Ps(1) = Ps0; Pf(1) = Ppf0; Pb(end) = PpbL;

Af = (gR/alphaP)*Ppf0 ;
Ab = (gR/alphaP)*PpbL ;
% Backward pump se calcula en ausencia de señales
for l=1:(length(Z)-1)
    z = Z(l);
    Pb(end-l) = Pb(end)*(exp(-alphaP*z) );
end
for l=1:(length(Z)-1)
    %Pf(l+1) = Pf(1)*(exp(-alphaP*z) );
    Pf(l+1) = Pf(l) + deltaZ*( -alphaP*Pf(l) - (lambdaS/lambdaP)*gR*Ps(l)*Pf(l) ) ; 
    Ps(l+1) = Ps(l) + deltaZ*( -alphaS*Ps(l) + gR*( Pf(l)+Pb(l) )*Ps(l) );
        %-alphaS*z + Af*( 1-exp(1-exp(-alphaP*z)) ) + Ab*exp(-alphaP*z)*( exp(alphaP*z)-1 ) 
    %if l>1
        %   undepleted pump approximation
        %Ps(l) = Ps(1)*exp(-alphaS*z + Af*( 1-exp(-alphaP*z) ) + Ab*exp(-alphaP*L)*( exp(alphaP*z)-1 ) ) ;
        %   EQ CFO amplificacion raman
        %Ps(l) = Ps(1)*exp( gR*Ppf0*Leff_P(l)/alphaP - alphaS*L ) ;
    %end
end

Pout(:,1) = Ps';
Pout(:,2) = Pf';
Pout(:,3) = Pb';

    
%    for m = 1:length(Modos) % Recorrido por modo 
%        for ch = 1:Nchan
    
      
        
    Raman.z = Z';
    Raman.Pump.forward = Pout(:,2);   %pump;
    Raman.Pump.backward = Pout(:,3);
    Raman.Sig = Pout(:,1);             %sig ;

end




