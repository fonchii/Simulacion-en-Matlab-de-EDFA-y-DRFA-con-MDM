% evaluates the Raman gain for given signal structure assuming a pumpt at 1445 nm
function [Raman] = RamanSMFWDM(In)

% % Load Raman responce for SMF 28
A                   = load('RamanResponse.txt');
g_R                 = smooth(A(:,2),20);
% Max_g_R           = P.gR;                                         % max Raman gain coeffiencit (G_R/Aeff)
C_R                 = In.Fibra.C_Rmax*(g_R/max(g_R));
FF_gR               = A(:,1)*1e12;

% % Function
Gr_fun              = @(f)interp1(FF_gR,C_R,f);

% % Calculate power profile for each channel and pump wavelength
    % Datos de Fibra
deltaZ              = 0.01;                                         % Paso de iteracion en fibra [km]
L                   = In.Fibra.Length;
Z                   = (0 : deltaZ : L);

alpha_Pf            = (In.Pump.Alpha.*log(10)/10) ;                 % fibre loss at the pump frequency Foward propagation (np/km)
alpha_Pp_b          = alpha_Pf*-1;                                  % fibre loss at the pump frequency Backward propagation np/km)
alphaP              = alpha_Pf;
alpha_Ss            = (In.Signal.Alpha*log(10)/10);                 % fibre loss at the signal frequency (np/km)
alphaS              = alpha_Ss;

Leff_P              = [1-exp(-alphaP.*Z)]./alphaP;

    % Wavelenghts & Powers

Ppf0                = In.Pump.Powers;                               % Pump forward powers (W)
PpbL                = In.Pump.Powers;                               % Pump backward powers (W)
Ppb0                = PpbL.*exp(-alpha_Pf.*L);                      % backward pump power at z = 0 km
Ps0                 = 10.^( (In.Signal.Powers-30)/10 ) ;            % Total Signal input power (W)
PumpWavelengths     = In.Pump.Wavelengths;                          % Pump Wavelengths (nm)
lambdaP             = PumpWavelengths ;

ModoS               = fieldnames(In.Signal.Wavelenghts);
modo                = ModoS{1};
lambdaS             = In.Signal.Wavelenghts.(modo).*1e-9;


F_s                 = 3e8./( lambdaS);  
F_p                 = 3e8./(lambdaP*1e-9);
for p = 1:length(lambdaP)
    for s = 1:length(lambdaS)
        gR(p,s) = Gr_fun( abs(F_p(p)-F_s(s)) ) ;
    end
end
fmn = 1;


% % Calculo potencias

Pb = zeros(length(lambdaP),length(Z));
Pf = zeros(length(lambdaP),length(Z));
Ps = zeros(length(lambdaS),length(Z));
Ps(:,1) = Ps0; 
switch In.Fibra.RamanMethod
    case 'Forward'
        Pf(:,1) = Ppf0;
        Pb(:,end) = 0;
    case 'Backward'
        Pf(:,1) = 0;
        Pb(:,end) = PpbL;
    case 'Forward&Backward'
        Pf(:,1) = Ppf0;
        Pb(:,end) = PpbL;
end


% Backward pump se calcula en ausencia de señales
if sum( Pb(:,end) ) ~= 0
    for l=1:(length(Z)-1)
        z = Z(l);
        for wP = 1:length(lambdaP)
            Pb(wP,end-l) = Pb(wP,end)*(exp(-alphaP*z) );
        end
    end
end
% Forward y Signal se calculan tomando en cuenta Backward
for l = 1:(length(Z)-1)
    for wS = 1:length(lambdaS)
        p_sum = sum( fmn.*gR(:,wS).*( Pf(:,l)+Pb(:,l) ) );
        Ps(wS,l+1) = Ps(wS,l) + deltaZ*( -alphaS*Ps(wS,l) + p_sum*Ps(wS,l) );
    end
    for wP = 1:length(lambdaP)
        s_sum = sum( fmn*gR(wP,:)'.*lambdaS(:).*Ps(:,l) ) ;
        Pf(wP,l+1) = Pf(wP,l) + deltaZ*( -alphaP*Pf(wP,l) - (1/lambdaP(wP))*s_sum*Pf(wP,l) ) ; 
    end
    
end

 
Raman.z = Z;
Raman.Pump.forward = Pf;   %pump;
Raman.Pump.backward = Pb;
Raman.Sig = Ps;             %sig ;

end




