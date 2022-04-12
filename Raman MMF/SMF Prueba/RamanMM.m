% evaluates the Raman gain for given signal structure assuming a pumpt at 1445 nm
function [Raman] = RamanMM(In,f)


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
ModoS                           = fieldnames(In.Signal);
ModoP                           = fieldnames(In.Pump);
deltaZ                          = 0.01;                                                     % Paso de iteracion en fibra [km]
L                               = In.Fibra.Length;
Z                               = (0 : deltaZ : L);
% --- Bombeos --- %
for i=1:length(ModoP)
  % Alphas
    %alphaP.(ModoP{i})         = In.Pump.(ModoP{i}).Alpha;                                 % fibre loss at the Pump frequency (np/km)
    alphaP.(ModoP{i})         = In.Pump.(ModoP{i}).Alpha .* log(10)/10;                    % fibre loss at the Pump frequency (np/km)
    alphaPb.(ModoP{i})         = alphaP.(ModoP{i})*-1;                                   % fibre loss at the pump frequency Backward propagation np/km)
  % Frecuencias y Wavelengths
    PumpWavelengths.(ModoP{i})  = In.Pump.(ModoP{i}).Wavelengths;
    lambdaP.(ModoP{i})          = ( PumpWavelengths.(ModoP{i}) ) .*1e-9 ;
    F_p.(ModoP{i})              = 3e8./(lambdaP.(ModoP{i}));
  % Potencias
    Ppf0.(ModoP{i})             = In.Pump.(ModoP{i}).Powers;                                % Pump forward powers (W)
    PpbL.(ModoP{i})             = In.Pump.(ModoP{i}).Powers;                                % Pump backward powers (W)
    Ppb0.(ModoP{i})             = (PpbL.(ModoP{i})) .*exp((-alphaP.(ModoP{i})) .*L);       % backward pump power at z = 0 km
    Pb.(ModoP{i})               = zeros( length(lambdaP.(ModoP{i})), length(Z) );
    Pf.(ModoP{i})               = zeros( length(lambdaP.(ModoP{i})), length(Z) );
end                                  

% --- Señales --- %
for i=1:length(ModoS)
  % Alphas
    alphaS.(ModoS{i})           = (In.Signal.(ModoS{i}).Alpha);                             % fibre loss at the signal frequency (np/km)
    %alphaS.(ModoS{i})           = In.Signal.(ModoS{i}).Alpha .* log(10)/10;                             % fibre loss at the signal frequency (np/km)
  % Frecuencias y Wavelengths
    lambdaS.(ModoS{i})          = (In.Signal.(ModoS{i}).Wavelengths) .*1e-9;
    F_s.(ModoS{i})              = 3e8./( lambdaS.(ModoS{i})); 
  % Potencias
    Ps0.(ModoS{i})              = 10.^( ((In.Signal.(ModoS{i}).Powers) -30)/10 ) ;          % Total Signal input power (W)
    Ps.(ModoS{i}) = zeros( length(lambdaS.(ModoS{i})), length(Z) );
    Ps.(ModoS{i})(:,1) = Ps0.(ModoS{i}); 
end


% Estructura de gR:
%   gR.Pump.j.i -> guarda valores gR del modo j de bombeo con modo i de señal 
%               -> tiene length(F_p(j)) filas x length(F_s(i)) columnas
%   
%   gR.Signal.a.b -> guarda valores gR del modo a de señal con modo b de bombeo 
%                 -> tiene length(F_s(a)) filas x length(F_p(b)) columnas     
%   
for mp = 1:length(ModoP)
    for ms = 1:length(ModoS)
        for p = 1:length( F_p.(ModoP{mp}) )
            for s = 1:length( F_s.(ModoS{ms}) )
                gR.Pump.(ModoP{mp}).(ModoS{ms})(p,s) = Gr_fun( abs(F_p.(ModoP{mp})(p)- F_s.(ModoS{ms})(s)) ) ;
                gR.Signal.(ModoS{ms}).(ModoP{mp})(p,s) = Gr_fun( abs(F_p.(ModoP{mp})(p)-F_s.(ModoS{ms})(s)) ) ;
            end
        end
    end
end

fmn = 1;


% % Calculo potencias

switch In.Fibra.RamanMethod
    case 'Forward'
        for i = 1:length(lambdaP)
            Pf.(ModoP{i})(:,1) = Ppf0.(ModoP{i});
            Pb.(ModoP{i})(:,end) = 0;
        end
    case 'Backward'
        for i = 1:length(lambdaP)
            Pf.(ModoP{i})(:,1) = 0;
            Pb.(ModoP{i})(:,end) = PpbL.(ModoP{i});
        end
    case 'Forward&Backward'
        for i = 1:length(lambdaP)
            Pf.(ModoP{i})(:,1) = Ppf0.(ModoP{i});
            Pb.(ModoP{i})(:,end) = PpbL.(ModoP{i});
        end
end


% Backward pump se calcula en ausencia de señales
for i = length(ModoP)
    if sum( Pb.(ModoP{i})(:,end) ) ~= 0
        for l=1:(length(Z)-1)
            z = Z(l);
            for wP = 1:length(lambdaP.(ModoP{i}))
                Pb.(ModoP{i})(wP,end-l) = Pb.(ModoP{i})(wP,end)*(exp(-alphaP.(ModoP{i})*z) );
            end
        end
    end
end
%% Forward y Signal se calculan tomando en cuenta Backward

for l = 1:(length(Z)-1)
    
    for ms = 1:length(ModoS)
        p_sum = 0;
        for wS = 1:length(lambdaS.(ModoS{ms}))
            for mp = 1:length(ModoP)
                p_sum = p_sum + sum( fmn .* gR.Signal.(ModoS{ms}).(ModoP{mp})(:,wS) .* ( Pf.(ModoP{mp})(:,l)+Pb.(ModoP{mp})(:,l) ) );
            end
            Ps.(ModoS{ms})(wS,l+1) = Ps.(ModoS{ms})(wS,l) + deltaZ*( -alphaS.(ModoS{ms})*Ps.(ModoS{ms})(wS,l) + p_sum*Ps.(ModoS{ms})(wS,l) );
        end
    end
    for mp = 1:length(ModoP)
        s_sum = 0;
        for wP = 1:length(lambdaP.(ModoP{mp}))
            for ms = 1:length(ModoS)
                s_sum = sum( fmn * gR.Pump.(ModoP{mp}).(ModoS{ms})(wP,:)' .* lambdaS.(ModoS{ms})(:).*Ps.(ModoS{ms})(:,l) ) ;
            end
            Pf.(ModoP{mp})(wP,l+1) = Pf.(ModoP{mp})(wP,l) + deltaZ*( -alphaP.(ModoP{mp})*Pf.(ModoP{mp})(wP,l) - (1/lambdaP.(ModoP{mp})(wP))*s_sum*Pf.(ModoP{mp})(wP,l) ) ; 
        end
    end
    
end

 
Raman.z = Z;
Raman.Pump.forward = Pf;   %pump;
Raman.Pump.backward = Pb;
Raman.Sig = Ps;             %sig ;
Raman.ModoS = ModoS; Raman.ModoP = ModoP;

end



