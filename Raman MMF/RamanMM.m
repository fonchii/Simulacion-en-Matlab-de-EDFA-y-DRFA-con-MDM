% evaluates the Raman gain for given signal structure assuming a pumpt at 1445 nm
function [Raman] = RamanMM(In)


% % Load Raman and Rayleigh responce 
A                   = load('RamanResponse.txt');
g_R                 = smooth(A(:,2),20);
% Max_g_R           = P.gR;                                         % max Raman gain coeffiencit (G_R/Aeff)
C_R                 = In.Fibra.PolarizationFactor*(g_R/max(g_R));
FF_gR               = A(:,1)*1e12;
eta = load('RADynamic_Rayleigh.dat');                               % Rayleigh coeficient
freq = [100 : 10 : 1390];
mag = zeros(1,length(freq));
freq2 = [1650 : 10 : 1800];
mag2 =  zeros(1,length(freq2));
f = [freq eta(:,1)' freq2]*1e-9;
magg = [mag eta(:,2)' mag2];

% % Raman Responce Function
Gr_fun              = @(f)interp1(FF_gR,C_R,f);
eta_fun = @(lambda)interp1(f,magg,lambda);
%figure(3) ; plot( FF_gR,Gr_fun(FF_gR) ) , xlabel("Δ F [Hz]") ; ylabel("Magnitude") ; title("RamanResponce")
%figure(3) ; plot(eta(:,1),eta(:,2)) , xlabel("λ [nm]") ; ylabel("Magnitude") ; title("Rayleigh Coeficient")


    % Datos de Fibra y Constantes
h = 6.6260*10^(-34) ;           % Constante de Planck
k = 1.380649*10^(-23) ;         % Constante de Boltzmann
T = In.Fibra.T + 273.15 ;       % Temperatura absoluta de la fibra       
Gamma = 1/In.Fibra.PolarizationFactor;

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

    Pf.Rayleigh.(ModoP{i})      = zeros( length(lambdaP.(ModoP{i})), length(Z) );
end                                  

% --- Señales --- %
for i=1:length(ModoS)
  % Alphas
    %alphaS.(ModoS{i})           = (In.Signal.(ModoS{i}).Alpha);                             % fibre loss at the signal frequency (np/km)
    alphaS.(ModoS{i})           = In.Signal.(ModoS{i}).Alpha .* log(10)/10;                             % fibre loss at the signal frequency (np/km)
  % Frecuencias y Wavelengths
    lambdaS.(ModoS{i})          = (In.Signal.(ModoS{i}).Wavelengths) .*1e-9;
    F_s.(ModoS{i})              = 3e8./( lambdaS.(ModoS{i})); 
  % Potencias
    Ps0.(ModoS{i})              = 10.^( ((In.Signal.(ModoS{i}).Powers) -30)/10 ) ;          % Total Signal input power (W)
    Ps.(ModoS{i}) = zeros( length(lambdaS.(ModoS{i})), length(Z) );
    Ps.(ModoS{i})(:,1) = Ps0.(ModoS{i}); 

    Ps.Rayleigh.(ModoS{i})      = zeros( length(lambdaS.(ModoS{i})), length(Z) );
    Ps.ASE.(ModoS{i})           = zeros( length(lambdaS.(ModoS{i})), length(Z) );
end


% Estructura de gR:
%   Tanto para gr.pump como para gr.signal se guardan en filas los modos de bombeo
%   y en columnas los modos de señal 
%
%   gR.Pump.j.i -> guarda valores gR del modo j de bombeo con modo i de señal 
%               -> tiene length(F_p(j)) filas x length(F_s(i)) columnas
%                   λp1 ...
%                   .
%                   .
%                   .
%                   λpm
%   gR.Signal.a.b -> guarda valores gR del modo a de señal con modo b de bombeo 
%                 -> tiene length(F_s(b)) filas x length(F_p(a)) columnas
%                   λs1 ... λsn

for mp = 1:length(ModoP)
    for ms = 1:length(ModoS)
        for p = 1:length( F_p.(ModoP{mp}) )
            for s = 1:length( F_s.(ModoS{ms}) )
                lambdas = 3e8 /(F_s.(ModoS{ms})(s)) ; lambdap = 3e8/(F_p.(ModoP{mp})(p)) ;
                gR.Pump.(ModoP{mp}).(ModoS{ms})(p,s) = Gr_fun( abs(F_p.(ModoP{mp})(p)- F_s.(ModoS{ms})(s)) ) ;
                gR.Signal.(ModoS{ms}).(ModoP{mp})(p,s) = Gr_fun( abs(F_p.(ModoP{mp})(p)-F_s.(ModoS{ms})(s)) ) ;
                fmn.Signal.(ModoS{ms}).(ModoP{mp})(p,s) = int_overlap_numerical(In.Fibra , ModoP{mp} , lambdap , ModoS{ms} , lambdas ) ;
                fmn.Pump.(ModoP{mp}).(ModoS{ms})(p,s) = int_overlap_numerical(In.Fibra , ModoS{ms} , lambdas , ModoP{mp} , lambdap ) ;
            end
        end
    end
end

%fmn = 1;

%% % Calculo potencias

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
    % Señal
    for ms = 1:length(ModoS)
        for wS = 1:length(lambdaS.(ModoS{ms}))
            p_sum = 0;
            eta_sum = 0;
            ase_sum = 0;
            for mp = 1:length(ModoP)
                p_sum = p_sum + sum( fmn.Signal.(ModoS{ms}).(ModoP{mp})(:,wS) .* gR.Signal.(ModoS{ms}).(ModoP{mp})(:,wS)...
                                    .* ( Pf.(ModoP{mp})(:,l)+Pb.(ModoP{mp})(:,l) ) );
                % Contribución ASE
                deltaV = abs( 3e8/lambdaP.(ModoP{mp})(:) - 3e8/lambdaS.(ModoS{ms})(wS) ) ;
                ase_sum = ase_sum + sum( ( fmn.Signal.(ModoS{ms}).(ModoP{mp})(:,wS) .* gR.Signal.(ModoS{ms}).(ModoP{mp})(:,wS)...
                                    .* ( Pf.(ModoP{mp})(:,l)+Pb.(ModoP{mp})(:,l) ) ) .* deltaV ) ;
            end
            % Aporte Rayleigh Scattering
                % Revisar si algún canal Pump en otro modo tiene igual frecuencia
            for mms = 1:length(ModoS)
                if string(ModoS{ms}) ~= string(ModoS{mms}) 
                    for wwS = 1:length(lambdaS.(ModoS{mms}))
                        if lambdaS.(ModoS{mms})(wwS) == lambdaS.(ModoS{ms})(wS)
                            eta_sum = eta_sum + eta_fun(lambdaS.(ModoS{mms})(wwS)) * Ps.(ModoS{mms})(wwS,l);
                        end
                    end
                end
            end

            Ps.(ModoS{ms})(wS,l+1) = Ps.(ModoS{ms})(wS,l) + deltaZ*( -alphaS.(ModoS{ms})*Ps.(ModoS{ms})(wS,l) + p_sum*Ps.(ModoS{ms})(wS,l) ) + ...
                                        ase_sum * h * Gamma * ( 3e8/lambdaS.(ModoS{ms})(wS) ) + eta_sum;
            Ps.Rayleigh.(ModoS{ms})(wS,l+1) = eta_sum;
            Ps.ASE.(ModoS{ms})(wS,l+1) = ase_sum * h * Gamma * ( 3e8/lambdaS.(ModoS{ms})(wS) );
        end
    end

    % Bombeo Forward
    for mp = 1:length(ModoP)
        for wP = 1:length(lambdaP.(ModoP{mp}))
            s_sum = 0;
            eta_sum = 0;
            for ms = 1:length(ModoS)
                s_sum = sum( fmn.Pump.(ModoP{mp}).(ModoS{ms})(wP,:) * gR.Pump.(ModoP{mp}).(ModoS{ms})(wP,:)'...
                                .* lambdaS.(ModoS{ms})(:).*Ps.(ModoS{ms})(:,l) ) ;
            end

            % Aporte Rayleigh Scattering
                % Revisar si algún canal Pump en otro modo tiene igual frecuencia
            for mmp = 1:length(ModoP)
                if string(ModoP{mp}) ~= string(ModoP{mmp}) 
                    for wwP = 1:length(lambdaP.(ModoP{mmp}))
                        if lambdaP.(ModoP{mmp})(wwP) == lambdaP.(ModoP{mp})(wP)
                            eta_sum = eta_sum + eta_fun(lambdaP.(ModoP{mp})(wP)) * ( Pf.(ModoP{mmp})(wwP,l) + Pb.(ModoP{mmp})(wwP,l) );
                        end
                    end
                end
            end

            Pf.(ModoP{mp})(wP,l+1) = Pf.(ModoP{mp})(wP,l) + deltaZ*( -alphaP.(ModoP{mp})*Pf.(ModoP{mp})(wP,l) - (1/lambdaP.(ModoP{mp})(wP))*s_sum*Pf.(ModoP{mp})(wP,l) )...
                                    + eta_sum; 
            Pf.Rayleigh.(ModoP{mp})(wP,l+1) = eta_sum;
        end
    end
    
end
for ms = 1:length(ModoS)
    Gain.(ModoS{ms}) =  10*log10( Ps.(ModoS{ms})(:,end) ./ Ps.(ModoS{ms})(:,1) );
end
 
Raman.z = Z;
Raman.Pump.forward = Pf;            %pump;
Raman.Pump.backward = Pb;
Raman.Sig.Power = Ps;               %sig ;
Raman.Sig.Gain = Gain ;
Raman.ModoS = ModoS; Raman.ModoP = ModoP;
Raman.fmn = fmn;
Raman.gR = gR;

end




