% evaluates the Raman gain for given signal structure assuming a pumpt at 1445 nm
function [Raman] = RamanMMv3(In)


% % % Load Raman and Rayleigh responce 
% % Implementación profe, g_R normalizado a 0.5
% A                   = load('RamanResponse.txt');
% g_R                 = smooth(A(:,2),20);
% %%%Max_g_R           = P.gR;                                         % max Raman gain coeffiencit (G_R/Aeff)
% C_R                 = In.Fibra.PolarizationFactor*(g_R/max(g_R));
% FF_gR               = A(:,1)*1e12;

% % % Implementación sin normalizar
% A                   = load('RamanGainEfficiency_SMF28.dat');
% g_R                 = smooth(A(:,2),1); %10
% C_R                 = In.Fibra.PolarizationFactor*(g_R)*1000; % Ajuste de polarización y paso a km
% FF_gR               = A(:,1);

A                   =load('RamanGain3.dat');
g_R                 = A(:,2);
C_R                 = In.Fibra.PolarizationFactor*(g_R);
FF_gR               = A(:,1).*1e12;

eta = load('RADynamic_Rayleigh.dat');                               % Rayleigh coeficient
freq = 100 : 10 : 1390;
mag = zeros(1,length(freq));
freq2 = 1650 : 10 : 1800;
mag2 =  zeros(1,length(freq2));
f = [freq eta(:,1)' freq2]*1e-9;
magg = [mag eta(:,2)' mag2];

% % Raman Responce Function
Gr_fun              = @(f)interp1(FF_gR,C_R,f);
eta_fun = @(lambda)interp1(f,magg,lambda);
%figure(1) ; plot( FF_gR,Gr_fun(FF_gR) ) , xlabel("Δ F [Hz]") ; ylabel("Magnitude") ; title("RamanResponce")
%figure(1) ; plot( FF_gR.*1e-12,Gr_fun(FF_gR)/1000 ) , set(gca,'ColorOrderIndex',1,'FontSize',13), xlabel("Corrimiento en Frecuencia [THz]","FontSize",16) ; ylabel("Ganancia Raman (m/W)","FontSize",16) ; title("Coeficiente de ganancia Raman","FontSize",18)
%figure(2) ; plot(eta(:,1),eta(:,2)) , xlabel("λ [nm]") ; ylabel("Magnitude") ; title("Rayleigh Coeficient")


    % Datos de Fibra y Constantes
h = 6.6260*10^(-34) ;                                                       % Constante de Planck
k = 1.380649*10^(-23) ;                                                     % Constante de Boltzmann
T = In.Fibra.T + 273.15 ;                                                   % Temperatura absoluta de la fibra       
c = 299792458;
Gamma = 1/In.Fibra.PolarizationFactor;

ModoS                           = fieldnames(In.Signal);
ModoP                           = fieldnames(In.Pump);
deltaZ                          = 0.04;                                                     % Paso de iteracion en fibra [km]
L                               = In.Fibra.Length;
Z                               = (0 : deltaZ : L);
% --- Bombeos --- %
for i=1:length(ModoP)
  % Frecuencias y Wavelengths
    PumpWavelengths.(ModoP{i})  = In.Pump.(ModoP{i}).Wavelengths;
    lambdaP.(ModoP{i})          = ( PumpWavelengths.(ModoP{i}) ).*1e-9  ;
    F_p.(ModoP{i})              = c./(lambdaP.(ModoP{i}));
  % Potencias
    Ppf0.(ModoP{i})             = In.Pump.(ModoP{i}).Powers;                                % Pump forward powers (W)
    PpbL.(ModoP{i})             = In.Pump.(ModoP{i}).Powers;                                % Pump backward powers (W)
    %Ppb0.(ModoP{i})             = (PpbL.(ModoP{i})) .*exp((-alphaP.(ModoP{i})) .*L);       % backward pump power at z = 0 km
    Pb.(ModoP{i})               = zeros( length(lambdaP.(ModoP{i})), length(Z) );
    Pf.(ModoP{i})               = zeros( length(lambdaP.(ModoP{i})), length(Z) );

    Pf.Rayleigh.(ModoP{i})      = zeros( length(lambdaP.(ModoP{i})), length(Z) );
end                                  

% --- Señales --- %
for i=1:length(ModoS)
  % Frecuencias y Wavelengths
    lambdaS.(ModoS{i})          = (In.Signal.(ModoS{i}).Wavelengths) .*1e-9;
    F_s.(ModoS{i})              = c./( lambdaS.(ModoS{i})); 
  % Potencias
    Ps0.(ModoS{i})              = 10.^( ((In.Signal.(ModoS{i}).Powers) -30)/10 ) ;          % Total Signal input power (W)
    Ps.(ModoS{i})               = zeros( length(lambdaS.(ModoS{i})), length(Z) );
    Ps.(ModoS{i})(:,1)          = Ps0.(ModoS{i}); 
    Psoff.(ModoS{i})            = zeros( length(lambdaS.(ModoS{i})), length(Z) );
    Psoff.(ModoS{i})(:,1)       = Ps0.(ModoS{i}); 

    Ps.Rayleigh.(ModoS{i})      = zeros( length(lambdaS.(ModoS{i})), length(Z) );
    Ps.ASE.(ModoS{i})           = zeros( length(lambdaS.(ModoS{i})), length(Z) );   % ASE sum for each signal chanel
    Pase.(ModoS{i})             = zeros( length(lambdaS.(ModoS{i})), length(Z) );   % ASE propagation vector forward direction
    Pase0.(ModoS{i})            = 10.^( ((In.ASE.(ModoS{i})) -30)/10 ) ;
    Pase.(ModoS{i})(:,1)        = Pase0.(ModoS{i}); 
    
end

switch In.Fibra.AttenuationMethod
    case 'Dynamic'
        alp = load('RADynamic_Attenuation.dat');
        alp1 = alp(:,1).*1e-9 ; alp2 = alp(:,2).* log(10)/10;
        for ms=1:length(ModoS)
            alphaS.(ModoS{ms}) = @(f) interp1(alp1,alp2,f);
        end
        for mp=1:length(ModoP)
            alphaP.(ModoP{mp}) = @(f) interp1(alp1,alp2,f);
        end
    case 'Static'
        %alp = load('RADynamic_Attenuation.dat');
        alp = (1300:1:1700)';
        alp1 = alp(:,1).*1e-9 ;  
        for ms=1:length(ModoS)
            alp2 = ones( 1,length(alp1 ) ).* ( In.Signal.(ModoS{ms}).Alpha .* log(10)/10) ;
            alphaS.(ModoS{ms}) = @(f) interp1(alp1,alp2,f);
        end
        for mp=1:length(ModoP)
            alp2 = ones( 1,length(alp1 ) ).* ( In.Pump.(ModoP{mp}).Alpha.* log(10)/10) ;
            alphaP.(ModoP{mp}) = @(f) interp1(alp1,alp2,f);
        end
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
% Forma Numerica
% norm =  int_overlap_numerical(In.Fibra , '01' , 1450e-9 , '01' , 1550e-9);
% for mp = 1:length(ModoP)
%     for ms = 1:length(ModoS)
%         for p = 1:length( F_p.(ModoP{mp}) )
%             for s = 1:length( F_s.(ModoS{ms}) )
%                 lambdas = c /(F_s.(ModoS{ms})(s)) ; lambdap = c/(F_p.(ModoP{mp})(p)) ;
%                 gR.Pump.(ModoP{mp}).(ModoS{ms})(p,s) = Gr_fun( abs(F_p.(ModoP{mp})(p)- F_s.(ModoS{ms})(s)) ) ;
%                 gR.Signal.(ModoS{ms}).(ModoP{mp})(p,s) = Gr_fun( abs(F_p.(ModoP{mp})(p)-F_s.(ModoS{ms})(s)) ) ;
%                 fmn.Signal.(ModoS{ms}).(ModoP{mp})(p,s) = int_overlap_numerical(In.Fibra , ModoP{mp} , lambdap , ModoS{ms} , lambdas ) / norm ;
%                 fmn.Pump.(ModoP{mp}).(ModoS{ms})(p,s) = int_overlap_numerical(In.Fibra , ModoS{ms} , lambdas , ModoP{mp} , lambdap ) / norm ;
%             end
%         end
%     end
% end

norm =  int_overlapv2(In.Fibra , '01' , 1450e-9 , '01' , 1550e-9); % Superposiciones de modos se normalizan respecto a superposicion de LP01 con LP01
for mp = 1:length(ModoP)
    for ms = 1:length(ModoS)
        for p = 1:length( F_p.(ModoP{mp}) )
            parfor s = 1:length( F_s.(ModoS{ms}) )
                lambdas = c /(F_s.(ModoS{ms})(s)) ; lambdap = c/(F_p.(ModoP{mp})(p)) ;

                gRPumpAux(p,s) = Gr_fun( abs(F_p.(ModoP{mp})(p)- F_s.(ModoS{ms})(s)) ) ;
                if isnan(gRPumpAux(p,s))
                    gRPumpAux(p,s) = 0;
                end
                gRSignalAux(p,s) = Gr_fun( abs(F_p.(ModoP{mp})(p)-F_s.(ModoS{ms})(s)) ) ;
                if isnan(gRSignalAux(p,s))
                    gRSignalAux(p,s) = 0;
                end
                fmnSignalAux(p,s) = int_overlapv2(In.Fibra , ModoP{mp} , lambdap , ModoS{ms} , lambdas ) / norm ;
                fmnPumpAux(p,s) = int_overlapv2(In.Fibra , ModoS{ms} , lambdas , ModoP{mp} , lambdap ) / norm ;
            end
            gR.Pump.(ModoP{mp}).(ModoS{ms})(p,:) = gRPumpAux(p,:);
            gR.Signal.(ModoS{ms}).(ModoP{mp})(p,:) = gRSignalAux(p,:);
            fmn.Signal.(ModoS{ms}).(ModoP{mp})(p,:) = fmnSignalAux(p,:);
            fmn.Pump.(ModoP{mp}).(ModoS{ms})(p,:) = fmnPumpAux(p,:);
        end
    end
end
ASE_BW = 1e12;  %6*10^12; %THz "paper: RAMAN amplifier gain dynamics with ASE"

%% % Calculo potencias
% Inicialización de variables
for i = 1:length(ModoP)
    Pf.(ModoP{i}) = zeros(1, length(Z) );
    Pb.(ModoP{i}) = zeros(1,length(Z) );
end
switch In.Fibra.RamanMethod
    case 'Forward'
        for i = 1:length(ModoP)
            Pf.(ModoP{i})(:,1) = Ppf0.(ModoP{i});
            Pb.(ModoP{i})(:,end) = 0;
            Ppbcalc = 1;
            AdjIter = false;
        end
    case 'Backward'
        for i = 1:length(ModoP)
            Pf.(ModoP{i})(:,1) = 0;
            Pb.(ModoP{i})(:,end) = PpbL.(ModoP{i});
            Ppbcalc = 0;
            AdjIter = true;
        end
    case 'Forward&Backward'
        for i = 1:length(ModoP)
            Pf.(ModoP{i})(:,1) = Ppf0.(ModoP{i});

            Pb.(ModoP{i})(:,end) = PpbL.(ModoP{i});
            Ppbcalc = 1;
            AdjIter = true;
        end
end


% Calculo Psoff - Potencia de señal sin amplificacion

for l = 1:(length(Z)-1)
    for ms = 1:length(ModoS)
        for wS = 1:length(lambdaS.(ModoS{ms}))
            lambda = lambdaS.(ModoS{ms})(wS); alpS = alphaS.(ModoS{ms})(lambda);
            Psoff.(ModoS{ms})(wS,l+1) =  Psoff.(ModoS{ms})(wS,l) + deltaZ*( -alpS*Psoff.(ModoS{ms})(wS,l) ) ;
        end
    end
end

for mp = 1:length(ModoP) % Método 3 - Pb Tomando en cuenta Psoff
    if sum( Pb.(ModoP{mp})(:,end) ) ~= 0
        for l=(length(Z)-1):-1:1
            for wP = 1:length(lambdaP.(ModoP{mp}))
                s_sumOff = 0;
                for ms = 1:length(ModoS)
                    s_sumOff = s_sumOff + sum( fmn.Pump.(ModoP{mp}).(ModoS{ms})(wP,:) * gR.Pump.(ModoP{mp}).(ModoS{ms})(wP,:)'...
                                    .* lambdaS.(ModoS{ms})(:).*Psoff.(ModoS{ms})(:,l) ) ;
                end

                lambda = lambdaP.(ModoP{mp})(wP); alpP = alphaP.(ModoP{mp})(lambda);
                Pb.(ModoP{mp})(wP,l) = Pb.(ModoP{mp})(wP,l+1) + deltaZ*( -alpP*Pb.(ModoP{mp})(wP,l+1) - (1/lambda)*s_sumOff*Pb.(ModoP{mp})(wP,l+1) ) ;

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
                deltaV = abs( c/lambdaP.(ModoP{mp})(:) - c/lambdaS.(ModoS{ms})(wS) ) ; BW = ASE_BW;
                ase_sum = ase_sum + sum(    sum( ( fmn.Signal.(ModoS{ms}).(ModoP{mp})(:,wS) .* gR.Signal.(ModoS{ms}).(ModoP{mp})(:,wS)...
                                    .* ( Pf.(ModoP{mp})(:,l)+Pb.(ModoP{mp})(:,l) ) ) .* (1+(exp(h*deltaV/(k*T))-1).^-1 ) .* BW ) ) ;
            end

            lambda = lambdaS.(ModoS{ms})(wS); alpS = alphaS.(ModoS{ms})(lambda);
            Ps.(ModoS{ms})(wS,l+1) = Ps.(ModoS{ms})(wS,l) + deltaZ*( -alpS*Ps.(ModoS{ms})(wS,l) + p_sum*Ps.(ModoS{ms})(wS,l)  + ...
                                        ase_sum * h  * ( c/lambdaS.(ModoS{ms})(wS) ) + eta_sum);
            Ps.Rayleigh.(ModoS{ms})(wS,l+1) = eta_sum;
            Ps.ASE.(ModoS{ms})(wS,l+1) = ase_sum * h  * ( c/lambdaS.(ModoS{ms})(wS) );

            % ASE PROPAGATION
            Pase.(ModoS{ms})(wS,l+1) = Pase.(ModoS{ms})(wS,l) + deltaZ*( -alpS*Pase.(ModoS{ms})(wS,l) + p_sum*Pase.(ModoS{ms})(wS,l)  + ...
                                        ase_sum * h  * ( c/lambdaS.(ModoS{ms})(wS) ));
            
        end
    end

    % Bombeo Forward
    if Ppbcalc
        for mp = 1:length(ModoP)
            for wP = 1:length(lambdaP.(ModoP{mp}))
                s_sum = 0;
                eta_sum = 0;
                ase_sum = 0;
                for ms = 1:length(ModoS)
                    s_sum = s_sum + sum( fmn.Pump.(ModoP{mp}).(ModoS{ms})(wP,:) * gR.Pump.(ModoP{mp}).(ModoS{ms})(wP,:)'...
                                    .* lambdaS.(ModoS{ms})(:).*Ps.(ModoS{ms})(:,l) ) ;

                    deltaV = abs( c/lambdaP.(ModoP{mp})(wP) - c/lambdaS.(ModoS{ms})(:) ) ; BW = ASE_BW;
                    lam = lambdaS.(ModoS{ms})(:)' ;
                    ase_sum = ase_sum + sum(    sum( ( fmn.Pump.(ModoP{mp}).(ModoS{ms})(:,wP) .* gR.Pump.(ModoP{mp}).(ModoS{ms})(:,wP)...
                                        .* ( Pf.(ModoP{mp})(wP,l) )  )  .* (1+(exp(h*deltaV/(k*T))-1).^-1 ) .* lam .* BW ) ) ;
                end

                lambda = lambdaP.(ModoP{mp})(wP); alpP = alphaP.(ModoP{mp})(lambda);
                
                % Aporte Rayleigh BackScattering , cuando hay bombeo bidireccional
                if Ppbcalc == 1 && AdjIter == true
                    eta_sum = eta_fun(lambda) * ( Pb.(ModoP{mp})(wP,l) );
                end

                Pf.(ModoP{mp})(wP,l+1) = Pf.(ModoP{mp})(wP,l) + deltaZ*( -alpP*Pf.(ModoP{mp})(wP,l) - (1/lambdaP.(ModoP{mp})(wP))*s_sum*Pf.(ModoP{mp})(wP,l) ...
                                        + eta_sum - ase_sum*(1/lambdaP.(ModoP{mp})(wP)) * h  * ( c/lambdaP.(ModoP{mp})(wP) ));  
                Pf.Rayleigh.(ModoP{mp})(wP,l+1) = eta_sum;
                
            end
        end
    end    
end

%% Segunda Iteración - Ajuste de Pbackward

if AdjIter
    Raman.Sig.PowerOneIter = Ps;
    for ms = 1:length(ModoS)
        Gain.(ModoS{ms}) =  10*log10( Ps.(ModoS{ms})(:,end) ./ Ps.(ModoS{ms})(:,1) );
    end
    Raman.Sig.GainOneIter = Gain;
    
    for iter=1:2 % Iteraciones para ajuste de potencia backward
        for mp = 1:length(ModoP) %Pbackward Tomando en cuenta Ps
            if sum( Pb.(ModoP{mp})(:,end) ) ~= 0
                for l=(length(Z)-1):-1:1
                    for wP = 1:length(lambdaP.(ModoP{mp}))
                        s_sumOn = 0;
                        eta_sum = 0;
                        for ms = 1:length(ModoS)
                            s_sumOn = s_sumOn + sum( fmn.Pump.(ModoP{mp}).(ModoS{ms})(wP,:) * gR.Pump.(ModoP{mp}).(ModoS{ms})(wP,:)'...
                                            .* lambdaS.(ModoS{ms})(:).*Ps.(ModoS{ms})(:,l) ) ;
                        end
                        lambda = lambdaP.(ModoP{mp})(wP); alpP = alphaP.(ModoP{mp})(lambda);
                        % Aporte Rayleigh BackScattering , cuando hay bombeo bidireccional
                        if Ppbcalc == 1 && AdjIter == true
                            eta_sum = eta_fun(lambda) * ( Pf.(ModoP{mp})(wP,l) );
                        end
                        Pb.(ModoP{mp})(wP,l) = Pb.(ModoP{mp})(wP,l+1) + deltaZ*( -alpP*Pb.(ModoP{mp})(wP,l+1) - (1/lambda)*s_sumOn*Pb.(ModoP{mp})(wP,l+1) ) ...
                            + eta_sum;
                    end
                end
            end
        end
        
        % Forward y Signal con Backward ajustado
        
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
                        deltaV = abs( c/lambdaP.(ModoP{mp})(:) - c/lambdaS.(ModoS{ms})(wS) ) ; BW = ASE_BW;
                        ase_sum = ase_sum + sum(    sum( ( fmn.Signal.(ModoS{ms}).(ModoP{mp})(:,wS) .* gR.Signal.(ModoS{ms}).(ModoP{mp})(:,wS)...
                                            .* ( Pf.(ModoP{mp})(:,l)+Pb.(ModoP{mp})(:,l) ) ) .* (1+(exp(h*deltaV/(k*T))-1).^-1 ) .* BW ) ) ;
                    end

                    lambda = lambdaS.(ModoS{ms})(wS); alpS = alphaS.(ModoS{ms})(lambda);
                    Ps.(ModoS{ms})(wS,l+1) = Ps.(ModoS{ms})(wS,l) + deltaZ*( -alpS*Ps.(ModoS{ms})(wS,l) + p_sum*Ps.(ModoS{ms})(wS,l)  + ...
                                                ase_sum * h  * ( c/lambdaS.(ModoS{ms})(wS) ) + eta_sum);
                    Ps.Rayleigh.(ModoS{ms})(wS,l+1) = eta_sum;
                    Ps.ASE.(ModoS{ms})(wS,l+1) = ase_sum * h  * ( c/lambdaS.(ModoS{ms})(wS) );

                    % ASE PROPAGATION
                    Pase.(ModoS{ms})(wS,l+1) = Pase.(ModoS{ms})(wS,l) + deltaZ*( -alpS*Pase.(ModoS{ms})(wS,l) + p_sum*Pase.(ModoS{ms})(wS,l)  + ...
                                        ase_sum * h  * ( c/lambdaS.(ModoS{ms})(wS) ));
                end
            end
        
            % Bombeo Forward
            if Ppbcalc
                for mp = 1:length(ModoP)
                    for wP = 1:length(lambdaP.(ModoP{mp}))
                        s_sum = 0;
                        eta_sum = 0;
                        ase_sum = 0;
                        for ms = 1:length(ModoS)
                            s_sum = s_sum + sum( fmn.Pump.(ModoP{mp}).(ModoS{ms})(wP,:) * gR.Pump.(ModoP{mp}).(ModoS{ms})(wP,:)'...
                                            .* lambdaS.(ModoS{ms})(:).*Ps.(ModoS{ms})(:,l) ) ;
                            
                            deltaV = abs( c/lambdaP.(ModoP{mp})(wP) - c/lambdaS.(ModoS{ms})(:) ) ; BW = ASE_BW;
                            lam = lambdaS.(ModoS{ms})(:)';
                            ase_sum = ase_sum + sum(    sum( ( fmn.Pump.(ModoP{mp}).(ModoS{ms})(:,wP) .* gR.Pump.(ModoP{mp}).(ModoS{ms})(:,wP)...
                                            .* ( Pf.(ModoP{mp})(wP,l) ) ) .* lam .* (1+(exp(h*deltaV/(k*T))-1).^-1 ) .* BW ) ) ;

                        end
            
                        lambda = lambdaP.(ModoP{mp})(wP); alpP = alphaP.(ModoP{mp})(lambda);
                        % Aporte Rayleigh BackScattering , cuando hay bombeo bidireccional
                        if Ppbcalc == 1 && AdjIter == true
                            eta_sum = eta_fun(lambda) * ( Pb.(ModoP{mp})(wP,l) );
                        end
                        

                        Pf.(ModoP{mp})(wP,l+1) = Pf.(ModoP{mp})(wP,l) + deltaZ*( -alpP*Pf.(ModoP{mp})(wP,l) - (1/lambdaP.(ModoP{mp})(wP))*s_sum*Pf.(ModoP{mp})(wP,l) ...
                                                + eta_sum - ase_sum*(1/lambdaP.(ModoP{mp})(wP)) * h  * ( c/lambdaP.(ModoP{mp})(wP) )); 
                        Pf.Rayleigh.(ModoP{mp})(wP,l+1) = eta_sum;
                    end
                end
            end

        end
    end
end


%% Ganancias
for ms = 1:length(ModoS)
    Gain.(ModoS{ms}) =  10*log10( Ps.(ModoS{ms})(:,end) ./ Ps.(ModoS{ms})(:,1) );
    GainOnOFF.(ModoS{ms}) =  10*log10( Ps.(ModoS{ms})(:,end) ./ Psoff.(ModoS{ms})(:,end) );
    %TempOSNR.(ModoS{ms}) = 10*log10( Ps.(ModoS{ms})(:,:)./Ps.ASE.(ModoS{ms})(:,:) );
    TempOSNR.(ModoS{ms}) = 10*log10( Ps.(ModoS{ms})(:,:)./Pase.(ModoS{ms})(:,:) );
    %OSNR.(ModoS{ms}) = 10*log10( Ps.(ModoS{ms})(:,end)./Ps.ASE.(ModoS{ms})(:,end) );
    OSNR.(ModoS{ms}) = 10*log10( Ps.(ModoS{ms})(:,end)./Pase.(ModoS{ms})(:,end) );
    NF.(ModoS{ms}) = TempOSNR.(ModoS{ms})(:,2)-TempOSNR.(ModoS{ms})(:,end);
end

Raman.z = Z;
Raman.Pump.forward = Pf;            
Raman.Pump.backward = Pb;
Raman.Sig.Power = Ps;               
Raman.Sig.Gain = Gain;
Raman.Sig.GainOnOFF = GainOnOFF;
Raman.Sig.Power.Off = Psoff;
Raman.ModoS = ModoS; Raman.ModoP = ModoP;
Raman.fmn = fmn;
Raman.gR = gR;
Raman.OSNR = OSNR;
Raman.OSNR_Z = TempOSNR;
Raman.NF = NF;
Raman.functions.gr = Gr_fun ;
Raman.functions.eta = eta_fun ;
Raman.LambdasS = lambdaS ;
Raman.LambdasP = lambdaP ;
Raman.Atenuation.Signal = alphaS;
Raman.Atenuation.Pump = alphaP;
Raman.ASE = Pase;

end




