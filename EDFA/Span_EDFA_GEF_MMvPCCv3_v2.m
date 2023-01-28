function edfa = Span_EDFA_GEF_MMvPCCv3_v2(Fibra,Signal,Pump,ASE)
% Datos de entrada:
% signal (struct)
%       signal.lambda -> int [nm]
%       signal.P0 -> struct con ( int [mW] ) por modo y frecuencia
%       signal.modos -> str
% pump (struct)
%       pump.lambda -> int [nm]
%       pump.P0 -> struct con ( int [mW] ) por modo y frecuencia
%       pump.modos -> str
% ASE -> potencia -> int [mW]
% largo -> int [m]
% canales espaciales -> int [numero canales espaciales]
% fibra (struct)
%       fibra.largo -> int
%       fibra.n1/n2 -> int (indices de refraccion)
%       fibra.N -> int ( Concentracion total de iones de Erbio)
%

% % Parametros de la fibra y constantes
c=299.792458e6;                         % Velocidad de la luz en el vacio m/s
h = 6.626*10^(-34);                     % constante de planck
L = Fibra.largo ;                       % largo
del_z = 1;                              % Tamaño de paso de iteraciones
Z = linspace(0,L,L*(1/del_z));          % se toman segmentos de del_z m
if length(Z)<20                         % Correccion para EDFA muy corto
    Z=linspace(0,L,20);
    del_z = (Z(2)-Z(1));
end
if L>20                                 % Correccion para EDFAs Largos
    Z=linspace(0,L,50);
    del_z = (Z(2)-Z(1));
end
if L>10
    Z=linspace(0,L,30);
    del_z = (Z(2)-Z(1));
end

Nz = length(Z);

d_vk = Fibra.dvk; %125*10^9;                        % 1 nm - Noise Bandwidth
tau = 10e-3;
nucleos = Fibra.nucleos;
Sch = nucleos;                          % N° de nucleos
N = Fibra.N ;

% ----- Datos de Bombeo y Señal ----- %

P_s0 = Signal.P0;
P_p0 = Pump.P0;

ModoS = strcat("LP_" , Signal.modos(:));
ModoP = strcat("LP_" , Pump.modos(:));

Smod = length(ModoS);                   % N° de modos de señal
Pmod = length(ModoP);                   % N° de modos de bombeo

P_ase0 = ASE;             % Potencia ASE entrada en Watts
P_ase0_Nulo = 1e-3*10.^(( -200 )/10); %-200

F_ase1= c / Signal.lambda.(ModoS(1))(end);  %189.25e12; % (Ex DifPump:191.19421875e12)
F_ase_end= c / Signal.lambda.(ModoS(1))(1);  %198.85e12; % (Ex DifPump: 193.64421875e12)
delta_noise=12.5e9;

% ----- Espectros emision y absorción ----- %

% --- Datos obtenidos de VPI --- %
if Fibra.SigmaSpectrum == "VPI"
    VPI = load('Erbium_VPI.dat');
    Sa = VPI(:,3); Se = VPI(:,2);
    lambda_cross = VPI(:,1).*1e-9;
end

% --- Datos OptiSystem --- %
if Fibra.SigmaSpectrum == "OptiSystem"
    OptiSystem = load('Erbium_OptiSystem.dat');
    Sa = OptiSystem(:,2); Se = OptiSystem(:,3);
    lambda_cross = OptiSystem(:,1).*1e-9;
end


sigma_abs = fit(lambda_cross,Sa,'linearinterp');
sigma_ems = fit(lambda_cross,Se,'linearinterp');


for i=1:length(ModoS)
    A_s.(ModoS(i)) = pi*Fibra.radio^2;                                % Area efectiva para señal
end
for i=1:length(ModoP)
    A_p.(ModoP(i)) = pi*Fibra.radio^2 ;                               % Area efectiva para bombeo
end


sdm = struct;


% ----- Calculo de Gammas ----- %
warning('off')

for p = 1:1:Pmod % Mode overlap factor for pump; entre modo y perfil de dopaje (uniforme)
    Nwlp = length(Pump.lambda.(ModoP(p)));
    for i=1:1:Nwlp % Cada longitud de onda del modo p
        lambda_p = Pump.lambda.(ModoP(p));
        [gamma_p.(ModoP(p)){i},beta0_p.(ModoP(p)){i}] = norm_intensity2(Fibra,Pump.modos(p),lambda_p(i));
        gamma_p.(ModoP(p)){i} = 0.77;
    end
end

for s = 1:1:Smod % Mode overlap factor for signal entre modo y perfil de dopaje
    Nwl = length(Signal.lambda.(ModoS(s)));
    for i=1:1:Nwl % Cada longitud de onda del modo s
        lambda_s = Signal.lambda.(ModoS(s));
        [gamma_s.(ModoS(s)){i},beta0_s.(ModoS(s)){i}] = norm_intensity2(Fibra,Signal.modos(s),lambda_s(i)) ; % *CoupCoef ;
    end
end
warning('on')

% ----- Power Coupling Coeficient ----- %
% Power coupling coeficient pump
h_pccP = pcc_calc(beta0_p,ModoP); 
% Power coupling coeficient signal
h_pccS = pcc_calc(beta0_s,ModoS);

% ----- ASE Spectrum ----- %
% Crear vector de longitudes de onda para ruido ase
allwavelengths = [];
for s = 1:1:Smod
    for w=1:1:length(Signal.lambda.(ModoS(s)))
        if ~(ismember(Signal.lambda.(ModoS(s))(w) , allwavelengths))
            allwavelengths = [allwavelengths Signal.lambda.(ModoS(s))(w)];
        end
    end
end
allwavelengths = sort(allwavelengths);

if Fibra.ASEFlag == 1  % retorna [0,0,0], evita calcular espectro ASE
    [lambda_ase,frequency_ase]= ase_freqVPI(F_ase1,F_ase_end,delta_noise,c);

else
    [lambda_ase,frequency_ase]= ase_freqVPI(F_ase1,F_ase_end,delta_noise,c);
end


%%
fprintf('Iniciando cálculo...\n')
for n = 1:1:Sch     % Iteración en nucleos
    ch = strcat('Nucleo',int2str(n));
    % Variables para cálculo de N2 y N1
    N2 = zeros(Nz,1);                    % Densidad de poblacion en estado excitado
    N1 = zeros(Nz,1);                    % Densidad de pobiacion en estado basal
    Nt = zeros(Nz,1);
    Psp = [];                           % Potencia de señal en dirección +z
    Ppp = [];                           % Potencia de bombeo en dirección +z
    Pap = [];                           % Potencia ASE en dirección +z
    Pap_amp = [];                       % Potencia ASE en dirección +z , Solo componente amplificado
    Pap_tot = [];                       % Potencia ASE en dirección +z , Ambas componentes
    %Pan = [];                           % Potencia ASE en dirección -z
    Pase = [];                          % Potencia total ASE en el EDFA
    OSNR = [];                          % Relacion señal a ruido

    for s = 1:1:Smod
        Nwl = length(Signal.lambda.(ModoS(s)));
        Pan.(ModoS(s)) = zeros(Nwl,Nz);
    end

    gain = [];

    Pap_sp = [];
    Pan_sp = [];
    ASE_Spectrum = [];

    % Calculo de densidades de iones en estado basal y excitado (N1 y N2)

    z_waitbar = waitbar(0 ,'','Name',"Cálculo a lo largo del EDFA: " , 'Visible', 'off') ;

    QQ = 3; % Iteraciones para aumentar precisión
    for Q = 1:QQ
        if (Q == 1) % Primera iteracion no calcula ASE en direccion -z
            for z = 1:1:Nz  % Iteraciones a lo largo del EDFA
                % % Eq19 Giles 1991 para N(z)
                % Usando los cambios de variable : sigma_a * gamma = alpha/nt ; A = pi*b_eff

                % Primera Iteración
                sig_xx = 0;
                sig_yy = 0;                    % Inicialización de variables
                ase_xx = 0;
                ase_yy = 0;
                pmp_xx = 0;
                pmp_yy = 0;
                if(z == 1)
                    % Potencia Bombeo
                    parfor p = 1:Pmod
                        Nwlp = length(Pump.lambda.(ModoP(p)));
                        lambda_p = Pump.lambda.(ModoP(p)); v_p = c./lambda_p;
                        for i=1:Nwlp
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            Pp0 = P_p0.(ModoP(p))(i);

                            pmp_xx = pmp_xx + (sigma_abs(lambda_p(i)))*(Pp0*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));                                   % Término en numerador
                            pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Pp0*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));          % Término en denominador
                        end
                    end
                    % Potencia de señal
                    parfor s = 1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        for i = 1:1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Ps0 = P_s0.(ModoS(s))(i);

                            sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(Ps0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));                        % Término en numerador
                            sig_yy = sig_yy+(sigma_abs(lambda_s(i))+sigma_ems(lambda_s(i)))*(Ps0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Término en denominador

                        end
                    end

                    % Potencia ASE
                    parfor s=1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        for i = 1:1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pase0 = P_ase0.(ModoS(s))(i);

                            ase_xx = ase_xx+ sigma_abs(lambda_s(i))*(Pase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));                           % Término en numerador
                            ase_yy = ase_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(Pase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Término en denominador
                            
                        end
                    end

                else
                    % Iteraciones Restantes z>1
                    % Bombeo
                    parfor p = 1:1:Pmod
                        Nwlp = length(Pump.lambda.(ModoP(p)));
                        lambda_p = Pump.lambda.(ModoP(p)); v_p = c./lambda_p;
                        for i = 1:1:Nwlp
                            Gamma_p = gamma_p.(ModoP(p)){i};

                            pmp_xx = pmp_xx + sigma_abs(lambda_p(i))*(Ppp.(ModoP(p))(i,z-1)*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));                        % Termino en numerador
                            pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Ppp.(ModoP(p))(i,z-1)*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));           % Termino en denominador
                        end
                    end

                    % Señal
                    parfor s = 1:1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        for i = 1:1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};

                            sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(Psp.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));
                            sig_yy = sig_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i))) * (Psp.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));
                        end
                    end

                    % ASE
                    parfor s = 1:1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        for i = 1:1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};

                            ase_xx = ase_xx+(sigma_abs(lambda_s(i))*(Pase.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));
                            ase_yy = ase_yy+( sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(Pase.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));
                        end
                    end
                end

                N2(z) = ( (pmp_xx + sig_xx + ase_xx) / (pmp_yy + sig_yy + ase_yy + (1/tau)) )*N;        % Densidad de iones de Erbio en estado excitado
                N1(z) = N-N2(z);                                                                        % Densidad de iones de Erbio en estado basal
                Nt(z) = N1(z) + N2(z);
%                 if z==1
%                     N1(z) = N-N2(z);                                                                        % Densidad de iones de Erbio en estado basal
%                 else
%                     N1(z) = N1(z-1);
%                 end
%                 Nt(z) = N;%N1(z) + N2(z);



                %% Ecuaciones de Potencias

                % Ecuacion diferencial para bombeo en direccion +z
                for p = 1:1:Pmod      % Iteracion en cada modo de bombeo
                    Nwlp = length(Pump.lambda.(ModoP(p)));
                    lambda_p = Pump.lambda.(ModoP(p));
                    if(z == 1)
                        parfor i = 1:Nwlp
                            Pp0 = P_p0.(ModoP(p))(i);
                            PppAux(i) = Pp0 ;
                        end
                        Ppp.(ModoP(p))(:,z) = PppAux(:);

                    else
                        parfor i = 1:Nwlp
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            PppAux(i) = Ppp.(ModoP(p))(i,z-1) + ((N2(z))*sigma_ems(lambda_p(i)) - (N1(z))*sigma_abs(lambda_p(i)))*Gamma_p*Ppp.(ModoP(p))(i,z-1)*del_z;
                        end
                        Ppp.(ModoP(p))(:,z) = PppAux(:);
                    end
                end
                % Añadiendo PCC en cada modo de bombeo 
                Ppp_aux=Ppp;
                for pu = 1:1:Pmod
                    for pv=1:Pmod
                        if (pu~=pv)
                            Ppp.(ModoP(pu))=Ppp_aux.(ModoP(pu))-h_pccP(pu,pv).*(Ppp_aux.(ModoP(pu))-Ppp_aux.(ModoP(pv)));
                        end
                    end
                end

                % Ecuacion diferencial para señal en direccion +z
                for s = 1:1:Smod
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    lambda_s = Signal.lambda.(ModoS(s));
                    if(z == 1)
                        parfor i = 1:Nwl
                            Ps0 = P_s0.(ModoS(s))(i);
                            PspAux(i) = Ps0
                        end
                        Psp.(ModoS(s))(:,z) = PspAux(:);
                    else
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PspAux(i) = Psp.(ModoS(s))(i,z-1) + ((N2(z))*sigma_ems(lambda_s(i))-(N1(z))*sigma_abs(lambda_s(i)))*Gamma_s*Psp.(ModoS(s))(i,z-1)*del_z;
                        end
                        Psp.(ModoS(s))(:,z) = PspAux(:);
                    end
                end
                % Añadiendo PCC en cada modo de señal
                Psp_aux=Psp;
                for su = 1:1:Smod
                   for sv=1:Smod
                        if (su~=sv)
                            Psp.(ModoS(su))=Psp_aux.(ModoS(su))-h_pccS(su,sv).*(Psp_aux.(ModoS(su))-Psp_aux.(ModoS(sv)));
                        end
                   end
                end

                % Ecuacion diferencial para ASE en direccion +z
                % Ruido ASE generado en amplificador (respecto a -200dBm)
                for s = 1:1:Smod
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    lambda_s = Signal.lambda.(ModoS(s));
                    v_s = c./lambda_s;
                    if(z == 1)
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux(i) = P_ase0_Nulo;%P_ase0.(ModoS(s))(i); 
                        end
                        Pap.(ModoS(s))(:,z) = PapAux(:);
                    else
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux(i) = Pap.(ModoS(s))(i,z-1) + (((N2(z))*sigma_ems(lambda_s(i))-(N1(z))*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(ModoS(s))(i,z-1)  + 2*(N2(z))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                        end
                        Pap.(ModoS(s))(:,z) = PapAux(:);
                    end
                end

                % Ecuacion diferencial para ASE Amplificado en direccion +z - Toma como valor de entrada la salida del EDFA anterior; No usa término con m*sigma*h*f*deltaF....
                % Amplifica ASE de entrada como si fuera señal
                for s = 1:1:Smod
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    if(z == 1)
                        for i = 1:1:Nwl
                            Pap_amp.(ModoS(s))(i,z) = P_ase0.(ModoS(s))(i);
                        end
                    else
                        for i = 1:1:Nwl
                            lambda_s = Signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux_amp(i) = Pap_amp.(ModoS(s))(i,z-1) + ((N2(z,:)*sigma_ems(lambda_s(i))-N1(z,:)*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(ModoS(s))(i,z-1)  )*del_z;
                        end
                        Pap_amp.(ModoS(s))(:,z) = PapAux_amp(:);
                    end
                end

                % Ecuacion diferencial para ASE en direccion -z
                % No se calcula en 1era iteracion en fibra

                % Calculo de Potencia ASE
                for s = 1:1:Smod
                    Pase.(ModoS(s))(:,z) = Pap.(ModoS(s))(:,z)+Pan.(ModoS(s))(:,z);
                    Pap_tot.(ModoS(s))(:,z) = Pap.(ModoS(s))(:,z) + Pap_amp.(ModoS(s))(:,z) ;
                end
                % Mostrar avance como prints en pantalla:
                if Fibra.Avance
                    if Fibra.ASEFlag == 0
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  ((z/Nz)*(1/QQ) + ((Q-1)/QQ))/2  ) % Mostrar % de avance del cálculo
                    else
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ) + ((Q-1)/QQ)  ) % Mostrar % de avance del cálculo
                    end
                end
                % Mostrar avance como WaitBar:
                if Fibra.WaitBar
                    if Fibra.ASEFlag == 0
                        avance = ((z/Nz)*(1/QQ) + ((Q-1)/QQ))/2 ; waitbar(avance, z_waitbar , sprintf("Cálculo a lo largo del EDFA: \n %.2f %%",avance*100 ) ) ; set(z_waitbar,'Visible', 'on');
                    else
                        avance = (z/Nz)*(1/QQ) + ((Q-1)/QQ) ; waitbar(avance, z_waitbar , sprintf("Cálculo a lo largo del EDFA: \n %.2f %%",avance*100 ) ) ; set(z_waitbar,'Visible', 'on');
                    end
                end
            end   % Fin de calculo en fibra

        else  % Fin 1era iteración, usar estos valores de N para una 2da iteración mas precisa
            % Agregando ASE en dirección -z

            % ----- Guardado de datos sin GEF ----- %
                    
            for s = 1:1:Smod
                Nwl = length(Signal.lambda.(ModoS(s)));
                for i = 1:1:Nwl
                    OSNR_sinGEF.(ModoS(s))(i,:) = 10*log10(Psp.(ModoS(s))(i,:)/1e-3) - 10*log10(Pap_tot.(ModoS(s))(i,:)/1e-3); %Pase
                    gain_sinGEF.(ModoS(s))(1,i) = 10*log10(Psp.(ModoS(s))(i,end)/Psp.(ModoS(s))(i,1));
                end
                sdm.(ch).DatosSinGEF.NF.(ModoS(s)) = OSNR_sinGEF.(ModoS(s))(:,1)./OSNR_sinGEF.(ModoS(s))(:,end);
                sdm.(ch).DatosSinGEF.OSNR.(ModoS(s)) = OSNR_sinGEF.(ModoS(s)) ;
                sdm.(ch).DatosSinGEF.Gain.(ModoS(s)) = gain_sinGEF.(ModoS(s)) ;
                sdm.(ch).DatosSinGEF.Ps.(ModoS(s)) = Psp.(ModoS(s));
                sdm.(ch).DatosSinGEF.Pase.(ModoS(s)) = Pase.(ModoS(s));
            end

            for z = 1:1:Nz %Iteraciones a lo largo del EDFA
                sig_xx = 0;
                sig_yy = 0;                    % Inicialización de variables
                ase_xx = 0;
                ase_yy = 0;
                pmp_xx = 0;
                pmp_yy = 0;
                % Siguientes Primeras Iteraciones

                if(z == 1)
                    % Potencia Bombeo
                    parfor p = 1:Pmod
                        Nwlp = length(Pump.lambda.(ModoP(p)));
                        lambda_p = Pump.lambda.(ModoP(p)); v_p = c./lambda_p;
                        for i=1:1:Nwlp
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            Pp0 = P_p0.(ModoP(p))(i);

                            pmp_xx = pmp_xx + (sigma_abs(lambda_p(i)))*(Pp0*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));                                   % Termino en numerador
                            pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Pp0*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));          % Termino en denominador
                        end
                    end
                    % Potencia de señal
                    parfor s = 1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        for i = 1:1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Ps0 = P_s0.(ModoS(s))(i);

                            sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(Ps0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));                        % Termino en numerador
                            sig_yy = sig_yy+(sigma_abs(lambda_s(i))+sigma_ems(lambda_s(i)))*(Ps0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Termino en denominador
                        end
                    end
                    % Potencia ASE
                    parfor s=1:1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        for i = 1:1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pase0 = P_ase0.(ModoS(s))(i);

                            ase_xx = ase_xx+ sigma_abs(lambda_s(i))*(Pase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));                           % Termino en numerador
                            ase_yy = ase_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(Pase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Termino en denominador
                        end
                    end

                else % Iteraciones Restantes z>1

                    if (z == round(Nz/2,0) && Q == QQ) % Aplicar Filtro Equalizador
                        %load('GEF_Filters/EDFA_5m_300mw_OptiSystem') ;
                        Filtro = Fibra.GEF.Filtro;
                        lambda = Signal.lambda.(ModoS(1));
                        GEF = fit(lambda',Filtro,'linearinterp');
    
                        % Aplicar filtro
                        for s = 1:1:Smod  
                            for f = 1:length(lambda)
                                Psp.(ModoS(s))(f,z-1) = Psp.(ModoS(s))(f,z-1) * ( 1 - GEF(lambda(f)) );
                                Pap.(ModoS(s))(f,z-1) = Pap.(ModoS(s))(f,z-1) * ( 1 - GEF(lambda(f)) );
                            end
                        end
                    end
                    
                    % Bombeo
                    parfor p = 1:Pmod
                        Nwlp = length(Pump.lambda.(ModoP(p)));
                        lambda_p = Pump.lambda.(ModoP(p)); v_p = c./lambda_p;
                        for i = 1:1:Nwlp
                            Gamma_p = gamma_p.(ModoP(p)){i};

                            pmp_xx = pmp_xx + sigma_abs(lambda_p(i))*(Ppp.(ModoP(p))(i,z-1)*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));                        % Termino en numerador
                            pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Ppp.(ModoP(p))(i,z-1)*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));           % Termino en denominador
                        end
                    end

                    % Señal
                    parfor s = 1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        for i = 1:1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};

                            sig_xx = sig_xx + (sigma_abs(lambda_s(i))*(Psp.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));
                            sig_yy = sig_yy + (sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i))) * (Psp.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));
                        end
                    end

                    % ASE
                    parfor s = 1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        for i = 1:1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};

                            ase_xx = ase_xx + (sigma_abs(lambda_s(i))*(Pase.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));
                            ase_yy = ase_yy + ( sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(Pase.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));
                        end
                    end
                end

                N2(z) = ( (pmp_xx + sig_xx + ase_xx)/(pmp_yy + sig_yy + ase_yy + (1/tau)) )*N;             % Densidad de iones de Erbio en estado excitado
                N1(z) = N-N2(z);                                                                        % Densidad de iones de Erbio en estado basal
                Nt(z) = N1(z) + N2(z);
%                 if z==1
%                     N1(z) = N-N2(z);                                                                        
%                 else
%                     N1(z) = N1(z-1);
%                 end
%                 Nt(z) = N;%N1(z) + N2(z);

                % % Ecuaciones de Potencias

                % Ecuacion diferencial para bombeo en direccion +z
                for p = 1:1:Pmod      % Iteracion en cada modo de bombeo
                    Nwlp = length(Pump.lambda.(ModoP(p)));
                    lambda_p = Pump.lambda.(ModoP(p));
                    if(z == 1)
                        parfor i = 1:Nwlp
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            Pp0 = P_p0.(ModoP(p))(i);
                            PppAux(i) = Pp0
                        end
                        Ppp.(ModoP(p))(:,z) = PppAux(:);
                    else
                        parfor i = 1:Nwlp
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            PppAux(i) = Ppp.(ModoP(p))(i,z-1) + ((N2(z))*sigma_ems(lambda_p(i)) - (N1(z))*sigma_abs(lambda_p(i)))*Gamma_p*Ppp.(ModoP(p))(i,z-1)*del_z;
                        end
                        Ppp.(ModoP(p))(:,z) = PppAux(:);
                    end
                end

                % Añadiendo PCC en cada modo de bombeo
                Ppp_aux=Ppp;
                for pu = 1:1:Pmod
                    for pv=1:Pmod
                        if (pu~=pv)
                                Ppp.(ModoP(pu))=Ppp_aux.(ModoP(pu))-h_pccP(pu,pv).*(Ppp_aux.(ModoP(pu))-Ppp_aux.(ModoP(pv)));
                        end
                    end
                end

                % Ecuacion diferencial para señal en direccion +z
                for s = 1:1:Smod
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    lambda_s = Signal.lambda.(ModoS(s));
                    if(z == 1)
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Ps0 = P_s0.(ModoS(s))(i);
                            PspAux(i) = Ps0 
                        end
                        Psp.(ModoS(s))(:,z) = PspAux(:);
                    else
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PspAux(i) = Psp.(ModoS(s))(i,z-1) + ((N2(z))*sigma_ems(lambda_s(i))-(N1(z))*sigma_abs(lambda_s(i)))*Gamma_s*Psp.(ModoS(s))(i,z-1)*del_z;
                        end
                        Psp.(ModoS(s))(:,z) = PspAux(:);
                    end
                end
                Psp_aux=Psp;
                 % Añadiendo PCC en cada modo de señal
                for su = 1:1:Smod
                   for sv=1:Smod
                        if (su~=sv)
                            Psp.(ModoS(su))=Psp_aux.(ModoS(su))-h_pccS(su,sv).*(Psp_aux.(ModoS(su))-Psp_aux.(ModoS(sv)));
                        end
                   end
                end

                % Ecuacion diferencial para ASE en direccion +z
                % Siempre respecto a -200 dBm entrada
                for s = 1:1:Smod
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                    if(z == 1)
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux(i) = P_ase0_Nulo;
                        end
                        Pap.(ModoS(s))(:,z) = PapAux(:);
                    else
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux(i) = Pap.(ModoS(s))(i,z-1) + (((N2(z))*sigma_ems(lambda_s(i)) - (N1(z))*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(ModoS(s))(i,z-1)  + 2*(N2(z))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                        end
                        Pap.(ModoS(s))(:,z) = PapAux(:);
                    end
                end

                % Ecuacion diferencial para ASE Amplificado en direccion +z - Toma como valor de entrada la salida del EDFA anterior; No usa término con m*sigma*h*f*deltaF....
                % Amplifica ruido ASE de entrada como si fuera señal
                for s = 1:1:Smod
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    if(z == 1)
                        for i = 1:1:Nwl
                            Pap_amp.(ModoS(s))(i,z) = P_ase0.(ModoS(s))(i);
                        end
                    else
                        for i = 1:1:Nwl
                            lambda_s = Signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux_amp(i) = Pap_amp.(ModoS(s))(i,z-1) + ((N2(z,:)*sigma_ems(lambda_s(i))-N1(z,:)*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(ModoS(s))(i,z-1) )*del_z; 
                        end
                        Pap_amp.(ModoS(s))(:,z) = PapAux_amp(:);
                    end
                end

                % Ecuacion diferencial para ASE en direccion -z
                for s = 1:1:Smod
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                    if(z == 1)
                        for i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pase0 = P_ase0_Nulo; 

                            Pan.(ModoS(s))(i,Nz-z+1) = Pase0; 
                        end
                        %Pan.(ModoS(s))(:,Nz-z+1) = PanAux(:);
                    else
                        for i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pan.(ModoS(s))(i,Nz-z+1) = Pan.(ModoS(s))(i,Nz-z+1+1) + (((N2(Nz-z+1))*sigma_ems(lambda_s(i)) - (N1(Nz-z+1))*sigma_abs(lambda_s(i)))*Gamma_s*Pan.(ModoS(s))(i,Nz-z+1+1) + 2*(N2(Nz-z+1))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                        end
                        %Pan.(ModoS(s))(:,Nz-z+1) = PanAux(:);
                    end
                end
                
                for s = 1:1:Smod
                    Pase.(ModoS(s))(:,z) = Pap.(ModoS(s))(:,z) + Pan.(ModoS(s))(:,z);
                    Pap_tot.(ModoS(s))(:,z) = Pap.(ModoS(s))(:,z) + Pap_amp.(ModoS(s))(:,z) ;
                end
                % Mostrar % de avance como prints en pantalla:
                if Fibra.Avance
                    if Fibra.ASEFlag == 0
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  ((z/Nz)*(1/QQ) + ((Q-1)/QQ))/2 ) % Mostrar % de avance del cálculo
                    else
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ) + ((Q-1)/QQ) ) % Mostrar % de avance del cálculo
                    end
                    
                end
                % Mostrar % de avance como WaitBar:
                if Fibra.WaitBar
                    if Fibra.ASEFlag == 0
                        avance = ((z/Nz)*(1/QQ) + ((Q-1)/QQ))/2 ; waitbar(avance , z_waitbar ,  sprintf("Cálculo a lo largo del EDFA: \n %.2f %%",avance*100 ) ) ; set(z_waitbar,'Visible', 'on');
                    else
                        avance = (z/Nz)*(1/QQ) + ((Q-1)/QQ) ; waitbar(avance , z_waitbar ,  sprintf("Cálculo a lo largo del EDFA: \n %.2f %%",avance*100 ) ) ; set(z_waitbar,'Visible', 'on');
                    end
                end
            end % fin iteraciones en largo de fibra

        end
    end % Fin iteraciones para estabilizar ganancias

    %% ASE SPECTRUM :
    % Ecuacion diferencial para ASE en direccion +z
    v_s_sp = c./lambda_ase;
    Nch_ase = length(lambda_ase);
    d_vk_sp = abs(c/lambda_ase(1) - c/lambda_ase(2));
    if ~(Fibra.ASEFlag == 1)
        for z = 1:Nz
            for s = 1:1:Smod
                P_ase0_sp = mean(P_ase0.(ModoS(s)));
                Gamma_s = sum(cell2mat(gamma_s.(ModoS(s))))/length(cell2mat(gamma_s.(ModoS(s))));
                if(z == 1)
                    parfor i = 1:Nch_ase
                        AuxMatrix(i) = P_ase0_sp;
                    end
                    Pap_sp.(ModoS(s))(:,z) = AuxMatrix(:) ;
        
                else
                    if (z == round(Nz/2,0) && Q == QQ) % Aplicar Filtro Equalizador
                        Filtro = Fibra.GEF.Filtro;
                        lambda = Signal.lambda.(ModoS(1));
                        GEF = fit(lambda',Filtro,'linearinterp');
    
                        % Aplicar filtro
                        for f = 1:Nch_ase
                            Pap_sp.(ModoS(s))(f,z-1) = Pap_sp.(ModoS(s))(f,z-1) * ( 1 - GEF(lambda_ase(f)) );
                        end
                    end
                    parfor i = 1:Nch_ase
                        if (i==1)
                            d_vk_sp=abs(v_s_sp(i)-v_s_sp(i+1));
                        else
                            d_vk_sp=abs(v_s_sp(i-1)-v_s_sp(i));
                        end
                        AuxMatrix(i) = Pap_sp.(ModoS(s))(i,z-1) + (((N2(z))*sigma_ems(lambda_ase(i)) - (N1(z))*sigma_abs(lambda_ase(i)))*Gamma_s*Pap_sp.(ModoS(s))(i,z-1) + 2*(N2(z))*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk_sp)*del_z;
                    end
                    Pap_sp.(ModoS(s))(:,z) = AuxMatrix(:) ;
                end
            end
        
            if Fibra.Avance && Fibra.ASEFlag == 0
                clc ; fprintf("Cálculo Espectro ASE: \n") ; fprintf('%.2f %% \n' ,  ( 50 + (((z-0.5)/Nz)*50 )) ) % Mostrar % de avance delase cálculo
            end
            % Mostrar % de avance como WaitBar:
            if Fibra.WaitBar && Fibra.ASEFlag == 0
                avance = 0.50 + (((z-0.5)/Nz))/2 ; waitbar(avance , z_waitbar ,  sprintf("Cálculo Espectro ASE: \n %.2f %%",avance*100 ) ) ; set(z_waitbar,'Visible', 'on');
            end
        
            % Ecuacion diferencial para ASE en direccion -z
            for s = 1:1:Smod
                Gamma_s = sum(cell2mat(gamma_s.(ModoS(s))))/length(cell2mat(gamma_s.(ModoS(s))));
                if(z == 1)
                    parfor i = 1:Nch_ase
                        AuxMatrix(i) = P_ase0_Nulo;%P_ase0_sp
                    end
                    Pan_sp.(ModoS(s))(:,Nz-z+1) = AuxMatrix(:) ;
                else
                    
                    parfor i = 1:Nch_ase
                        if (i==1)
                            d_vk_sp=abs(v_s_sp(i)-v_s_sp(i+1));
                        else
                            d_vk_sp=abs(v_s_sp(i-1)-v_s_sp(i));
                        end
                        AuxMatrix(i) = Pan_sp.(ModoS(s))(i,Nz-z+1+1) + (((N2(Nz-z+1))*sigma_ems(lambda_ase(i)) - (N1(Nz-z+1))*sigma_abs(lambda_ase(i)))*Gamma_s*Pan_sp.(ModoS(s))(i,Nz-z+1+1) + 2*(N2(Nz-z+1))*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk_sp)*del_z;
                    end
                    Pan_sp.(ModoS(s))(:,Nz-z+1) = AuxMatrix(:) ;
                end
            end
            % Mostrar % de avance como prints en pantalla:
            if Fibra.Avance && Fibra.ASEFlag == 0
                clc ; fprintf("Cálculo Espectro ASE: \n") ; fprintf('%.2f %% \n' ,  ( 50 + ((z/Nz)*50 )) ) % Mostrar % de avance del cálculo
            end
            % Mostrar % de avance como WaitBar:
            if Fibra.WaitBar && Fibra.ASEFlag == 0
                avance = 0.50 + ((z/Nz))/2 ; waitbar(avance , z_waitbar ,  sprintf("Cálculo Espectro ASE: \n %.2f %%",avance*100 ) ) ; set(z_waitbar,'Visible', 'on');
            end
        end % Fin calculo ASE en posicion z
        
        % ASE Spectrum es el valor del ruido ASE al final del amplificador para
        % cada frecuencia
        for s = 1:1:Smod
            for i = 1:1:Nch_ase
                ASE_Spectrum.(ModoS(s))(i) = 10*log10( (Pap_sp.(ModoS(s))(i,end) + Pan_sp.(ModoS(s))(i,end))/1e-3 );
            end
        end
    end

    if Fibra.WaitBar
        close(z_waitbar)
    end

    % ---- Representation using Noise Bins ---- %
    TotalNoise_bandwidth=F_ase_end-F_ase1;
    NoiseBin_num=floor(TotalNoise_bandwidth/d_vk) + 1; % Number of Noise Bins
    NoiseBin_aux=d_vk/delta_noise;

    for s = 1:Smod
        Pap_spectrum.(ModoS(s))(:) = Pap_sp.(ModoS(s))(:,end);
        aux_f=1;
        for i=1:NoiseBin_num % Noise Bins power 
            if(i==1)
                f_left(i) = frequency_ase(aux_f)-d_vk/2;
                f_right(i) = frequency_ase(aux_f+NoiseBin_aux/2);
                NoiseBin_power.(ModoS(s))(i) = mean(Pap_spectrum.(ModoS(s))(aux_f:NoiseBin_aux/2+1));
                NoiseBin_power_z.(ModoS(s))(i,:) = mean(Pap_sp.(ModoS(s))(aux_f:NoiseBin_aux/2+1,:));
                aux_f = NoiseBin_aux/2+1;
            else if (i<NoiseBin_num)
                f_left(i) = frequency_ase(aux_f);
                f_right(i) = frequency_ase(aux_f+NoiseBin_aux);
                NoiseBin_power.(ModoS(s))(i) = mean(Pap_spectrum.(ModoS(s))(aux_f:aux_f+NoiseBin_aux));
                NoiseBin_power_z.(ModoS(s))(i,:) = mean(Pap_sp.(ModoS(s))(aux_f:aux_f+NoiseBin_aux,:));
                aux_f = aux_f+NoiseBin_aux;
            else % i==NoiseBin_num
                f_left(i) = frequency_ase(aux_f);
                f_right(i) = frequency_ase(aux_f+NoiseBin_aux/2);
                NoiseBin_power.(ModoS(s))(i) = mean(Pap_spectrum.(ModoS(s))(aux_f:aux_f+NoiseBin_aux/2));
                NoiseBin_power_z.(ModoS(s))(i,:) = mean(Pap_sp.(ModoS(s))(aux_f:aux_f+NoiseBin_aux/2,:));
            end
            end
        end
        Frequency_gridS = c./Signal.lambda.(ModoS(s)) ;
        for j = 1:length(Frequency_gridS)
            for jj = 1:NoiseBin_num
                if(f_left(jj)<Frequency_gridS(j))&&(Frequency_gridS(j)<f_right(jj))
                    Pap_forsignal.(ModoS(s))(j) = NoiseBin_power.(ModoS(s))(jj);
                    Pap_forsignal_z.(ModoS(s))(length(Frequency_gridS)+1-j,:) = NoiseBin_power_z.(ModoS(s))(jj,:);
                elseif jj == NoiseBin_num
                    if(f_left(end)<=Frequency_gridS(j))&&(Frequency_gridS(j)<=frequency_ase(end))
                        Pap_forsignal.(ModoS(s))(j) = NoiseBin_power.(ModoS(s))(jj); 
                        Pap_forsignal_z.(ModoS(s))(length(Frequency_gridS)+1-j,:) = NoiseBin_power_z.(ModoS(s))(jj,:); 
                    end
                end
            end
        end
        Pap_forsignal.(ModoS(s)) = fliplr(Pap_forsignal.(ModoS(s))); % ASE Noise a la salida por canal por modo
    end


    %%

    % Calculo OSNR
    for s = 1:1:Smod
        Nwl = length(Signal.lambda.(ModoS(s)));
        for i = 1:1:Nwl
            OSNR.(ModoS(s))(i,:) = 10*log10(Psp.(ModoS(s))(i,:)/1e-3) - 10*log10(Pap_tot.(ModoS(s))(i,:)/1e-3); %Pase
            OSNRForward.(ModoS(s))(i) = 10*log10(Psp.(ModoS(s))(i,end)/1e-3) - 10*log10(Pap_forsignal.(ModoS(s))(i)/1e-3);
            OSNRForward_z.(ModoS(s))(i,:) = 10*log10(Psp.(ModoS(s))(i,:)/1e-3) - 10*log10(Pap_forsignal_z.(ModoS(s))(i,:)/1e-3);
            gain.(ModoS(s))(1,i) = 10*log10(Psp.(ModoS(s))(i,end)/Psp.(ModoS(s))(i,1));
            freq = c/Signal.lambda.(ModoS(s))(i);
            sdm.(ch).NF_v2.(ModoS(s))(i) = 10*log10( (( (Pap_forsignal_z.(ModoS(s))(i,end))/(h*freq*d_vk) ) + 1 )* (Psp.(ModoS(s))(i,1)/Psp.(ModoS(s))(i,end))   );
            sdm.(ch).OSNR_v2.(ModoS(s))(i) = OSNR.(ModoS(s))(i,2) - sdm.(ch).NF_v2.(ModoS(s))(i);
        end
    end


% CALCULAR GRAFICO ASE SPECTRUM CON GANANCIAS DE SEÑAL v2
    if ~(Fibra.ASEFlag == 1)
        for s = 1:1:Smod
            ASESpectFun.(ModoS(s)) = fit(lambda_ase',(ASE_Spectrum.(ModoS(s))'),'linearinterp');
            
            ase_lam = ase_lambdas(allwavelengths);
            lambda_s = Signal.lambda.(ModoS(s));
            ASE_SP.(ModoS(s)) = (ASESpectFun.(ModoS(s))(ase_lam))' ; 
            
            for w = 1:length(lambda_s)
                for i = 1:length(ase_lam)
                    if abs(ase_lam(i)*1e9 - lambda_s(w)*1e9) < 1/11
                        ASE_SP.(ModoS(s))(i) = 10*log10( Psp.(ModoS(s))(w,end) / 1e-3 );
                    end
                end
            end
        end
    else
        ASESpectFun = "No calculado";
        ASE_SP = "No calculado";
        ase_lam = "No calculado";
    end

    % Guardar los datos calculados
    sdm.(ch).N2 = N2 ; sdm.(ch).N1 = N1 ; sdm.(ch).Nt = Nt;
    for p=1:Pmod
        sdm.(ch).pump.Potencia_dBm.(ModoP(p)) = 10*log10(Ppp.(ModoP(p))./1e-3) ;
        sdm.(ch).salida.pump.potencia_dBm.(ModoP(p)) = 10*log10(Ppp.(ModoP(p))(:,end)./1e-3);
    end
    for s=1:Smod
        sdm.(ch).signal.Potencia_dBm.(ModoS(s)) = 10*log10(Psp.(ModoS(s))./1e-3) ;
        sdm.(ch).salida.signal.potencia_dBm.(ModoS(s)) = 10*log10(Psp.(ModoS(s))(:,end)./1e-3);

        sdm.(ch).Pase.(ModoS(s)) = 10*log10(Pase.(ModoS(s))./1e-3);
        sdm.(ch).Pap.(ModoS(s)) = 10*log10(Pap.(ModoS(s))./1e-3);
        sdm.(ch).Pap_amp.(ModoS(s)) = 10*log10(Pap_amp.(ModoS(s))./1e-3);
        sdm.(ch).Pap_tot.(ModoS(s)) = 10*log10(Pap_tot.(ModoS(s))./1e-3);
        sdm.(ch).Pan.(ModoS(s)) = 10*log10(Pan.(ModoS(s))./1e-3);
        sdm.(ch).salida.ASE.potencia_dBm.(ModoS(s)) = 10*log10(Pase.(ModoS(s))(:,end)./1e-3);

        sdm.(ch).salida.ganancias.(ModoS(s)) = gain.(ModoS(s));
        sdm.(ch).NF.(ModoS(s)) = OSNR.(ModoS(s))(:,2) - OSNR.(ModoS(s))(:,end);
        sdm.(ch).NF_NoiseBins.(ModoS(s)) = OSNRForward_z.(ModoS(s))(:,2) - OSNRForward_z.(ModoS(s))(:,end);

    end

    sdm.(ch).signal.lambdas = lambda_s ; sdm.(ch).pump.lambdas = lambda_p;
    sdm.(ch).OSNR = OSNR ;
    sdm.(ch).salida.OSNR = OSNR(:,end);
    sdm.(ch).OSNRF= OSNRForward;
    sdm.(ch).OSNRF_z = OSNRForward_z ;
    sdm.(ch).z = Z;
    sdm.(ch).ASE_Spectrum.mag = ASE_SP;
    sdm.(ch).ASE_Spectrum.fun = ASESpectFun;
    sdm.(ch).ASE_Spectrum.lambdas = ase_lam;
    sdm.(ch).Gamma.Gamma_p = gamma_p;
    sdm.(ch).Gamma.Gamma_s = gamma_s;
    sdm.(ch).ASE_Spectrum.Pap_sp = Pap_sp;
    sdm.(ch).ASE_Spectrum.Pan_sp = Pan_sp;

end % Fin cálculo en todos los nucleos

edfa = sdm;

end % Fin funcion
