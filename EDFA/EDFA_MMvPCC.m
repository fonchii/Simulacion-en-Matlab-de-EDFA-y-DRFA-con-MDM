function edfa = EDFA_MMvPCC(fibra,signal,pump,ASE)
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
c = 3e8;                                % Velocidad de la luz en el vacio m/s
h = 6.626*10^(-34);                     % constante de planck
L = fibra.largo ;                       % largo
del_z = 1;                              % Tamaño de paso de iteraciones
Z = linspace(0,L,L*(1/del_z));          % se toman segmentos de del_z m
if length(Z)<20                         % Correccion para EDFA muy corto
    Z=linspace(0,L,20);
    del_z = (Z(2)-Z(1));
end
Nz = length(Z);

d_vk = 125*10^9;                        % 1 nm - Noise Bandwidth
tau = 10e-3;
nucleos = fibra.nucleos;
Sch = nucleos;                          % N° de nucleos
N = fibra.N ;

% Datos de Bombeo y Señal

P_s0 = signal.P0;
P_p0 = pump.P0;

ModoS = strcat("LP_" , signal.modos(:));
ModoP = strcat("LP_" , pump.modos(:));

Smod = length(ModoS);                   % N° de modos de señal
Pmod = length(ModoP);                   % N° de modos de bombeo

P_ase0 = 1e-3*10.^(ASE/10);             % Potencia ASE entrada en Watts

Psat = 10;                              % Potencia de saturación (dBm)
Psat = 1e-3*10^(Psat/10);


% %   Espectros emision y absorción

% Datos obtenidos de OptiSystem.
% E = load('Erbium.dat');
% lambda_cross = E(:,1).*1e-9;
% Sa = E(:,2); Se = E(:,3);
% sigma_abs = fit(lambda_cross,Sa,'linearinterp');
% sigma_ems = fit(lambda_cross,Se,'linearinterp');

% Datos obtenidos de VPI
VPI = load('Erbium_VPI.dat');
lambda_cross = VPI(:,1).*1e-9;
Sa = VPI(:,3); Se = VPI(:,2);
sigma_abs = fit(lambda_cross,Sa,'linearinterp');
sigma_ems = fit(lambda_cross,Se,'linearinterp');


for i=1:length(ModoS)
    A_s.(ModoS(i)) = pi*fibra.radio^2;                                % Area efectiva para señal
end
for i=1:length(ModoP)
    A_p.(ModoP(i)) = pi*fibra.radio^2 ;                               % Area efectiva para bombeo
end

% De Paper ó VPIphotonics 
% Para radio=10e-6 y AN=0.2
%A_s.LP_01=80e-12; A_s.LP_11_a=76e-12; A_s.LP_11_b=76e-12; A_s.LP_02=83e-12; A_s.LP_21_a=86e-12; A_s.LP_21_b=86e-12;
%A_p.LP_01=80e-12; A_p.LP_11_a=76e-12; A_p.LP_11_b=76e-12; A_p.LP_02=83e-12; A_p.LP_21=86e-12;
% Para radio=6e-6 y AN=0.2
% A_s.LP_01=80e-12; A_s.LP_11_a=76e-12; A_s.LP_11_b=76e-12; A_s.LP_02=83e-12; A_s.LP_21_a=86e-12; A_s.LP_21_b=86e-12;
% A_p.LP_01=80e-12; A_p.LP_11_a=76e-12; A_p.LP_11_b=76e-12; A_p.LP_02=83e-12; A_p.LP_21=86e-12;
sdm = struct;

% Calculo de Gammas
warning('off')

for p = 1:1:length(pump.modos) % Mode overlap factor for pump; entre modo y perfil de dopaje (uniforme)
    Nwlp = length(pump.lambda.(ModoP(p)));
    for i=1:1:Nwlp % Cada longitud de onda del modo p
        lambda_p = pump.lambda.(ModoP(p));
        [gamma_p.(ModoP(p)){i},beta0_p.(ModoP(p)){i}] = norm_intensity(fibra,pump.modos(p),lambda_p(i));
    end
end

for s = 1:1:length(signal.modos) % Mode overlap factor for signal entre modo y perfil de dopaje
    Nwl = length(signal.lambda.(ModoS(s)));
    for i=1:1:Nwl % Cada longitud de onda del modo s
        lambda_s = signal.lambda.(ModoS(s));
        [gamma_s.(ModoS(s)){i},beta0_s.(ModoS(s)){i}] = norm_intensity(fibra,signal.modos(s),lambda_s(i)) ; 
    end
end
warning('on')

% Power coupling coeficient pump
h_pccP = pcc_calc(beta0_p,ModoP); 
% Power coupling coeficient signal
h_pccS = pcc_calc(beta0_s,ModoS);

% ASE Spectrum
% Crear vector de longitudes de onda para ruido ase
allwavelengths = [];
for s = 1:1:Smod
    for w=1:1:length(signal.lambda.(ModoS(s)))
        if ~(ismember(signal.lambda.(ModoS(s))(w) , allwavelengths))
            allwavelengths = [allwavelengths signal.lambda.(ModoS(s))(w)];
        end
    end
end
allwavelengths = sort(allwavelengths);

if fibra.ASEFlag == 1  % retorna [0,0,0], evita calcular espectro ASE
    lambda_ase = ase_lambdas(allwavelengths,1);
else
    lambda_ase = ase_lambdas(allwavelengths);
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
    %Pan = [];                           % Potencia ASE en dirección -z
    Pase = [];                          % Potencia total ASE en el EDFA
    OSNR = [];                          % Relacion señal a ruido

    for s = 1:1:Smod
        Nwl = length(signal.lambda.(ModoS(s)));
        Pan.(ModoS(s)) = zeros(Nwl,Nz);
    end

    gain = [];

    Pap_sp = [];
    Pan_sp = [];
    ASE_Spectrum = [];

    % Calculo de densidades de iones en estado basal y excitado (N1 y N2)
    sig_xx = 0;
    sig_yy = 0;                    % Inicialización de variables
    ase_xx = 0;
    ase_yy = 0;
    pmp_xx = 0;
    pmp_yy = 0;

    z_waitbar = waitbar(0 ,'','Name',"Cálculo a lo largo del EDFA: " , 'Visible', 'off') ;
    QQ = 3; % Iteraciones para aumentar precisión
    for Q = 1:QQ
        if (Q == 1) % Primera iteracion no calcula ASE en direccion -z
            for z = 1:1:Nz  % Iteraciones a lo largo del EDFA
                %                 sig_xx = 0;
                %                 sig_yy = 0;                   % Inicialización de variables
                %                 ase_xx = 0;
                %                 ase_yy = 0;
                %                 pmp_xx = 0;
                %                 pmp_yy = 0;

                % % Eq19 Giles 1991 para N(z)
                % Usando los cambios de variable : sigma_a * gamma = alpha/nt ; A = pi*b_eff

                % Primera Iteración
                if(z == 1)
                    % Potencia Bombeo
                    parfor p = 1:Pmod
                        Nwlp = length(pump.lambda.(ModoP(p)));
                        for i=1:Nwlp
                            lambda_p = pump.lambda.(ModoP(p));
                            Gamma_p = gamma_p.(ModoP(p)){i};

                            v_p = c./lambda_p; Pp0 = P_p0.(ModoP(p))(i);

                            pmp_xx = pmp_xx + (sigma_abs(lambda_p(i)))*(Pp0*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));                                   % Término en numerador
                            pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Pp0*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));          % Término en denominador
                        end
                    end
                    % Potencia de señal
                    parfor s = 1:1:Smod
                        Nwl = length(signal.lambda.(ModoS(s)));
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            v_s = c./lambda_s;  Ps0 = P_s0.(ModoS(s))(i);

                            sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(Ps0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));                        % Término en numerador
                            sig_yy = sig_yy+(sigma_abs(lambda_s(i))+sigma_ems(lambda_s(i)))*(Ps0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Término en denominador

                        end
                    end

                    % Potencia ASE
                    parfor s=1:1:Smod
                        Nwl = length(signal.lambda.(ModoS(s)));
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            v_s = c./lambda_s;

                            ase_xx = ase_xx+ sigma_abs(lambda_s(i))*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));                           % Término en numerador
                            ase_yy = ase_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Término en denominador
                            
                        end
                    end

                else
                    % Iteraciones Restantes z>1
                    % Bombeo
                    parfor p = 1:1:Pmod
                        Nwlp = length(pump.lambda.(ModoP(p)));
                        for i = 1:1:Nwlp
                            lambda_p = pump.lambda.(ModoP(p));
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            v_p = c./lambda_p;

                            pmp_xx = pmp_xx + sigma_abs(lambda_p(i))*(Ppp.(ModoP(p))(i,z-1)*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));                        % Termino en numerador
                            pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Ppp.(ModoP(p))(i,z-1)*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));           % Termino en denominador
                        end
                    end

                    % Señal
                    parfor s = 1:1:Smod
                        Nwl = length(signal.lambda.(ModoS(s)));
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            v_s = c./lambda_s;

                            sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(Psp.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));
                            sig_yy = sig_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i))) * (Psp.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));
                        end
                    end

                    % ASE
                    parfor s = 1:1:Smod
                        Nwl = length(signal.lambda.(ModoS(s)));
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            v_s = c./lambda_s;

                            ase_xx = ase_xx+(sigma_abs(lambda_s(i))*(Pase.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));
                            ase_yy = ase_yy+( sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(Pase.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));
                        end
                    end
                end

                N2(z) = ( (pmp_xx + sig_xx + ase_xx) / (pmp_yy + sig_yy + ase_yy + (1/tau)) )*N;        % Densidad de iones de Erbio en estado excitado
                N1(z) = N-N2(z);                                                                        % Densidad de iones de Erbio en estado basal
                Nt(z) = N1(z) + N2(z);
                %         for i = 1:1:Nwl
                %             G(i,z) = G(i,z) + del_z*( sigma_ems(lambda_s(i))*N2(z) - sigma_abs(lambda_s(i))*N1(z) );
                %         end


                %% Ecuaciones de Potencias

                % Ecuacion diferencial para bombeo en direccion +z
                for p = 1:1:Pmod      % Iteracion en cada modo de bombeo
                    Nwlp = length(pump.lambda.(ModoP(p)));
                    if(z == 1)
                        for i = 1:1:Nwlp
                            lambda_p = pump.lambda.(ModoP(p));
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            Pp0 = P_p0.(ModoP(p))(i);
                            Ppp.(ModoP(p))(i,z) = Pp0 + (sum(N2(z,:))*sigma_ems(lambda_p(i)) - sum(N1(z,:))*sigma_abs(lambda_p(i)))*Gamma_p*Pp0*del_z;
                        end

                    else
                        for i = 1:1:Nwlp
                            lambda_p = pump.lambda.(ModoP(p));
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            Ppp.(ModoP(p))(i,z) = Ppp.(ModoP(p))(i,z-1) + (sum(N2(z,:))*sigma_ems(lambda_p(i)) - sum(N1(z,:))*sigma_abs(lambda_p(i)))*Gamma_p*Ppp.(ModoP(p))(i,z-1)*del_z;
                        end
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
                    Nwl = length(signal.lambda.(ModoS(s)));
                    if(z == 1)
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Ps0 = P_s0.(ModoS(s))(i);
                            Psp.(ModoS(s))(i,z) = Ps0 + (sum(N2(z,:))*sigma_ems(lambda_s(i))- sum(N1(z,:))*sigma_abs(lambda_s(i)))*Gamma_s*Ps0*del_z/(1+Ps0/Psat);
                        end
                    else
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Psp.(ModoS(s))(i,z) = Psp.(ModoS(s))(i,z-1) + (sum(N2(z,:))*sigma_ems(lambda_s(i))-sum(N1(z,:))*sigma_abs(lambda_s(i)))*Gamma_s*Psp.(ModoS(s))(i,z-1)*del_z/(1+Psp.(ModoS(s))(i,z-1)/Psat);
                        end
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
                for s = 1:1:Smod
                    Nwl = length(signal.lambda.(ModoS(s)));
                    if(z == 1)
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            v_s = c./lambda_s;
                            Pap.(ModoS(s))(i,z) = P_ase0 + ((sum(N2(z,:))*sigma_ems(lambda_s(i))-sum(N1(z,:))*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0 + 2*sum(N2(z,:))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+P_ase0/Psat);
                        end
                    else
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            v_s = c./lambda_s;
                            Pap.(ModoS(s))(i,z) = Pap.(ModoS(s))(i,z-1) + ((sum(N2(z,:))*sigma_ems(lambda_s(i))-sum(N1(z,:))*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(ModoS(s))(i,z-1) + 2*sum(N2(z,:))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pap.(ModoS(s))(i,z-1)/Psat);
                        end
                    end
                end

                % Ecuacion diferencial para ASE en direccion -z
                % No se calcula en 1era iteracion en fibra

                % Calculo de Potencia ASE
                for s = 1:1:Smod
                    Pase.(ModoS(s))(:,z) = Pap.(ModoS(s))(:,z)+Pan.(ModoS(s))(:,z);
                end
                % Mostrar avance como prints en pantalla:
                if fibra.Avance
                    if fibra.ASEFlag == 0
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 50  ) % Mostrar % de avance del cálculo
                    else
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 100  ) % Mostrar % de avance del cálculo
                    end
                end
                % Mostrar avance como WaitBar:
                if fibra.WaitBar
                    if fibra.ASEFlag == 0
                        avance = ((z/Nz)*(1/QQ) + ((Q-1)/QQ))/2 ; waitbar(avance, z_waitbar ,  sprintf("Cálculo a lo largo del EDFA: \n %.2f %%",avance*100 ) ) ; set(z_waitbar,'Visible', 'on');
                    else
                        avance = (z/Nz)*(1/QQ) + ((Q-1)/QQ) ; waitbar(avance, z_waitbar ,  sprintf("Cálculo a lo largo del EDFA: \n %.2f %%",avance*100 ) ) ; set(z_waitbar,'Visible', 'on');
                    end
                end
            end   % Fin de calculo en fibra

        else  % Fin 1era iteración, usar estos valores de N para una 2da iteración mas precisa
            % Agregando ASE en dirección -z

            for z = 1:1:Nz %Iteraciones a lo largo del EDFA
                % Siguientes Primeras Iteraciones
                if(z == 1)
                    % Potencia Bombeo
                    parfor p = 1:1:Pmod
                        Nwlp = length(pump.lambda.(ModoP(p)));
                        for i=1:1:Nwlp
                            lambda_p = pump.lambda.(ModoP(p));
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            v_p = c./lambda_p; Pp0 = P_p0.(ModoP(p))(i);

                            pmp_xx = pmp_xx + (sigma_abs(lambda_p(i)))*(Pp0*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));                                   % Termino en numerador
                            pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Pp0*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));          % Termino en denominador
                        end
                    end
                    % Potencia de señal
                    parfor s = 1:1:Smod
                        Nwl = length(signal.lambda.(ModoS(s)));
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            v_s = c./lambda_s;  Ps0 = P_s0.(ModoS(s))(i);

                            sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(Ps0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));                        % Termino en numerador
                            sig_yy = sig_yy+(sigma_abs(lambda_s(i))+sigma_ems(lambda_s(i)))*(Ps0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Termino en denominador
                        end
                    end
                    % Potencia ASE
                    parfor s=1:1:Smod
                        Nwl = length(signal.lambda.(ModoS(s)));
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            v_s = c./lambda_s;
                            Gamma_s = gamma_s.(ModoS(s)){i};

                            ase_xx = ase_xx+ sigma_abs(lambda_s(i))*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));                           % Termino en numerador
                            ase_yy = ase_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Termino en denominador
                        end
                    end

                else
                    % Iteraciones Restantes z>1
                    % Bombeo
                    parfor p = 1:1:Pmod
                        Nwlp = length(pump.lambda.(ModoP(p)));
                        for i = 1:1:Nwlp
                            lambda_p = pump.lambda.(ModoP(p));
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            v_p = c./lambda_p;

                            pmp_xx = pmp_xx + sigma_abs(lambda_p(i))*(Ppp.(ModoP(p))(i,z-1)*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));                        % Termino en numerador
                            pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Ppp.(ModoP(p))(i,z-1)*Gamma_p/A_p.(ModoP(p)))/(h*v_p(i));           % Termino en denominador
                        end
                    end

                    % Señal
                    parfor s = 1:1:Smod
                        Nwl = length(signal.lambda.(ModoS(s)));
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            v_s = c./lambda_s;
                            Gamma_s = gamma_s.(ModoS(s)){i};

                            sig_xx = sig_xx + (sigma_abs(lambda_s(i))*(Psp.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));
                            sig_yy = sig_yy + (sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i))) * (Psp.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));
                        end
                    end

                    % ASE
                    parfor s = 1:1:Smod
                        Nwl = length(signal.lambda.(ModoS(s)));
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            v_s = c./lambda_s;
                            Gamma_s = gamma_s.(ModoS(s)){i};

                            ase_xx = ase_xx + (sigma_abs(lambda_s(i))*(Pase.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)));
                            ase_yy = ase_yy + ( sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(Pase.(ModoS(s))(i,z-1)*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));
                        end
                    end
                end

                N2(z) = ( (pmp_xx + sig_xx + ase_xx)/(pmp_yy + sig_yy + ase_yy + (1/tau)) )*N;             % Densidad de iones de Erbio en estado excitado
                N1(z) = N-N2(z);                                                               % Densidad de iones de Erbio en estado basal
                Nt(z) = N1(z) + N2(z);

                % % Ecuaciones de Potencias

                % Ecuacion diferencial para bombeo en direccion +z
                for p = 1:1:Pmod      % Iteracion en cada modo de bombeo
                    Nwlp = length(pump.lambda.(ModoP(p)));
                    if(z == 1)
                        for i = 1:1:Nwlp
                            lambda_p = pump.lambda.(ModoP(p));
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            Pp0 = P_p0.(ModoP(p))(i);
                            Ppp.(ModoP(p))(i,z) = Pp0 + (sum(N2(z,:))*sigma_ems(lambda_p(i)) - sum(N1(z,:))*sigma_abs(lambda_p(i)))*Gamma_p*Pp0*del_z;
                        end
                    else
                        for i = 1:1:Nwlp
                            lambda_p = pump.lambda.(ModoP(p));
                            Gamma_p = gamma_p.(ModoP(p)){i};
                            Ppp.(ModoP(p))(i,z) = Ppp.(ModoP(p))(i,z-1) + (sum(N2(z,:))*sigma_ems(lambda_p(i)) - sum(N1(z,:))*sigma_abs(lambda_p(i)))*Gamma_p*Ppp.(ModoP(p))(i,z-1)*del_z;
                        end
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
                    Nwl = length(signal.lambda.(ModoS(s)));
                    if(z == 1)
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Ps0 = P_s0.(ModoS(s))(i);
                            Psp.(ModoS(s))(i,z) = Ps0 + (sum(N2(z,:))*sigma_ems(lambda_s(i)) - sum(N1(z,:))*sigma_abs(lambda_s(i)))*Gamma_s*Ps0*del_z/(1+Ps0/Psat);
                        end
                    else
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Psp.(ModoS(s))(i,z) = Psp.(ModoS(s))(i,z-1) + (sum(N2(z,:))*sigma_ems(lambda_s(i))-sum(N1(z,:))*sigma_abs(lambda_s(i)))*Gamma_s*Psp.(ModoS(s))(i,z-1)*del_z/(1+Psp.(ModoS(s))(i,z-1)/Psat);
                        end
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
                for s = 1:1:Smod
                    Nwl = length(signal.lambda.(ModoS(s)));
                    if(z == 1)
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            v_s = c./lambda_s;
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pap.(ModoS(s))(i,z) = P_ase0 + ((sum(N2(z,:))*sigma_ems(lambda_s(i)) - sum(N1(z,:))*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0 + 2*sum(N2(z,:))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+P_ase0/Psat);
                        end
                    else
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            v_s = c./lambda_s;
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pap.(ModoS(s))(i,z) = Pap.(ModoS(s))(i,z-1) + ((sum(N2(z,:))*sigma_ems(lambda_s(i)) - sum(N1(z,:))*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(ModoS(s))(i,z-1) + 2*sum(N2(z,:))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pap.(ModoS(s))(i,z-1)/Psat);
                        end
                    end
                end

                % Ecuacion diferencial para ASE en direccion -z
                for s = 1:1:Smod
                    Nwl = length(signal.lambda.(ModoS(s)));
                    if(z == 1)
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            v_s = c./lambda_s;
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pan.(ModoS(s))(i,Nz-z+1) = 0+((sum(N2(Nz-z+1,:))*sigma_ems(lambda_s(i))-sum(N1(Nz-z+1,:))*sigma_abs(lambda_s(i)))*Gamma_s*0 + 2*sum(N2(Nz-z+1,:))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                        end
                    else
                        for i = 1:1:Nwl
                            lambda_s = signal.lambda.(ModoS(s));
                            v_s = c./lambda_s;
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pan.(ModoS(s))(i,Nz-z+1) = Pan.(ModoS(s))(i,Nz-z+1+1) + ((sum(N2(Nz-z+1,:))*sigma_ems(lambda_s(i)) - sum(N1(Nz-z+1,:))*sigma_abs(lambda_s(i)))*Gamma_s*Pan.(ModoS(s))(i,Nz-z+1+1) + 2*sum(N2(Nz-z+1,:))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pan.(ModoS(s))(i,Nz-z+1+1)/Psat);
                        end
                    end
                end
                % Mostrar % de avance como prints en pantalla:
                if fibra.Avance
                    if fibra.ASEFlag == 0
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ)* 50 + ((Q-1)/QQ) * 50 ) % Mostrar % de avance del cálculo
                    else
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 100 ) % Mostrar % de avance del cálculo
                    end
                    
                end
                % Mostrar % de avance como WaitBar:
                if fibra.WaitBar
                    if fibra.ASEFlag == 0
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

    for z = 1:Nz
        for s = 1:1:Smod
            Nwl = length(signal.lambda.(ModoS(s)));
            for w = 1:1:Nwl
                Nch_ase = length(lambda_ase);
                P_ase0_sp = P_ase0;
                v_s_sp = c./lambda_ase;
                Gamma_s = gamma_s.(ModoS(s)){w};

                if(z == 1)
                    for i = 1:1:Nch_ase
                        Pap_sp.(ModoS(s))(i,z) = P_ase0_sp + ((sum(N2(z,:))*sigma_ems(lambda_ase(i)) - sum(N1(z,:))*sigma_abs(lambda_ase(i)))*Gamma_s*P_ase0_sp + 2*sum(N2(z,:))*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+P_ase0_sp/Psat);
                    end
                else
                    for i = 1:1:Nch_ase
                        Pap_sp.(ModoS(s))(i,z) = Pap_sp.(ModoS(s))(i,z-1) + ((sum(N2(z,:))*sigma_ems(lambda_ase(i)) - sum(N1(z,:))*sigma_abs(lambda_ase(i)))*Gamma_s*Pap_sp.(ModoS(s))(i,z-1) + 2*sum(N2(z,:))*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+Pap_sp.(ModoS(s))(i,z-1)/Psat);
                    end
                end
            end
        end

        % Ecuacion diferencial para ASE en direccion -z
        for s = 1:1:Smod
            if(z == 1)
                for i = 1:1:Nch_ase
                    Pan_sp.(ModoS(s))(i,Nz-z+1) = 0+((sum(N2(Nz-z+1,:))*sigma_ems(lambda_ase(i)) - sum(N1(Nz-z+1,:))*sigma_abs(lambda_ase(i)))*Gamma_s*0 + 2*sum(N2(Nz-z+1,:))*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z;
                end
            else
                for i = 1:1:Nch_ase
                    Pan_sp.(ModoS(s))(i,Nz-z+1) = Pan_sp.(ModoS(s))(i,Nz-z+1+1) + ((sum(N2(Nz-z+1,:))*sigma_ems(lambda_ase(i)) - sum(N1(Nz-z+1,:))*sigma_abs(lambda_ase(i)))*Gamma_s*Pan_sp.(ModoS(s))(i,Nz-z+1+1) + 2*sum(N2(Nz-z+1,:))*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+Pan_sp.(ModoS(s))(i,Nz-z+1+1)/Psat);
                end
            end
        end
        % Mostrar % de avance como prints en pantalla:
        if fibra.Avance && fibra.ASEFlag == 0
            clc ; fprintf("Cálculo Espectro ASE: \n") ; fprintf('%.2f %% \n' ,  ( 50 + ((z/Nz)*50 )) ) % Mostrar % de avance del cálculo
        end
        % Mostrar % de avance como WaitBar:
        if fibra.WaitBar && fibra.ASEFlag == 0
            avance = 0.50 + ((z/Nz))/2 ; waitbar(avance , z_waitbar ,  sprintf("Cálculo Espectro ASE: \n %.2f %%",avance*100 ) ) ; set(z_waitbar,'Visible', 'on');
        end
    end % Fin calculo ASE en posicion z


    for s = 1:1:Smod
        for i = 1:1:Nch_ase
            ASE_Spectrum.(ModoS(s))(i,1) = 10*log10( (Pap_sp.(ModoS(s))(i,end) + Pan_sp.(ModoS(s))(i,end)));
        end
    end

    if fibra.WaitBar
        close(z_waitbar)
    end



    %%

    % Calculo OSNR
    for s = 1:1:Smod
        Nwl = length(signal.lambda.(ModoS(s)));
        for i = 1:1:Nwl
            OSNR.(ModoS(s))(i,:) = 10*log10(Psp.(ModoS(s))(i,:)/1e-3) - 10*log10(Pase.(ModoS(s))(i,:)/1e-3);
            gain.(ModoS(s))(1,i) = 10*log10(Psp.(ModoS(s))(i,end)/Psp.(ModoS(s))(i,1));
        end
    end

    % CALCULAR GRAFICO ASE SPECTRUM CON GANANCIAS DE SEÑAL
    for s = 1:1:Smod
        ASE_Spectrum.(ModoS(s))(:,2) = ASE_Spectrum.(ModoS(s))(:,1);

        lambda_s = signal.lambda.(ModoS(s));
        for w = 1:length(lambda_s)
            for i = 1:length(lambda_ase)
                if abs(lambda_ase(i)*1e9 - lambda_s(w)*1e9) < 1/11
                    %ASE_Spectrum(i,2,s) = gain(1,w,s);
                    ASE_Spectrum.(ModoS(s))(i,2) = 10*log10( (1e-3*10.^(ASE_Spectrum.(ModoS(s))(i,1)/10) + Psp.(ModoS(s))(w,end)) ./1e-3);
                    %ASE_Spectrum(i,2,s) = 10*log10( (( Psp(w,end,s)) ./1e-3) );
                end
            end
        end
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
        sdm.(ch).salida.ASE.potencia_dBm.(ModoS(s)) = 10*log10(Pase.(ModoS(s))(:,end)./1e-3);

        sdm.(ch).salida.ganancias.(ModoS(s)) = gain.(ModoS(s));
        sdm.(ch).NF.(ModoS(s)) = OSNR.(ModoS(s))(:,1)./OSNR.(ModoS(s))(:,end);
    end

    sdm.(ch).signal.lambdas = lambda_s ; sdm.(ch).pump.lambdas = lambda_p;
    sdm.(ch).OSNR = OSNR ;
    sdm.(ch).salida.OSNR = OSNR(:,end);
    sdm.(ch).z = Z;
    sdm.(ch).ASE_Spectrum.mag = ASE_Spectrum;
    sdm.(ch).ASE_Spectrum.lambdas = lambda_ase;
    sdm.(ch).Gamma.Gamma_p = gamma_p;
    sdm.(ch).Gamma.Gamma_s = gamma_s;
    

end % Fin cálculo en todos los nucleos

edfa = sdm;

end % Fin funcion

%% Pendiente
% Revisar potencia de saturacion y Area efectiva
% Superposición entre señal/bombeo y perfil de dopaje es muy bajo
