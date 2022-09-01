function edfa = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE)
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

d_vk = 125*10^9;                        % 1 nm - Noise Bandwidth
tau = 10e-3;
nucleos = Fibra.nucleos;
Sch = nucleos;                          % N° de nucleos
N = Fibra.N ;

% Datos de Bombeo y Señal

P_s0 = Signal.P0;
P_p0 = Pump.P0;

ModoS = strcat("LP_" , Signal.modos(:));
ModoP = strcat("LP_" , Pump.modos(:));

Smod = length(ModoS);                   % N° de modos de señal
Pmod = length(ModoP);                   % N° de modos de bombeo

P_ase0 = 1e-3*10.^(ASE/10);             % Potencia ASE entrada en Watts

% %   Espectros emision y absorción

%       % % Datos obtenidos de VPI
% VPI = load('Erbium_VPI.dat');
% Sa = VPI(:,3); Se = VPI(:,2);
% lambda_cross = VPI(:,1).*1e-9;

% %       % % Datos OptiSystem
OptiSystem = load('Erbium_OptiSystem.dat');
Sa = OptiSystem(:,2); Se = OptiSystem(:,3);
lambda_cross = OptiSystem(:,1).*1e-9;


sigma_abs = fit(lambda_cross,Sa,'linearinterp');
sigma_ems = fit(lambda_cross,Se,'linearinterp');


for i=1:length(ModoS)
    A_s.(ModoS(i)) = pi*Fibra.radio^2;                                % Area efectiva para señal
end
for i=1:length(ModoP)
    A_p.(ModoP(i)) = pi*Fibra.radio^2 ;                               % Area efectiva para bombeo
end

% De Paper ó VPIphotonics 
% Para radio=10e-6 y AN=0.2
%A_s.LP_01=80e-12; A_s.LP_11_a=76e-12; A_s.LP_11_b=76e-12; A_s.LP_02=83e-12; A_s.LP_21_a=86e-12; A_s.LP_21_b=86e-12;
%A_p.LP_01=80e-12; A_p.LP_11_a=76e-12; A_p.LP_11_b=76e-12; A_p.LP_02=83e-12; A_p.LP_21=86e-12;
% Para radio=6e-6 y AN=0.2
% A_s.LP_01=80e-12; A_s.LP_11_a=76e-12; A_s.LP_11_b=76e-12; A_s.LP_02=83e-12; A_s.LP_21_a=86e-12; A_s.LP_21_b=86e-12;
% A_p.LP_01=80e-12; A_p.LP_11_a=76e-12; A_p.LP_11_b=76e-12; A_p.LP_02=83e-12; A_p.LP_21=86e-12;
sdm = struct;

% % Calculo de Gammas
% warning('off')
% 
% for p = 1:1:length(Pump.modos) % Mode overlap factor for pump; entre modo y perfil de dopaje (uniforme)
%     Nwlp = length(Pump.lambda.(ModoP(p)));
%     for i=1:1:Nwlp % Cada longitud de onda del modo p
%         lambda_p = Pump.lambda.(ModoP(p));
%         [gamma_p.(ModoP(p)){i},beta0_p.(ModoP(p)){i}] = norm_intensity(Fibra,Pump.modos(p),lambda_p(i));
%     end
% end
% 
% for s = 1:1:length(Signal.modos) % Mode overlap factor for signal entre modo y perfil de dopaje
%     Nwl = length(Signal.lambda.(ModoS(s)));
%     for i=1:1:Nwl % Cada longitud de onda del modo s
%         lambda_s = Signal.lambda.(ModoS(s));
%         [gamma_s.(ModoS(s)){i},beta0_s.(ModoS(s)){i}] = norm_intensity(Fibra,Signal.modos(s),lambda_s(i)) ; 
%     end
% end
% warning('on')

% Calculo de Gammas
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


% Power coupling coeficient pump
h_pccP = pcc_calc(beta0_p,ModoP); 
% Power coupling coeficient signal
h_pccS = pcc_calc(beta0_s,ModoS);

% ASE Spectrum
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
    lambda_ase = ase_lambdas(allwavelengths,1);
else
    %lambda_ase = ase_lambdas(allwavelengths);
    lambda_ase = 1520e-9 : 1*1e-9: 1600e-9;
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
        Nwl = length(Signal.lambda.(ModoS(s)));
        Pan.(ModoS(s)) = zeros(Nwl,Nz);
    end

    gain = [];

    Pap_sp = [];
    Pan_sp = [];
    ASE_Spectrum = [];

    % Calculo de densidades de iones en estado basal y excitado (N1 y N2)

    z_waitbar = waitbar(0 ,'','Name',"Cálculo a lo largo del EDFA: " , 'Visible', 'off') ;
    QQ = 4; % Iteraciones para aumentar precisión
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

                            ase_xx = ase_xx+ sigma_abs(lambda_s(i))*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));                           % Término en numerador
                            ase_yy = ase_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Término en denominador
                            
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



                %% Ecuaciones de Potencias

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
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    lambda_s = Signal.lambda.(ModoS(s));
                    v_s = c./lambda_s;
                    if(z == 1)
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux(i) = P_ase0 
                        end
                        Pap.(ModoS(s))(:,z) = PapAux(:);
                    else
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux(i) = Pap.(ModoS(s))(i,z-1) + (((N2(z))*sigma_ems(lambda_s(i))-(N1(z))*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(ModoS(s))(i,z-1) + 2*(N2(z))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                        end
                        Pap.(ModoS(s))(:,z) = PapAux(:);
                    end
                end

                % Ecuacion diferencial para ASE en direccion -z
                % No se calcula en 1era iteracion en fibra

                % Calculo de Potencia ASE
                for s = 1:1:Smod
                    Pase.(ModoS(s))(:,z) = Pap.(ModoS(s))(:,z);%+Pan.(ModoS(s))(:,z);
                end
                % Mostrar avance como prints en pantalla:
                if Fibra.Avance
                    if Fibra.ASEFlag == 0
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 50  ) % Mostrar % de avance del cálculo
                    else
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 100  ) % Mostrar % de avance del cálculo
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

        elseif Q < QQ  % Fin 1era iteración, usar estos valores de N para una 2da iteración mas precisa
            % Agregando ASE en dirección -z

            for z = 1:1:Nz %Iteraciones a lo largo del EDFA
                sig_xx = 0;
                sig_yy = 0;                    % Inicialización de variables
                ase_xx = 0;
                ase_yy = 0;
                pmp_xx = 0;
                pmp_yy = 0;

%                 % GEF - Filtro de equalización
% 
%                 if (z == round(Nz/2,0) && Q == QQ) % Centro de la fibra
%                     for s = 1:1:Smod
%                         offsetPot.(ModoS(s)) = Psp.(ModoS(s))(:,z-1) - 0.8.*min( Psp.(ModoS(s))(:,z-1) );    % Curva trasladada a cero
%                         maxDiffPot.(ModoS(s)) = max( Psp.(ModoS(s))(:,z-1) ) - 0.8.*min( Psp.(ModoS(s))(:,z-1) ); % "Amplitud"
%                         normPot.(ModoS(s)) = offsetPot.(ModoS(s))./maxDiffPot.(ModoS(s)) ; % curva en cero normalizada
%                         
%                         weigth = 0.4;
%                         Psp.(ModoS(s))(:,z-1) = Psp.(ModoS(s))(:,z-1) .* (1-weigth*normPot.(ModoS(s)));
%                         sdm.(ch).GEF.(ModoS(s)) = 1-weigth*normPot.(ModoS(s)) ;
%                     end
%                 end

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

                            ase_xx = ase_xx+ sigma_abs(lambda_s(i))*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));                           % Termino en numerador
                            ase_yy = ase_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Termino en denominador
                        end
                    end

                else
                    % Iteraciones Restantes z>1
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
                N1(z) = N-N2(z);                                                               % Densidad de iones de Erbio en estado basal
                Nt(z) = N1(z) + N2(z);

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
                for s = 1:1:Smod
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                    if(z == 1)
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux(i) = P_ase0 %+ (((N2(z))*sigma_ems(lambda_s(i)) - (N1(z))*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0 + 2*(N2(z))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+P_ase0/Psat);
                        end
                        Pap.(ModoS(s))(:,z) = PapAux(:);
                    else
                        parfor i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            PapAux(i) = Pap.(ModoS(s))(i,z-1) + (((N2(z))*sigma_ems(lambda_s(i)) - (N1(z))*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(ModoS(s))(i,z-1) + 2*(N2(z))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                        end
                        Pap.(ModoS(s))(:,z) = PapAux(:);
                    end
                end

                % Ecuacion diferencial para ASE en direccion -z
                for s = 1:1:Smod
                    Nwl = length(Signal.lambda.(ModoS(s)));
                    lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                    if(z == 1)
                        for i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pan.(ModoS(s))(i,Nz-z+1) = P_ase0+(((N2(Nz-z+1))*sigma_ems(lambda_s(i))-(N1(Nz-z+1))*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0 + 2*(N2(Nz-z+1))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                            %PanAux(i) = P_ase0+(((N2(Nz-z+1))*sigma_ems(lambda_s(i))-(N1(Nz-z+1))*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0 + 2*(N2(Nz-z+1))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                        end
                        %Pan.(ModoS(s))(:,Nz-z+1) = PanAux(:);
                    else
                        for i = 1:Nwl
                            Gamma_s = gamma_s.(ModoS(s)){i};
                            Pan.(ModoS(s))(i,Nz-z+1) = Pan.(ModoS(s))(i,Nz-z+1+1) + (((N2(Nz-z+1))*sigma_ems(lambda_s(i)) - (N1(Nz-z+1))*sigma_abs(lambda_s(i)))*Gamma_s*Pan.(ModoS(s))(i,Nz-z+1+1) + 2*(N2(Nz-z+1))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                            %PanAux(i) = Pan.(ModoS(s))(i,Nz-z+1+1) + (((N2(Nz-z+1))*sigma_ems(lambda_s(i)) - (N1(Nz-z+1))*sigma_abs(lambda_s(i)))*Gamma_s*Pan.(ModoS(s))(i,Nz-z+1+1) + 2*(N2(Nz-z+1))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pan.(ModoS(s))(i,Nz-z+1+1)/Psat);
                        end
                        %Pan.(ModoS(s))(:,Nz-z+1) = PanAux(:);
                    end
                end

                for s = 1:1:Smod
                    Pase.(ModoS(s))(:,z) = Pap.(ModoS(s))(:,z);%+Pan.(ModoS(s))(:,z);
                end

                % Mostrar % de avance como prints en pantalla:
                if Fibra.Avance
                    if Fibra.ASEFlag == 0
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ)* 50 + ((Q-1)/QQ) * 50 ) % Mostrar % de avance del cálculo
                    else
                        clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% \n' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 100 ) % Mostrar % de avance del cálculo
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

        % Guardado de ganancias sin GEF
        for s = 1:1:Smod
            Nwl = length(Signal.lambda.(ModoS(s)));
            for i = 1:1:Nwl
                gain_sinGEF.(ModoS(s))(1,i) = 10*log10(Psp.(ModoS(s))(i,end)/Psp.(ModoS(s))(i,1));
            end
        end
        
        sdm.(ch).salida.ganancias_sinGEF.(ModoS(s)) = gain_sinGEF.(ModoS(s));
        
        % Ultima iteracion - GEF
        InitialPsp = Psp;
        InitialPap = Pap;

        elseif Q == QQ  
            initial_flag = 0;
            if initial_flag == 0
                InitialDiffPot = max( Psp.(ModoS(s))(:,Nz ) ) - min( Psp.(ModoS(s))(:,Nz ) );
                ActualDiffPot = max( Psp.(ModoS(s))(:,Nz ) ) - min( Psp.(ModoS(s))(:,Nz ) );
                minGain = 0.8*min( Psp.(ModoS(s))(:,Nz ) );
                Logro = [ActualDiffPot/InitialDiffPot , 1 ] ; 
            end 
            
            while_count = 0; weight_change_flag = 0; break_flag = 0;

            while ActualDiffPot>0.15*InitialDiffPot
                Psp = InitialPsp; Pap = InitialPap;
                while_count = while_count +1;
                
                Logro(2) = Logro(1) ; Logro(1) = ActualDiffPot/InitialDiffPot  ; 
                
                if Logro(1) > Logro(2) % Empeora el resultado
                    normPot = best_normPot;
                    weight = weight+0.1 ; weight_change_flag = 1;
                    if weight >=0.95
                        sprintf('mejor resultado posible...')
                        normPot = best_normPot; weight = best_weight;
                        break_flag = 1;
                        %break
                    end
                    %sprintf('No hay mejora, modificando weight...')
                
                else % mejora el resultado
                    weight_change_flag = 0;
                    if while_count>1
                        best_normPot = normPot; best_weight = weight;
                    end
                end

                for z = 1:1:Nz %Iteraciones a lo largo del EDFA  
                    sig_xx = 0;
                    sig_yy = 0;                    % Inicialización de variables
                    ase_xx = 0;
                    ase_yy = 0;
                    pmp_xx = 0;
                    pmp_yy = 0;
    
                    % GEF - Filtro de equalización
                    
                    if (z == round(Nz/2,0) && Q == QQ) % Centro de la fibra
                        if while_count == 1
                            offsetPot = Psp.(ModoS(1))(:,Nz) - min( Psp.(ModoS(1))(:,Nz) );    % Curva trasladada a cero
                            maxDiffPot = max( Psp.(ModoS(1))(:,Nz) ) - min( Psp.(ModoS(1))(:,Nz) ); % "Amplitud"
                            normPot = offsetPot./maxDiffPot ; % curva en cero normalizada
                            for s = 1:1:Smod  
                                weight = 0.25; % 0.35 para OptiSystem ; 0.85 VPI / 5m largo 250mw
                                Psp.(ModoS(s))(:,z-1) = Psp.(ModoS(s))(:,z-1) .* ( 1-weight*normPot );
                                Pap.(ModoS(s))(:,z-1) = Pap.(ModoS(s))(:,z-1) .* ( 1-weight*normPot );
                            end
                            sdm.(ch).GEF.Iteracion1 = weight*normPot ;
                            best_normPot = normPot; best_weight = weight;

                        else
                            ajuste = 0.03;
                            if while_count > 3
                                ajuste = 0.01;
                            elseif while_count > 20
                                ajuste = 0.005;
%                             elseif while_count > 30
%                                 ajuste = 0.0001;
                            end

                            Nwl = length(Signal.lambda.(ModoS(1)));

                            % Ajustar Valores del polinomio (filtro)
                            if break_flag==0 % Si no se encontró mejor resultado usa el mejor polinomio guardado
                                for f=1:Nwl
                                    if Psp.(ModoS(1))(f,Nz) > minGain
                                        normPot(f) = normPot(f)-ajuste;
                                        if normPot(f)>1
                                            normPot(f) = 1;
                                        elseif normPot(f)<0
                                            normPot(f) = 0;
                                        end
                                    elseif Psp.(ModoS(1))(f,Nz) < minGain
                                        normPot(f) = normPot(f)+ajuste;
                                        if normPot(f)>1
                                            normPot(f)=1;
                                        elseif normPot(f)<0
                                            normPot(f) = 0;
                                        end
                                    end
                                end
                            end
                            sdm.(ch).GEF.(strcat('Iteracion',num2str(while_count))) = weight*normPot ;
                            % Aplicar filtro
                            for s = 1:1:Smod  
                                Psp.(ModoS(s))(:,z-1) = Psp.(ModoS(s))(:,z-1) .* ( 1-weight*normPot );
                                Pap.(ModoS(s))(:,z-1) = Pap.(ModoS(s))(:,z-1) .* ( 1-weight*normPot );
                            end
                        end
                    end
    
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
    
                                ase_xx = ase_xx+ sigma_abs(lambda_s(i))*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i));                           % Termino en numerador
                                ase_yy = ase_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(P_ase0*Gamma_s/A_s.(ModoS(s)))/(h*v_s(i)); % Termino en denominador
                            end
                        end
    
                    else
                        % Iteraciones Restantes z>1
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
                    N1(z) = N-N2(z);                                                               % Densidad de iones de Erbio en estado basal
                    Nt(z) = N1(z) + N2(z);
    
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
                    for s = 1:1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        if(z == 1)
                            parfor i = 1:Nwl
                                Gamma_s = gamma_s.(ModoS(s)){i};
                                PapAux(i) = P_ase0 %+ (((N2(z))*sigma_ems(lambda_s(i)) - (N1(z))*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0 + 2*(N2(z))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+P_ase0/Psat);
                            end
                            Pap.(ModoS(s))(:,z) = PapAux(:);
                        else
                            parfor i = 1:Nwl
                                Gamma_s = gamma_s.(ModoS(s)){i};
                                PapAux(i) = Pap.(ModoS(s))(i,z-1) + (((N2(z))*sigma_ems(lambda_s(i)) - (N1(z))*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(ModoS(s))(i,z-1) + 2*(N2(z))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                            end
                            Pap.(ModoS(s))(:,z) = PapAux(:);
                        end
                    end
    
                    % Ecuacion diferencial para ASE en direccion -z
                    for s = 1:1:Smod
                        Nwl = length(Signal.lambda.(ModoS(s)));
                        lambda_s = Signal.lambda.(ModoS(s)); v_s = c./lambda_s;
                        if(z == 1)
                            for i = 1:Nwl
                                Gamma_s = gamma_s.(ModoS(s)){i};
                                Pan.(ModoS(s))(i,Nz-z+1) = P_ase0+(((N2(Nz-z+1))*sigma_ems(lambda_s(i))-(N1(Nz-z+1))*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0 + 2*(N2(Nz-z+1))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                                %PanAux(i) = P_ase0+(((N2(Nz-z+1))*sigma_ems(lambda_s(i))-(N1(Nz-z+1))*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0 + 2*(N2(Nz-z+1))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                            end
                            %Pan.(ModoS(s))(:,Nz-z+1) = PanAux(:);
                        else
                            for i = 1:Nwl
                                Gamma_s = gamma_s.(ModoS(s)){i};
                                Pan.(ModoS(s))(i,Nz-z+1) = Pan.(ModoS(s))(i,Nz-z+1+1) + (((N2(Nz-z+1))*sigma_ems(lambda_s(i)) - (N1(Nz-z+1))*sigma_abs(lambda_s(i)))*Gamma_s*Pan.(ModoS(s))(i,Nz-z+1+1) + 2*(N2(Nz-z+1))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                                %PanAux(i) = Pan.(ModoS(s))(i,Nz-z+1+1) + (((N2(Nz-z+1))*sigma_ems(lambda_s(i)) - (N1(Nz-z+1))*sigma_abs(lambda_s(i)))*Gamma_s*Pan.(ModoS(s))(i,Nz-z+1+1) + 2*(N2(Nz-z+1))*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pan.(ModoS(s))(i,Nz-z+1+1)/Psat);
                            end
                            %Pan.(ModoS(s))(:,Nz-z+1) = PanAux(:);
                        end
                    end

                    for s = 1:1:Smod
                        Pase.(ModoS(s))(:,z) = Pap.(ModoS(s))(:,z);%+Pan.(ModoS(s))(:,z);
                    end

                    % Mostrar % de avance como prints en pantalla:
                    if Fibra.Avance
                        if Fibra.ASEFlag == 0
                            if weight_change_flag==1 && break_flag==0
                                clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% , Iteracion GEF: %.0f \n Actual/Initial = %.2f \n No hay mejora, modificando weight... \n' ,  (z/Nz)*(1/QQ)* 50 + ((Q-1)/QQ) * 50 , while_count , ActualDiffPot/InitialDiffPot) % Mostrar % de avance del cálculo
                            elseif break_flag==1
                                clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% , Iteracion GEF: %.0f \n  Actual/Initial = %.2f \n Utilizando mejor filtro encontrado... \n' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 100 , while_count , ActualDiffPot/InitialDiffPot) % Mostrar % de avance del cálculo
                            else
                                clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% , Iteracion GEF: %.0f \n Actual/Initial = %.2f \n ' ,  (z/Nz)*(1/QQ)* 50 + ((Q-1)/QQ) * 50 , while_count , ActualDiffPot/InitialDiffPot) % Mostrar % de avance del cálculo
                            end
                        else
                            if weight_change_flag==1 && break_flag==0
                                clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% , Iteracion GEF: %.0f \n  Actual/Initial = %.2f \n No hay mejora, modificando weight... \n' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 100 , while_count , ActualDiffPot/InitialDiffPot) % Mostrar % de avance del cálculo
                            elseif break_flag==1
                                clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% , Iteracion GEF: %.0f \n  Actual/Initial = %.2f \n Utilizando mejor filtro encontrado... \n' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 100 , while_count , ActualDiffPot/InitialDiffPot) % Mostrar % de avance del cálculo
                            else
                                clc ; fprintf("Cálculo a lo largo del EDFA: \n") ; fprintf('%.2f %% , Iteracion GEF: %.0f \n  Actual/Initial = %.2f \n ' ,  (z/Nz)*(1/QQ)* 100 + ((Q-1)/QQ) * 100 , while_count , ActualDiffPot/InitialDiffPot) % Mostrar % de avance del cálculo
                            end
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
                ActualDiffPot = max( Psp.(ModoS(s))(:,Nz) ) - min( Psp.(ModoS(s))(:,Nz) );
                if break_flag==1
                    break
                end
            end % END WHILE GEF

        end
        
    end % Fin iteraciones para estabilizar ganancias

    %% ASE SPECTRUM :
    % Ecuacion diferencial para ASE en direccion +z
    v_s_sp = c./lambda_ase;
    Nch_ase = length(lambda_ase);
    d_vk = abs(c/lambda_ase(1) - c/lambda_ase(2));
    if ~(Fibra.ASEFlag == 1)
        for z = 1:Nz
            for s = 1:1:Smod
                P_ase0_sp = P_ase0;
                Gamma_s = sum(cell2mat(gamma_s.(ModoS(s))))/length(cell2mat(gamma_s.(ModoS(s))));
                if(z == 1)
                    parfor i = 1:Nch_ase
                        AuxMatrix(i) = P_ase0_sp;
                    end
                    Pap_sp.(ModoS(s))(:,z) = AuxMatrix(:) ;
        
                else
                    parfor i = 1:Nch_ase
                        if (i==1)
                            d_vk=abs(v_s_sp(i)-v_s_sp(i+1));
                        else
                            d_vk=abs(v_s_sp(i-1)-v_s_sp(i));
                        end
                        AuxMatrix(i) = Pap_sp.(ModoS(s))(i,z-1) + (((N2(z))*sigma_ems(lambda_ase(i)) - (N1(z))*sigma_abs(lambda_ase(i)))*Gamma_s*Pap_sp.(ModoS(s))(i,z-1) + 2*(N2(z))*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z;
                    end
                    Pap_sp.(ModoS(s))(:,z) = AuxMatrix(:) ;
                end
            end
        
            if Fibra.Avance && Fibra.ASEFlag == 0
                clc ; fprintf("Cálculo Espectro ASE: \n") ; fprintf('%.2f %% \n' ,  ( 50 + (((z-0.5)/Nz)*50 )) ) % Mostrar % de avance del cálculo
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
                        AuxMatrix(i) = P_ase0_sp%+(((N2(Nz-z+1))*sigma_ems(lambda_ase(i)) - (N1(Nz-z+1))*sigma_abs(lambda_ase(i)))*Gamma_s*P_ase0_sp + 2*(N2(Nz-z+1))*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z;
                    end
                    Pan_sp.(ModoS(s))(:,Nz-z+1) = AuxMatrix(:) ;
                else
                    parfor i = 1:Nch_ase
                        if (i==1)
                            d_vk=abs(v_s_sp(i)-v_s_sp(i+1));
                        else
                            d_vk=abs(v_s_sp(i-1)-v_s_sp(i));
                        end
                        AuxMatrix(i) = Pan_sp.(ModoS(s))(i,Nz-z+1+1) + (((N2(Nz-z+1))*sigma_ems(lambda_ase(i)) - (N1(Nz-z+1))*sigma_abs(lambda_ase(i)))*Gamma_s*Pan_sp.(ModoS(s))(i,Nz-z+1+1) + 2*(N2(Nz-z+1))*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z;
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



    %%

    % Calculo OSNR
    for s = 1:1:Smod
        Nwl = length(Signal.lambda.(ModoS(s)));
        for i = 1:1:Nwl
            OSNR.(ModoS(s))(i,:) = 10*log10(Psp.(ModoS(s))(i,:)/1e-3) - 10*log10(Pase.(ModoS(s))(i,:)/1e-3);
            gain.(ModoS(s))(1,i) = 10*log10(Psp.(ModoS(s))(i,end)/Psp.(ModoS(s))(i,1));
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
        sdm.(ch).salida.ASE.potencia_dBm.(ModoS(s)) = 10*log10(Pase.(ModoS(s))(:,end)./1e-3);

        sdm.(ch).salida.ganancias.(ModoS(s)) = gain.(ModoS(s));
        sdm.(ch).NF.(ModoS(s)) = OSNR.(ModoS(s))(:,1)./OSNR.(ModoS(s))(:,end);
    end
    sdm.(ch).GEF.best_weight = best_weight;

    sdm.(ch).signal.lambdas = lambda_s ; sdm.(ch).pump.lambdas = lambda_p;
    sdm.(ch).OSNR = OSNR ;
    sdm.(ch).salida.OSNR = OSNR(:,end);
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

