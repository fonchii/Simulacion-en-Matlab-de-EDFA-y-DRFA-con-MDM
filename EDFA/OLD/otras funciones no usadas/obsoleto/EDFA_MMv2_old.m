function edfa = EDFA_MMv2(fibra,signal,pump,ASE)  
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

L = fibra.largo ; N = fibra.N;

    % % Parametros de la fibra y constantes
c = 3e8;                                % Velocidad de la luz en el vacio m/s
m = 2;                                  % numero de polarizaciones
h = 6.626*10^(-34);                     % constante de planck
del_z = 1;                           % Tamaño de paso de iteraciones
Z = linspace(0,L,L*(1/del_z));          % se toman segmentos de 10cm
if length(Z)<30 % Correccion para EDFA muy corto
    Z=linspace(0,L,30);
    del_z = (Z(2)-Z(1));
end
Nz = length(Z);
d_vk = 125*10^9;                        % 1 nm - Noise Bandwidth                         
tau = 10e-3;
nucleos = fibra.nucleos; 
Sch = nucleos;                          % N° de nucleos

% Datos de Bombeo y Señal

P_s0 = signal.P0; 
P_p0 = pump.P0; 

Smod = length(signal.modos);            % N° de modos de señal
Pmod = length(pump.modos);              % N° de modos de bombeo

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


A_s = pi*fibra.radio^2;                                % Area efectiva para señal
A_p = pi*fibra.radio^2 ; %58.1e-12;                    % Area efectiva para bombeo

sdm = struct;

% Calculo de Gammas
warning('off')

for p = 1:1:length(pump.modos) % Mode overlap factor for pump entre modo y perfil de dopaje
    Nwlp = length(pump.lambda.(strcat("LP_",pump.modos(p))));
    for i=1:1:Nwlp
        lambda_p = pump.lambda.(strcat("LP_",pump.modos(p)));
        gamma_p.(strcat("LP_",pump.modos(p))){i} = norm_intensity(fibra,pump.modos(p),lambda_p(i));
    end
end 

for s = 1:1:length(signal.modos) % Mode overlap factor for signal entre modo y perfil de dopaje
    Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
    
    for i=1:1:Nwl
        lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
        CoupCoef = 0;
        
        for p = 1:1:length(pump.modos) % Agregar coeficiente de acoplamiento entre bombeo y señal
            Nwlp = length(pump.lambda.(strcat("LP_",pump.modos(p))));
             for pp=1:1:Nwlp
                 lambda_p = pump.lambda.(strcat("LP_",pump.modos(p)));
                 CoupCoef = CoupCoef + coupling_coef(fibra,signal.modos(s),pump.modos(p),lambda_s(i),lambda_p(pp));
             end
        end
        
        gamma_s.(strcat("LP_",signal.modos(s))){i} = norm_intensity(fibra,signal.modos(s),lambda_s(i)) * CoupCoef;
        %gamma_s.(strcat("LP_",signal.modos(s))){i} = CoupCoef;

    end
end
warning('on')

% ASE Spectrum 
allwavelengths = [];
for s = 1:1:Smod
    for w=1:1:length(signal.lambda.(strcat("LP_",signal.modos(s))))
        if ~(ismember(signal.lambda.(strcat("LP_",signal.modos(s)))(w) , allwavelengths))
            allwavelengths = [allwavelengths signal.lambda.(strcat("LP_",signal.modos(s)))(w)];
        end
    end
end
allwavelengths = sort(allwavelengths);
if(isfield(fibra,'ASEFlag'))
    lambda_ase = ase_lambdas(allwavelengths,1); % retorna [0,0,0] no calcula espectro ASE
else 
    lambda_ase = ase_lambdas(allwavelengths);
end

for n = 1:1:Sch     % Iteración en nucleos
    ch = strcat('Nucleo',int2str(n));
    % Variables para cálculo de N2 y N1
    N2 = zeros(1,Nz);                    % Densidad de poblacion en estado excitado
    N1 = zeros(1,Nz);                    % Densidad de pobiacion en estado basal
    Nt = zeros(1,Nz);
    Psp = [];                           % Potencia de señal en dirección +z
    Ppp = [];                           % Potencia de bombeo en dirección +z
    Pap = [];                           % Potencia ASE en dirección +z
    Pan = [];                           % Potencia ASE en dirección -z
    Pase = [];                          % Potencia total ASE en el EDFA
    OSNR = [];                          % Relacion señal a ruido
    
    G = [];                             % G_k coeficiente de ganancia
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
    
    for z = 1:1:Nz  % Iteraciones a lo largo del EDFA
              
        % % Eq19 Giles 1991 para N(z)
        % Usando los cambios de variable : sigma_a * gamma = alpha/nt ; A = pi*b_eff
               
        % Primera Iteración
        if(z == 1)
            % Potencia Bombeo
            for p = 1:1:Pmod
                Nwlp = length(pump.lambda.(strcat("LP_",pump.modos(p))));
                for i=1:1:Nwlp
                    lambda_p = pump.lambda.(strcat("LP_",pump.modos(p)));
                    Gamma_p = gamma_p.(strcat("LP_",pump.modos(p))){i};
                    v_p = c./lambda_p; Pp0 = P_p0.(strcat("LP_",pump.modos(p)))(i);
                    
                    pmp_xx = pmp_xx + (sigma_abs(lambda_p(i)))*(Pp0*Gamma_p/A_p)/(h*v_p(i));                                   % Termino en numerador
                    pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Pp0*Gamma_p/A_p)/(h*v_p(i));          % Termino en denominador
                end
            end
            % Potencia de señal
            for s = 1:1:Smod
                Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    v_s = c./lambda_s;  Ps0 = P_s0.(strcat("LP_",signal.modos(s)))(i);
                    
                    sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(Ps0*Gamma_s/A_s)/(h*v_s(i)));                        % Termino en numerador
                    sig_yy = sig_yy+(sigma_abs(lambda_s(i))+sigma_ems(lambda_s(i)))*(Ps0*Gamma_s/A_s)/(h*v_s(i)); % Termino en denominador
                end
            end
            % Potencia ASE
            for s=1:1:Smod
                Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    v_s = c./lambda_s;
                    
                    ase_xx = ase_xx+ sigma_abs(lambda_s(i))*(P_ase0*Gamma_s/A_s)/(h*v_s(i));                           % Termino en numerador
                    ase_yy = ase_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(P_ase0*Gamma_s/A_s)/(h*v_s(i)); % Termino en denominador
                end
            end
            
        else
            % Iteraciones Restantes
                % Bombeo
            for p = 1:1:Pmod  
                Nwlp = length(pump.lambda.(strcat("LP_",pump.modos(p))));
                for i = 1:1:Nwlp
                    lambda_p = pump.lambda.(strcat("LP_",pump.modos(p)));
                    Gamma_p = gamma_p.(strcat("LP_",pump.modos(p))){i};
                    v_p = c./lambda_p;
                                        
                    pmp_xx = pmp_xx + sigma_abs(lambda_p(i))*(Ppp.(strcat("LP_",pump.modos(p)))(i,z-1)*Gamma_p/A_p)/(h*v_p(i));                        % Termino en numerador
                    pmp_yy = pmp_yy + (sigma_abs(lambda_p(i)) + sigma_ems(lambda_p(i)))*(Ppp.(strcat("LP_",pump.modos(p)))(i,z-1)*Gamma_p/A_p)/(h*v_p(i));           % Termino en denominador
                end
            end
                % Señal
            for s = 1:1:Smod
                Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    v_s = c./lambda_s;
                    
                    sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(Psp.(strcat("LP_",signal.modos(s)))(i,z-1)*Gamma_s/A_s)/(h*v_s(i)));
                    sig_yy = sig_yy+(sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i))) * (Psp.(strcat("LP_",signal.modos(s)))(i,z-1)*Gamma_s/A_s)/(h*v_s(i));
                end
            end
                
                % ASE
            for s = 1:1:Smod
                Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    v_s = c./lambda_s;
                    
                    ase_xx = ase_xx+(sigma_abs(lambda_s(i))*(Pase.(strcat("LP_",signal.modos(s)))(i,z-1)*Gamma_s/A_s)/(h*v_s(i)));
                    ase_yy = ase_yy+( sigma_abs(lambda_s(i)) + sigma_ems(lambda_s(i)) )*(Pase.(strcat("LP_",signal.modos(s)))(i,z-1)*Gamma_s/A_s)/(h*v_s(i));
                end
            end
        end
        
        N2(z) = ((pmp_xx+sig_xx+ase_xx)/(pmp_yy+sig_yy+ase_yy+(1/tau)))*N;             % Densidad de iones de Erbio en estado excitado
        %((pmp_xx+sig_xx+ase_xx)/(pmp_yy+sig_yy+ase_yy+(1/tau)))
        N1(z) = N-N2(z);                                                               % Densidad de iones de Erbio en estado basal
        Nt(z) = N1(z) + N2(z);
%         for i = 1:1:Nwl
%             G(i,z) = G(i,z) + del_z*( sigma_ems(lambda_s(i))*N2(z) - sigma_abs(lambda_s(i))*N1(z) );
%         end
        
        %% Ecuaciones de Potencias

        % Ecuacion diferencial para bombeo en direccion +z
        for p = 1:1:Pmod      % Iteracion en cada modo de bombeo
            Nwlp = length(pump.lambda.(strcat("LP_",pump.modos(p))));
            if(z == 1)
                for i = 1:1:Nwlp
                    lambda_p = pump.lambda.(strcat("LP_",pump.modos(p)));
                    Gamma_p = gamma_p.(strcat("LP_",pump.modos(p))){i};
                    Pp0 = P_p0.(strcat("LP_",pump.modos(p)))(i);
                    
                    Ppp.(strcat("LP_",pump.modos(p)))(i,z) = Pp0+(N2(z)*sigma_ems(lambda_p(i)) - N1(z)*sigma_abs(lambda_p(i)))*Gamma_p*Pp0*del_z;
                end
            else
                for i = 1:1:Nwlp
                    lambda_p = pump.lambda.(strcat("LP_",pump.modos(p)));
                    Gamma_p = gamma_p.(strcat("LP_",pump.modos(p))){i};
                    
                    Ppp.(strcat("LP_",pump.modos(p)))(i,z) = Ppp.(strcat("LP_",pump.modos(p)))(i,z-1)+(N2(z)*sigma_ems(lambda_p(i)) - N1(z)*sigma_abs(lambda_p(i)))*Gamma_p*Ppp.(strcat("LP_",pump.modos(p)))(i,z-1)*del_z;
                end
            end
        end
        
        % Ecuacion diferencial para señal en direccion +z
        for s = 1:1:Smod
            Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
            if(z == 1)
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    Ps0 = P_s0.(strcat("LP_",signal.modos(s)))(i);
                    
                    Psp.(strcat("LP_",signal.modos(s)))(i,z) = Ps0+(N2(z)*sigma_ems(lambda_s(i))-N1(z)*sigma_abs(lambda_s(i)))*Gamma_s*Ps0*del_z/(1+Ps0/Psat);
                end
            else
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    
                    Psp.(strcat("LP_",signal.modos(s)))(i,z) = Psp.(strcat("LP_",signal.modos(s)))(i,z-1)+(N2(z)*sigma_ems(lambda_s(i))-N1(z)*sigma_abs(lambda_s(i)))*Gamma_s*Psp.(strcat("LP_",signal.modos(s)))(i,z-1)*del_z/(1+Psp.(strcat("LP_",signal.modos(s)))(i,z-1)/Psat);
                end
            end
        end
        
        % Ecuacion diferencial para ASE en direccion +z
        for s = 1:1:Smod
            Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
            if(z == 1)
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    
                    Pap.(strcat("LP_",signal.modos(s)))(i,z) = P_ase0+((N2(z)*sigma_ems(lambda_s(i))-N1(z)*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0 + m*N2(z)*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+P_ase0/Psat);
                end
            else
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    
                    Pap.(strcat("LP_",signal.modos(s)))(i,z) = Pap.(strcat("LP_",signal.modos(s)))(i,z-1)+((N2(z)*sigma_ems(lambda_s(i))-N1(z)*sigma_abs(lambda_s(i)))*Gamma_s*Pap.(strcat("LP_",signal.modos(s)))(i,z-1) + m*N2(z)*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pap.(strcat("LP_",signal.modos(s)))(i,z-1)/Psat);
                end
            end
        end
        
        % Ecuacion diferencial para ASE en direccion -z
        for s = 1:1:Smod
            Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
            if(z == 1)
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    
                    Pan.(strcat("LP_",signal.modos(s)))(i,Nz-z+1) = 0+((N2(Nz-z+1)*sigma_ems(lambda_s(i))-N1(Nz-z+1)*sigma_abs(lambda_s(i)))*Gamma_s*0 + m*N2(Nz-z+1)*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                end
            else
                for i = 1:1:Nwl
                    lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
                    Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){i};
                    
                    Pan.(strcat("LP_",signal.modos(s)))(i,Nz-z+1) = Pan.(strcat("LP_",signal.modos(s)))(i,Nz-z+1+1)+((N2(Nz-z+1)*sigma_ems(lambda_s(i))-N1(Nz-z+1)*sigma_abs(lambda_s(i)))*Gamma_s*Pan.(strcat("LP_",signal.modos(s)))(i,Nz-z+1+1) + 2*N2(Nz-z+1)*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pan.(strcat("LP_",signal.modos(s)))(i,Nz-z+1+1)/Psat);

                end
            end
        end
        
        % ASE SPECTRUM :
        % Ecuacion diferencial para ASE en direccion +z
        for s = 1:1:Smod
            Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
            for w = 1:1:Nwl
                Nch_ase = length(lambda_ase);
                P_ase0_sp = P_ase0;
                v_s_sp = c./lambda_ase;
                Gamma_s = gamma_s.(strcat("LP_",signal.modos(s))){w};
                if(z == 1)
                    for i = 1:1:Nch_ase
                        Pap_sp.(strcat("LP_",signal.modos(s)))(i,z) = P_ase0_sp + ((N2(z)*sigma_ems(lambda_ase(i))-N1(z)*sigma_abs(lambda_ase(i)))*Gamma_s*P_ase0_sp + 2*N2(z)*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+P_ase0_sp/Psat);
                    end
                else
                    for i = 1:1:Nch_ase
                        Pap_sp.(strcat("LP_",signal.modos(s)))(i,z) = Pap_sp.(strcat("LP_",signal.modos(s)))(i,z-1)+((N2(z)*sigma_ems(lambda_ase(i))-N1(z)*sigma_abs(lambda_ase(i)))*Gamma_s*Pap_sp.(strcat("LP_",signal.modos(s)))(i,z-1) + 2*N2(z)*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+Pap_sp.(strcat("LP_",signal.modos(s)))(i,z-1)/Psat);
                    end
                end
            end
        end
        
        % Ecuacion diferencial para ASE en direccion -z
        for s = 1:1:Smod
            if(z == 1)
                for i = 1:1:Nch_ase
                    Pan_sp.(strcat("LP_",signal.modos(s)))(i,Nz-z+1) = 0+((N2(Nz-z+1)*sigma_ems(lambda_ase(i))-N1(Nz-z+1)*sigma_abs(lambda_ase(i)))*Gamma_s*0 + 2*N2(Nz-z+1)*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z;
                end
            else
                for i = 1:1:Nch_ase
                    Pan_sp.(strcat("LP_",signal.modos(s)))(i,Nz-z+1) = Pan_sp.(strcat("LP_",signal.modos(s)))(i,Nz-z+1+1)+((N2(Nz-z+1)*sigma_ems(lambda_ase(i))-N1(Nz-z+1)*sigma_abs(lambda_ase(i)))*Gamma_s*Pan_sp.(strcat("LP_",signal.modos(s)))(i,Nz-z+1+1) + 2*N2(Nz-z+1)*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+Pan_sp.(strcat("LP_",signal.modos(s)))(i,Nz-z+1+1)/Psat);
                end
            end
        end
        
        % Calculo de Potencia ASE
        for s = 1:1:Smod
            Pase.(strcat("LP_",signal.modos(s)))(:,z) = Pap.(strcat("LP_",signal.modos(s)))(:,z)+Pan.(strcat("LP_",signal.modos(s)))(:,z);
            for i = 1:1:Nch_ase
                ASE_Spectrum.(strcat("LP_",signal.modos(s)))(i,1) = 10*log10( (Pap_sp.(strcat("LP_",signal.modos(s)))(i,end) + Pan_sp.(strcat("LP_",signal.modos(s)))(i,end)));
            end
        end
        
    end   % Fin de calculo en fibra
    
    % Calculo OSNR
    for s = 1:1:Smod
        Nwl = length(signal.lambda.(strcat("LP_",signal.modos(s))));
        for i = 1:1:Nwl
            OSNR.(strcat("LP_",signal.modos(s)))(i,:) = 10*log10(Psp.(strcat("LP_",signal.modos(s)))(i,:)/1e-3) - 10*log10(Pase.(strcat("LP_",signal.modos(s)))(i,:)/1e-3);
            gain.(strcat("LP_",signal.modos(s)))(1,i) = 10*log10(Psp.(strcat("LP_",signal.modos(s)))(i,end)/Psp.(strcat("LP_",signal.modos(s)))(i,1));
        end
    end
    
    % CALCULAR GRAFICO ASE SPECTRUM CON GANANCIAS DE SEÑAL
    for s = 1:1:Smod
        ASE_Spectrum.(strcat("LP_",signal.modos(s)))(:,2) = ASE_Spectrum.(strcat("LP_",signal.modos(s)))(:,1);
        
        lambda_s = signal.lambda.(strcat("LP_",signal.modos(s)));
        for w = 1:length(lambda_s)
            for i = 1:length(lambda_ase)
                if abs(lambda_ase(i)*1e9 - lambda_s(w)*1e9) < 1/11
                    %ASE_Spectrum(i,2,s) = gain(1,w,s);
                    ASE_Spectrum.(strcat("LP_",signal.modos(s)))(i,2) = 10*log10( (1e-3*10.^(ASE_Spectrum.(strcat("LP_",signal.modos(s)))(i,1)/10) + Psp.(strcat("LP_",signal.modos(s)))(w,end)) ./1e-3);
                    %ASE_Spectrum(i,2,s) = 10*log10( (( Psp(w,end,s)) ./1e-3) );
                end
            end
        end
    end
    % Guardar los datos calculados
    sdm.(ch).N2 = N2 ; sdm.(ch).N1 = N1 ; sdm.(ch).Nt = Nt;
    for p=1:Pmod
        sdm.(ch).pump.Potencia_dBm.(strcat("LP_",pump.modos(p))) = 10*log10(Ppp.(strcat("LP_",pump.modos(p)))./1e-3) ;
        sdm.(ch).salida.pump.potencia_dBm.(strcat("LP_",pump.modos(p))) = 10*log10(Ppp.(strcat("LP_",pump.modos(p)))(:,end)./1e-3);
    end
    for s=1:Smod
        sdm.(ch).signal.Potencia_dBm.(strcat("LP_",signal.modos(s))) = 10*log10(Psp.(strcat("LP_",signal.modos(s)))./1e-3) ;
        sdm.(ch).salida.signal.potencia_dBm.(strcat("LP_",signal.modos(s))) = 10*log10(Psp.(strcat("LP_",signal.modos(s)))(:,end)./1e-3);
        
        sdm.(ch).Pase.(strcat("LP_",signal.modos(s))) = 10*log10(Pase.(strcat("LP_",signal.modos(s)))./1e-3);
        sdm.(ch).salida.ASE.potencia_dBm.(strcat("LP_",signal.modos(s))) = 10*log10(Pase.(strcat("LP_",signal.modos(s)))(:,end)./1e-3);
        
        sdm.(ch).salida.ganancias.(strcat("LP_",signal.modos(s))) = gain.(strcat("LP_",signal.modos(s)));
        sdm.(ch).NF.(strcat("LP_",signal.modos(s))) = OSNR.(strcat("LP_",signal.modos(s)))(:,1)./OSNR.(strcat("LP_",signal.modos(s)))(:,end);
    end

    sdm.(ch).signal.lambdas = lambda_s ; sdm.(ch).pump.lambdas = lambda_p;
    sdm.(ch).OSNR = OSNR ;
    sdm.(ch).salida.OSNR = OSNR(:,end);
    sdm.(ch).z = Z;
    sdm.(ch).ASE_Spectrum.mag = ASE_Spectrum;
    sdm.(ch).ASE_Spectrum.lambdas = lambda_ase;
    sdm.(ch).Gamma.Gamma_p = gamma_p;
    sdm.(ch).Gamma.Gamma_s = gamma_s;
    sdm.(ch).Gamma.Gamma_s.info = "Este coeficiente es la multiplicación entre la superposición del modo con el dopaje multiplicado y la suma de los acoplamientos con cada bombeo";
end % Fin cálculo en todos los nucleos

    edfa = sdm;
    
end % Fin funcion
    
%% Pendiente
% Revisar potencia de saturacion y Area efectiva
% Superposición entre señal/bombeo y perfil de dopaje es muy bajo
