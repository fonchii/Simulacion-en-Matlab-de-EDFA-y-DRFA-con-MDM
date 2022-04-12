function edfa = EDFA_MM(fibra,signal,pump,ASE)  
    % Datos de entradas:
    % signal (struct) 
    %       signal.lambda -> int [nm]
    %       signal.P0 -> int [mW]
    %       signal.modos -> str
    % pump (struct)
    %       pump.lambda -> int [nm]
    %       pump.P0 -> int [mW]
    %       pump.modos -> str
    % ASE -> potencia -> int [mW]
    % largo -> int [m]
    % canales espaciales -> int [numero canales espaciales]
    % fibra (struct)
    %       fibra.largo -> int
    %       fibra.n1/n2 -> int (indices de refraccion)
    %       fibra.N -> int ( Concentracion total de iones de Erbio)
    %       

lambda_s = signal.lambda;   P_s0 = signal.P0; 
lambda_p = pump.lambda; P_p0 = pump.P0; 
modos_s = signal.modos ; modos_p = pump.modos;
Pase0_dBm = ASE*ones(1,length(P_s0));
nucleos = fibra.nucleos;

AN = fibra.AN ; L = fibra.largo ; N = fibra.N;
a = fibra.radio;

% Confirmar que Lambdas esten en nm
for i = 1:1:length(lambda_s)
    if lambda_s(i) > 1
        lambda_s(i) = lambda_s(i)*1e-9;
    end
end
for i = 1:1:length(lambda_p)
    if lambda_p(i) > 1
        lambda_p(i) = lambda_p(i)*1e-9;
    end
end
    
    
    % % Parametros de la fibra y constantes
c = 3e8;                                % Velocidad de la luz en el vacio m/s
m = 2;                                  % numero de polarizaciones
h = 6.626*10^(-34);                     % constante de planck
del_z = 1/10;                          % Tamaño de paso de iteraciones
z = linspace(0,L,L*(1/del_z));          % se toman segmentos de 10mm
Nz = length(z);
d_vk = 125*10^9;                        % 1 nm - Noise Bandwidth
%N = 3e24;                               
tau = 10e-3;


% Datos de Bombeo y Señal

v_p = 3*10^8./lambda_p;                  % frecuencia de bombeo
v_s = c./lambda_s;                      % frecuencias de señal
Nch = length(lambda_s);                 % N° total de wavelengths
Nchp = length(lambda_p);
Sch = nucleos;                          % N° de nucleos
P_ase0 = 1e-3*10.^(Pase0_dBm/10);       % Potencia ASE entrada en Watts

Psat = 3;                               % Potencia de saturación (dBm)
Psat = 1e-3*10^(Psat/10);

V_s = (2*pi ./ lambda_s)*AN*a;          % frecuencias normalizadas
V_p = (2*pi ./ lambda_p)*AN*a;
M_s = floor(V_s.^2 / 2) ;               % Cantidad de modos
M_p = floor(V_p.^2 / 2) ; 

% %   Espectros emision y absorción
 
E = load('D:\Documentos\UTFSM\Elo\2021-2\Proyecto de titulacion\OptiSystem\Data Files\EDFAs\Erbium.dat');
lambda_cross = E(:,1).*1e-9;
Gamma_s = 1;                    % mode overlap factor for signal & ASE entre modo y perfil de dopaje
Gamma_p = 1;                    % mode overlap factor for pump entre modo y perfil de dopaje
alpha_cross = E(:,2); alpha_cross = alpha_cross./(Gamma_s);
g_cross = E(:,3); g_cross = g_cross./(Gamma_s);
sigma_abs = fit(lambda_cross,alpha_cross,'linearinterp'); 
sigma_ems = fit(lambda_cross,g_cross,'linearinterp');

A_s = 58.1e-12;                    % Area efectiva para señal
%A_p = 34.6e-12;                   % Area efectiva para bombeo
A_p = 58.1e-12;                    % Area efectiva para bombeo

sigma_p_a = sigma_abs(lambda_p);
sigma_p_e = sigma_ems(lambda_p);

sdm = struct;

% % Resolución de ecuaciones

for s = 1:1:Sch
    ch = strcat('Nucleo',int2str(s));
    % Variables para cálculo de N2 y N1
    N2 = zeros(1,Nz);                    % Densidad de poblacion en estado excitado
    N1 = zeros(1,Nz);                    % Densidad de pobiacion en estado basal
    Nt = zeros(1,Nz);
    Psp = zeros(Nch,Nz);                 % Potencia de señal en dirección +z
    Ppp = zeros(Nchp,Nz);                % Potencia de bombeo en dirección +z
    Pap = zeros(Nch,Nz);                 % Potencia ASE en dirección +z
    Pan = zeros(Nch,Nz);                 % Potencia ASE en dirección -z
    Pase = zeros(Nch,Nz);                % Potencia total ASE en el EDFA
    OSNR = zeros(Nch,Nz);                % Relacion señal a ruido 

    G = zeros(Nch,Nz);                   % G_k coeficiente de ganancia 
    gain = zeros(1,Nch);
    
    %lambda_ase = 1520e-9:0.5*1e-9:1570e-9;
    lambda_ase = ase_lambdas(lambda_s);
    Nch_ase = length(lambda_ase);
    ASE_Spectrum = zeros(Nch_ase,2);
    Pap_sp = zeros(Nch_ase,Nz);                 % Potencia ASE en dirección +z
    Pan_sp = zeros(Nch_ase,Nz);                 % Potencia ASE en dirección -z
    P_ase0_sp = 1e-3*10.^(ASE*ones(1,length(lambda_ase))/10);
    v_s_sp = c./lambda_ase;
    
    for j = 1:1:Nz  
 
    
% Calculo de densidades de iones en estado basal y excitado (N1 y N2)
        sig_xx = 0;                 
        sig_yy = 0;                    % Inicialización de variables 
        ase_xx = 0;                     
        ase_yy = 0;
        pmp_xx = 0;
        pmp_yy = 0;
         
% % eq19 Giles 1991 para N(z)
% usando los cambios de variable : sigma_a * gamma = alpha/nt ; A = pi*b_eff

    % Primera Iteración
            if(j == 1)
    % Potencia Bombeo
                for i=1:1:Nchp
                    pmp_xx = pmp_xx + (sigma_p_a(i))*(P_p0(i)*Gamma_p/A_p)/(h*v_p(i));                    % Termino en numerador
                    pmp_yy = pmp_yy + (sigma_p_a(i)+sigma_p_e(i))*(P_p0(i)*Gamma_p/A_p)/(h*v_p(i));          % Termino en denominador
                end
	% Potencia de señal
                for i = 1:1:Nch                                                           
                    sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(P_s0(i)*Gamma_s/A_s)/(h*v_s(i)));                        % Termino en numerador               
                    sig_yy = sig_yy+(sigma_abs(lambda_s(i))+sigma_ems(lambda_s(i)))*(P_s0(i)*Gamma_s/A_s)/(h*v_s(i)); % Termino en denominador
                end
    % Potencia ASE            
                for i = 1:1:Nch                                                           
                    ase_xx = ase_xx+(sigma_abs(lambda_s(i))*(P_ase0(i)*Gamma_s/A_s)/(h*v_s(i)));                        % Termino en numerador
                    ase_yy = ase_yy+(sigma_abs(lambda_s(i)+sigma_ems(lambda_s(i)))*(P_ase0(i)*Gamma_s/A_s)/(h*v_s(i))); % Termino en denominador
                end

            else        
    % Iteraciones Restantes
                for i = 1:1:Nchp
                    pmp_xx = pmp_xx + (sigma_p_a(i))*(Ppp(i,j-1)*Gamma_p/A_p)/(h*v_p(i));                     % Termino en numerador
                    pmp_yy = pmp_yy + (sigma_p_a(i)+sigma_p_e(i))*(Ppp(i,j-1)*Gamma_p/A_p)/(h*v_p(i));           % Termino en denominador
                end
                for i = 1:1:Nch                                                           
                    sig_xx = sig_xx+(sigma_abs(lambda_s(i))*(Psp(i,j-1)*Gamma_s/A_s)/(h*v_s(i)));                                                       
                    sig_yy = sig_yy+(sigma_abs(lambda_s(i)+sigma_ems(lambda_s(i)))*(Psp(i,j-1)*Gamma_s/A_s)/(h*v_s(i)));
                end
                for i = 1:1:Nch                                                           
                    ase_xx = ase_xx+(sigma_abs(lambda_s(i))*(Pase(i,j-1)*Gamma_s/A_s)/(h*v_s(i)));                                                         
                    ase_yy = ase_yy+(sigma_abs(lambda_s(i)+sigma_ems(lambda_s(i)))*(Pase(i,j-1)*Gamma_s/A_s)/(h*v_s(i)));
                end
            end

            N2(j) = ((pmp_xx+sig_xx+ase_xx)/(pmp_yy+sig_yy+ase_yy+(1/tau)))*N;             % Densidad de iones de Erbio en estado excitado
            N1(j) = N-N2(j);                                                               % Densidad de iones de Erbio en estado basal
            Nt(j) = N1(j) + N2(j);
            for i = 1:1:Nch
                G(i,j) = G(i,j) + del_z*( sigma_ems(lambda_s(i))*N2(j) - sigma_abs(lambda_s(i))*N1(j) );
            end

    %% Ecuaciones diferenciales 
    % Ecuacion diferencial para bombeo en direccion +z
            if(j == 1)
                for i = 1:1:Nchp 
                    Ppp(i,j) = P_p0(i)+(N2(j)*sigma_p_e(i) - N1(j)*sigma_p_a(i))*Gamma_p*P_p0(i)*del_z;
                end
            else
                for i = 1:1:Nchp
                    Ppp(i,j) = Ppp(i,j-1)+(N2(j)*sigma_p_e(i)-N1(j)*sigma_p_a(i))*Gamma_p*Ppp(i,j-1)*del_z;
                end
            end

    % Ecuacion diferencial para señal en direccion +z
            if(j == 1)
                for i = 1:1:Nch
                    Psp(i,j) = P_s0(i)+(N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*Gamma_s*P_s0(i)*del_z/(1+P_s0(i)/Psat);
                    %Psp(i,j) = P_s0(i)+(N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*gamma_s*P_s0(i)*del_z;
                end
            else
                for i = 1:1:Nch
                    Psp(i,j) = Psp(i,j-1)+(N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*Gamma_s*Psp(i,j-1)*del_z/(1+Psp(i,j-1)/Psat);
                    %Psp(i,j) = Psp(i,j-1)+(N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*gamma_s*Psp(i,j-1)*del_z;
                end
            end

    % Ecuacion diferencial para ASE en direccion +z
            if(j == 1)
                for i = 1:1:Nch
                    Pap(i,j) = P_ase0(i)+((N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*Gamma_s*P_ase0(i) + m*N2(j)*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+P_ase0(i)/Psat);
                    %Pap(i,j) = P_ase0(i)+((N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*gamma_s*P_ase0(i) + m*N2(j)*sigma_ems(lambda_s(i))*gamma_s*h*v_s(i)*d_vk)*del_z;
                    %Pap(i,j) = P_ase0(i)+((N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*gamma_s*P_ase0(i) + M_s(i)*N2(j)*sigma_ems(lambda_s(i))*gamma_s*h*v_s(i)*d_vk)*del_z/(1+P_ase0(i)/Psat);
                end
            else
                for i = 1:1:Nch
                    Pap(i,j) = Pap(i,j-1)+((N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*Gamma_s*Pap(i,j-1) + m*N2(j)*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pap(i,j-1)/Psat);
                    %Pap(i,j) = Pap(i,j-1)+((N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*gamma_s*Pap(i,j-1) + m*N2(j)*sigma_ems(lambda_s(i))*gamma_s*h*v_s(i)*d_vk)*del_z;
                    %Pap(i,j) = Pap(i,j-1)+((N2(j)*sigma_ems(lambda_s(i))-N1(j)*sigma_abs(lambda_s(i)))*gamma_s*Pap(i,j-1) + M_s(i)*N2(j)*sigma_ems(lambda_s(i))*gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pap(i,j-1)/Psat);
                end
            end

    % Ecuacion diferencial para bombeo en direccion -z
            if(j == 1)
                for i = 1:1:Nch
                    Pan(i,Nz-j+1) = 0+((N2(Nz-j+1)*sigma_ems(lambda_s(i))-N1(Nz-j+1)*sigma_abs(lambda_s(i)))*Gamma_s*0 + m*N2(Nz-j+1)*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z;
                    %Pan(i,Nz-j+1) = 0+((N2(Nz-j+1)*sigma_ems(lambda_s(i))-N1(Nz-j+1)*sigma_abs(lambda_s(i)))*gamma_s*0 + M_s(i)*N2(Nz-j+1)*sigma_ems(lambda_s(i))*gamma_s*h*v_s(i)*d_vk)*del_z;
                end
            else
                for i = 1:1:Nch
                    Pan(i,Nz-j+1) = Pan(i,Nz-j+1+1)+((N2(Nz-j+1)*sigma_ems(lambda_s(i))-N1(Nz-j+1)*sigma_abs(lambda_s(i)))*Gamma_s*Pan(i,Nz-j+1+1) + 2*N2(Nz-j+1)*sigma_ems(lambda_s(i))*Gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pan(i,Nz-j+1+1)/Psat);
                    %Pan(i,Nz-j+1) = Pan(i,Nz-j+1+1)+((N2(Nz-j+1)*sigma_ems(lambda_s(i))-N1(Nz-j+1)*sigma_abs(lambda_s(i)))*gamma_s*Pan(i,Nz-j+1+1) + 2*N2(Nz-j+1)*sigma_ems(lambda_s(i))*gamma_s*h*v_s(i)*d_vk)*del_z;
                    %Pan(i,Nz-j+1) = Pan(i,Nz-j+1+1)+((N2(Nz-j+1)*sigma_ems(lambda_s(i))-N1(Nz-j+1)*sigma_abs(lambda_s(i)))*gamma_s*Pan(i,Nz-j+1+1) + M_s(i)*N2(Nz-j+1)*sigma_ems(lambda_s(i))*gamma_s*h*v_s(i)*d_vk)*del_z/(1+Pan(i,Nz-j+1+1)/Psat);
                end
            end
% ASE SPECTRUM : 
    % Ecuacion diferencial para ASE en direccion +z
            if(j == 1)
                for i = 1:1:Nch_ase
                    Pap_sp(i,j) = P_ase0_sp(i)+((N2(j)*sigma_ems(lambda_ase(i))-N1(j)*sigma_abs(lambda_ase(i)))*Gamma_s*P_ase0_sp(i) + 2*N2(j)*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+P_ase0_sp(i)/Psat);
                    %Pap(i,j) = P_ase0_sp(i)+((N2(j)*sigma_ems(lambda_ase(i))-N1(j)*sigma_abs(lambda_ase(i)))*gamma_s*P_ase0_sp(i) + M_s(i)*N2(j)*sigma_ems(lambda_ase(i))*gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+P_ase0_sp(i)/Psat);
                end
            else
                for i = 1:1:Nch_ase
                    Pap_sp(i,j) = Pap_sp(i,j-1)+((N2(j)*sigma_ems(lambda_ase(i))-N1(j)*sigma_abs(lambda_ase(i)))*Gamma_s*Pap_sp(i,j-1) + 2*N2(j)*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+Pap_sp(i,j-1)/Psat);
                    %Pap(i,j) = Pap(i,j-1)+((N2(j)*sigma_ems(lambda_ase(i))-N1(j)*sigma_abs(lambda_ase(i)))*gamma_s*Pap(i,j-1) + M_s(i)*N2(j)*sigma_ems(lambda_ase(i))*gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+Pap(i,j-1)/Psat);
                end
            end

    % Ecuacion diferencial para bombeo en direccion -z
            if(j == 1)
                for i = 1:1:Nch_ase
                    Pan_sp(i,Nz-j+1) = 0+((N2(Nz-j+1)*sigma_ems(lambda_ase(i))-N1(Nz-j+1)*sigma_abs(lambda_ase(i)))*Gamma_s*0 + 2*N2(Nz-j+1)*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z;
                    %Pan(i,Nz-j+1) = 0+((N2(Nz-j+1)*sigma_ems(lambda_ase(i))-N1(Nz-j+1)*sigma_abs(lambda_ase(i)))*gamma_s*0 + M_s(i)*N2(Nz-j+1)*sigma_ems(lambda_ase(i))*gamma_s*h*v_s_sp(i)*d_vk)*del_z;
                end
            else
                for i = 1:1:Nch_ase
                    Pan_sp(i,Nz-j+1) = Pan_sp(i,Nz-j+1+1)+((N2(Nz-j+1)*sigma_ems(lambda_ase(i))-N1(Nz-j+1)*sigma_abs(lambda_ase(i)))*Gamma_s*Pan_sp(i,Nz-j+1+1) + 2*N2(Nz-j+1)*sigma_ems(lambda_ase(i))*Gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+Pan_sp(i,Nz-j+1+1)/Psat);
                    %Pan(i,Nz-j+1) = Pan(i,Nz-j+1+1)+((N2(Nz-j+1)*sigma_ems(lambda_ase(i))-N1(Nz-j+1)*sigma_abs(lambda_ase(i)))*gamma_s*Pan(i,Nz-j+1+1) + M_s(i)*N2(Nz-j+1)*sigma_ems(lambda_ase(i))*gamma_s*h*v_s_sp(i)*d_vk)*del_z/(1+Pan(i,Nz-j+1+1)/Psat);
                end
            end

    % Calculo de Potencia ASE 
            Pase(:,j) = Pap(:,j)+Pan(:,j); 
            for i = 1:1:Nch_ase
                %ASE_Spectrum(i,1) = 10*log10( (Pap_sp(i,end) + Pan_sp(i,end)) ./ (Pap_sp(i,1) + Pan_sp(i,1)) );
                ASE_Spectrum(i,1) = 10*log10( (Pap_sp(i,end) + Pan_sp(i,end)));
            end
            
    end   % Fin de calculo en fibra

    % Calculo OSNR
            for i = 1:1:Nch
                OSNR(i,:) = 10*log10(Psp(i,:)/1e-3) - 10*log10(Pase(i,:)/1e-3);
                gain(1,i) = 10*log10(Psp(i,end)/Psp(i,1));
            end
            
    % CALCULAR GRAFICO ASE SPECTRUM CON GANANCIAS DE SEÑAL
    ASE_Spectrum(:,2) = ASE_Spectrum(:,1);
    for j = 1:length(lambda_s)
        for i = 1:length(lambda_ase)
            if abs(lambda_ase(i)*1e9 - lambda_s(j)*1e9) < 1/11
                ASE_Spectrum(i,2) = gain(1,j);
            end
        end
    end
                
        % Guardar los datos calculados en la struct sdm
            sdm.(ch).N2 = N2 ; sdm.(ch).N1 = N1 ; sdm.(ch).Nt = Nt;
            %SEÑAL
            sdm.(ch).signal.Potencia_dBm = 10*log10(Psp./1e-3) ; sdm.(ch).pump.Potencia_dBm = 10*log10(Ppp./1e-3) ;
            sdm.(ch).signal.lambdas = lambda_s ; sdm.(ch).pump.lambdas = lambda_p;
            sdm.(ch).Pase = 10*log10(Pase./1e-3);
            sdm.(ch).OSNR = OSNR ;
            sdm.(ch).salida.signal.potencia_dBm = 10*log10(Psp(:,end)./1e-3);
            sdm.(ch).salida.pump.potencia_dBm = 10*log10(Ppp(:,end)./1e-3);
            sdm.(ch).salida.ASE.potencia_dBm = 10*log10(Pase(:,end)./1e-3);
            sdm.(ch).salida.OSNR = OSNR(:,end);
            sdm.(ch).salida.ganancias = gain;
            sdm.(ch).z = z;
            sdm.(ch).NF = OSNR(:,1)./OSNR(:,end);
            sdm.(ch).ASE_Spectrum.mag = ASE_Spectrum;
            sdm.(ch).ASE_Spectrum.lambdas = lambda_ase;

end % Fin cálculo en todos los canales espaciales

edfa = sdm;
end


% FALTA : 
%    IMPLEMENTAR NUMERO DE MODOS
%    IMPLEMENTAR INTEGRAL DE SUPERPOSICION PARA CADA MODO
