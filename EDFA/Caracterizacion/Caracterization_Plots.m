    %% Graficos

close all ; clear all

%% Cargar datos y guardar en struct

% ---------- Ganancias vs lamda ---------- %
% %%load(strcat("Data_GananciavsLargo_OptiSystem.mat")) ; 
% load(strcat("Data_GananciaEspectro.mat")) ; 
% for larg=[3 5 7 10 15 20] %1:length(fieldnames(Largos))
%     cont=1;
%     for j=["150mw" "200mw" "250mw" "300mw" "400mw" "500mw" "700mw" "1000mw"]
%         Ganancias.(strcat('L',num2str(larg),'m'))(cont,:) = GainSpectrum.(strcat("L",num2str(larg),'m')).(strcat("EDFA_",j)).Nucleo1.salida.ganancias.LP_01; 
%         cont=cont+1;
%     end
% end
% Ganancias.ejex = GainSpectrum.(strcat("L",num2str(larg),'m')).(strcat("EDFA_",j)).Nucleo1.signal.lambdas;



% ---------- Ganancias vs Largo ---------- %

% % %load(strcat("Largos_Potencias")) % Datos VPI
% load(strcat("Data_Largos_Potencias_OptiSystem")) % Datos OptiSystem
% for j=["150mw" "200mw" "250mw" "300mw" "400mw" "500mw" "700mw" "1000mw" "1500mw"]
%     cont = 1;
%     for i=[1,3,5,7,9,11,13,15,17,19]
%         GananciasvsLargo_Temp(cont) = Largos_v2.(strcat("EDFA_",num2str(i),'m')).(strcat('Pump',j)).Nucleo1.salida.ganancias.LP_01(16); % 1555 nm
%         cont = cont+1;
%     end
%     GananciasvsLargo.(strcat("P",j)) = GananciasvsLargo_Temp;
% end
% 



% ----------- Ganancias vs Potencia ---------- %
% %load("PotenciasvsLargo.mat")
% load("Data_GananciavsLargo_OptiSystem.mat")
% %for i=[3,5,7,9,11,13,15,17,19]
% for i=[3,5,7,10,15,20]
%     cont = 1;
%     for j=["150mw" "200mw" "250mw" "300mw" "350mw" "400mw" "450mw" "500mw" "550mw" "600mw" "650mw" "700mw" "750mw" "800mw" "850mw" "900mw" "950mw" "1000mw"]
%         GananciasvsPot_Temp(cont) = Potencias.(strcat('L',num2str(i),'m')).(strcat('EDFA_',j)).Nucleo1.salida.ganancias.LP_01(16); % 1555 nm
%         cont = cont+1;
%     end
%     GananciasvsPot.(strcat("L",num2str(i),'m')) = GananciasvsPot_Temp;
% end


%% GRAFICAR

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% SAVE:
% print -dpdf 'NAME'

% ---------- Ganancia vs Lambda ---------- %

load("Ganancias_OptiSystem.mat")
larg = ["3m" "5m" "7m" "10m" "15m" "20m"];
pots = ["150mw" "200mw" "250mw" "300mw" "400mw" "500mw" "700mw" "1000mw"];
cont = 1;

for l = "L10m" %larg
    figure(cont)
    for p = 1:length(Ganancias.L10m(:,1))
        leyenda = strcat('Pump ',pots(p));
        ejex = Ganancias.ejex.*1e9;
        plot(ejex,Ganancias.(l)(p,:) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
    end %; clear s ejex leyenda;
    set(gca,'FontSize',8)
    legend(Location="southoutside",FontSize=9,Box="off",Orientation="horizontal", NumColumns=4);  
    title('Distribución espectral de Ganancias para EDFA de largo 10m','FontSize',14) ; xlabel('Longitud de onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
    ylim([0 35])
    cont = cont + 1;
end


% ---------- GANANCIA VS PumpPower ---------- %
% load("GananciavsPumpPower_OptiSystem.mat")
% ejex = [150:50:1000];
% %for j=["3m" "5m" "7m" "9m" "11m" "13m" "15m" "17m" "19m"]
% for j=["3m" "5m" "7m" "10m" "15m" "20m"]
%     %plot(ejex , PotenciasvsLargo.(strcat("L",j)) , 'DisplayName' , strcat('Largo= ',j)  )  ; hold on;
%     plot(ejex , GananciasvsPot.(strcat("L",j)) , 'DisplayName' , strcat('Largo= ',j)  )  ; hold on;
% end ; clear s ejex leyenda;
% set(gca,'FontSize',8)
% legend(Location="southoutside",FontSize=9,Box="off",Orientation="horizontal",NumColumns=3)
% title('Ganancia vs Potencia de Bombeo para canal de 1555 nm','FontSize',14) ; xlabel('Potencia de Bombeo [mw]','FontSize',14) ; ylabel('Ganancia [dB]','FontSize',14)


% --------- Ganancia vs Largo --------- %

% xlargos=[1,3,5,7,9,11,13,15,17,19];
% for g = [150,200,250,300,400,500,700,1000,1500]
%     plot(xlargos , GananciasvsLargo.(strcat('P',int2str(g),'mw')) , 'DisplayName' , strcat('Pump= ',int2str(g),' mw') )  ; hold on;
%     
% end ; clear s ejex leyenda;
% set(gca,'FontSize',8)
% legend(Location="southoutside",FontSize=9,Box="off",Orientation="horizontal",NumColumns=5)
% title('Ganancia vs Largo del amplificador para canal de 1555 nm ','FontSize',14) ; xlabel('Largo del EDFA [m]','FontSize',14) ; ylabel('Ganancia [dB]','FontSize',14)


%% %% Parámetros de entrada
% 
%     % Señal : Modos y Canales
% NCh = 30;
% Signal.modos = ["01"] ;
% 
% Signal.lambda.LP_01     = linspace(1525e-9,1585e-9,NCh);              P0_signal.LP_01     = -15*ones(1,length(Signal.lambda.LP_01));
% 
%     % Bombeo : Modos y Canales
% Pump.modos = ["01" ]   ;
% 
% Pump.lambda.LP_01   = 980e-9;                                 P0_pump.LP_01   = [250e-3]  ;  
% 
%     % POTENCIAS
% 
% for i=1:length(Signal.modos)        % Potencia: dBm -> mW
%     for j=1:length( P0_signal.(strcat("LP_",Signal.modos(i))) )
%         P0_signal.( strcat("LP_",Signal.modos(i) ) )(j) = 1e-3*10^(P0_signal.( strcat("LP_",Signal.modos(i) ) )(j)/10);
%     end
% end ;clear i j;
% 
% Signal.P0 = P0_signal; 
% Pump.P0 = P0_pump;
% ASE = -200;                                                  %dBm  -50
% Signal.NumberOfChannels=NCh;
% 
%     % Datos de la fibra
% 
% Fibra.nucleos = 1;                                           % Numero de nucleos
% Fibra.largo = 1     ; Fibra.radio = 5e-6   ; Fibra.N = 7e24; % fibra.N = 3e24; 
% 
% Fibra.dvk=300e9;
% Signal.NumberOfChannels=30;
% 
% 
% Fibra.n1 = 1.45 ;   
% Fibra.n2 = 1.4354 ;
% %Fibra.dvk= P.OpticalBW; % diferencia : max_lambda - min_lambda 
% 
% % Despliegue de Info
% Fibra.WaitBar = 1; Fibra.Avance = 1;    
% 
% % Calculo ASE?
% Fibra.ASEFlag = 1;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)

%Largos.(strcat("EDFA_",num2str(largos),'m')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);       % Con efecto acomplamiento de Potencia intermodal
