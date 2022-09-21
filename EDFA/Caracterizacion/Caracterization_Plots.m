    %% Graficos

close all ; clear all

%% Cargar datos y guardar en struct

% ---------- Ganancias vs lamda ---------- %
% %%load(strcat("Data_GananciavsLargo_OptiSystem.mat")) ; 
% load(strcat("Data_GainSpectrum_5m-100mw.mat")) ; 
% for larg=[3 5 7 10 15 20] %1:length(fieldnames(Largos))
%     cont=1;
%     for j=["50mw" "100mw"  "150mw" "200mw" "300mw" "400mw" "500mw"]
%         Ganancias.(strcat('L',num2str(larg),'m'))(cont,:) = GainSpectrum.(strcat("L",num2str(larg),'m')).(strcat("EDFA_",j)).Nucleo1.salida.ganancias.LP_01; 
%         cont=cont+1;
%     end
% end
% Ganancias.ejex = GainSpectrum.(strcat("L",num2str(larg),'m')).(strcat("EDFA_",j)).Nucleo1.signal.lambdas;

% Data_GananciavsPumpPower_OptiSystem2 -> -35dBm entrada
% Data_GananciavsPumpPower_OptiSystem2_5m100mw -> -20dBm entrada
% OptiSystem -> comienza desde 150mw hasta 1000
% OptiSystem2 -> comienza desde 50mw hasta 500



% ---------- Ganancias vs Largo ---------- %

% load(strcat("Data_GananciavsLargo_OptiSystem2_5m100mw")) % Datos OptiSystem
% %for j=["150mw" "200mw" "250mw" "300mw" "400mw" "500mw" "700mw" "1000mw" "1500mw"]
% for j=["50mw" "100mw" "150mw" "200mw" "250mw" "300mw" "400mw" "500mw"]
%     cont = 1;
%     for i=[3,5,7,10,15,20]
%         GananciasvsLargo_Temp(cont) = Largos_v2.(strcat("EDFA_",num2str(i),'m')).(strcat('Pump',j)).Nucleo1.salida.ganancias.LP_01(16); % 1555 nm
%         cont = cont+1;
%     end
%     GananciasvsLargo.(strcat("P",j)) = GananciasvsLargo_Temp;
% end




% ----------- Ganancias vs Potencia ---------- %

% load("Data_GananciavsPumpPower_OptiSystem2_5m100mw.mat") 
%  
% for i=[3,5,7,10,15,20] 
%     cont = 1;
%     %for j=["150mw" "200mw" "250mw" "300mw" "350mw" "400mw" "450mw" "500mw" "550mw" "600mw" "650mw" "700mw" "750mw" "800mw" "850mw" "900mw" "950mw" "1000mw"]
%     for j=["50mw" "100mw" "150mw" "200mw" "250mw" "300mw" "350mw" "400mw" "450mw" "500mw"]
%         GananciasvsPot_Temp(cont) = Potencias.(strcat('L',num2str(i),'m')).(strcat('EDFA_',j)).Nucleo1.salida.ganancias.LP_01(16); % 1555 nm
%         cont = cont+1;
%     end
%     GananciasvsPot.(strcat("L",num2str(i),'m')) = GananciasvsPot_Temp;
% end


% % ----------- Ganancias vs Potencia Entrada ---------- %   Pout vs Pin
% load("Data_Gain-vs-PsIn_5m-100mw") % Data_Gain-vs-PsIn_10m-150mw.mat
% load("Data_Gain-vs-PsIn_10m-200mw.mat")
% %for i=[3,5,7,9,11,13,15,17,19] 
% row = 0;
% for Pin=10:2:50
%     row = row+1;
%     g = P_signalIn.(strcat('Pin',num2str(Pin),'dB')).Nucleo1.salida.ganancias.LP_01;
%     pout = P_signalIn.(strcat('Pin',num2str(Pin),'dB')).Nucleo1.salida.signal.potencia_dBm.LP_01;
%     GananciasvsPin.(strcat("Pin",num2str(Pin),'dB')) = g;%[g(1) g(11) g(21) g(31)]; % 1530 , 1540, 1550, 1560 nm
%     g1(row) = pout(1); g2(row) = pout(11) ; g3(row) = pout(21) ; g4(row) = pout(31);
% end
% GananciasvsPin.Potencias_ejex =  P_signalIn.(strcat('Pin',num2str(Pin),'dB')).Nucleo1.signal.lambdas;
% GananciasvsPin.Lambda1530nm = g1; GananciasvsPin.Lambda1540nm = g2; GananciasvsPin.Lambda1550nm = g3; GananciasvsPin.Lambda1560nm = g4;
% GananciasvsPin.Lambda_ejex = [-10:-2:-50];




%% GRAFICAR

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% SAVE:
% print -dpdf 'NAME'

% ---------- Ganancia vs Lambda ---------- %

load("Struct_SpectralGain_5m-100mw.mat")
larg = ["3m" "5m" "7m" "10m" "15m" "20m"];
pots = ["50mw" "100mw" "150mw" "200mw"  "300mw"  "400mw"  "500mw" ];
cont = 1;

% for l = "L5m" %larg
%     figure(cont)
%     for p = 1:length(Ganancias.L5m(:,1))
%         leyenda = strcat('Pump ',pots(p));
%         ejex = Ganancias.ejex.*1e9;
%         plot(ejex,Ganancias.(l)(p,:) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
%     end %; clear s ejex leyenda;
%     set(gca,'FontSize',8)
%     legend(Location="southoutside",FontSize=9,Box="off",Orientation="horizontal", NumColumns=4);  
%     title('Distribución espectral de Ganancias para EDFA de largo 5m','FontSize',14) ; xlabel('Longitud de onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)
%     ylim([5 30])
%     cont = cont + 1;
% end
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0]) ; print -dpdf 'GananciaEspectral_vsPump_5m' 

% ---------- GANANCIA VS PumpPower ---------- %
% load("Struct_GananciavsPumpPower_OptiSystem2_5m100mw.mat")
% ejex = [50:50:500]; %ejex = [150:50:1000];
% for j=["3m" "5m" "7m" "10m" "15m" ]%for j=["3m" "5m" "7m" "10m" "15m" "20m"]
%     %plot(ejex , PotenciasvsLargo.(strcat("L",j)) , 'DisplayName' , strcat('Largo= ',j)  )  ; hold on;
%     plot(ejex , GananciasvsPot.(strcat("L",j))(1:end) , 'DisplayName' , strcat('Largo= ',j)  )  ; hold on;
% end ; clear s ejex leyenda;
% set(gca,'FontSize',8)
% legend(Location="southoutside",FontSize=9,Box="off",Orientation="horizontal",NumColumns=3)
% title('Ganancia vs Potencia de Bombeo para canal de 1555 nm','FontSize',14) ; xlabel('Potencia de Bombeo [mw]','FontSize',14) ; ylabel('Ganancia [dB]','FontSize',14)
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% ylim([5 30])
% % % SAVE:
% print -dpdf 'GananciavsPumpPower'

% --------- Ganancia vs Largo --------- %
% load("Struct_GananciasvsLargo2_5m100mw.mat")
% xlargos=[3,5,7,10,15]; %xlargos=[3,5,7,10,15,20];
% %for g = [150,200,250,300,400,500,700,1000,1500]
% for g = [100,150,200,250,300,400,500]
%     plot(xlargos , GananciasvsLargo.(strcat('P',int2str(g),'mw'))(1:end-1) , 'DisplayName' , strcat('Pump= ',int2str(g),' mw') )  ; hold on;
%     
% end ; clear s ejex leyenda;
% set(gca,'FontSize',8)
% legend(Location="southoutside",FontSize=9,Box="off",Orientation="horizontal",NumColumns=5)
% title('Ganancia vs Largo del amplificador para canal de 1555 nm ','FontSize',14) ; xlabel('Largo del EDFA [m]','FontSize',14) ; ylabel('Ganancia [dB]','FontSize',14)
%  set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% % % SAVE:
% print -dpdf 'GananciavsLargo'


% % ----------- Ganancias vs Potencia Entrada ---------- %
% load("Struct_Ganancias_vs_Pin_5m-100mw.mat") %  Struct_Ganancias_vs_Pin_10m-150mw.mat
% figure(1)
% for Pin = [10 16 20 26 34 40]%10:6:50
%     plot(GananciasvsPin.Potencias_ejex.*1e9 , GananciasvsPin.(strcat("Pin",num2str(Pin),'dB')) , "DisplayName" , strcat("Pin " , num2str(-Pin) , " dB") ) ; hold on ; 
% end
% title("Distribución espectral de ganancia vs Potencia de entrada","FontSize",14) ; ylabel("Ganancia [dB]","FontSize",14) ; xlabel("Longitud de onda [nm]","FontSize",14)
% legend('Location','southoutside','Box','off','Orientation','horizontal','FontSize',9,'NumColumns',6)
% 
% set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% % % % SAVE:
% print -dpdf 'GananciaEspectral_vsPin' ;% close all;

% %% ----- Pout vs Pin ----- %%

% figure(2)
% for i = [1530,1540,1550,1560]
%     plot( GananciasvsPin.Lambda_ejex  , GananciasvsPin.(strcat("Lambda",num2str(i),"nm")) , "DisplayName",strcat(num2str(i)," nm")) ; hold on
% end
% xline(-20,'--',"DisplayName","Punto de Operación")
% set(gca,  'fontSize',8); % 'xdir','reverse',
% title("Potencia de Salida vs Potencia de Entrada","FontSize",14) ; ylabel("Potencia de Salida [dBm]","FontSize",14) ; xlabel("Potencia de Entrada [dBm]","FontSize",14)
% legend('Location','southoutside','Box','off','Orientation','horizontal','FontSize',9)

%set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% % % %SAVE:
%print -dpdf 'Pout_vs_Pin'


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
