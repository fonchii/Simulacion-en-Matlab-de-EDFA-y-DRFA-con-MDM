    %% Graficos

close all ; clear all


set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% SAVE:
% print -dpdf 'NAME'





% % Ganancias vs lamda
% Cargar datos y guardar en struct
% for j=["150mw" "200mw" "250mw" "300mw" "400mw" "500mw" "700mw" "1000mw" "1500mw"]
%     load(strcat("Largo_Caracterizacion_",j,".mat"))
%     for i=1:length(fieldnames(Largos))
%         Ganancias.(strcat('g_',num2str(i),'m')) = Largos.(strcat("EDFA_",num2str(i),'m')).Nucleo1.salida.ganancias.LP_01(11:end-10); %
%         GananciasvsLargo_Temp(i) = Ganancias.(strcat('g_',num2str(i),'m'))(6); %1555 nm
%     end
%     GananciasvsLargo.(strcat("P",j)) = GananciasvsLargo_Temp;
% end


% for g = [1,2,3,5,7,10,15]
%     leyenda = strcat('Largo: ',int2str(g)," m");
%     ejex = Largos.EDFA_1m.Nucleo1.signal.lambdas(11:end-10).*1e9;
%     plot(ejex,Ganancias.(strcat('g_',num2str(g),'m')) , '-o' , 'DisplayName',leyenda ) ; hold on ; 
% end ; clear s ejex leyenda;
% set(gca,'FontSize',8)
% legend(Location="southoutside",FontSize=9,Box="off",Orientation="horizontal");  
% title('Distribución espectras de Ganancias','FontSize',14) ; xlabel('Longitud de onda [nm]','FontSize',14) ; ylabel('Magnitud [dB]','FontSize',14)


% Ganancias vs Largo
% load("GananciavsLargo2.mat") ; load("Largo_Caracterizacion_1500mw.mat")
% xlargos=[1:length(fieldnames(Largos))];
% for g = [150,200,250,300,400,500,700,1000,1500]
%     plot(xlargos(1:end-15) , GananciasvsLargo.(strcat('P',int2str(g),'mw'))(1:end-15) , 'DisplayName' , strcat('Pump= ',int2str(g),' mw') )  ; hold on;
%     
% end ; clear s ejex leyenda;
% set(gca,'FontSize',8)
% legend(Location="southoutside",FontSize=9,Box="off",Orientation="horizontal",NumColumns=5)
% title('Ganancia vs Largo del amplificador','FontSize',14) ; xlabel('Largo del EDFA [m]','FontSize',14) ; ylabel('Ganancia [dB]','FontSize',14)



% ----------- Ganancias vs Potencia
% % % Cargar datos y guardar en struct
% for j=["3m" "5m" "7m" "10m" "15m" "20m"]
%     load(strcat("Potencia_Caracterizacion_",j,".mat"))
%     %load(strcat("Potencia_Caracterizacion_10m",".mat"))
%     for i=1:length(fieldnames(Potencias))
%         Ganancias.(strcat('g_',num2str(i),'m')) = Potencias.(strcat("EDFA_",num2str(100+50*i),'mw')).Nucleo1.salida.ganancias.LP_01;%;(11:end-10); %
%         GananciasvsLargos_Temp(i) = Ganancias.(strcat('g_',num2str(i),'m'))(21); %1560 nm
%     end
%     PotenciasvsLargo.(strcat("L",j)) = GananciasvsLargos_Temp;
% end


load("PotenciasvsLargo.mat")
ejex = [150:50:1000];
for j=["3m" "5m" "7m" "10m" "15m" "20m"]
    plot(ejex , PotenciasvsLargo.(strcat("L",j)) , 'DisplayName' , strcat('Largo= ',j)  )  ; hold on;
end ; clear s ejex leyenda;
set(gca,'FontSize',8)
legend(Location="southoutside",FontSize=9,Box="off",Orientation="horizontal",NumColumns=5)
title('Ganancia vs Potencia de Bombeo del amplificador','FontSize',14) ; xlabel('Potencia de Bombeo [mw]','FontSize',14) ; ylabel('Ganancia [dB]','FontSize',14)




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
