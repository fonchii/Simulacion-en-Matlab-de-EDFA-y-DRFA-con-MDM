% Simulador EDFA MM
close all; 
clear all; clc

%% Parámetros de entrada

    % Señal : Modos y Canales
NCh = 31;
Signal.modos = ["01" "11a"] ;
Pin = -20; %-35 ; %-15

Signal.lambda.LP_01     = linspace(1530e-9,1560e-9,NCh);              P0_signal.LP_01     = Pin*ones(1,length(Signal.lambda.LP_01));
Signal.lambda.LP_11a     = linspace(1530e-9,1560e-9,NCh);              P0_signal.LP_11a     = Pin*ones(1,length(Signal.lambda.LP_11a));

    % Bombeo : Modos y Canales
Pump.modos = ["01" ]   ;

Pump.lambda.LP_01   = 980e-9;                                 P0_pump.LP_01   = 100e-3  ;  

    % POTENCIAS

for i=1:length(Signal.modos)        % Potencia: dBm -> mW
    for j=1:length( P0_signal.(strcat("LP_",Signal.modos(i))) )
        P0_signal.( strcat("LP_",Signal.modos(i) ) )(j) = 1e-3*10^(P0_signal.( strcat("LP_",Signal.modos(i) ) )(j)/10);
    end
end ;clear i j;

Signal.P0 = P0_signal; 
Pump.P0 = P0_pump;
ASE = -200;                                                  %dBm  -50
Signal.NumberOfChannels=NCh;

    % Datos de la fibra

Fibra.nucleos = 1;                                           % Numero de nucleos
Fibra.largo = 5     ; Fibra.radio = 5e-6   ; Fibra.N = 7e24; % fibra.N = 3e24; 

Fibra.dvk=300e9;



Fibra.n1 = 1.45 ;   
Fibra.n2 = 1.4354 ;
Fibra.CrossSectionParameters = "OptiSystem" ;
%Fibra.dvk= P.OpticalBW; % diferencia : max_lambda - min_lambda 

% Despliegue de Info
Fibra.WaitBar = 1; Fibra.Avance = 1;    

% Calculo ASE?
Fibra.ASEFlag = 1;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)


%% 


tic;
% % %% --------- Iterar por Largos --------- %%
% l_cont = 1;
% for largos=[3 5 7 10 15 20]%largos=1:2:20
%         for pump = [50,100,150,200,250,300,400,500]
%             
%             P0_pump.LP_01   = pump*1e-3  ; Pump.P0 = P0_pump;
%             Fibra.largo = largos;
%             fprintf('Iniciando Fibra %.0f de %.0f\n', l_cont,length([50,100,150,200,250,300,400,500])*length([3 5 7 10 15 20]));
%             l_cont = l_cont +1;
% 
%             Largos_v2.(strcat("EDFA_",num2str(largos),'m')).(strcat('Pump',num2str(pump),'mw')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);
%         end
%     
% end



% --------- Iterar por Potencias --------- %

% cont = 1;
% for largos=[3 5 7 10 15 20]
%     for potencias=50:50:500
%         P0_pump.LP_01 = potencias.*1e-3; Pump.P0 = P0_pump;
%         fprintf('Iniciando Fibra %.0f de %.0f\n', cont,(length(50:50:500)*length([3 5 7 10 15 20]))); cont = cont+1;
%         Fibra.largo = largos; % 3,5,7,10
%         Potencias.(strcat('L',num2str(largos),'m')).(strcat("EDFA_",num2str(potencias),'mw')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);       % Con efecto acomplamiento de Potencia intermodal
%     
%     end
% end


% --------- Ganancias Espectral --------- %
% % Cambiar arriba a estos datos:
% %Signal.lambda.LP_01 ;    P0_signal.LP_01   ;

cont = 1;
for largos=[3 5 7 10 15 20]
    for potencias=50:50:500
        P0_pump.LP_01 = potencias.*1e-3; Pump.P0 = P0_pump;
        fprintf('Iniciando Fibra %.0f de %.0f\n', cont,(length(50:50:500)*length([3 5 7 10 15 20]))); cont = cont+1;
        Fibra.largo = largos; % 3,5,7,10
        GainSpectrum.(strcat('L',num2str(largos),'m')).(strcat("EDFA_",num2str(potencias),'mw')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);       % Con efecto acomplamiento de Potencia intermodal
    end
end



% % --------- Potencia saturacion a la salida --------- % (Ganancia vs Pin)
% cont = 1;
% for Pin=10:2:50
%     Signal.modos = ["01" "11a"] ;
%     Signal.lambda.LP_01     = linspace(1530e-9,1560e-9,NCh);              P0_signal.LP_01     = -1*Pin*ones(1,length(Signal.lambda.LP_01));
%     Signal.lambda.LP_11a    = linspace(1530e-9,1560e-9,NCh);              P0_signal.LP_11a    = -1*Pin*ones(1,length(Signal.lambda.LP_11a));
%     for i=1:length(Signal.modos)        % Potencia: dBm -> mW
%         for j=1:length( P0_signal.(strcat("LP_",Signal.modos(i))) )
%             P0_signal.( strcat("LP_",Signal.modos(i) ) )(j) = 1e-3*10^(P0_signal.( strcat("LP_",Signal.modos(i) ) )(j)/10);
%         end
%     end ;clear i j;
%     Signal.P0 = P0_signal; 
%     
%     Fibra.largo = 10;  % 5
%     Pump.lambda.LP_01   = 980e-9;                                 P0_pump.LP_01   = 150e-3  ; %100  
%     Pump.P0 = P0_pump;
% 
%     fprintf('Iniciando Fibra %.0f de %.0f\n', cont,(length(10:2:50))); cont = cont+1;
%     
%     P_signalIn.(strcat('Pin',num2str(Pin),'dB')) = EDFA_MMvPCCv3(Fibra,Signal,Pump,ASE);       % Con efecto acomplamiento de Potencia intermodal
%     %P_signalInGEF.(strcat('Pin',num2str(Pin),'dB')) = EDFA_MMvPCCv3_GEF(Fibra,Signal,Pump,ASE);
% end




Iterations_time=toc; fprintf('Tiempo Total de cómputo: %.2f segundos\n', Iterations_time )




