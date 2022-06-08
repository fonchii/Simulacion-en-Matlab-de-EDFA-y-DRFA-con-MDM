% close all; 
clear all; clc

% Parámetros de entrada

    % Modos y canales de señal y bombeo
signal.NumberOfChannels=20;
signal.modos = ["01" "11_a"] ;
Frequency_gridS=linspace(191.07234e12,196.7723e12,signal.NumberOfChannels);
c=299.792458e6; % [m/s]
Wavelength_gridS=c./Frequency_gridS;

Pin=-10; %[dBm]

signal.lambda.LP_01     = Wavelength_gridS;                                  P0_signal.LP_01     = Pin*ones(1,length(signal.lambda.LP_01));
signal.lambda.LP_11_a   = Wavelength_gridS;                                  P0_signal.LP_11_a   = Pin*ones(1,length(signal.lambda.LP_11_a));


pump.modos = "01" ;
Wavelength_gridP=980e-9;
Ppump= 250e-3; %[W]

pump.lambda.LP_01   = Wavelength_gridP;                         P0_pump.LP_01   = Ppump  ;  

ModoS=strcat("LP_",signal.modos(:));
ModoP=strcat("LP_",pump.modos(:));

    % POTENCIAS
for i=1:length(signal.modos)        % Potencia de señal a W
    for j=1:length(P0_signal.(ModoS(i)))
        P0_signal.(ModoS(i))(j) = 1e-3*10^(P0_signal.(ModoS(i))(j)/10);
    end
end ;clear i j;
signal.P0 = P0_signal; 
pump.P0 = P0_pump;
h=6.62607015*10^(-34);
P.Np=2; 
P.Fc=c/Wavelength_gridS(ceil(length(Wavelength_gridS)/2)); P.Fb = 300e9; 
ASE= -200;

    % Datos de la fibra
fibra.nucleos = 1;                                           % Numero de nucleos
fibra.largo = 3; fibra.radio = 5e-6 ; fibra.N = 7e24; % fibra.N = 3e24; 
%fibra.AN = 0.2;   %fibra.AN = 0.2 ; %fibra.AN=sqrt(fibra.n1^2-fibra.n2^2);
fibra.n1 = 1.45 ;   fibra.IndexContrast=0.01;
fibra.AN=fibra.n1*sqrt(2*fibra.IndexContrast);
fibra.n2 =sqrt((fibra.n1^2-fibra.AN^2));
fibra.M = 10; fibra.Nalpha = inf;
fibra.dvk=P.Fb;

fibra.WaitBar = 1; fibra.Avance = 1;    % Despliegue de info
fibra.ASEFlag = 0;                      % 1 : Evita Calculo Espectro ASE ; 0 : Lo Calcula (lento)

%%
tic;
EDFA = EDFA_MMvpi2(fibra,signal,pump,ASE);%EDFA_MMvPCCv3(fibra,signal,pump,ASE);
% EDFA = EDFA_MM(fibra,signal,pump,ASE); %197.82 segundos
t_end = toc; fprintf('Tiempo de cómputo: %.2f segundos\n', t_end);


    %% Graficos

% close all
graf.Nc = length(fieldnames(EDFA)); 
graf.z = EDFA.(strcat('Nucleo',int2str(1))).z;
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
for n = 1:graf.Nc
% colorlist= [0 0.4470 0.7410; 0.8500 0.325 0.098];
% n=1;
% figure(n)
%%%%%%%%%%%%%%%%%%%      solo para gráficos similares al VPI - Getting Started %%%%%%%%%%%%%%%%%%%
% for ase = 1:length(signal.modos)
%     for aux=1:length(EDFA.Nucleo1.ASE_Spectrum.lambdas)
%         if (graf.ASE_Spectrum.(ModoS(ase))(aux,2)>0 && ase==1)
%             aux1=stem(graf.freq_ase(aux),graf.ASE_Spectrum.(ModoS(ase))(aux,2)','^','BaseValue',-10,'Color',colorlist(ase,:)) ; hold on;
%     
%         else 
%             if (graf.ASE_Spectrum.(ModoS(ase))(aux,2)>0 && ase==2)
%                 aux2=stem(graf.freq_ase(aux),graf.ASE_Spectrum.(ModoS(ase))(aux,2)','^','BaseValue',-10,'Color',colorlist(ase,:)) ; hold on;
%             end
%         end
%         %xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Espectro ASE') ; %legend()
%     end 
%    
% end  %clear ase xlab ylab t_end ;
% 
% figure(n)
% for ase = 1:length(signal.modos)
%     area( graf.freq_ase, graf.ASE_Spectrum.(ModoS(ase))(:,1) ,'BaseValue',-60,'FaceColor',colorlist(ase,:),'DisplayName',strcat( "LP",signal.modos(ase) )) ; hold on;
% end 
% xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Noise and Signals Amplified by EDFA MM') ; 
% legend([aux1 aux2],{'{LP}_{01}','{LP}_{11_a}'})
% ylim([-60 12])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Señal
%     for s = 1:length(signal.modos) % Grafico
%         graf.signal.(ModoS(s)) = EDFA.(strcat("Nucleo",int2str(n))).signal.Potencia_dBm.(ModoS(s));
%         subplot 221
%         for f = 1:length(signal.lambda.(ModoS(s)))
%             plot(graf.z,graf.signal.(ModoS(s))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",signal.modos(s))," @"),strcat(int2str(signal.lambda.(strcat("LP_",signal.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Signal}') ; legend(); grid on ; legend('location', 'best');
%         end
%     end; clear s f ;
%     
% %     % Bombeo
% %     for s = 1:length(pump.modos) % Grafico
% %         graf.pump.(strcat("LP_",pump.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).pump.Potencia_dBm(:,:,s);
% %         figure(n)
% %         subplot 222
% %         for f = 1:length(pump.lambda.LP_01)
% %             plot(z,graf.pump.(strcat("LP_",pump.modos(s)))(f,:) , 'DisplayName',strcat(strcat(strcat("LP",pump.modos(s))," @"),strcat(int2str(pump.lambda.(strcat("LP_",pump.modos(s)))(f)*1e9) ,' nm')) ) ; hold on ; xlabel(xlab) ; ylabel(ylab); title('P_{Pump}') ; legend(); grid on
% %         end
% %     end
% 
%     % Ganancias
%     for s = 1:length(signal.modos)
%         graf.ganancias.(ModoS(s)) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(ModoS(s));
%         leyenda = strcat((strcat("LP",signal.modos(s))));
%         ax = subplot (2,2,2);
%         
%         ejex = signal.lambda.(ModoS(s)).*1e9;
%         plot(ejex,graf.ganancias.(ModoS(s)) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%         %ax.ColorOrderIndex = s;
%         xlabel('λ [nm]') ; ylabel(ylab)
%         
%     end ; clear s ejex leyenda;
%     
    % Espectro ASE
    if fibra.ASEFlag == 0
        for ase = 1:length(signal.modos)
            graf.freq_ase = EDFA.(strcat('Nucleo',int2str(n))).ASE_Spectrum.lambdas*1e9;
            graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase))) = EDFA.(strcat("Nucleo",int2str(n))).ASE_Spectrum.mag.(strcat("LP_",signal.modos(ase)));
            %subplot(2,2,[3,4])
            plot( graf.freq_ase, graf.ASE_Spectrum.(strcat("LP_",signal.modos(ase))) , 'DisplayName', strcat( "LP",signal.modos(ase) ) ) ; hold on;
            xlabel('Longitud de onda (\lambda) [nm]') ; ylabel(ylab); title('Potencia a la salida del EDFA') ; legend()
        end ; clear ase xlab ylab t_end ;
    end
    
end
% 
% %% Graficos Modos
% 
% % Señal
% n=n+1 ; figure(n)
% graficar_modos(fibra,signal.modos,signal.lambda.LP_01)
% sgtitle('Modos de Señal') 
% 
% % Bombeo
% n=n+1 ; figure(n)
% graficar_modos(fibra,pump.modos,pump.lambda.LP_01)
% sgtitle('Modos de Bombeo') ; %clear n ;

%% Graficar potencia vs frecuencia
% n=n+1; figure(n)
figure(2)
xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
graf.Nc = length(fieldnames(EDFA)); 
for n = 1:graf.Nc
    for s = 1:length(signal.modos)
        graf.ganancias.(ModoS(s)) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(ModoS(s));
        leyenda = strcat((ModoS(s)));
        ejex = fliplr((3*10^8./(signal.lambda.(ModoS(s)).*1e9)));
        plot(ejex/1000,fliplr(graf.ganancias.(ModoS(s))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
        %ax.ColorOrderIndex = s;
        xlabel('Frecuencia [THz]') ; ylabel(ylab)

    end ; clear s ejex leyenda;
end


%% Graficar potencia vs wavelenght
%n=n+1; figure(n)
%xlab = 'Posición en fibra [m]'; ylab = 'Potencia [dBm]';
%graf.Nc = length(fieldnames(EDFA)); 
%for n = 1:graf.Nc
%    for s = 1:length(signal.modos)
%        graf.ganancias.(strcat("LP_",signal.modos(s))) = EDFA.(strcat("Nucleo",int2str(n))).salida.ganancias.(strcat("LP_",signal.modos(s)));
%        leyenda = strcat((strcat("LP",signal.modos(s))));
%        ejex = fliplr((3*10^8./(signal.lambda.(strcat("LP_",signal.modos(s))).*1e9)));
%        plot(ejex/1000,fliplr(graf.ganancias.(strcat("LP_",signal.modos(s)))) , '-o' , 'DisplayName',leyenda ) ; hold on ; title('Ganancias') ; legend(); grid on ;
%        %ax.ColorOrderIndex = s;
%        xlabel('Frecuencia [THz]') ; ylabel(ylab)
%        set(gca,'xdir','reverse')
%    end ; clear s ejex leyenda;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% archivos de salida
% aux_out1=[graf.freq_ase' graf.ASE_Spectrum.LP_01(:,1) graf.ASE_Spectrum.LP_01(:,2) graf.ASE_Spectrum.LP_11_a(:,1) graf.ASE_Spectrum.LP_11_a(:,2)]; 
% aux_out1=[Frequency_gridS' Wavelength_gridS' EDFA.Nucleo1.salida.signal.potencia_dBm.LP_01 EDFA.Nucleo1.salida.signal.potencia_dBm.LP_11_a];
% aux_out2= [graf.ganancias.LP_01' graf.ganancias.LP_11_a'];
% 
% save 'C:\Users\HP -\Documents\UTFSM\Investigación Doctorado\Simulaciones Matlab\registro de resultados\Datos importados del simulador FMF\noise_simul.dat' aux_out1 -ascii
% save 'C:\Users\HP -\Documents\UTFSM\Investigación Doctorado\Simulaciones Matlab\registro de resultados\Datos importados del simulador FMF\signal_simul.dat' aux_out1 -ascii
% save 'C:\Users\HP -\Documents\UTFSM\Investigación Doctorado\Simulaciones Matlab\registro de resultados\Datos importados del simulador FMF\gain_simul.dat' aux_out2 -ascii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Datos VPI
% p_pump_vpi=vpiedfamm.axialdist.Pump;
% p_sig_LP01_vpi=vpiedfamm.axialdist.Sig.LP_01;
% p_sig_LP11a_vpi=vpiedfamm.axialdist.Sig.LP_11a;
% p_ase_LP01_vpi=vpiedfamm.axialdist.ASE.LP_01;
% p_ase_LP11a_vpi=vpiedfamm.axialdist.ASE.LP_11a;
% 
% Z=linspace(0,3,101);
% colorlist=[1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0 0 0; 0.85 0.3250 0.980; 0.4940 0.184 0.5560; 0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840];
% 
% % Figuras
% % Pump
% figure(1)
% plot(Z,EDFA.Nucleo1.pump.Potencia_dBm.LP_01,'DisplayName','Matlab'); hold on
% plot(Z,p_pump_vpi,'DisplayName','VPI')
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend()
% title('Power axial distibution Pump (LP01,980 [nm], 250 [mW] ')
% 
% %Signal primeras 10 frecuencias LP01
% figure(2)
% for aux_i=1:length(colorlist)-1
%     plot(Z,EDFA.Nucleo1.signal.Potencia_dBm.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(Frequency_gridS(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i,:));hold on
%     plot(Z,p_sig_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Signal (LP01,-10[dBm])')
% 
% %Signal últimas 10 frecuencias LP01
% figure(3)
% for aux_i=1+10:length(colorlist)-1+10
%     plot(Z,EDFA.Nucleo1.signal.Potencia_dBm.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(Frequency_gridS(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-10,:));hold on
%     plot(Z,p_sig_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-10,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Signal (LP01,-10[dBm])')
% 
% % ASE forward LP01
% figure(4)
% for aux_i=1:length(colorlist)
%     plot(Z,EDFA.Nucleo1.PaseF.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i,:));hold on
%     plot(Z,p_ase_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP01,-200 [dBm])')
% 
% figure(5)
% for aux_i=1+11:length(colorlist)+11
%     plot(Z,EDFA.Nucleo1.PaseF.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-11,:));hold on
%     plot(Z,p_ase_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-11,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP01,-200 [dBm])')
% 
% figure(6)
% for aux_i=1+11*2:length(colorlist)+11*2
%     plot(Z,EDFA.Nucleo1.PaseF.LP_01(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-11*2,:));hold on
%     plot(Z,p_ase_LP01_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-11*2,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP01,-200 [dBm])')
% 
% % LP_11a
% %Signal LP01a
% figure(7)
% for aux_i=1:length(colorlist)-1
%     plot(Z,EDFA.Nucleo1.signal.Potencia_dBm.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(Frequency_gridS(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i,:));hold on
%     plot(Z,p_sig_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Signal (LP11a,-10[dBm])')
% 
% %Signal últimas 10 frecuencias LP01
% figure(8)
% for aux_i=1+10:length(colorlist)-1+10
%     plot(Z,EDFA.Nucleo1.signal.Potencia_dBm.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(Frequency_gridS(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-10,:));hold on
%     plot(Z,p_sig_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-10,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Signal (LP11a,-10[dBm])')
% 
% % ASE forward LP01
% figure(9)
% for aux_i=1:length(colorlist)
%     plot(Z,EDFA.Nucleo1.PaseF.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i,:));hold on
%     plot(Z,p_ase_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP11a,-200 [dBm])')
% 
% figure(10)
% for aux_i=1+11:length(colorlist)+11
%     plot(Z,EDFA.Nucleo1.PaseF.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-11,:));hold on
%     plot(Z,p_ase_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-11,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP11a,-200 [dBm])')
% 
% figure(11)
% for aux_i=1+11*2:length(colorlist)+11*2
%     plot(Z,EDFA.Nucleo1.PaseF.LP_11_a(aux_i,:),'DisplayName',['Matlab',num2str(EDFA.Nucleo1.Pase.frequency_vector(aux_i)/1e12),'[THz]'],'Color',colorlist(aux_i-11*2,:));hold on
%     plot(Z,p_ase_LP11a_vpi(:,aux_i),'DisplayName','VPI','Color',colorlist(aux_i-11*2,:)); hold on
% end
% xlabel('Axial Position [m]') ; ylabel('Power [dBm]')
% legend('Location','best')
% title('Power axial distibution Forward ASE (LP11a,-200 [dBm])')
% 


%% GUARDAR datos

S01 = EDFA.(strcat("Nucleo",int2str(1))).signal.Potencia_dBm.(ModoS(1));
S11_a = EDFA.(strcat("Nucleo",int2str(1))).signal.Potencia_dBm.(ModoS(2));

frecuencias = Frequency_gridS';
lambdas = Wavelength_gridS';
G1 = graf.ganancias.LP_01';
G2 = graf.ganancias.LP_11_a';

DAT = [lambdas G1 lambdas G2 frecuencias];

