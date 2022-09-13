% Graficos de Spans
close all; clear all; clc
load("Datos_Span5m_100mw.mat")

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% SAVE:
% print -dpdf 'NAME' % Name: Cascada_1555nm


figure(1)
t = tiledlayout(1,7,'TileSpacing','none'); % compact - espacio pequeño ; none - sin espacio
bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
bgAx.Layout.TileSpan = [1 2];

GrafSpans.z_total = [Span.EDFA1.Nucleo1.z]; 
GrafSpans.Signal_total = Span.EDFA1.Nucleo1.signal.Potencia_dBm;
GrafSpans.ASE_total = Span.EDFA1.Nucleo1.Pap;

z = Span.EDFA1.Nucleo1.z; Nspans = 3; LargoFibra = 90; LargoEDFA = 5;

ymin = -30 ; ymax = 0;

ModoS = ["01" "11_a"];
modos_LP = ["LP_01" "LP_11_a"];

% for span = 1:Nspans
%     GrafSpans.z_total = [GrafSpans.z_total (z+(0.1*LargoFibra+GrafSpans.z_total(end)))];
%     for i=1:2
%         GrafSpans.Signal_total.(strcat("LP_",ModoS(i))) = [GrafSpans.Signal_total.(strcat("LP_",ModoS(i))) Span.(strcat("EDFA",num2str(span+1))).Nucleo1.signal.Potencia_dBm.(strcat("LP_",ModoS(i)))];
%         GrafSpans.ASE_total.(strcat("LP_",ModoS(i))) = [GrafSpans.ASE_total.(strcat("LP_",ModoS(i))) Span.(strcat("EDFA",num2str(span+1))).Nucleo1.Pap.(strcat("LP_",ModoS(i)))];
%     end ; clear i;
% end

Amp1 = Span.EDFA1.Nucleo1.signal.Potencia_dBm;
Amp2 = Span.EDFA2.Nucleo1.signal.Potencia_dBm;
Amp3 = Span.EDFA3.Nucleo1.signal.Potencia_dBm;
Amp4 = Span.EDFA4.Nucleo1.signal.Potencia_dBm;

% LP01
fibra1_01 = [Amp1.(modos_LP(1))(:,end) Amp2.(modos_LP(1))(:,1)];
fibra2_01 = [Amp2.(modos_LP(1))(:,end) Amp3.(modos_LP(1))(:,1)];
fibra3_01 = [Amp3.(modos_LP(1))(:,end) Amp4.(modos_LP(1))(:,1)];

% LP11a
fibra1_21 = [Amp1.(modos_LP(2))(:,end) Amp2.(modos_LP(2))(:,1)];
fibra2_21 = [Amp2.(modos_LP(2))(:,end) Amp3.(modos_LP(2))(:,1)];
fibra3_21 = [Amp3.(modos_LP(2))(:,end) Amp4.(modos_LP(2))(:,1)];

% % LP21a
% fibra1_02 = [Amp1.(modos_LP(3))(:,end) Amp2.(modos_LP(3))(:,1)];
% fibra2_02 = [Amp2.(modos_LP(3))(:,end) Amp3.(modos_LP(3))(:,1)];
% fibra3_02 = [Amp3.(modos_LP(3))(:,end) Amp4.(modos_LP(3))(:,1)];
z_fibra = [0 LargoFibra]; lambdas = Span.EDFA4.Nucleo1.signal.lambdas.*1e9;


% Grafico 1 - EDFA 1
ax1 = axes(t);
plot(ax1,z,Amp1.(modos_LP(1))(26,:) , "DisplayName","Modo LP11a") ; hold on
plot(ax1,z,Amp1.(modos_LP(2))(26,:) , "DisplayName","Modo LP21a") ; 
%plot(ax1,z,Amp1.(modos_LP(3))(26,:) , "DisplayName","Modo LP02") ; 
xline(ax1,LargoEDFA,':');
ax1.Box = 'off';
xlim(ax1,[0 LargoEDFA])
xlabel(ax1, 'EDFA 1 [m]')
ylim(ax1,[ymin ymax])
%legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 3,"FontSize",9)

% Grafico 2 - FIBRA 1
ax2 = axes(t);
ax2.Layout.Tile = 2;
plot(ax2,z_fibra,fibra1_01(26,:)) ; hold on
plot(ax2,z_fibra,fibra1_21(26,:)) ;
%plot(ax2,z_fibra,fibra1_02(26,:)) ;
xline(ax2,LargoFibra,':');
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2,[0 LargoFibra])
xlabel(ax2,'Fibra [km]')
ylim(ax2,[ymin ymax])



% Grafico 3 - EDFA 2
ax3 = axes(t);
ax3.Layout.Tile = 3;
plot(ax3,z,Amp2.(modos_LP(1))(26,:)) ; hold on
plot(ax3,z,Amp2.(modos_LP(2))(26,:)) ;
%plot(ax3,z,Amp2.(modos_LP(3))(26,:)) ;
xline(ax3,0,':');
xline(ax3,LargoEDFA,':');
ax3.YAxis.Visible = 'off';
ax3.Box = 'off';
xlim(ax3,[0 LargoEDFA])
xlabel(ax3,'EDFA 2 [m]')
ylim(ax3,[ymin ymax])

% Grafico 4 - FIBRA 2
ax4 = axes(t);
ax4.Layout.Tile = 4;
plot(ax4,z_fibra,fibra2_01(26,:) , "DisplayName","Modo LP01") ; hold on
plot(ax4,z_fibra,fibra2_21(26,:) , "DisplayName","Modo LP11a") ;
%plot(ax4,z_fibra,fibra2_02(26,:) , "DisplayName","Modo LP02") ;
ax4.YAxis.Visible = 'off';
ax4.Box = 'off';
xlim(ax4,[0 LargoFibra])
xlabel(ax4,'Fibra [km]')
ylim(ax4,[ymin ymax])
legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 3,"FontSize",9)

% Grafico 5 - EDFA 3
ax5 = axes(t);
ax5.Layout.Tile = 5;
plot(ax5,z,Amp3.(modos_LP(1))(26,:)) ; hold on
plot(ax5,z,Amp3.(modos_LP(2))(26,:)) ;
%plot(ax5,z,Amp3.(modos_LP(3))(26,:)) ;
xline(ax5,0,':');
xline(ax5,LargoEDFA,':');
ax5.YAxis.Visible = 'off';
ax5.Box = 'off';
xlim(ax5,[0 LargoEDFA])
xlabel(ax5,'EDFA 3 [m]')
ylim(ax5,[ymin ymax])


% Grafico 6 - FIBRA 3
ax6 = axes(t);
ax6.Layout.Tile = 6;
plot(ax6,z_fibra,fibra3_01(26,:)) ; hold on
plot(ax6,z_fibra,fibra3_21(26,:)) ;
%plot(ax6,z_fibra,fibra3_02(26,:)) ;
xline(ax6,LargoFibra,':');
ax6.YAxis.Visible = 'off';
ax6.Box = 'off';
xlim(ax6,[0 LargoFibra])
xlabel(ax6,'Fibra [km]')
ylim(ax6,[ymin ymax])

% Grafico 7 - EDFA 4
ax7 = axes(t);
ax7.Layout.Tile = 7;
plot(ax7,z,Amp4.(modos_LP(1))(26,:)) ; hold on
plot(ax7,z,Amp4.(modos_LP(2))(26,:)) ;
%plot(ax7,z,Amp4.(modos_LP(3))(26,:)) ;
xline(ax7,0,':');
ax7.YAxis.Visible = 'off';
ax7.Box = 'off';
xlim(ax7,[0 LargoEDFA])
xlabel(ax7,'EDFA 4 [m]')
ylim(ax7,[ymin ymax])

title(t,'Propagación de Señal a 1555 nm',"FontSize",14)
ylabel(t,'Potencia [dBm]')

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
print -dpdf 'Cascada_1555nm'

