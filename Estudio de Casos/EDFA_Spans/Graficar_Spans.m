% Graficos de Spans
close all; clear all; clc
load("Span.mat")

set( gcf,'PaperSize',[29.7 21.0], 'PaperPosition',[0 0 29.7 21.0])
% SAVE:
% print -dpdf 'NAME'


figure(1)
t = tiledlayout(1,7,'TileSpacing','none'); % compact - espacio pequeño ; none - sin espacio
bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
bgAx.Layout.TileSpan = [1 2];

GrafSpans.z_total = [Span.EDFA1.Nucleo1.z]; 
GrafSpans.Signal_total = Span.EDFA1.Nucleo1.signal.Potencia_dBm;
GrafSpans.ASE_total = Span.EDFA1.Nucleo1.Pap;

z = Span.EDFA1.Nucleo1.z; Nspans = 3; LargoFibra = 100;

ModoS = ["01" "11_a"];

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

fibra1_01 = [Amp1.LP_01(:,end) Amp2.LP_01(:,1)];
fibra2_01 = [Amp2.LP_01(:,end) Amp3.LP_01(:,1)];
fibra3_01 = [Amp3.LP_01(:,end) Amp4.LP_01(:,1)];
z_fibra = [0 100]; lambdas = Span.EDFA4.Nucleo1.signal.lambdas.*1e9;

ymin = -15 ; ymax = 25;

% Grafico 1 - EDFA 1
ax1 = axes(t);
plot(ax1,z,Amp1.LP_01(1,:) , "DisplayName",strcat(num2str(lambdas(1))," nm")) ; hold on
plot(ax1,z,Amp1.LP_01(21,:) , "DisplayName",strcat(num2str(lambdas(21)), " nm"))
plot(ax1,z,Amp1.LP_01(31,:) , "DisplayName",strcat(num2str(lambdas(31)), " nm"))
xline(ax1,3,':');
ax1.Box = 'off';
xlim(ax1,[0 3])
xlabel(ax1, 'EDFA 1 [m]')
ylim(ax1,[ymin ymax])
%legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 3,"FontSize",9)

% Grafico 2 - FIBRA 1
ax2 = axes(t);
ax2.Layout.Tile = 2;
plot(ax2,z_fibra,fibra1_01(1,:)) ; hold on
plot(ax2,z_fibra,fibra1_01(21,:))
plot(ax2,z_fibra,fibra1_01(31,:))
xline(ax2,100,':');
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2,[0 100])
xlabel(ax2,'Fibra [km]')
ylim(ax2,[ymin ymax])



% Grafico 3 - EDFA 2
ax3 = axes(t);
ax3.Layout.Tile = 3;
plot(ax3,z,Amp2.LP_01(1,:)) ; hold on
plot(ax3,z,Amp2.LP_01(21,:))
plot(ax3,z,Amp2.LP_01(31,:))
xline(ax3,3,':');
ax3.YAxis.Visible = 'off';
ax3.Box = 'off';
xlim(ax3,[0 3])
xlabel(ax3,'EDFA 2 [m]')
ylim(ax3,[ymin ymax])

% Grafico 4 - FIBRA 2
ax4 = axes(t);
ax4.Layout.Tile = 4;
plot(ax4,z_fibra,fibra2_01(1,:) , "DisplayName",strcat(num2str(lambdas(1))," nm")) ; hold on
plot(ax4,z_fibra,fibra2_01(21,:) , "DisplayName",strcat(num2str(lambdas(21))," nm"))
plot(ax4,z_fibra,fibra2_01(31,:) , "DisplayName",strcat(num2str(lambdas(31))," nm"))
ax4.YAxis.Visible = 'off';
ax4.Box = 'off';
xlim(ax4,[0 100])
xlabel(ax4,'Fibra [km]')
ylim(ax4,[ymin ymax])
legend('Location', 'southoutside','Orientation','horizontal','Box','off', "NumColumns" , 3,"FontSize",9)

% Grafico 5 - EDFA 3
ax5 = axes(t);
ax5.Layout.Tile = 5;
plot(ax5,z,Amp3.LP_01(1,:)) ; hold on
plot(ax5,z,Amp3.LP_01(21,:))
plot(ax5,z,Amp3.LP_01(31,:))
xline(ax5,0,':');
xline(ax5,3,':');
ax5.YAxis.Visible = 'off';
ax5.Box = 'off';
xlim(ax5,[0 3])
xlabel(ax5,'EDFA 3 [m]')
ylim(ax5,[ymin ymax])


% Grafico 6 - FIBRA 3
ax6 = axes(t);
ax6.Layout.Tile = 6;
plot(ax6,z_fibra,fibra3_01(1,:)) ; hold on
plot(ax6,z_fibra,fibra3_01(21,:))
plot(ax6,z_fibra,fibra3_01(31,:))
xline(ax6,100,':');
ax6.YAxis.Visible = 'off';
ax6.Box = 'off';
xlim(ax6,[0 100])
xlabel(ax6,'Fibra [km]')
ylim(ax6,[ymin ymax])

% Grafico 7 - EDFA 4
ax7 = axes(t);
ax7.Layout.Tile = 7;
plot(ax7,z,Amp4.LP_01(1,:)) ; hold on
plot(ax7,z,Amp4.LP_01(21,:))
plot(ax7,z,Amp4.LP_01(31,:))
xline(ax7,45,':');
ax7.YAxis.Visible = 'off';
ax7.Box = 'off';
xlim(ax7,[0 3])
xlabel(ax7,'EDFA 4 [m]')
ylim(ax7,[ymin ymax])

title(t,'Propagación de Señal en Modo LP01',"FontSize",14)
ylabel(t,'Magnitud [dBm]')

