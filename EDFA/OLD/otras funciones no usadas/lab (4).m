A=1;T=1;
t = linspace(0,20,1000);
unistep(t>=T)=1;
f=60; c = 10*1e-9; l = 10*1e-6; 
C = 1/(1i*c*2*pi*f); L = 1i*l*2*pi*f;
Zserie = C+L ; Zpar = (1/C + 1/L)^-1;
Z0 = 50;
vs = A/2 - A/2*(unistep) + A*exp((-t+T)/(C*L)^0.5);
plot(t,(vs))
%%
clc
syms C L Z0 R w 
% ZL = w*L ; Zc = 1/(w*C)
%Zrcp = R / sqrt(1+w^2*R^2*C^2);     %rc paralelo
%Zrlp = w*L/sqrt(R^2+w^2*L^2);       %rl paralelo
Zrcp = (1/R + 1i*w*C)^-1;
Zrlp = (1/R + 1/(1i*w*L));
Zrcs = R+1/(1i*w*C);                   %rc serie
Zrls = R+1i*w*L;                       %rl serie

gamma1 = simplify((Zrcp - Z0) / (Zrcp + Z0)); %rc paralelo
gamma2 = simplify((Zrlp - Z0) / (Zrlp + Z0)); %rl paralelo
gamma3 = simplify((Zrcs - Z0) / (Zrcs + Z0)); %rc serie
gamma4 = simplify((Zrls - Z0) / (Zrls + Z0)); %rl serie

%%
clc; close all; clear all
w = 2*pi*linspace(100e3,1.5e9,100000);
%C = 10e-12; L = 10e-6; R = 50; Z0 = 50;
R=80;  C=100e-12;   L=100e-6; Z0 = 50;
f = linspace(100e3,1.5e9,100000);

%Zrcp = R ./ sqrt(1+w.^2.*R^2*C^2);     %rc paralelo
%Zrlp = w.*L./sqrt(R^2+w.^2.*L^2);       %rl paralelo
Zrcp = ((1/R + (1i.*w.*C))).^-1;
Zrlp = ((1/R + 1./1i.*w.*L)).^-1;
Zrcs = R+1./(1i.*w.*C);                   %rc serie
Zrls = R+1i.*w.*L;                       %rl serie

gamma1 = (Zrcp/Z0 - Z0) ./ (Zrcp/Z0 + Z0); %rc paralelo
gamma2 = (Zrlp/Z0 - Z0) ./ (Zrlp/Z0 + Z0); %rl paralelo
gamma3 = (Zrcs/Z0 - Z0) ./ (Zrcs/Z0 + Z0); %rc serie
gamma4 = (Zrls/Z0 - Z0) ./ (Zrls/Z0 + Z0); %rl serie

figure(1)
smithplot(f,gamma1) ; title('\gamma RC Paralelo')

figure(2)
smithplot(f,gamma2) ; title('\gamma RL Paralelo')

figure(3)
smithplot(f,gamma3) ; title('\gamma RC Serie')

figure(4)
smithplot(f,gamma4) ; title('\gamma RL Serie')


