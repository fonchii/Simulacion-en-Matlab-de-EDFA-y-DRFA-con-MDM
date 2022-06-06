clear all;clc;close all

A                   = load('RamanGainEfficiency_SMF28.dat');
g_R                 = smooth(A(:,2),1); %10
C_R                 = 0.5*(g_R)*1000; % Ajuste de polarización y paso a km
FF_gR               = A(:,1);
% % Raman Responce Function
Gr_fun              = @(f)interp1(FF_gR,C_R,f);

Gr = Gr_fun(FF_gR);

pump = linspace(1400,1500,12).*1e-9; % nm
pump = 3e8./pump ; % Hz

sign = pump - FF_gR; % Hz

lambdas=3e8./(linspace(1400,1700,1000).*1e-9);

%DeltaF = pump+sign;
for i = 1:length(sign(1,:))
    DeltaFun.(strcat("F",string(i))) = @(f)interp1(sign(:,i),Gr,f);
end
SumGain = zeros(1,length(lambdas));

magn = [0.7 0.2*ones(1,length(sign(1,:))-2) 1];
%magn = 2.5.*magn./sum(magn);

for i =1:length(lambdas)
    for j = 1:length(sign(1,:))
        if ~isnan(DeltaFun.(strcat("F",string(j)))(lambdas(i)))
            SumGain(i) = SumGain(i) + magn(j)*DeltaFun.(strcat("F",string(j)))(lambdas(i));
        end
    end
end

plot(fliplr(3e8./sign).*1e9 ,fliplr(Gr)) ; xlabel("Espectro de Ganancia [nm]") ; ylabel("gR [1/W*km]") ; hold on;
plot(fliplr(3e8./lambdas).*1e9, SumGain)
