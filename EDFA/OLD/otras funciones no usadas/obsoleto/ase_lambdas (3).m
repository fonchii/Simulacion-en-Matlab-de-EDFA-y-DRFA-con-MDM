%% Vector de lambdas para ruido ASE
% reescalado para mayor presicion en torno a longitudes de onda de señal

function ase_wl = ase_lambdas(lambda_s)  

% if length(lambda_s)>1
%     lambda_ase = 1520e-9 : 0.5*1e-9: (lambda_s(1) - 0.5*1e-9);
%     for i = 1:1:(length(lambda_s)-1)
%         ase = linspace(max(lambda_ase) , lambda_s(i) + 0.5*1e-9, 10);
%         ase2 = max(ase) : 0.5*1e-9: (lambda_s(i+1) - 0.5*1e-9);
%         lambda_ase = [lambda_ase ase ase2];
%     end
% else 
%     lambda_ase = 1520e-9 : 0.5*1e-9: (lambda_s - 0.5*1e-9);
%     ase = linspace(max(lambda_ase) , lambda_s + 0.5*1e-9, 10);
% end
% if length(lambda_s)>1
%     ase = linspace(max(lambda_ase) , lambda_s(end) + 0.5*1e-9, 10);
% end
%     
% ase = max(lambda_ase) : 0.5*1e-9 : 1570e-9;
% lambda_ase = [lambda_ase ase];  

if length(lambda_s)>1
%     %lambda_ase = 1490e-9 : 0.5*1e-9: (lambda_s(1) - 0.5*1e-9);     % paso 0.5 nm
%     lambda_ase = 1490e-9 : 1*1e-9: (lambda_s(1) - 0.5*1e-9);        % paso 1 nm
%     ase = lambda_s(1)-0.5*1e-9 : 0.1*1e-9 : lambda_s(end)+0.5*1e-9; % paso 0.1 nm entre señales
%     %ase2 = ase(end):0.5*1e-9:1600e-9;                              % paso 0.5 nm
%     ase2 = ase(end):1*1e-9:1600e-9;                                 % paso 1 nm
%     ase_wl =  [lambda_ase ase ase2];
    lambdas = [];
    lambda_ase = 1520e-9 : 1*1e-9: (lambda_s(1) - 0.5*1e-9);
    for i=1:length(lambda_s)-1 % entre longitudes de onda entregadas
        ase = lambda_s(i)-0.5*1e-9 : 0.1*1e-9 : lambda_s(i)+0.5*1e-9;
        ase2 = ase(end) : 1*1e-9 : lambda_s(i+1)-0.5*1e-9;
        lambdas = [lambdas ase ase2];
    end
    ase = lambda_s(end)-0.5*1e-9 : 0.1*1e-9 : lambda_s(end)+0.5*1e-9;
    ase2 = ase(end):1*1e-9:1600e-9;
    ase_wl =  [lambda_ase lambdas ase ase2];

else
    %lambda_ase = 1520e-9 : 0.5*1e-9: (lambda_s(1) - 0.5*1e-9);
    lambda_ase = 1520e-9 : 1*1e-9: (lambda_s(1) - 0.5*1e-9);
    ase = lambda_s-0.5*1e-9 : 0.1*1e-9 : lambda_s+0.5*1e-9;
    %ase2 = ase(end):0.5*1e-9:1570e-9;
    ase2 = ase(end):1*1e-9:1570e-9;
    ase_wl =  [lambda_ase ase ase2];


end
