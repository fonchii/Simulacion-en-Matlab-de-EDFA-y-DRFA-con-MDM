%% Vector de lambdas para ruido ASE
% reescalado para mayor presicion en torno a longitudes de onda de señal

function ase_wl = ase_lambdas(lambda_s,flag)  
if nargin<2
    if length(lambda_s)>1
        
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
        
        lambda_ase = 1520e-9 : 1*1e-9: (lambda_s(1) - 0.5*1e-9);
        ase = lambda_s-0.5*1e-9 : 0.1*1e-9 : lambda_s+0.5*1e-9;
        ase2 = ase(end):1*1e-9:1570e-9;
        ase_wl =  [lambda_ase ase ase2];
    end

else 
    ase_wl = [0,0,0];

end
