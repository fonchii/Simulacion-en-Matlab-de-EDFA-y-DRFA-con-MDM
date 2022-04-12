%% Frecuencias de señal
% entregar cantidad de frecuencias a multiplexar y centro en nm
function freqs = wavelength(cantidad,centro)  
if(centro<1) centro = centro * 1e9; end     % escrito en nm?
for i = 1:cantidad
    if(i == 1)
        freqs(i) = centro - floor(cantidad/2);
    else
        freqs(i) = freqs(i-1) + 1;
    end
end
freqs = freqs.*1e-9;
    