function [Gamma,Beta_0] = norm_intensity2(fibra,modos,lambda_s)

% Find all the LP modes
param = modesparam(modos,lambda_s,fibra) ;

if isstruct(param) == 0
    error('No es posible transmitir el modo en la fibra')
end

Beta_0=param.beta*1e6; %[1/m]

% Define the plot area
maxPlotRadius = (fibra.radio*1e6) * 1.2;  % Lets us see the power in the cladding
gridSize = 300;                  % pixels



% Calculate the 2D field amplitudes of the LP modes
[modeSin, ~] = plot_all_LP_modes(param.allsolution, fibra.radio*1e6, maxPlotRadius, gridSize);

% Calculate power in core for each mode
Gamma = power_in_core(modeSin, fibra.radio*1e6, maxPlotRadius);

end


