% -------------------------------------------------------------------------
% example_general
% Fibre Optic LP Mode Solver and Simulator
% For Release 2 June 2020
% ------------------------------------------------------------------------
% Michael Hughes   m.r.hughes@kent.ac.uk
% Applied Optics Group, University of Kent
%
% License: BSD [https://opensource.org/licenses/BSD-3-Clause]
% -------------------------------------------------------------------------
% Demonstrates use of find_LP_modes and plot_all_LP_modes to find and plot
% LP modes of a step index fibre.
% -------------------------------------------------------------------------

clear
addpath('lib');

% Define the fibre characteristics and wavelength
nCore = 1.43;
nCladding = 1.42;
wavelength = 0.5;    % microns
coreRadius = 15;     % microns

% Define the plot area
maxPlotRadius = coreRadius * 1.2;  % Lets us see the power in the cladding
gridSize = 300;                  % pixels

% Find all the LP modes
[modes, nModes] = find_LP_modes(coreRadius, nCore, nCladding, wavelength);

% Report how many modes are found
fprintf('Estimating %d LP modes from V number. \n', round(est_num_modes(wavelength, coreRadius, nCore, nCladding) /2));
fprintf('Found %d LP modes (%d including rotations). \n', length(modes), nModes);

% Calculate the 2D field amplitudes of the LP modes
[modeSin, modeCos] = plot_all_LP_modes(modes, coreRadius, maxPlotRadius, gridSize);

% Calculate power in core for each mode
powerInCore = power_in_core(modeSin, coreRadius, maxPlotRadius);

% Couple a Gaussian beam into fibre
pixelSize = (maxPlotRadius * 2) / gridSize;
imSize = gridSize * pixelSize;
beamRadius = 5;
inIntensity = fspecial('Gaussian', gridSize, beamRadius /2 / pixelSize) ;
inField = sqrt(inIntensity);
[modeCouplingSin, modeCouplingCos, modeCouplingIntensity] = couple_beam(inField, modeSin, modeCos);
coupledPower = sum(modeCouplingIntensity);

% Select which mode to display
displayMode = 12;

% Display a mode amplitude
zLim = max(abs(modeSin(:)));
figure(1); imagesc(modeSin(:,:,displayMode),[-zLim, zLim])
axis equal;
colormap(amplitude_colour_map);
title(['Amplitude plot for mode ', num2str(displayMode)]); 


% Display a mode intensity
figure(2); imagesc(modeSin(:,:,displayMode).^2)
axis equal;
colormap('hot')
title(['Intensity plot for mode ', num2str(displayMode)]); 


% Display a rotated mode intensity
figure(3); imagesc((modeSin(:,:,displayMode) + modeCos(:,:,displayMode)).^2)
axis equal;
colormap('hot')
title(['Intensity plot for rotated mode ', num2str(displayMode)]); 

% Display a mode radial profile
figure(4); [profile, intProfile, rPoints] = plot_LP_mode_profile(modes(displayMode).l, modes(displayMode).u, modes(displayMode).w, coreRadius, maxPlotRadius, gridSize);
plot(rPoints, profile); hold on; plot(rPoints, intProfile); hold off;
title(['Radial profile for mode ', num2str(displayMode)]); xlabel('Radial position (microns)'); ylabel('Relative power/amplitude');
legend('Amplitude', 'Power');

% Display power in core for each mode
figure(5); plot(powerInCore); title('Fraction of power in core for each mode'); xlabel('Mode No.'); ylabel('Fraction of power in core');

% Display mode coupling
figure(6); plot(modeCouplingIntensity); title('Distribution of coupled power across modes'); xlabel('Mode No.'); ylabel('Fraction of power in mode');


