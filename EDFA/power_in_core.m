% ------------------------------------------------------------------------
% power_in_core
% version 1.0
% ------------------------------------------------------------------------
% Michael Hughes    m.r.hughes@kent.ac.uk
% Applied Optics Group, University of Kent
%
% License: BSD [https://opensource.org/licenses/BSD-3-Clause]
% ------------------------------------------------------------------------
% Computes the power in the core from an intensity plot of power 
% distribution, given a core radius. Eessentially integrates power
% inside a circle and expresses as fraction of total power.
% WARNING: Will be innaccuate if maxPlotRadius is similar to or less than 
% coreRadius.
% ------------------------------------------------------------------------
% Usage:
%  powerInCore = power_in_core(modePlot, coreRadius, maxPlotRadius)
%
% Parameters:
%   modePlot:      2D array, map of intensity
%   coreRadius:    radius of core in physical units
%   maxPlotRadius: semi-diameter of modePlot in physical unnits
%
% Returns
%   powerInCore:   scalar, should be < 1
% ------------------------------------------------------------------------

function powerInCore = power_in_core(modePlot, coreRadius, maxPlotRadius)

    gridSize = size(modePlot,1);
    gridCentre = gridSize / 2;

    nModes = size(modePlot,3);
    powerInCore = zeros(nModes,1);
    
    % Binary mask of all pixels that are inside core radius
    gridPoints = (1:gridSize) - gridCentre;
    [xMesh, yMesh] = meshgrid(gridPoints, gridPoints);
    [~, rad] = cart2pol(xMesh,yMesh);
    rad = rad ./ (gridSize / 2) * maxPlotRadius;
    inCore = (rad <= coreRadius) & (rad <= maxPlotRadius);
   
    % Calculate power inside core
    for i = 1: nModes    
        cMode = modePlot(:,:,i);
        powerInCore(i) = sum(sum(abs(cMode(inCore)).^2))./ sum(sum(abs(cMode).^2));
    end
    
end