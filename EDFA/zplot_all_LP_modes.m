% ------------------------------------------------------------------------
% plot_all_LP_modes
% ------------------------------------------------------------------------
% Michael Hughes 
% Applied Optics Group, University of Kent
%
% License: BSD [https://opensource.org/licenses/BSD-3-Clause]
% ------------------------------------------------------------------------
% Plots a series of fibre LP modes. The Sin and Cos rotations of each
% mode are provided in separate arrays. 
% 
% Usage:
%     [modeSin, modeCos] = plot_all_LP_modes(modes, coreRadius, ...
%                                 maxPlotRadius, gridSize)
%
% Parameters:
%   modes          : struct returned by find_LP_modes
%   coreRadius     : radius of fibre core (micfons)
%   maxPlotRadius  : radius at edge of grid (microns)
%   gridSize       : number of pixels in grid
% 
% Returns:
%   modeSin  : 3D array of dimensions (gridSize, gridSize, numModes)
%   modeCos  : 3D array of dimensions (gridSize, gridSize, numModes)
% ------------------------------------------------------------------------
function [modeSin, modeCos] = plot_all_LP_modes(modes, coreRadius, maxPlotRadius, gridSize)
   
    % Initialise output arrays
    modeSin = zeros(gridSize, gridSize, length(modes));
    modeCos = zeros(gridSize, gridSize, length(modes));
    
    for i = 1: length(modes)
       
        [mS, mC] = plot_LP_mode(modes(i).l, modes(i).u, modes(i).w, coreRadius, maxPlotRadius, gridSize);
        modeSin(:,:,i) = mS;
        modeCos(:,:,i) = mC;
        
    end        
            
end