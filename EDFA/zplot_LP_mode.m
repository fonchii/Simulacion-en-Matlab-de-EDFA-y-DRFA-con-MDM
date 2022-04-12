% ------------------------------------------------------------------------
% plot_LP_mode
% ------------------------------------------------------------------------
% Michael Hughes    m.r.hughes@kent.ac.uk
% Applied Optics Group, University of Kent
%
% License: BSD [https://opensource.org/licenses/BSD-3-Clause]
% ------------------------------------------------------------------------
% Plots a normalised fibre LP mode (sin and cos version).
% 
% Usage:
%     [modeCos, modeSin] = plot_LP_mode(order, u, w, coreRadius, maxPlotRadius, gridSize)
%
% Parameters:
%   coreRadius     : radius of the core (microns)
%   order          : mode order
%   u              : mode core u term (microns)
%   w              : mode core w term (microns)
%   maxPlotRadius  : radius at edge of grid (i.e. grid is 2x this)
%   gridSize       : number of pixels in grid
%
% Returns:
%   modeSin        : 2D array of mode amplitude (sin)
%   modeCos        : 2D array of mode amplitude (cos)
% ------------------------------------------------------------------------
function [modeCos, modeSin, profile, pixelSize , num , den] = zplot_LP_mode(order, u, w, coreRadius, maxPlotRadius, gridSize)
    if nargin<1
        fibra.radio=25e-6; fibra.n1=1.46; fibra.n2=1.455;
        modos="21"; lambda_s = 1550e-9 ; 
        param = modesparam(modos,lambda_s,fibra) ;
        order = str2num(extract(modos,1)) ; u = param.u ; w = param.w ; coreRadius = fibra.radio*1e6 ; 
        maxPlotRadius = coreRadius*1.5 ; gridSize = 300;
    end
    
    % Find centre of grid
    gridCentre = gridSize / 2;
    
    %Initialise output array
    modeCos = zeros(gridSize, gridSize);
    modeSin = zeros(gridSize, gridSize);
    profile = zeros(gridSize,1);
    
    % Calculate grid points
    gridPoints = (1:gridSize) - gridCentre;
    [xMesh, yMesh] = meshgrid(gridPoints, gridPoints);
    
    % Convert  to polar co-ordinates
    [angle, rad] = cart2pol(xMesh,yMesh);
    rad = rad ./ (gridSize / 2) * maxPlotRadius;
    
    % Return the size of each pixel in the plot
    pixelSize = maxPlotRadius * 2 ./ gridPoints;
    
    % Calculate the cos term
    cosTerm = cos(order .* angle);
    sinTerm = sin(order .* angle);
         
      
    % Calculate the core and cladding 2D functions
    coreBessel = besselj(order,u ./ coreRadius .* rad)/besselj(order,u);
    claddingBessel = besselk(order, w ./ coreRadius .* rad)/besselk(order,w);
    
    % Work out which grid points are in core and which in cladding
    inCore = (rad <= coreRadius) & (rad <= maxPlotRadius);
    inCladding = (rad > coreRadius) & (rad <= maxPlotRadius);
    
    % Calculate sin and cos versions of core field
    modeCos(inCore) = coreBessel(inCore) .* cosTerm(inCore);
    modeSin(inCore) = coreBessel(inCore) .* sinTerm(inCore); 
    
    % Calculate sin and cos versions of cladding field
    modeCos(inCladding) = claddingBessel(inCladding) .* cosTerm(inCladding);
    modeSin(inCladding) = claddingBessel(inCladding) .* sinTerm(inCladding);    
    
    % Normalise to total of 1.
    modeCos = modeCos / sqrt(sum(sum(modeCos.^2)));
    modeSin = modeSin / sqrt(sum(sum(modeSin.^2)));

    modeCos(isnan(modeCos)) = 0;
    modeSin(isnan(modeSin)) = 0;    
    
    num = sum ( sqrt ( ( modeCos(inCore).^2 + modeSin(inCore).^2 ) /2 ) ) ;
    den = num + sum (sqrt( (modeCos(inCladding).^2 + modeSin(inCladding).^2)/2 ) ) ;
    close all;
    %imagesc(sqrt(modeCos.^2 + modeSin.^2) ) ; hold on
    subplot(2,1,1) ; imagesc(sqrt(modeCos.^2 ) ); title('Coseno') ;
    viscircles( [gridSize/2,gridSize/2] ,gridSize*1/3)
    subplot(2,1,2) ; imagesc(sqrt(modeSin.^2 ) ) ; title('Seno') ;
    viscircles( [gridSize/2,gridSize/2] ,gridSize*1/3)
end