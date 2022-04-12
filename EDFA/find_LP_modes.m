% ------------------------------------------------------------------------
% find_LP_modes
% ------------------------------------------------------------------------
% Michael Hughes 
% Applied Optics Group, University of Kent
% research.kent.ac.uk/applied-optics
%
% License: BSD [https://opensource.org/licenses/BSD-3-Clause]
%
% ------------------------------------------------------------------------
% Finds all LP modes for a step index fibre by finding solutions to the
% equation:
%
% besselj(l, u) / (u * besselj(l - 1, u))
%           = besselk(l, w) / (w * besselk(l - 1, w))
%
% where w = sqrt(v^2-u^2) and 
%       v = (2*pi*coreRadius/wavelength)*sqrt(nCore^2-nCladding^2)  
%
% Usage:
%     [solution] = find_LP_modes(coreRadius, nCore, nCladding, wavelength)
%
% Parameters:
%   coreRadius    : radius of the core in microns
%   nCore         : core refractive index
%   ncladding     : cladding refractive index
%   wavelength    : wavlength of light in microns
%
% Returns
%   solutions     : 1D array struct with parameters:
%       .u            : argument for core bessel function (microns)
%       .w            : argument for cladding bessel function (microns)
%       .l            : l number of mode 
%       .m            : m number of mode
%       .beta         : effective wavenumber of mode
%   totalNumModes : The number of modes including orientations (i.e. the
%                   number of modes with m > 1. 
% ------------------------------------------------------------------------
function [solution, totalNumModes] = find_LP_modes(coreRadius, nCore, nCladding, wavelength)

    % Calculate fibre V number
    v = (2*pi*coreRadius/wavelength)*sqrt(nCore^2-nCladding^2);

    % Calculate wavenumber 
    k = 2 * pi / wavelength;
          
    % Initialise return
    solution.u = [];
    solution.l = [];
    solution.w = [];
    solution.beta = [];
    solution.m = [];
        
    % Sets the coarse search parameters
    fineNPoints = 1000;                     % May need to be larger for very large V fibres.
    fineURange = linspace(0,v -.1, fineNPoints);
    
    % Initialise parameters
    iL = 0;
    residual = zeros(1, length(fineURange));
    signChange = zeros(1, length(fineURange));
   
    while 1==1            
    
        % Search for approximate solutions where there is a sign change and
        % the second derivative is positive
        for i = 1:length(fineURange)
            residual(i) = calculate_LP_mismatch(coreRadius, nCore, nCladding, wavelength, iL, fineURange(i));
            if i > 1
                signChange(i) = (sign(residual(i)) ~= sign(residual(i-1))) && (residual(i) > residual(i-1));
            end 
        end
        
        % Pull out all the sign changes we found
        coarseSolutions = find(signChange);
        
        % Search for exact solutions around each sign change
        for iSolution = 1:length(coarseSolutions)
                  
            % Search where we know there is a change of sign
            minRange = fineURange(coarseSolutions(iSolution) - 1);
            maxRange = fineURange(coarseSolutions(iSolution));
            
            % Find the exact point of the change of sign
            solution(end+1).u = fzero(@(u)calculate_LP_mismatch(coreRadius, nCore, nCladding, wavelength, iL, u),[minRange, maxRange]);
                        
            % Record the solution
            solution(end).l = iL;
            solution(end).beta = sqrt(k^2 * nCore^2 - (solution(end).u/coreRadius)^2);
            solution(end).w = sqrt(v.^2-solution(end).u.^2);        
            solution(end).m = iSolution;
        
        end
        
        % If there are no solutions for this l then we are done
        if isempty(coarseSolutions)
            break;
        end

        % Move to the next l index
        iL = iL + 1;
  
    end
   
    solution(1) = [];
    
    % Calculate total number of modes, taking acount of rotational
    % variance of modes with m > 1
    totalNumModes = length(solution) + sum(cell2mat({solution.m}) > 1);
    
end