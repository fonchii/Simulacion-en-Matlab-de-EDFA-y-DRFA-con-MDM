% ------------------------------------------------------------------------
% calculate_LP_mismatch
% ------------------------------------------------------------------------
% Michael Hughes  m.r.hughes@kent.ac.uk
% Applied Optics Group, University of Kent
% research.kent.ac.uk/applied-optics
% ------------------------------------------------------------------------
% License: BSD [https://opensource.org/licenses/BSD-3-Clause]
% ------------------------------------------------------------------------
% Calculates:
%
% besselj(order, u) / (u * besselj(order - 1, u)) - 
%            besselk(order, w) / (w * besselk(order - 1, w))
%
% for the specified value of u, where :
% w = sqrt(v^2-u^2) and 
%       v = (2*pi*coreRadius/wavelength)*sqrt(nCore^2-nCladding^2)  
%
% Usage:
%     [residual, coreTerm, claddingTerm] = 
%      = calculate_LP_mismatch(coreRadius, nCore, nCladding, wavelength, order, u)
%
% Parameters:
%   coreRadius    : radius of the core in microns
%   nCore         : core refractive index
%   nCladding     : cladding refractive index
%   wavelength    : wavlength of light in microns
%   order         : mode l
%   u             : mode trial u 
%
% Returns:
%   residual      : difference between core and cladding terms, when
%                   this is zero a mode has been found
%   coreTerm      : the LHS of the equation
%   claddingTerm  : the RHS of the equation
% ------------------------------------------------------------------------
function [residual, coreTerm, claddingTerm] = calculate_LP_mismatch(coreRadius, nCore, nCladding, wavelength, order, u)

    % Wavenumber
    k = 2 * pi / wavelength;

    % Beta value
    beta = sqrt(k^2 * nCore^2 - (u/coreRadius)^2);

    % w value (cladding)
    w = coreRadius * sqrt(beta^2 - k^2 * nCladding^2);

    % core function
    coreTerm = besselj(order, u) / (u * besselj(order - 1, u));

    % cladding function
    claddingTerm = besselk(order, w) / (w * besselk(order - 1, w));

    % Calculate mismatch
    residual = coreTerm + claddingTerm;

end
