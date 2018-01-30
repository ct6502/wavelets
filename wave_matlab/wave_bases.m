%WAVE_BASES  1D Wavelet functions Morlet, Paul, or DOG
%
%  [DAUGHTER,FOURIER_FACTOR,COI,DOFMIN] = ...
%      wave_bases(MOTHER,K,SCALE,PARAM);
%
%   Computes the wavelet function as a function of Fourier frequency,
%   used for the wavelet transform in Fourier space.
%   (This program is called automatically by WAVELET)
%
% INPUTS:
%
%    MOTHER = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
%    K = a vector, the Fourier frequencies at which to calculate the wavelet
%    SCALE = a number, the wavelet scale
%    PARAM = the nondimensional parameter for the wavelet function
%
% OUTPUTS:
%
%    DAUGHTER = a vector, the wavelet function
%    FOURIER_FACTOR = the ratio of Fourier period to scale
%    COI = a number, the cone-of-influence size at the scale
%    DOFMIN = a number, degrees of freedom for each point in the wavelet power
%             (either 2 for Morlet and Paul, or 1 for the DOG)
%
%----------------------------------------------------------------------------
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%----------------------------------------------------------------------------
function [daughter,fourier_factor,coi,dofmin] = ...
	wave_bases(mother,k,scale,param);

mother = upper(mother);
n = length(k);

if (strcmp(mother,'MORLET'))  %-----------------------------------  Morlet
	if (param == -1), param = 6.;, end
	k0 = param;
	expnt = -(scale.*k - k0).^2/2.*(k > 0.);
	norm = sqrt(scale*k(2))*(pi^(-0.25))*sqrt(n);    % total energy=N   [Eqn(7)]
	daughter = norm*exp(expnt);
	daughter = daughter.*(k > 0.);     % Heaviside step function
	fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
	coi = fourier_factor/sqrt(2);                  % Cone-of-influence [Sec.3g]
	dofmin = 2;                                    % Degrees of freedom
elseif (strcmp(mother,'PAUL'))  %--------------------------------  Paul
	if (param == -1), param = 4.;, end
	m = param;
	expnt = -(scale.*k).*(k > 0.);
	norm = sqrt(scale*k(2))*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n);
	daughter = norm*((scale.*k).^m).*exp(expnt);
	daughter = daughter.*(k > 0.);     % Heaviside step function
	fourier_factor = 4*pi/(2*m+1);
	coi = fourier_factor*sqrt(2);
	dofmin = 2;
elseif (strcmp(mother,'DOG'))  %--------------------------------  DOG
	if (param == -1), param = 2.;, end
	m = param;
	expnt = -(scale.*k).^2 ./ 2.0;
	norm = sqrt(scale*k(2)/gamma(m+0.5))*sqrt(n);
	daughter = -norm*(i^m)*((scale.*k).^m).*exp(expnt);
	fourier_factor = 2*pi*sqrt(2./(2*m+1));
	coi = fourier_factor/sqrt(2);
	dofmin = 1;
else
	error('Mother must be one of MORLET,PAUL,DOG')
end

return

% end of code

