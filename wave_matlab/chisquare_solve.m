function PDIFF = chisquare_solve(XGUESS,P,V);
%CHISQUARE_SOLVE  Internal function used by CHISQUARE_INV
%
%   PDIFF=chisquare_solve(XGUESS,P,V)  Given XGUESS, a percentile P,
%   and degrees-of-freedom V, return the difference between
%   calculated percentile and P.

% Uses GAMMAINC
%
% Written January 1998 by C. Torrence

% extra factor of V is necessary because X is Normalized
	PGUESS = gammainc(V*XGUESS/2,V/2);  % incomplete Gamma function
	
	PDIFF = abs(PGUESS - P);            % error in calculated P
	
	TOL = 1E-4;
	if (PGUESS >= 1-TOL)  % if P is very close to 1 (i.e. a bad guess)
		PDIFF = XGUESS;   % then just assign some big number like XGUESS
	end

	return

% end of code

