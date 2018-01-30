function X = chisquare_inv(P,V);
%CHISQUARE_INV  Inverse of chi-square cumulative distribution function (cdf).
%
%   X = chisquare_inv(P,V) returns the inverse of chi-square cdf with V
%   degrees of freedom at fraction P.
%   This means that P*100 percent of the distribution lies between 0 and X.
%
%   To check, the answer should satisfy:   P==gammainc(X/2,V/2)

% Uses FMIN and CHISQUARE_SOLVE
%
% Written January 1998 by C. Torrence

	if (nargin < 2), error('Must input both P and V');, end
	if ((1-P) < 1E-4), error('P must be < 0.9999');, end
	
	if ((P==0.95) & (V==2)) % this is a no-brainer
		X = 5.9915;
		return
	end
	
	MINN = 0.01;         % hopefully this is small enough
	MAXX = 1;            % actually starts at 10 (see while loop below)
	X = 1;
	TOLERANCE = 1E-4;    % this should be accurate enough
    vers = version;
    vers = str2num(vers(1));

	while ((X+TOLERANCE) >= MAXX)  % should only need to loop thru once
		MAXX = MAXX*10.;
% this calculates value for X, NORMALIZED by V
% Note: We need two different versions, depending upon the version of Matlab.
        if (vers >= 6)
            X = fminbnd('chisquare_solve',MINN,MAXX,optimset('TolX',TOLERANCE),P,V);
        else
    		X = fmin('chisquare_solve',MINN,MAXX,[0,TOLERANCE],P,V);
        end
		MINN = MAXX;
	end
	
	X = X*V;  % put back in the goofy V factor

	return

% end of code

