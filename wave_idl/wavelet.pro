;******************************************************************* WAVELET
;+
; NAME:   WAVELET
;
; PURPOSE:   Compute the WAVELET transform of a 1D time series.
;       
;
; CALLING SEQUENCE:
;
;      wave = WAVELET(Y,DT)
;
;
; INPUTS:
;
;    Y = the time series of length N.
;
;    DT = amount of time between each Y value, i.e. the sampling time.
;
;
; OUTPUTS:
;
;    WAVE is the WAVELET transform of Y. This is a complex array
;    of dimensions (N,J+1). FLOAT(WAVE) gives the WAVELET amplitude,
;    ATAN(IMAGINARY(WAVE),FLOAT(WAVE)) gives the WAVELET phase.
;    The WAVELET power spectrum is ABS(WAVE)^2.
;
;
; OPTIONAL KEYWORD INPUTS:
;
;    S0 = the smallest scale of the wavelet.  Default is 2*DT.
;
;    DJ = the spacing between discrete scales. Default is 0.125.
;         A smaller # will give better scale resolution, but be slower to plot.
;
;    J = the # of scales minus one. Scales range from S0 up to S0*2^(J*DJ),
;        to give a total of (J+1) scales. Default is J = (LOG2(N DT/S0))/DJ.
;
;    MOTHER = A string giving the mother wavelet to use.
;            Currently, 'Morlet','Paul','DOG' (derivative of Gaussian)
;            are available. Default is 'Morlet'.
;
;    PARAM = optional mother wavelet parameter.
;            For 'Morlet' this is k0 (wavenumber), default is 6.
;            For 'Paul' this is m (order), default is 4.
;            For 'DOG' this is m (m-th derivative), default is 2.
;
;    PAD = if set, then pad the time series with enough zeroes to get
;         N up to the next higher power of 2. This prevents wraparound
;         from the end of the time series to the beginning, and also
;         speeds up the FFT's used to do the wavelet transform.
;         This will not eliminate all edge effects (see COI below).
;
;    LAG1 = LAG 1 Autocorrelation, used for SIGNIF levels. Default is 0.0
;
;    SIGLVL = significance level to use. Default is 0.95
;
;    VERBOSE = if set, then print out info for each analyzed scale.
;
;    RECON = if set, then reconstruct the time series, and store in Y.
;            Note that this will destroy the original time series,
;            so be sure to input a dummy copy of Y.
;
;    FFT_THEOR = theoretical background spectrum as a function of
;                Fourier frequency. This will be smoothed by the
;                wavelet function and returned as a function of PERIOD.
;
;
; OPTIONAL KEYWORD OUTPUTS:
;
;    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
;           to the SCALEs.
;
;    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J
;            where J+1 is the total # of scales.
;
;    COI = if specified, then return the Cone-of-Influence, which is a vector
;        of N points that contains the maximum period of useful information
;        at that particular time.
;        Periods greater than this are subject to edge effects.
;        This can be used to plot COI lines on a contour plot by doing:
;            IDL>  CONTOUR,wavelet,time,period
;            IDL>  PLOTS,time,coi,NOCLIP=0
;
;    YPAD = returns the padded time series that was actually used in the
;         wavelet transform.
;
;    DAUGHTER = if initially set to 1, then return the daughter wavelets.
;         This is a complex array of the same size as WAVELET. At each scale
;         the daughter wavelet is located in the center of the array.
;
;    SIGNIF = output significance levels as a function of PERIOD
;
;    FFT_THEOR = output theoretical background spectrum (smoothed by the
;                wavelet function), as a function of PERIOD.
;
;
; [ Defunct INPUTS:
; [   OCT = the # of octaves to analyze over.           ]
; [         Largest scale will be S0*2^OCT.             ]
; [         Default is (LOG2(N) - 1).                   ]
; [   VOICE = # of voices in each octave. Default is 8. ]
; [          Higher # gives better scale resolution,    ]
; [          but is slower to plot.                     ]
; ]
;
;----------------------------------------------------------------------------
;
; EXAMPLE:
;
;    IDL> ntime = 256
;    IDL> y = RANDOMN(s,ntime)       ;*** create a random time series
;    IDL> dt = 0.25
;    IDL> time = FINDGEN(ntime)*dt   ;*** create the time index
;    IDL> 
;    IDL> wave = WAVELET(y,dt,PERIOD=period,COI=coi,/PAD,SIGNIF=signif)
;    IDL> nscale = N_ELEMENTS(period)
;    IDL> LOADCT,39
;    IDL> CONTOUR,ABS(wave)^2,time,period, $
;       XSTYLE=1,XTITLE='Time',YTITLE='Period',TITLE='Noise Wavelet', $
;       YRANGE=[MAX(period),MIN(period)], $   ;*** Large-->Small period
;       /YTYPE, $                             ;*** make y-axis logarithmic
;       NLEVELS=25,/FILL
;    IDL>
;    IDL> signif = REBIN(TRANSPOSE(signif),ntime,nscale)
;    IDL> CONTOUR,ABS(wave)^2/signif,time,period, $
;           /OVERPLOT,LEVEL=1.0,C_ANNOT='95%'
;    IDL> PLOTS,time,coi,NOCLIP=0   ;*** anything "below" this line is dubious
;
;
;----------------------------------------------------------------------------
; Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
;
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.
; This routine is provided as is without any express or implied warranties
; whatsoever.
;
; Notice: Please acknowledge the use of the above software in any publications:
;    ``Wavelet software was provided by C. Torrence and G. Compo,
;      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
;
; Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
;            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
;
; Please send a copy of such publications to either C. Torrence or G. Compo:
;  Dr. Christopher Torrence               Dr. Gilbert P. Compo
;  Research Systems, Inc.                 Climate Diagnostics Center
;  4990 Pearl East Circle                 325 Broadway R/CDC1
;  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
;  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
;----------------------------------------------------------------------------
;-

FUNCTION morlet, $ ;*********************************************** MORLET
	k0,scale,k,period,coi,dofmin,Cdelta,psi0

	IF (k0 EQ -1) THEN k0 = 6d
	n = N_ELEMENTS(k)
	expnt = -(scale*k - k0)^2/2d*(k GT 0.)
	dt = 2*!PI/(n*k(1))
	norm = SQRT(2*!PI*scale/dt)*(!PI^(-0.25))   ; total energy=N   [Eqn(7)]
	morlet = norm*EXP(expnt > (-100d))
	morlet = morlet*(expnt GT -100)  ; avoid underflow errors
	morlet = morlet*(k GT 0.)  ; Heaviside step function (Morlet is complex)
	fourier_factor = (4*!PI)/(k0 + SQRT(2+k0^2)) ; Scale-->Fourier [Sec.3h]
	period = scale*fourier_factor
	coi = fourier_factor/SQRT(2)   ; Cone-of-influence [Sec.3g]
	dofmin = 2   ; Degrees of freedom with no smoothing
	Cdelta = -1
	IF (k0 EQ 6) THEN Cdelta = 0.776 ; reconstruction factor
	psi0 = !PI^(-0.25)
;	PRINT,scale,n,SQRT(TOTAL(ABS(morlet)^2,/DOUBLE))
	RETURN,morlet
END

FUNCTION paul, $ ;************************************************** PAUL
	m,scale,k,period,coi,dofmin,Cdelta,psi0

	IF (m EQ -1) THEN m = 4d
	n = N_ELEMENTS(k)
	expnt = -(scale*k)*(k GT 0.)
	dt = 2d*!PI/(n*k(1))
	norm = SQRT(2*!PI*scale/dt)*(2^m/SQRT(m*FACTORIAL(2*m-1)))
	paul = norm*((scale*k)^m)*EXP(expnt > (-100d))*(expnt GT -100)
	paul = paul*(k GT 0.)
	fourier_factor = 4*!PI/(2*m+1)
	period = scale*fourier_factor
	coi = fourier_factor*SQRT(2)
	dofmin = 2   ; Degrees of freedom with no smoothing
	Cdelta = -1
	IF (m EQ 4) THEN Cdelta = 1.132 ; reconstruction factor
	psi0 = 2.^m*FACTORIAL(m)/SQRT(!PI*FACTORIAL(2*m))
;	PRINT,scale,n,norm,SQRT(TOTAL(paul^2,/DOUBLE))*SQRT(n)
	RETURN,paul
END

FUNCTION dog, $ ;*************************************************** DOG
	m,scale,k,period,coi,dofmin,Cdelta,psi0

	IF (m EQ -1) THEN m = 2
	n = N_ELEMENTS(k)
	expnt = -(scale*k)^2/2d
	dt = 2d*!PI/(n*k(1))
	norm = SQRT(2*!PI*scale/dt)*SQRT(1d/GAMMA(m+0.5))
	I = DCOMPLEX(0,1)
	gauss = -norm*(I^m)*(scale*k)^m*EXP(expnt > (-100d))*(expnt GT -100)
	fourier_factor = 2*!PI*SQRT(2./(2*m+1))
	period = scale*fourier_factor
	coi = fourier_factor/SQRT(2)
	dofmin = 1   ; Degrees of freedom with no smoothing
	Cdelta = -1
	psi0 = -1
	IF (m EQ 2) THEN BEGIN
		Cdelta = 3.541 ; reconstruction factor
		psi0 = 0.867325
	ENDIF
	IF (m EQ 6) THEN BEGIN
		Cdelta = 1.966 ; reconstruction factor
		psi0 = 0.88406
	ENDIF
;	PRINT,scale,n,norm,SQRT(TOTAL(ABS(gauss)^2,/DOUBLE))*SQRT(n)
	RETURN,gauss
END


;****************************************************************** WAVELET
FUNCTION wavelet,y1,dt, $   ;*** required inputs
	S0=s0,DJ=dj,J=j, $   ;*** optional inputs
	PAD=pad,MOTHER=mother,PARAM=param, $
	VERBOSE=verbose,NO_WAVE=no_wave,RECON=recon, $
	LAG1=lag1,SIGLVL=siglvl,DOF=dof,GLOBAL=global, $   ;*** optional inputs
	SCALE=scale,PERIOD=period,YPAD=ypad, $  ;*** optional outputs
	DAUGHTER=daughter,COI=coi, $
	SIGNIF=signif,FFT_THEOR=fft_theor, $
	OCT=oct,VOICE=voice   ;*** defunct inputs
	
	ON_ERROR,2
	r = CHECK_MATH(0,1)
	n = N_ELEMENTS(y1)
	n1 = n
	base2 = FIX(ALOG(n)/ALOG(2) + 0.4999)   ; power of 2 nearest to N

;....check keywords & optional inputs
	IF (N_ELEMENTS(s0) LT 1) THEN s0 = 2.0*dt
	IF (N_ELEMENTS(voice) EQ 1) THEN dj = 1./voice
	IF (N_ELEMENTS(dj) LT 1) THEN dj = 1./8
	IF (N_ELEMENTS(oct) EQ 1) THEN J = FLOAT(oct)/dj
	IF (N_ELEMENTS(J) LT 1) THEN J=FIX((ALOG(FLOAT(n)*dt/s0)/ALOG(2))/dj)  ;[Eqn(10)]
	IF (N_ELEMENTS(mother) LT 1) THEN mother = 'MORLET'
	IF (N_ELEMENTS(param) LT 1) THEN param = -1
	IF (N_ELEMENTS(siglvl) LT 1) THEN siglvl = 0.95
	IF (N_ELEMENTS(lag1) LT 1) THEN lag1 = 0.0
	lag1 = lag1(0)
	verbose = KEYWORD_SET(verbose)
	do_daughter = KEYWORD_SET(daughter)
	do_wave = NOT KEYWORD_SET(no_wave)
	recon = KEYWORD_SET(recon)
	IF KEYWORD_SET(global) THEN MESSAGE, $
		'Please use WAVE_SIGNIF for global significance tests'
		
;....construct time series to analyze, pad if necessary
	ypad = y1 - TOTAL(y1)/n    ; remove mean
	IF KEYWORD_SET(pad) THEN BEGIN   ; pad with extra zeroes, up to power of 2
		ypad = [ypad,FLTARR(2L^(base2 + 1) - n)]
		n = N_ELEMENTS(ypad)
	ENDIF

;....construct SCALE array & empty PERIOD & WAVE arrays
	na = J + 1                  ; # of scales
	scale = DINDGEN(na)*dj      ; array of j-values
	scale = 2d0^(scale)*s0      ; array of scales  2^j   [Eqn(9)]
	period = FLTARR(na,/NOZERO) ; empty period array (filled in below)
	wave = COMPLEXARR(n,na,/NOZERO)  ; empty wavelet array
	IF (do_daughter) THEN daughter = wave   ; empty daughter array

;....construct wavenumber array used in transform [Eqn(5)]
	k = (DINDGEN(n/2) + 1)*(2*!PI)/(DOUBLE(n)*dt)
	k = [0d,k,-REVERSE(k(0:(n-1)/2 - 1))]

;....compute FFT of the (padded) time series
	yfft = FFT(ypad,-1,/DOUBLE)  ; [Eqn(3)]

	IF (verbose) THEN BEGIN  ;verbose
		PRINT
		PRINT,mother
		PRINT,'#points=',n1,'   s0=',s0,'   dj=',dj,'   J=',FIX(J)
		IF (n1 NE n) THEN PRINT,'(padded with ',n-n1,' zeroes)'
		PRINT,['j','scale','period','variance','mathflag'], $
			FORMAT='(/,A3,3A11,A10)'
	ENDIF  ;verbose
	IF (N_ELEMENTS(fft_theor) EQ n) THEN fft_theor_k = fft_theor ELSE $
		fft_theor_k = (1-lag1^2)/(1-2*lag1*COS(k*dt)+lag1^2)  ; [Eqn(16)]
	fft_theor = FLTARR(na)
	
;....loop thru each SCALE
	FOR a1=0,na-1 DO BEGIN  ;scale
		psi_fft=CALL_FUNCTION(mother, $
			param,scale(a1),k,period1,coi,dofmin,Cdelta,psi0)
		IF (do_wave) THEN $
			wave(*,a1) = FFT(yfft*psi_fft,1,/DOUBLE)  ;wavelet transform[Eqn(4)]
		period(a1) = period1   ; save period
		fft_theor(a1) = TOTAL((ABS(psi_fft)^2)*fft_theor_k)/n
		IF (do_daughter) THEN $
			daughter(*,a1) = FFT(psi_fft,1,/DOUBLE)   ; save daughter
		IF (verbose) THEN PRINT,a1,scale(a1),period(a1), $
				TOTAL(ABS(wave(*,a1))^2),CHECK_MATH(0), $
				FORMAT='(I3,3F11.3,I6)'
	ENDFOR  ;scale

	coi = coi*[FINDGEN((n1+1)/2),REVERSE(FINDGEN(n1/2))]*dt   ; COI [Sec.3g]
	
	IF (do_daughter) THEN $   ; shift so DAUGHTERs are in middle of array
		daughter = [daughter(n-n1/2:*,*),daughter(0:n1/2-1,*)]

;....significance levels [Sec.4]
	sdev = (MOMENT(y1))(1)
	fft_theor = sdev*fft_theor  ; include time-series variance
	dof = dofmin
	signif = fft_theor*CHISQR_CVF(1. - siglvl,dof)/dof   ; [Eqn(18)]

	IF (recon) THEN BEGIN  ; Reconstruction [Eqn(11)]
		IF (Cdelta EQ -1) THEN BEGIN
			y1 = -1
			MESSAGE,/INFO, $
				'Cdelta undefined, cannot reconstruct with this wavelet'
		ENDIF ELSE BEGIN
			y1=dj*SQRT(dt)/(Cdelta*psi0)*(FLOAT(wave) # (1./SQRT(scale)))
			y1 = y1[0:n1-1]
		ENDELSE
	ENDIF
	
	RETURN,wave(0:n1-1,*)    ; get rid of padding before returning

END
