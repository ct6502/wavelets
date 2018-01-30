;************************************************************** WAVE_SIGNIF
;+
; NAME:   WAVE_SIGNIF
;
; PURPOSE:   Compute the significance levels for a wavelet transform.
;       
;
; CALLING SEQUENCE:
;
;      result = WAVE_SIGNIF(y,dt,scale,sigtest)
;
;
; INPUTS:
;
;    Y = the time series, or, the VARIANCE of the time series.
;        (If this is a single number, it is assumed to be the variance...)
;
;    DT = amount of time between each Y value, i.e. the sampling time.
;
;    SCALE = the vector of scale indices, from previous call to WAVELET.
;
;    SIGTEST = 0, 1, or 2.    If omitted, then assume 0.
;
;          If 0 (the default), then just do a regular chi-square test,
;       		 i.e. Eqn (18) from Torrence & Compo.
;          If 1, then do a "time-average" test, i.e. Eqn (23).
;       		 In this case, DOF should be set to NA, the number
;       		 of local wavelet spectra that were averaged together.
;       		 For the Global Wavelet Spectrum, this would be NA=N,
;       		 where N is the number of points in your time series.
;          If 2, then do a "scale-average" test, i.e. Eqns (25)-(28).
;       		 In this case, DOF should be set to a
;       		 two-element vector [S1,S2], which gives the scale
;       		 range that was averaged together.
;       		 e.g. if one scale-averaged scales between 2 and 8,
;                     then DOF=[2,8].
;
;
; OUTPUTS:
;
;    result = significance levels as a function of SCALE,
;             or if /CONFIDENCE, then confidence intervals
;
;
; OPTIONAL KEYWORD INPUTS:
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
;    LAG1 = LAG 1 Autocorrelation, used for SIGNIF levels. Default is 0.0
;
;    SIGLVL = significance level to use. Default is 0.95
;
;    DOF = degrees-of-freedom for signif test.
;          IF SIGTEST=0, then (automatically) DOF = 2 (or 1 for MOTHER='DOG')
;          IF SIGTEST=1, then DOF = NA, the number of times averaged together.
;          IF SIGTEST=2, then DOF = [S1,S2], the range of scales averaged.
;
;   	 Note: IF SIGTEST=1, then DOF can be a vector (same length as SCALEs),
;   		   in which case NA is assumed to vary with SCALE.
;   		   This allows one to average different numbers of times
;   		   together at different scales, or to take into account
;   		   things like the Cone of Influence.
;   		   See discussion following Eqn (23) in Torrence & Compo.
;
;    GWS = global wavelet spectrum. If input then this is used
;          as the theoretical background spectrum,
;          rather than white or red noise.
;
;    CONFIDENCE = if set, then return a Confidence INTERVAL.
;                 For SIGTEST=0,2 this will be two numbers, the lower & upper.
;                 For SIGTEST=1, this will return an array (J+1)x2,
;                 where J+1 is the number of scales.
;
;
; OPTIONAL KEYWORD OUTPUTS:
;
;    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
;           to the SCALEs.
;
;    FFT_THEOR = output theoretical red-noise spectrum as fn of PERIOD.
;
;
;----------------------------------------------------------------------------
;
; EXAMPLE:
;
;    IDL> wave = WAVELET(y,dt,PERIOD=period,SCALE=scale)
;    IDL> signif = WAVE_SIGNIF(y,dt,scale)
;    IDL> signif = REBIN(TRANSPOSE(signif),ntime,nscale)
;    IDL> CONTOUR,ABS(wave)^2/signif,time,period, $
;           LEVEL=1.0,C_ANNOT='95%'
;
;
;----------------------------------------------------------------------------
; Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo,
; University of Colorado, Program in Atmospheric and Oceanic Sciences.
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties whatsoever.
;
; Notice: Please acknowledge the use of the above software in any publications:
;    ``Wavelet software was provided by C. Torrence and G. Compo,
;      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
;
;----------------------------------------------------------------------------
;-
;************************************************************* WAVE_SIGNIF
FUNCTION wave_signif,y,dt,scale,sigtest, $   ;*** required inputs
	MOTHER=mother,PARAM=param, $   ;*** optional inputs
	LAG1=lag1,SIGLVL=siglvl,DOF=dof, $   ;*** optional inputs
	GWS=gws,CONFIDENCE=confidence, $   ;*** optional inputs
	FFT_THEOR=fft_theor,PERIOD=period, $   ;*** optional outputs
	SAVG=Savg,SMID=Smid,CDELTA=CDelta,PSI0=psi0   ;*** optional outputs
	
	ON_ERROR,2
	IF (N_ELEMENTS(y) LT 1) THEN MESSAGE,'Time series Y must be input'
	IF (N_ELEMENTS(dt) LT 1) THEN MESSAGE,'DT must be input'
	IF (N_ELEMENTS(scale) LT 1) THEN MESSAGE,'Scales must be input'
	IF (N_PARAMS() LT 4) THEN sigtest = 0   ; the default
	IF (N_ELEMENTS(y) EQ 1) THEN variance=y ELSE variance=(MOMENT(y))(1)
	
;....check keywords & optional inputs
	IF (N_ELEMENTS(mother) LT 1) THEN mother = 'MORLET'
	IF (N_ELEMENTS(param) LT 1) THEN param = -1
	IF (N_ELEMENTS(siglvl) LT 1) THEN siglvl = 0.95
	IF (N_ELEMENTS(lag1) LT 1) THEN lag1 = 0.0
	confidence = KEYWORD_SET(confidence)
	
	lag1 = lag1(0)
	
	J = N_ELEMENTS(scale) - 1
	s0 = MIN(scale)
	dj = ALOG(scale(1)/scale(0))/ALOG(2)
	
	CASE (STRUPCASE(mother)) OF
	'MORLET': BEGIN
			IF (param EQ -1) THEN k0=6d ELSE k0=param
			fourier_factor = (4*!PI)/(k0 + SQRT(2+k0^2)) ; [Sec.3h]
			empir = [2.,-1,-1,-1]
			IF (k0 EQ 6) THEN empir(1:3)=[0.776,2.32,0.60]
		END
	'PAUL': BEGIN ;****************** PAUL
			IF (param EQ -1) THEN m=4d ELSE m=param
			fourier_factor = 4*!PI/(2*m+1)
			empir = [2.,-1,-1,-1]
			IF (m EQ 4) THEN empir(1:3)=[1.132,1.17,1.5]
		END
	'DOG': BEGIN ;******************* DOG
			IF (param EQ -1) THEN m=2 ELSE m=param
			fourier_factor = 2*!PI*SQRT(2./(2*m+1))
			empir = [1.,-1,-1,-1]
			IF (m EQ 2) THEN empir(1:3) = [3.541,1.43,1.4]
			IF (m EQ 6) THEN empir(1:3) = [1.966,1.37,0.97]
		END
	ENDCASE
	
	period = scale*fourier_factor
	dofmin = empir(0) ; Degrees of freedom with no smoothing
	Cdelta = empir(1) ; reconstruction factor
	gamma = empir(2)  ; time-decorrelation factor
	dj0 = empir(3)    ; scale-decorrelation factor

;....significance levels [Sec.4]
	freq = dt/period  ; normalized frequency
	fft_theor = (1-lag1^2)/(1-2*lag1*COS(freq*2*!PI)+lag1^2)  ; [Eqn(16)]
	fft_theor = variance*fft_theor  ; include time-series variance
	IF (N_ELEMENTS(gws) EQ (J+1)) THEN fft_theor = gws
	signif = fft_theor

	CASE (sigtest) OF

	0: BEGIN   ; no smoothing, DOF=dofmin
		dof = dofmin
		signif = fft_theor*CHISQR_CVF(1. - siglvl,dof)/dof   ; [Eqn(18)]
		IF confidence THEN BEGIN
			sig = (1. - siglvl)/2.
			chisqr = dof/[CHISQR_CVF(sig,dof),CHISQR_CVF(1.-sig,dof)]
			signif = fft_theor # chisqr
		ENDIF
		END

	1: BEGIN   ; time-averaged, DOFs depend upon scale [Sec.5a]
		IF (N_ELEMENTS(dof) LT 1) THEN dof = dofmin
		IF (gamma EQ -1) THEN MESSAGE, $
			'Gamma (decorrelation factor) not defined for '+mother+ $
			' with param='+STRTRIM(param,2)
		IF (N_ELEMENTS(dof) EQ 1) THEN dof = FLTARR(J+1) + dof
		dof = dof > 1
		dof = dofmin*SQRT( 1 + (dof*dt/gamma/scale)^2 ) ; [Eqn(23)]
		dof = dof > dofmin   ; minimum DOF is dofmin
		IF (NOT confidence) THEN BEGIN
			FOR a1=0,J DO BEGIN
				chisqr = CHISQR_CVF(1. - siglvl,dof(a1))/dof(a1)
				signif(a1) = fft_theor(a1)*chisqr
			ENDFOR
		ENDIF ELSE BEGIN
			signif = FLTARR(J+1,2)
			sig = (1. - siglvl)/2.
			FOR a1=0,J DO BEGIN
				chisqr = dof(a1)/ $
					[CHISQR_CVF(sig,dof(a1)),CHISQR_CVF(1.-sig,dof(a1))]
				signif(a1,*) = fft_theor(a1)*chisqr
			ENDFOR
		ENDELSE
		END

	2: BEGIN  ; scale-averaged, DOFs depend upon scale range [Sec.5b]
		IF (N_ELEMENTS(dof) NE 2) THEN $
			MESSAGE,'DOF must be set to [S1,S2], the range of scale-averages'
		IF (Cdelta EQ -1) THEN MESSAGE, $
			'Cdelta & dj0 not defined for '+mother+ $
			' with param='+STRTRIM(param,2)
		s1 = dof(0)
		s2 = dof(1)
		avg = WHERE((scale GE s1) AND (scale LE s2),navg)
		IF (navg LT 1) THEN MESSAGE,'No valid scales between ' + $
			STRTRIM(s1,2) + ' and ' + STRTRIM(s2,2)
		s1 = MIN(scale(avg))
		s2 = MAX(scale(avg))
		Savg = 1./TOTAL(1./scale(avg))       ; [Eqn(25)]
		Smid = EXP((ALOG(s1)+ALOG(s2))/2.)   ; power-of-two midpoint
		dof = (dofmin*navg*Savg/Smid)*SQRT(1 + (navg*dj/dj0)^2)  ; [Eqn(28)]
		fft_theor = Savg*TOTAL(fft_theor(avg)/scale(avg))  ; [Eqn(27)]
		chisqr = CHISQR_CVF(1. - siglvl,dof)/dof
		IF confidence THEN BEGIN
			sig = (1. - siglvl)/2.
			chisqr = dof/[CHISQR_CVF(sig,dof),CHISQR_CVF(1.-sig,dof)]
		ENDIF
		signif = (dj*dt/Cdelta/Savg)*fft_theor*chisqr  ; [Eqn(26)]
		END

	ENDCASE

	RETURN,signif

END
