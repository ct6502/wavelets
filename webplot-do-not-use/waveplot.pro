PRO waveplot,filename,time_str,TEST=test

	errors = ['<FONT SIZE=6><CENTER>&#160;<BR>*** ERROR ***']
	error_status = 0
	CATCH,error_status
	IF (error_status NE 0) THEN BEGIN
		errors = [errors,'[IDL] '+!ERR_STRING+'</CENTER></FONT>'] + '<P>'
		error_status = 0
		GOTO,oops
	ENDIF

;********************************************************** INITIALIZE
	dataset = 'nino3sst'
	mother = 'Morlet'
	mother_old = ''
	param = 6d
	pad = 0b
	coi = 0b
	gws = 0b
	tstart = -1E37
	tend = 1E37
	dj = 1./4
	oct = 10
;	lag1 = 0.0
	siglvl = 10
	do_global = 1
	title = ' '
	xtitle = 'Time'
	ytitle = 'Period'
	t_units = 'years'
	eqn_string = ''
	units = ''
	data = FLTARR(4096)
	dataset1 = ''
	source = ''
	printflag = 0
	do_color = 1b   ; color
	encapsulated = 0b  ; PS or EPSF
	window_size = [600,400]
	test = KEYWORD_SET(test)
	IF (test) THEN BEGIN
		filename = 'test'
		time_str = ''
	ENDIF
	
;********************************************************** READ DATAFILES

;---------------------------------------------------------- Read user files
	file_in = filename
	flag = 0
	str = ''
	number = ''
	ON_IOERROR,skip
	OPENW,2,'/tmp/wave'
	FOR i=0,1 DO BEGIN
		link = file_in + '.dat'
		readfile = ([file_in+'.dat',filename+'.opt'])(i)
		OPENR,1,readfile
		READF,1,str
		WHILE (str NE 'data') DO BEGIN
			READF,1,number
			IF (STRMID(str,0,1) EQ '$') THEN BEGIN
				temp = number
				str_exec = STRMID(str,1,79) + '=temp'
			ENDIF ELSE str_exec = str + '=' + STRING(number)
			r = EXECUTE(str_exec)
			PRINTF,2,readfile+' '+str_exec
			IF (EOF(1)) THEN GOTO,skip
			READF,1,str
		ENDWHILE
		IF (i EQ 0) THEN READF,1,data
skip:
		IF (STRTRIM(eqn_string,2) NE '') THEN BEGIN ; user-defined eqn
			n = 2048
			pi = !PI
			t = t0 + dt*FINDGEN(n)
; get rid of multiple statements
			cutamp = STRPOS(eqn_string,'&')
			IF (cutamp GT 0) THEN eqn_string = STRMID(eqn_string,0,cutamp-1)
			r = EXECUTE('data = ' + eqn_string)
		ENDIF
		n = MAX(WHERE(data NE 0)) + 1
		IF (n EQ 0) THEN BEGIN ;---------------------- Read default data files
			IF (flag) THEN BEGIN
				errors=[errors,'[IDL] '+!ERR_STRING+'</CENTER></FONT>'] + '<P>'
				ON_IOERROR, NULL
				GOTO,oops
			ENDIF
			flag = 1  ; should only reach this point once
			file_in = 'data/' + dataset
			i = -1  ; loop again to read in default.dat and user.opt 
		ENDIF
		CLOSE,1
	ENDFOR
ON_IOERROR,NULL
CLOSE,2

;********************************************************** CALCULATE WAVELET
    IF (N_ELEMENTS(s0) LT 1) THEN BEGIN
    	s0 = 2.
    	coi = 1b
    	pad = 1b
    ENDIF
	dj = 0.1 > dj < 1
	IF (dt LT 1E-9) THEN errors = [errors,'dt must be greater than 0.0']
	IF (n LT 5) THEN errors = [errors,'Dataset must have at least 5 points']
	
    s0dt = s0*dt   ; put back in the delta-t scaling
	data = data(0:n-1)
	data0 = data
	data = data - TOTAL(data)/n  ;*** remove the mean
    variance = STDEV(data)^2
;	data = data/SQRT(variance)  ;*** normalize by standard deviation
	time = t0 + dt*FINDGEN(n)
	coi1 = KEYWORD_SET(coi) ;*** save a copy
	IF (mother_old NE mother) THEN dummy = TEMPORARY(param) ; use default
	CASE (mother) OF
		'Morlet': BEGIN
			IF (N_ELEMENTS(param) LT 1) THEN param = 6d ; default
			param = 2d > param < 20d
			param_str = ROUND_ANY(param,SIG=4)
			END
		'Paul': BEGIN
			IF (N_ELEMENTS(param) LT 1) THEN param = 4d ; default
			param = 1d > FIX(param) < 20d
			param_str = STRTRIM(FIX(param),2)
			END
		'DOG': BEGIN
			IF (N_ELEMENTS(param) LT 1) THEN param = 2d ; default
			param = 1d > FIX(param) < 20d
			param_str = STRTRIM(FIX(param),2)
			END
	ENDCASE
	base2 = FIX(ALOG(n/s0)/ALOG(2) + 0.4999)
	oct = 1 > oct < base2

	IF (N_ELEMENTS(errors) GT 1) THEN BEGIN ;************ E R R O R **********
		errors = [errors,'</CENTER></FONT>'] + '<P>'
		GOTO,oops
	ENDIF
	
	wave = WAVELET(data,dt, $   ;*** required inputs
		S0=s0dt,OCT=oct,DJ=dj, $   ;*** optional inputs
		PAD=pad,MOTHER=mother,PARAM=param, $   ;*** optional inputs
		SCALE=scale,PERIOD=period,COI=coi) ;*** optional outputs

;----------------------------------------------------- truncate data
	time0 = time[0]
	time1 = time[n-1]
	IF (tstart EQ -1E37) THEN tstart = time0
	IF (tend EQ 1E37) THEN tend = time1
    tstart = tstart < time(n-5)
    tend = time(4) > tend
    keep = WHERE((time GE tstart) AND (time LE tend),n)
    time = time(keep)
    data = data(keep)
    data0 = data0(keep)
    wave = wave(keep,*)
    coi = coi(keep)

; construct power, global wavelet
	power = ABS(wave)^2
	nt = N_ELEMENTS(time)
	na = N_ELEMENTS(period)
	global_ws = TOTAL(power,1)/nt

; significance level
	IF (N_ELEMENTS(lag1) LT 1) THEN BEGIN
		lag1 = A_CORRELATE(data,[1,2])
		lag1 = (lag1(0) + SQRT(lag1(1)))/2.
		lag1 = FIX(lag1*100. + 0.5)/100.
	ENDIF
	lag1 = (-1) > lag1 < 1
	IF (do_global) THEN BEGIN  ; use GWS as background
		fft_theor = global_ws
	ENDIF ELSE BEGIN  ; use AR1 as background
		freq = dt/period  ; normalized frequency
		fft_theor = (1-lag1^2)/(1-2*lag1*COS(freq*2*!PI)+lag1^2)  ; [Eqn(16)]
		fft_theor = variance*fft_theor
	ENDELSE
	
	siglvl = 0 > siglvl < 50
	signif = fft_theor*CHISQR_CVF((0>siglvl<50)/100.,2)/2.

; check if normalized by GWS
	IF (gws) THEN BEGIN
		power = TEMPORARY(power)/REBIN(TRANSPOSE(global_ws),nt,na)
		signif = signif/global_ws
		levels = [0.5,1,2,3]
	ENDIF ELSE BEGIN
		p1 = power(SORT(power))
		n1 = N_ELEMENTS(p1)
;		boxsize = [0.5,0.75,0.95]
		boxsize = [0.25,0.5,0.75,0.95]
		levels = [p1(boxsize*n1)]
	ENDELSE

	labels = levels  ; make a copy for MAKE_KEY

;----------------------------------------------------- Global significance
	IF (siglvl GT 0) THEN BEGIN
		dof = n;-scale  ;*** the -scale corrects for the edge effects
		globalsignif = WAVE_SIGNIF(data,dt,scale,1, $
			DOF=dof,LAG1=lag1,SIGLVL=1.-(siglvl<50)/100., $
			MOTHER=mother,PARAM=param)
	ENDIF



;********************************************************** START PLOTTING
	postscript = (printflag GT 0)
	CASE (printflag) OF
	0: BEGIN
		SET_PLOT,'Z'
		DEVICE,Z_BUFFERING=0,SET_RESOLUTION=window_size
		END
	1: BEGIN
		SET_PLOT,'ps'
		DEVICE,FILE=filename + time_str + '.ps', $
			BITS=8,/COLOR, $
			ISOLATIN=0, $
			/LANDSCAPE, $
			/INCHES,XSIZE=9,XOFFSET=1,YSIZE=6.5,YOFFSET=10
		DEVICE,FONT_INDEX=5,/HELV,/BOLD
		END
	2: BEGIN
		SET_PLOT,'ps'
		DEVICE,FILE=filename + time_str + '.ps', $
			BITS=8,/COLOR, $
			ISOLATIN=0, $
			/ENCAPSULATED, $
			/INCHES,XSIZE=8,YSIZE=6
		DEVICE,FONT_INDEX=5,/HELV,/BOLD
		END
	ENDCASE
	
	LOADCT,42*do_color,/SILENT
;	LOADCT,49,/SILENT
	!P.BACKGROUND = 255b
	!P.COLOR = 0
	!X.STYLE = 1
	!Y.STYLE = 1
	!P.FONT = -1 + postscript
	!P.CHARSIZE = 0.8 + 0.2*postscript
	!P.CHARTHICK = 1 - postscript
	!X.TICKLEN = -0.03
	!Y.TICKLEN = -0.01

	grid = ([1,5])(postscript)
	yrange = [MAX(period),MIN(period)]
	ytickv = 2.^(INDGEN(oct) + FIX(ALOG(MIN(period))/ALOG(2.)))
	xtitle1 = xtitle
	ytitle1 = ytitle
	ytitle2 = 'Scale'
	IF (N_ELEMENTS(dataset1) LT 1) THEN dataset1 = dataset
	dataset2 = dataset1
	IF (STRLEN(t_units) GT 0) THEN BEGIN  ;*** if "units", then tack onto labels
		xtitle1 = xtitle1 + ' (' + t_units + ')'
		ytitle1 = ytitle1 + ' (' + t_units + ')'
		ytitle2 = ytitle2 + ' (' + t_units + ')'
	ENDIF
	IF (STRLEN(units) GT 0) THEN BEGIN  ;*** if "units", then tack onto labels
		dataset2 = dataset1 + ' (' + units + ')'
	ENDIF

;---------------------------------------------------------- Time Series Plot
	!P.MULTI = 0
	pcolor = 16*do_color
	pos1 = [0.1,0.7,0.7,0.945]
	PLOT,time,data0, $
		XSTYLE=9,YSTYLE=9, $
		XRANGE=[tstart,tend], $
		XTITLE=xtitle1,YTITLE=dataset2,/NODATA,YMINOR=2,POSITION=pos1
	OPLOT,time,data0,THICK=1+postscript,COLOR=pcolor
	dummy = PLOT_LABEL(0,0.75,/AFTER,NUMBER=1,TITLE=title)
	
;---------------------------------------------------------- Wavelet Contour Plot
	colors = BYTSCL(FINDGEN(N_ELEMENTS(levels)),TOP=200) + 48
	colors = [208,176,128,80]
	IF (do_color) THEN colors = [64,144,180,240]
	pos2 = [pos1(0),0.21,pos1(2),0.55]
	
	CONTOUR,power,time,period,/NOERASE, $
		XRANGE=[tstart,tend], $
		POSITION=pos2, $
		YRANGE=yrange,YSTYLE=1, $
		/YTYPE,/FILL,LEVELS=levels,C_COLORS=colors, $
		YTICKS=N_ELEMENTS(ytickv)-1,YTICKV=ytickv,YMINOR=1, $
		XTITLE=xtitle1,YTITLE=ytitle1
;		XGRIDS=grid,YGRIDS=grid,XTICKLEN=0.5,YTICKLEN=1
;	FIT_IN_WINDOW,postscript,ROTATE(SQRT(power),7),TOP=254,BOTTOM=48

;---------------------------------------------------------- signif level
	IF (siglvl GT 0) THEN BEGIN
		signif = power/REBIN(TRANSPOSE(signif),n,na)
		CONTOUR,signif,time,period,/OVERPLOT,LEVEL=1.0,THICK=2+3*postscript
	ENDIF
	
;---------------------------------------------------------- COI
	IF (coi1) THEN BEGIN
		color = ([100,4])(do_color)
		maxx = MAX(coi) > MAX(period)
		POLYFILL,[time(0),time,MAX(time)],[maxx,coi,maxx], $
			ORIEN=+45,SPACING=0.5,COLOR=color,NOCLIP=0,THICK=1
		POLYFILL,[time(0),time,MAX(time)],[maxx,coi,maxx], $
			ORIEN=-45,SPACING=0.5,COLOR=color,NOCLIP=0,THICK=1
		PLOTS,time,coi,COLOR=color,NOCLIP=0,THICK=1
		xx=0.33 & yy=-0.07
   	ENDIF

;---------------------------------------------------------- WPS title
	units_sqr = ''
	IF (units NE '') THEN units_sqr = ' (' + units + ')!U2!N'
;	variance_str = ROUND_ANY(variance,SIG=3) + units_sqr
	IF (gws) THEN title1=' (normalized by global)' ELSE title1=''
	dummy = PLOT_LABEL(0,1,/AFTER,NUMBER=2, $
		TITLE='Wavelet Power Spectrum'+title1)
	
;---------------------------------------------------------- Make color key
	pos3 = [pos2[0],pos2(1),pos2[2]-(pos2[2]-pos2[0])/3.,pos2(3)]
	PLOT,[0,1],/NODATA,XSTYLE=5,YSTYLE=5,POSITION=pos3,/NOERASE
	nlab = N_ELEMENTS(labels)
	labels = [0.,labels]
;	labels = ROUND_ANY(labels,SIG=2)
;	labels(nlab-1) = ' ' + labels(nlab-1) + '!U !Ns!U2!N'
	labels = STRING(DOUBLE(labels),FORMAT='(G9.2)')
	FOR i=0,N_ELEMENTS(labels)-1 DO BEGIN
		expo = STRPOS(labels[i],'E')
		IF (expo GE 1) THEN BEGIN
			num = STRMID(labels[i],0,expo)
			expnt = FIX(STRMID(labels[i],expo+1,255))
			labels[i] = num + 'x10!U' + STRING(expnt) + '!N'
		ENDIF
	ENDFOR
	labels = STRCOMPRESS(labels,/REMOVE_ALL)
	colors = [!P.BACKGROUND,colors]
;	IF (N_ELEMENTS(boxsize) GT 0) THEN extra = {BOXSIZE:[boxsize,1]-[0,boxsize]}
	IF (gws) THEN units2=' (relative to global)' ELSE units2=units_sqr
	MAKE_KEY,LABELS=labels,COLORS=colors,_EXTRA=extra, $
		TITLE='Power'+units2

;---------------------------------------------------------- Title
	XYOUTS,0.84,0.9, $
		'!5INTERACTIVE WAVELET PLOT!X'+$
		'!C!CC. Torrence & G.P. Compo'+$
		'!C!Chttp://paos.colorado.edu/!Cresearch/wavelets/', $
		/NORMAL,ALIGN=0.5,CHARSIZE=0.8


;********************************************************** PLOTTING GLOBAL
	pos4 = [0.74,pos2(1),0.95,pos2(3)]
	blank = REPLICATE(' ',29)
	
	PLOT,global_ws,period,/NOERASE, $
		XSTYLE=4,YSTYLE=9, $
		/XTYPE, $
		XTICK_GET=xtick_get, $
		YRANGE=yrange,/YTYPE, $
		YTICKS=N_ELEMENTS(ytickv)-1,YTICKV=ytickv,YTICKNAME=blank, $
		YTICKLEN=-0.02, $
		POSITION=pos4, $
		THICK=2+2*postscript
;	xtickname = ROUND_ANY(DOUBLE(xtick_get),SIG=2)
	xtickname = STRING(DOUBLE(xtick_get),FORMAT='(G9.2)')
	FOR i=0,N_ELEMENTS(xtickname)-1 DO BEGIN
		expo = STRPOS(xtickname[i],'E')
		IF (expo GE 1) THEN BEGIN
			num = STRMID(xtickname[i],0,expo)
			expnt = FIX(STRMID(xtickname[i],expo+1,255))
			xtickname[i] = '10!U' + STRING(expnt) + '!N'
		ENDIF
	ENDFOR
	xtickname = STRCOMPRESS(xtickname,/REMOVE_ALL)
	AXIS,XAXIS=0,XTICKNAME=xtickname, $
		XTITLE='Variance'+units_sqr
	IF ((siglvl GT 0) AND (N_ELEMENTS(globalsignif) GT 1)) THEN BEGIN
		OPLOT,globalsignif,period, $
			COLOR=208*do_color, $
			LINES=2-postscript, $
			NOCLIP=0
	ENDIF
	dummy = PLOT_LABEL(0,1,/AFTER,NUMBER=3,TITLE='Global Wavelet')

;********************************************************** PLOTTING MOTHER
	I = COMPLEX(0,1)
	t0 = (DINDGEN(101)-50)/10.
	pname = 'm'
	CASE (mother) OF
		'Morlet': BEGIN
			psi = EXP(I*param*t0)*EXP(-t0^2/2)
			pname = 'k!D0!N'
			END
		'Paul': BEGIN
			psi = (1-I*t0)^(-param+1)
			END
		'DOG': BEGIN
			sw = 8*2*!DPI*DINDGEN(101)/100d
			psi = FLOAT(FFT(-I^param*(sw^param)*EXP(-sw^2/2),1,/DOUBLE))
			psi = [psi[50:*],psi[0:49]]
			END
	ENDCASE
	keep = WHERE(ABS(psi) GE 0.001*MAX(ABS(psi)))
	t0 = t0[keep]
	psi = psi[keep]
	psiReal = FLOAT(psi)
	psiImag = IMAGINARY(COMPLEX(psi))
	mn = MIN([psiReal,psiImag],MAX=mx)

	PLOT,t0,psiImag,/NOERASE,LINES=1, $
		XSTYLE=5,YSTYLE=5, $
		/NOCLIP, $
		YRANGE=[mn,mx], $
		POSITION=[pos2[2]-(pos2[2]-pos2[0])/4.,0.02,pos2[2],pos2[1]-0.06]
;	OPLOT,t0,psiImag,LINES=1,COLOR=!P.BACKGROUND
	OPLOT,t0,psiReal
	title1 = ''
	IF (MAX(ABS(psiImag)) GT 0) THEN title1='!C !Dsolid=Real, dashed=Imaginary'
	XYOUTS,MAX(t0),0,' ' + mother + ' ' + $
		pname + '=' + param_str + title1,/NOCLIP
	

oops:
	CATCH,error_status
	IF (error_status NE 0) THEN EXIT

;********************************************************** OUTPUT OPTIONS .html
	OPENR,1,'waveopt.html'
	outf = 2
	OPENW,outf,filename + '.html'
	WHILE NOT EOF(1) DO BEGIN
		READF,1,str
		CASE (STRMID(str,0,1)) OF
			'#': r = EXECUTE(STRMID(str,1,255))
			ELSE: PRINTF,outf,str
		ENDCASE
	ENDWHILE
	CLOSE,/ALL
	
;********************************************************** OUTPUT .gif FILE
	IF (N_ELEMENTS(errors) LT 2) THEN BEGIN
		IF (postscript) THEN BEGIN
			DEVICE,/CLOSE
		ENDIF ELSE BEGIN
			image1 = TVRD()
			WRITE_GIF,filename + time_str + '.gif',image1
		ENDELSE
	ENDIF
END
