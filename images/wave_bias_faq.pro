
	dt = 1
	n = 128
	t = FINDGEN(n)*dt
	t1 = 8
	t2 = 32
	x = SIN(2*!PI*t/t1) + SIN(2*!PI*t/t2)
;	x = RANDOMN(s,n)
	
	f = SPECTRUM(x,dt,PERIOD=fperiod)
	w = WAVELET(x,dt,PERIOD=wperiod,DJ=0.125)
	gws = TOTAL(ABS(w)^2,1)/n^2

	WINDOW,XSIZE=240,YSIZE=350
	LOADCT,39
	!P.COLOR = 192
	!P.MULTI = [0,1,3]
	!X.STYLE = 9
	!Y.STYLE = 9
	!X.MARGIN = [8,2]
	!Y.MARGIN = [4,3]
	!Y.MINOR = 1
	!X.TICKLEN = -0.04
	!P.CHARSIZE = 1.75
	xticks = FIX(ALOG(n)/ALOG(2))-1
	xtickv = 2L^(LINDGEN(xticks+1)+1)
	xrange = [MAX(xtickv),MIN(xtickv)]
	
	PLOT,t,x, $
		XRANGE=[0,n], $
		XMINOR=1,XTICKS=xticks-2, $
		XTITLE='Time (days)', $
		TITLE='Time series with two sine waves'
	OPLOT,t,x,COLOR=254

	PLOT,fperiod,f,/XTYPE, $
		XRANGE=xrange,YRANGE=[0,0.5],YTICKS=2, $
		XTICKS=xticks,XTICKV=xtickv, $
		XTITLE='Period (days)',YTITLE='Power', $
		TITLE='Fourier power spectrum'
	OPLOT,fperiod,f,PSYM=-4,/NOCLIP,COLOR=254,SYMSIZE=0.75

	PLOT,wperiod,gws,/XTYPE, $
		XRANGE=xrange,YRANGE=[0,0.25],YTICKS=2, $
		XTICKS=xticks,XTICKV=xtickv, $
		XTITLE='Period (days)',YTITLE='Power', $
		TITLE='Global wavelet spectrum'
	OPLOT,wperiod,gws,PSYM=-4,/NOCLIP,COLOR=254,SYMSIZE=0.75
	
	TVLCT,r,g,b,/GET
	WRITE_GIF,'wave_bias_faq.gif',TVRD(),r,g,b
END
