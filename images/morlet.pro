

FUNCTION psi,k,t,WAVE_PART=wave_part,AMP_PART=amp_part

	wave_part = COMPLEX(COS(k*t),SIN(k*t))   ; exp(ikt) = cos(kt) + i*sin(kt)
	amp_part = EXP( -ABS(t)^2/2.)   ; exp(-|t|^2/2)
	RETURN,wave_part*amp_part   ; morlet

END


	k0 = 6.0
	time = FINDGEN(401)/20. - 10.
	
	a = 2
	b = 0
	
	morlet = (1./SQRT(a))*PSI(k0,(time - b)/a,WAVE=wave,AMP=amp)
	
	WINDOW,0,COLOR=256,XSIZE=512,YSIZE=128
	LOADCT,42
	!P.BACKGROUND = 255
	!P.COLOR = 0
	!P.FONT = -1
	!P.CHARSIZE = 1.5
	!X.STYLE = 5
	!Y.STYLE = 5
	!P.MULTI = [0,2,1]
	!X.MARGIN = [2,2]
	!Y.MARGIN = [0,0]
	!Y.RANGE = [MIN(FLOAT(wave))*0.95, MAX(amp)*1.05]
	dy = (!Y.CRANGE(1) - !Y.CRANGE(0))/5.
	
	PLOT,time,FLOAT(morlet), $
		THICK=2,COLOR=50
	XYOUTS,!X.CRANGE(0),!Y.CRANGE(1) - dy,'(a)',CHARSIZE=1.2,ALIGN=0.5

	PLOT,time,FLOAT(wave)/SQRT(2), $
		COLOR=128,THICK=2
	OPLOT,time,amp,COLOR=208,THICK=2
	OPLOT,time,FLOAT(morlet),THICK=2,LINES=2,COLOR=50
	XYOUTS,!X.CRANGE(0),!Y.CRANGE(1) - dy,'(b)',CHARSIZE=1.2,ALIGN=0.5



END
