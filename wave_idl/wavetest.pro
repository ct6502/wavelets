;************************************************************** WAVETEST
;+
; NAME:   WAVETEST
;
; PURPOSE:   Example IDL program for WAVELET, using ENSO SST dataset
;
; EXECUTION:
; 
;            IDL> .run wavetest
;
;
; See "http://paos.colorado.edu/research/wavelets/"
; Written January 1998 by C. Torrence
;
;-
;**************************************************************

  compile_opt idl2

  if (0) then begin  ; original dataset
    n = 504
    sst = FLTARR(n)
    OPENR,1,'sst_nino3.dat'   ; input SST time series
    READF,1,sst
    CLOSE,1
    dt = 0.25
    dtName = 'seasonal'
    years = [1871, 1996]
    contourLevels = [0.5,1,2,4]
    dataset = 'NINO3'
  endif else begin
    years = intarr(2)
    openr, 1, 'nino34.long.data'
    readf, 1, years
    data = fltarr(13, years[1] - years[0] + 1)
    readf, 1, data
    close, 1
    data = data[1:12,*]
    fullYears = total(data gt 0, 1) eq 12
    data = data[*, where(fullYears, nyears)]
    years[1] = years[0] + nyears - 1
    annualCycle = total(data, 2)/nyears
    sst = data - rebin(annualCycle, 12, nyears)
    sst = sst[*]
    dt = 1d/12
    dtName = 'monthly'
    n = n_elements(sst)
    contourLevels = [1, 2, 4, 8]
    dataset = 'NINO3.4'
  endelse

;------------------------------------------------------ Computation

  sst = sst - MEAN(sst)

  time = FINDGEN(n)*dt + years[0]  ; construct time array
  xrange = [years[0]/10*10,(years[1] + 9)/10*10]  ; plotting range
  pad = 1
  s0 = dt    ; this says start at a scale of 3 months
  dj = 0.25  ; this will do 4 sub-octaves per octave
  j1 = 9./dj  ; this says do 9 powers-of-two with dj sub-octaves each
  mother = 'Morlet'
  recon_sst = sst   ; save an extra copy, so we don't erase original sst

; estimate lag-1 autocorrelation, for red-noise significance tests
  lag1 = (A_CORRELATE(sst,1) + SQRT(A_CORRELATE(sst,2)))/2.
  
; Wavelet transform:
  wave = WAVELET(recon_sst,dt,PERIOD=period,SCALE=scale,S0=s0, $
    PAD=pad,COI=coi,DJ=dj,J=j1,MOTHER=mother,/RECON)
  power = (ABS(wave))^2  ; compute wavelet power spectrum
  global_ws = TOTAL(power,1)/n   ; global wavelet spectrum (GWS)
  J = N_ELEMENTS(scale) - 1

; Significance levels:
  signif = WAVE_SIGNIF(sst,dt,scale,0, $
    LAG1=lag1,SIGLVL=0.95,MOTHER=mother)
  signif = REBIN(TRANSPOSE(signif),n,J+1)  ; expand signif --> (J+1)x(N) array
  signif = power/signif   ; where ratio > 1, power is significant

; GWS significance levels:
  dof = n - scale   ; the -scale corrects for padding at edges
  global_signif = WAVE_SIGNIF(sst,dt,scale,1, $
    LAG1=lag1,DOF=dof,MOTHER=mother,CDELTA=Cdelta,PSI0=psi0)

; check total variance (Parseval's theorem) [Eqn(14)]
  scale_avg = REBIN(TRANSPOSE(scale),n,J+1)  ; expand scale-->(J+1)x(N) array
  power_norm = power/scale_avg
  variance = (MOMENT(sst))[1]
  recon_variance = dj*dt/(Cdelta*n)*TOTAL(power_norm)  ; [Eqn(14)]
  
  IF (N_ELEMENTS(recon_sst) GT 1) THEN BEGIN
    recon_variance = (MOMENT(recon_sst))[1]
; RMS of Reconstruction [Eqn(11)]
    rms_error = SQRT(TOTAL((sst - recon_sst)^2)/n)
    PRINT
    PRINT,'        ******** RECONSTRUCTION ********'
    PRINT,'original variance =',variance,' degC^2'
    PRINT,'reconstructed var =',FLOAT(recon_variance),' degC^2'
    PRINT,'Ratio = ',recon_variance/variance
    PRINT,'root-mean-square error of reconstructed sst = ',rms_error,' degC'
    PRINT
    IF (mother EQ 'DOG') THEN BEGIN
      PRINT,'Note: for better reconstruction with the DOG, you need'
      PRINT,'      to use a very small s0.'
    ENDIF
    PRINT
  ENDIF
  
; Scale-average between El Nino periods of 2--8 years
  avg = WHERE((scale GE 2) AND (scale LT 8))
  scale_avg = dj*dt/Cdelta*TOTAL(power_norm[*,avg],2)  ; [Eqn(24)]
  scaleavg_signif = WAVE_SIGNIF(sst,dt,scale,2, $
    LAG1=lag1,SIGLVL=0.95,DOF=[2,7.9],MOTHER=mother)


;------------------------------------------------------ Plotting

  w = WINDOW(DIMENSIONS=[900,1000], LOCATION=[100, 0])

;--- Plot time series
  pos1 = [0.1,0.75,0.7,0.95]
  p = PLOT(time, sst, '2', XRANGE=xrange, /CURRENT, $
    XTITLE='Time (year)',YTITLE=dataset + ' SST (!Uo!NC)', $
    TITLE='a) ' + dataset + ' Sea Surface Temperature (' + dtName + ')', $
    LAYOUT=[1,3,1], MARGIN=[0.1,0.15,0.3,0.15], YMINOR=1)
  IF (N_ELEMENTS(recon_sst) GT 1) THEN begin
;    p = PLOT(time,recon_sst,COLOR='red', /OVERPLOT)
  endif
  t = TEXT(0.85,0.85,/NORM,ALIGN=0.5,FONT_SIZE=12, $
    'Wavelet Analysis'+$
    '!C!CC. Torrence & G.P. Compo'+$
    '!C!Chttp://paos.colorado.edu/!Cresearch/wavelets/')

;--- Contour plot wavelet power spectrum
  yrange = [64,0.5]   ; years
  colors = [80,120,160,200]
  colors = ['bisque','orange','orange_red','dark_red']
  c_colors = (OrderedHash(!color))[colors.ToUpper()]
  c_colors = (c_colors.Values()).ToArray(/TRANSPOSE)
  period2 = FIX(ALOG(period)/ALOG(2))   ; integer powers of 2 in period
  ytickv = 2.^(period2[UNIQ(period2)])  ; unique powers of 2
  ytickv = ytickv[WHERE(ytickv ge MIN(yrange) and ytickv le MAX(yrange))]
  pos2 = [pos1[0],0.35,pos1[2],0.65]

  cLabel = STRING(contourLevels, format='(10(g0,:,","))')

  c = CONTOUR(power,time,period, /CURRENT, $
    LAYOUT=[1,3,2], MARGIN=[0.1,0.15,0.3,0.15], $
    AXIS_STYLE=1, XRANGE=xrange,YRANGE=yrange, /YLOG, $
    YTICKV=ytickv, YMINOR=0, $
    C_VALUE=contourLevels,C_COLOR=c_colors,/FILL, $
    XTITLE='Time (year)',YTITLE='Period (years)', $
    TITLE='b) Wavelet Power Spectrum (contours at ' + cLabel + '!Uo!NC!U2!N)')

; significance contour
  c = CONTOUR(signif,time,period,/OVERPLOT,C_VALUE=1,C_THICK=2, $
    C_LABEL_SHOW=0, COLOR='black')

; cone-of-influence, anything "below" is dubious
  x = [time[0],time,time[-1]]
  y = [MAX(period),coi > 0.00001,MAX(period)]
  for angle=-45,45,90 do begin
    p = POLYGON(x, y, /DATA, TARGET=c, FILL_COLOR='gray', $
      /FILL_BACKGROUND, PATTERN_ORIENT=angle, PATTERN_SPACING=15)
  endfor
  p = PLOT(time, coi, /OVERPLOT, COLOR='gray')

;--- Plot global wavelet spectrum
  pos3 = [0.74,pos2[1],0.95,pos2[3]]
  blank = REPLICATE(' ',29)
  p = PLOT(global_ws,period, /CURRENT, $
    LAYOUT=[2,3,4], MARGIN=[0.49,0.15,0.1,0.15], $
    AXIS_STYLE=1, XRANGE=[0,MAX(global_ws)], YRANGE=yrange, /YLOG, $
    YTICKV=ytickv, YTICKNAME=blank, YMINOR=0, YTICKLEN=0.08, $
    XTITLE='Power (!Uo!NC!U2!N)',TITLE='c) Global Wavelet Spectrum')
  p = PLOT(global_signif,period,'--', /OVERPLOT)

;--- Plot 2--8 yr scale-average time series
  pos4 = [pos1[0],0.05,pos1[2],0.25]
  p = PLOT(time, scale_avg, /CURRENT, $
    LAYOUT=[1,3,3], MARGIN=[0.1,0.15,0.3,0.15], $
    XRANGE=xrange,YRANGE=[0,MAX(scale_avg)*1.25],THICK=2, $
    XTITLE='Time (year)',YTITLE='Avg variance (!Uo!NC!U2!N)', $
    TITLE='d) 2-8 yr Scale-average Time Series')
  p = PLOT(xrange, scaleavg_signif+[0,0], '--', /OVERPLOT)

END
