C****************************************************************************
C WAVETEST: Example Fortran program for WAVELET, using NINO3 SST dataset
C
C COMPILE:   f77 chisqr.f cfftpack.f wavelet.f wavetest.f
C
C See "http://paos.colorado.edu/research/wavelets/"
C
C  Copyright (C) 1998, Christopher Torrence and Gilbert P. Compo
C This software may be used, copied, or redistributed as long as it is not
C sold and this copyright notice is reproduced on each copy made.  This
C routine is provided as is without any express or implied warranties
C whatsoever.
C
C Modified: November 1999 by Arjan van Dijk to include IMPLICIT NONE and
C           to convert all routines to DOUBLE precision.
C****************************************************************************

      PROGRAM wavetest

      IMPLICIT none

      INTEGER n,subscale,jtot
      DOUBLE PRECISION dt,s0,dj

C these parameters depend on the particular time series
      PARAMETER (n=504,dt=0.25D0,s0=dt)
      PARAMETER (subscale=4)
      PARAMETER (dj=1.D0/subscale,jtot=11*subscale)
C Note: for accurate reconstruction and wavelet-derived variance
C     do not pad with zeroes, set s0=dt (for Paul set s0=dt/4), and use
C     a large "jtot" (even though the extra scales will be within
C     the cone of influence).
C     For plotting purposes, it is only necessary to use
C     s0=2dt (for Morlet) and "jtot" from Eqn(10) Torrence&Compo(1998).

      INTEGER mother,ibase2,npad
      DOUBLE PRECISION sst(n),recon_sst(n),param,pi
      DOUBLE PRECISION scale(jtot),period(jtot),coi(n)
      DOUBLE COMPLEX wave(n,jtot)

      INTEGER i,j,isigtest,javg1,javg2
      DOUBLE PRECISION lag1,siglvl,dof(jtot)
      DOUBLE PRECISION fft_theor(jtot),signif(jtot),ymean,variance
      DOUBLE PRECISION recon_mean,recon_vari
      DOUBLE PRECISION Cdelta,psi0
      DOUBLE PRECISION global_ws(jtot),global_signif(jtot)
      DOUBLE PRECISION savg_dof(jtot),savg_signif(jtot),sstENSO(n)

      pi = 4.D0*ATAN(1.D0)
      ibase2 = NINT(LOG(DBLE(n))/LOG(2.D0))+1
      npad = INT(2.D0**ibase2)
C      npad = n  ! this is for no padding with zeroes

C*************************************************** Wavelet transform

C** let the WAVELET subroutine choose the defaults for these:
      mother = 0
      param = 6.D0

C** read in the NINO3 SST data
      OPEN(UNIT=11,FILE='sst_nino3.dat',STATUS='old')
      READ(11,*) sst
      CLOSE(11)
      PRINT'(/,"sst(1)=",F6.2,"  sst(n) = ",F6.2,/)',sst(1),sst(n)

C** get the wavelet transform
      CALL WAVELET(n,sst,dt,mother,param,s0,dj,jtot,npad,
     &             wave,scale,period,coi)


C*************************************************** Significance testing

C** local significance test
      isigtest = 0
      lag1 = 0.72D0
      siglvl = 0.05D0
      CALL WAVE_SIGNIF (isigtest,n,sst,dt,mother,param,dj,jtot,
     &       scale,period,lag1,siglvl,dof,fft_theor,signif,
     &       ymean,variance,Cdelta,psi0)


C** global wavelet spectrum & significance test
      isigtest = 1
      lag1 = 0.72D0
      siglvl = 0.05D0
      DO 10 j=1,jtot
        DO 20 i=1,n
          global_ws(j) = global_ws(j) + ABS(wave(i,j))**2
20      CONTINUE
        global_ws(j) = global_ws(j)/n
        dof(j) = n - scale(j)
10    CONTINUE

      CALL WAVE_SIGNIF (isigtest,n,sst,dt,mother,param,dj,jtot,
     &       scale,period,lag1,siglvl,dof,fft_theor,global_signif,
     &       ymean,variance,Cdelta,psi0)


C** scale-average time series & significance test
      isigtest = 2
      lag1 = 0.72D0
      siglvl = 0.05D0
C    scale average between 2 and 7.9 years
      savg_dof(1) = 2.0D0
      savg_dof(2) = 7.9D0
C    find the "j"-values that correspond to savg_dof(1) & savg_dof(2)
      javg1 = 0
      javg2 = 0
      DO 30 j=1,jtot
        IF ((scale(j).GE.savg_dof(1)).AND.(javg1.EQ.0)) javg1 = j
        IF (scale(j).LE.savg_dof(2)) javg2 = j
30    CONTINUE
C   call wave_signif first, to get the value of "Cdelta"
      CALL WAVE_SIGNIF (isigtest,n,sst,dt,mother,param,dj,jtot,
     &     scale,period,lag1,siglvl,savg_dof,fft_theor,savg_signif,
     &     ymean,variance,Cdelta,psi0)
C   construct the scale-averaged time series [Eqn(24)]
      DO 50 i=1,n
        sstENSO(i) = 0.D0
        DO 60 j=javg1,javg2
          sstENSO(i) = sstENSO(i) + (ABS(wave(i,j))**2)/scale(j)
60      CONTINUE
        sstENSO(i) = dj*dt*sstENSO(i)/Cdelta
50    CONTINUE


C************************************************************* print results
      PRINT*,' n=',n
      PRINT*,' dt=',dt
      PRINT*,' mother=',mother
      PRINT*,' param=',param
      PRINT*,' s0=',s0
      PRINT*,' dj=',dj
      PRINT*,' jtot=',jtot
      PRINT*,' npad=',npad
      PRINT'(/,"Let w = wave(n/2,j)",/)'
      PRINT'(A4,7A10)',"j","Scale","Period","ABS(w)^2","phase(w)",
     &  "5%signif","Global","GWS5%sig"
      PRINT'(I4,7F10.3)',(j,scale(j),period(j),
     &   ABS(wave(n/2,j))**2,
     &   ATAN2(DIMAG(wave(n/2,j)),DBLE(wave(n/2,j)))*180.D0/pi,
     &   signif(j),global_ws(j),global_signif(j),j=1,jtot)
      PRINT'(/,A,F10.3)',
     &    ' Scale-average degrees of freedom = ',savg_dof(1)
      PRINT'(A,F10.3,/)',
     &    ' Scale-avg 5% significance level  = ',savg_signif(1)


C************************************************************ Reconstruction

C** construct the wavelet derived variance (Parseval's theorem)  [Eqn(14)]
C   Cdelta & psi0 are returned from WAVE_SIGNIF
      recon_vari = 0.D0
      DO 900 i=1,n
        DO 1000 j=1,jtot
          recon_vari = recon_vari + (ABS(wave(i,j))**2)/scale(j)
1000    CONTINUE
900   CONTINUE
      recon_vari = dj*dt*recon_vari/(Cdelta*n)
      PRINT'(A,F14.5)',' Reconstructed variance=',recon_vari
      PRINT'(A,F14.5)',' Original variance   =',variance
      PRINT'(A,F14.5,A,/)',' Ratio = ',recon_vari/variance,
     &     ' (this is low due to padding with zeroes)'

C** reconstruct the time series [Eqn(11)]
C   check mean and RMS difference of reconstructed time series
      recon_mean=0.D0
      recon_vari = 0.D0
      DO 1100 i=1,n
        recon_sst(i)=0.D0
        DO 1200 j=1,jtot
          recon_sst(i) = recon_sst(i)+(DBLE(wave(i,j)))/SQRT(scale(j))
1200    CONTINUE
        recon_sst(i) = dj*SQRT(dt)*recon_sst(i)/(Cdelta*psi0)
        recon_vari = recon_vari+(sst(i)-ymean-recon_sst(i))**2
        recon_mean = recon_mean + recon_sst(i)
1100  CONTINUE
      recon_mean = recon_mean/n
      recon_vari = SQRT(recon_vari/n)

      PRINT'(A,F14.6)',' Reconstructed mean=',recon_mean
      PRINT'(A,F14.6)',' Original mean   =',ymean
      PRINT'(A,F14.6,/)',' Root-mean-square difference of time series=',
     &      recon_vari

      END
