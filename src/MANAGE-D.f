      SUBROUTINE MANAGE
C
C     Changes since version published in RIWC 44 (1994) 158-160 :
C     9/10/95  IS replaced by ILAST (typo in published version)
C     9/10/95  Added check that abundance estimates are in year order
C     7/1/98   Changed to double precision
C
C **** This program sets a catch limit for a specified species and area. ****
C **** If catch capping and/or catch cascading is required, they must    ****
C **** be applied to this catch limit afterwards.                        ****
C
C ****
C **** INPUT FILES:  CLC.PAR   containing the input parameters
C                    CLC.DAT   containing the catch and abundance estimates
C      OUTPUT FILE:  CLC.OUT   
C
C ****
C **** GLOBAL INPUT PARAMETERS:
C
      COMMON /MANPAR/ PPROB, PYMAX, PNYSTP, PKSTEP, PDSTEP,PNBSTP,PBMIN,
     1                PBMAX, PSCALE,PHASET, PHASEP, PCYCLE,PLEVEL,PSLOPE
      DOUBLE PRECISION PPROB,PYMAX, PNYSTP, PKSTEP, PDSTEP,PNBSTP,PBMIN, 
     1                PBMAX, PSCALE,PHASET, PHASEP, PCYCLE,PLEVEL,PSLOPE
C     PPROB   Probability level
C     PYMAX   Maximum value of the productivity parameter (Y) for integration 
C     PNYSTP  Number of steps in integration over Y  
C     PKSTEP  Maximum relative step size in integration over K
C     PDSTEP  Target step size for integration over depletion
C     PNBSTP  Number of steps for integration over the bias parameter  
C     PBMIN   Minimum multiplicative bias eg 0.5 means 50% downward bias  
C     PBMAX   Maximum multiplicative bias eg 1.5 means 50% upward bias    
C     PSCALE  Raw deviance scaling factor S = 1/PSCALE**2
C     PHASET  Number of years without surveys before phaseout invoked 
C     PHASEP  Phaseout: annual reduction proportion
C     PCYCLE  Maximum number of years before a new CLC 
C             = number of years a CLC is valid 
C     PLEVEL  Internal protection level     
C     PSLOPE  `Slope' of catch control law    
C
C
C     MANDAT: All years are scaled so 0 = a year prior to the 1st input data
      COMMON   /MANDAT/  ISTART,IYEAR,NS,NZ,RKLO,RKHI
      INTEGER            ISTART,IYEAR,NS,NZ
      DOUBLE PRECISION   RKLO, RKHI
C     ISTART Year of first input data (catch or abundance data)
C     IYEAR  Year for which to set first catch limit
C     NS     Number of non-zero estimates
C     NZ     Number of zero estimates
C     RKLO   Lower bound used in integration over K
C     RKHI   Upper bound used in integration over K
C
C
      INTEGER MAXYR,MAXEST,MSIZE
      PARAMETER (MAXYR=200, MAXEST=100)
      PARAMETER (MSIZE=(MAXEST*(MAXEST+1)/2))
      DOUBLE PRECISION CATCH(0:MAXYR), SIGHT(0:MAXEST), FMATRX(0:MSIZE),
     +     ZMULT(0:MAXEST),POP(0:MAXYR), G(0:MAXYR), RAWCL, CL, CLIMIT
      INTEGER ISYR(0:MAXEST), IZYR(0:MAXEST), POUT
C
C     CATCH  Catch array, indexed by year
C     SIGHT  Abundance estimate array, indexed by estimate number
C     FMATRX Information matrix (H) of the log sightings estimates (excluding 
C            zero estimates, stored as a lower triangle
C            FMATRX((i*(i+1))/2 + j = element H(i,j)   NB i & j start at 0
C     ISYR(N) Year of Nth abundance estimate SIGHT(N),  where N=0,NS-1
C     ZMULT(N) Poisson multiplier for Nth zero estimate
C     IZYR(N) Year of Nth zero estimate,  where N=0,NZ-1
C     POP(I)  Modelled population size in year I (set in STKSIM)
C     G(I)    (set & used in DEVIAN)
C     RAWCL   Nominal catch limit i.e. output from the Catch Limit Algorithm
C     CL      The catch limit for the area considered
C     POUT    Parameter determining whether phaseout may be applied if 
C             necessary.  Test for phaseout if POUT=1.  (Phaseout is
C             not applied to Medium or Large area nominal catch limits) 
C
C ****
C **** LOCAL VARIABLES
C
      CHARACTER STOCK*30, FORMT*50
      DOUBLE PRECISION TOTALC,C
      INTEGER IY, INITYR, I, IYRCL, ILAST, N, IN, IOUT,N1
C
C     INITYR Year in which first premanagement catch taken 
C     ILAST  Year of the most recent abundance estimate
C
C ****
C *** Read in the data required for application of the CLA

      DATA IN/2/, IOUT/3/ 
      OPEN (IN, FILE='CLC.DAT')
      OPEN (IOUT, FILE='CLC.OUT')
C
C     Read in the parameter values
      CALL RDPARS
C
C     Read management unit name (species and area)
      READ (IN,'(A30/)') STOCK
      WRITE (IOUT,'(A30/)') STOCK
C
C     Read the year of the first catch, the first year for which catch 
C     limits are to be set & the phaseout option
      READ (IN,'(T30,I4)') INITYR
      WRITE (IOUT,'(A,T30,I4)') 'Year of first input data',INITYR
      READ (IN,'(T30,I4)') IYRCL
      WRITE (IOUT,'(A,T30,I4)') 'Year of first catch limit',IYRCL
      READ (IN,'(T30,I4)') POUT
      IF (POUT.EQ.1) THEN
        WRITE (IOUT,'(A,T30,A4)') 'Apply phaseout if necessary','Yes'
      ELSE
        WRITE (IOUT,'(A,T30,A4)') 'Apply phaseout','No'
      ENDIF
C
C     Re-scale IYRCL such that 0 is the year prior to the first input data
      ISTART = 0      
      IYEAR = IYRCL - INITYR
      IF (IYEAR.LE.0 .OR. IYEAR.GT.MAXYR) STOP 'INVALID YEAR'
C
C     Initialise the catch array
      DO 10 I=0,MAXYR
        CATCH(I) = 0.D0
   10 CONTINUE
C
C     Read in the catch data, scaling each year to the initial year
C     A negative year indicates the end of the catch data
C     Check that there is at least one catch.
      READ (IN,'(// A)') FORMT
      WRITE (IOUT,'(/A)') 'Historic catches:'
      TOTALC = 0.D0
      DO 20 I = 0,MAXYR
        READ (IN,FORMT) IY,C
        IF (IY.LT.0) GO TO 25
        WRITE (IOUT,FORMT) IY,C
        IY = IY-INITYR
        IF (IY.LT.0 .OR. IY.GE.IYEAR) STOP 
     +                              ' ERROR: Catch year out of range'
        CATCH(IY) = C
        TOTALC = TOTALC + C
   20 CONTINUE
C
   25 IF (TOTALC.LT.1.D0) STOP ' ERROR: No historic catch input'
C
C     Read in the non-zero sightings estimates and information matrix
      READ (IN,'(//T30,I4 / A)') NS,FORMT
      WRITE (IOUT,'(/A)') 'Abundance estimates:'
      IF (NS.GT.MAXEST) STOP ' ERROR: Abundance year out of range'
      DO 30 N=0,NS-1
        N1 = (N*(N+1))/2
c        READ (IN,FORMT) ISYR(N),SIGHT(N), (FMATRX(N1+J),J=0,N)
c        WRITE (IOUT,FORMT) ISYR(N),SIGHT(N), (FMATRX(N1+J),J=0,N)
        READ (IN,FORMT) ISYR(N),SIGHT(N), FMATRX(N1+N)
        WRITE (IOUT,FORMT) ISYR(N),SIGHT(N), FMATRX(N1+N)
        ISYR(N) = ISYR(N)-INITYR
        IF (ISYR(N).LT.0 .OR. ISYR(N).GE.IYEAR) STOP 
     +                              ' ERROR: Sight year out of range'
        IF (SIGHT(N).LE.0.D0) STOP ' ERROR: Estimate not positive'
        IF (N.GT.0 .AND. ISYR(N).LT.ISYR(N-1)) STOP 
     +                              ' ERROR: Sight year out of order'
   30 CONTINUE
C
C     Read in the Poisson multipliers for any zero sightings estimates
      READ (IN,'(/T30,I3/A)') NZ,FORMT
      IF (NZ.GT.0) WRITE (IOUT,'(/A)') 'Zero Abundance estimates:'
      IF (NZ.GT.MAXEST) STOP ' ERROR: Zero estimate array too small'
      DO 40 N=0,NZ-1
        N1 = (N*(N+1))/2
        READ (IN,FORMT) IZYR(N),ZMULT(N)
        WRITE (IOUT,FORMT) IZYR(N),ZMULT(N)
        IZYR(N) = IZYR(N)-INITYR
        IF (IZYR(N).LT.0 .OR. IZYR(N).GE.IYEAR) STOP 
     +                                       ' Sight year out of range'
        IF (ZMULT(N).LE.0.D0) STOP ' ERROR: Multiplier not positive'
        IF (N.GT.0 .AND. IZYR(N).LT.IZYR(N-1)) STOP 
     +                              ' ERROR: Sight year out of order'
   40 CONTINUE
      WRITE (IOUT,'()')
C
C     Bound the range for the integration over K
      RKHI = 1.D7
      RKLO = 0.D0
C
C ****
C **** Run the CLA to obtain the nominal catch limit
C ****
      RAWCL = CLIMIT(CATCH,SIGHT,FMATRX,ISYR,IZYR,ZMULT,POP,G)
      OPEN(UNIT=99,FILE='RES.OUT')
      WRITE(99,*) RAWCL
      CLOSE(99)
C
C
C     Set the catch limits for PCYCLE years. If the catch limit may be 
C     subject to phaseout, call PHOUT to apply the phaseout rule.

C     First set ILAST = year of the most recent abundance estimate
      IF (NS.GT.0) ILAST = ISYR(NS-1)
      IF (NZ.GT.0) ILAST = MAX(ILAST, IZYR(NZ-1))
C
      DO 100 IY = IYEAR,IYEAR+(INT(PCYCLE+0.0001))-1
        IF (POUT.EQ.1) THEN 
          CALL PHOUT (CL,RAWCL,ILAST,IY)
        ELSE
          CL = RAWCL
        ENDIF
        WRITE (IOUT,'(A6,I5,A17,I6)')'Year:',IY+INITYR,'Catch limit:',
     +         NINT(CL)
  100 CONTINUE
C
      CLOSE(IN)
      CLOSE(IOUT)
      RETURN
      END
C     INCLUDE 'CLC-D.FOR'