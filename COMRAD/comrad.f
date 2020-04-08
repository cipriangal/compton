C
      PROGRAM COMRAD
C
C **********************************************************************
C * Author: M. Swartz                              Date: November 1997 *
C * Reference:                                                         *
C *   SLAC-PUB-7701; hep-ph/9711447; Phys. Rev. D58, 014010 (1998)     *
C *                                                                    *
C * Monte Carlo for e-gamma=>e-gamma(gamma) and e-gamma=>e-e+e-        *
C * to full order-alpha**3 including initial state polarization        *
C * (circular polarization for gamma, general spin direction for e-).  *
C * Uses Stuart and Gongora, Z.Phys. C42,617 (1989) for e+2gam final   *
C * state; Tsai, DeRaad, and Milton, Phys Rev D6, 1411 (1972) for      *
C * virtual and soft-photon corrections; and home-grown calculation    *
C * (based upon the Gongora-Stuart spinor technique) for the 3e final  *
C * state.  Includes explicit unpol xsection calcs as diagnostics.     *
C * Consists of three weighted-generators: COMTN2 which the does the   *
C * e-gamma final state, COMEGG which does the e-2gamma final state,   *
C * and COMEEE which does the 3e final state.  The user interface to   *
C * all 3 generators is via the routine WGTHST (see comments there).   *
C *                                                                    *
C * Requires two external CERNLIB routines:                            *
C *   RANMAR random number generator from MATHLIB (V113)               *
C *   DDILOG double precision dilogarithm from MATHLIB (C332)          *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  The block /CONTRL/ contains the adjustable parameters needed for the
C  calculation and some useful physical constants: EB - e- beam energy in
C  GeV, EPHOT - photon energy in GeV, XME - electron mass (in GeV),
C  XMG - the fictitious photon mass (GeV), KGMIN - the minimum resolvable
C  photon energy (GeV), ALPHA - the fine structure constant, PI - 3.1416...
C  ROOT2 - SQRT(2), BARN - (hbar*c)**2 in units of mb-GeV**(-2),
C  SPIN(3) - the initial state electron spin (in the me rest frame),
C  LDIAG - logical flag to turn on diagnostic T-DR-M unpol xsection, 
C  LBF - logical flag to turn on diagnostic Brown+Feynman unpol xsection, 
C  NTRY - the number of event trials to generate with COMTN2,COMEGG,COMEEE
C         (COMEEE generates NTRY/20)
C
      COMMON /CONTRL/ EB,EPHOT,XME,XMG,KGMIN,ALPHA,PI,ROOT2,BARN,
     >                SPIN(3),LDIAG,LBF,NTRY,WGT(4)
      LOGICAL LDIAG,LBF
      REAL*8 KGMIN
C
C  Define a bunch of physical constants
C
      XME=0.511D-3
      ALPHA=1.d0/137.0359895D0
      PI=ACOS(-1.d0)
      ROOT2=SQRT(2.d0)
C
C  XMG is a small mass for the photon to regulate the infrared divergence
C
      XMG=0.00001D0*XME
C
C  KGMIN is the minimum detectable energy of the photons in the electron
C  rest frame (not the cm frame)
C
      KGMIN=0.002D0*XME
C
C  Logical flags for virtual order-alpha**3 corrections to unpol xsection
C  (useful diagnostics to polarized calculation)
C  
      LDIAG=.TRUE.
      LBF=.TRUE.
C
C  The value of (hbar*c)**2 in units of mb-GeV**2
C
      BARN=3.8937966D-1
C
C  The electron and photon energies in GeV
C
C  The incident e- is moving in the +z direction
C
C  (45.65 is SLC/LEP, 500. is NLC, 27.5 for HERA)
C
      EB=18.0d0
C
C  The incident photon is moving in the -z direction
C
C  (2.3305d-9 is frequency-doubled Nd:YAG, 2.41d-9 for Ar-ion laser)
C
      EPHOT=2.3305d-9
C
C  The initial state electron spin direction (in its rest frame)
C
      SPIN(1)=0.d0
      SPIN(2)=1.d0
      SPIN(3)=0.d0
C
C  Insure proper normalization
C
      SNORM=SQRT(SPIN(1)**2+SPIN(2)**2+SPIN(3)**2)
      SPIN(1)=SPIN(1)/SNORM
      SPIN(2)=SPIN(2)/SNORM
      SPIN(3)=SPIN(3)/SNORM
C
C  Initialize the random number generation
C
      NTOTIN=0
      NTO2IN=0
C      CALL RMARIN(123456789,NTOTIN,NTO2IN)
      CALL RMARIN(54217137,NTOTIN,NTO2IN)
C
C  Initialize event weight accounting system
C
      CALL WGTHST(0,NEM,PP,NGAM,QP,NEP,PPOS,WGT)
C
C  Print adjustable parameters
C
      print 2000, eb,ephot*1.d9
2000  format(1x,'Beam energy = ',F7.3,' GeV, photon energy = ',F7.4,
     >       ' eV')
      print 2100, spin
2100  format(1x,'electron spin direction = ',2(f7.3,','),f7.3)
      print 3000, xmg*1000.,KGMIN*1000.
3000  format(1x,'photon mass = ',e11.4,' MeV, photon energy cutoff = ',
     >       e11.4,' MeV')
C
C  Set the number of trials
C  (4X10**6 trials takes 3-4 CPU hrs on an IBM RS6000 workstation)
C
      NTRY=500000
C
C  Generate the 2-body final states
C
      CALL COMTN2
C
C  Generate the e+gamma+gamma final states
C
C      CALL COMEGG
C
C  Generate the e-e+e- final states (returns immediately if E_cm<3m_e)
C
C      CALL COMEEE
C
C  Print the answers
C
C      CALL WGTHST(2,NEM,PP,NGAM,QP,NEP,PPOS,WGT)
C        print *,PP(1,1),PP(2,1)
C        print *, WGT(1),WGT(2),WGT(3),WGT(4)
C
      STOP
      END
C
      SUBROUTINE COMTN2
C
C **********************************************************************
C * This routine calculates the laboratory cross sections for 2=>2     *
C * Compton scattering to order alpha3 using the expressions given in  *
C * Stuart and Gongora, Z.Phys. C42,617 (1989) and  Tsai, DeRaad, and  *
C * Milton Phys Rev D6, 1411 (1972).                                   *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  The block /CONTRL/ contains the adjustable parameters needed for the
C  calculation and some useful physical constants: EB - e- beam energy in
C  GeV, EPHOT - photon energy in GeV, XME - electron mass (in GeV),
C  XMG - the fictitious photon mass (GeV), KGMIN - the minimum resolvable
C  photon energy (GeV), ALPHA - the fine structure constant, PI - 3.1416...
C  ROOT2 - SQRT(2), BARN - (hbar*c)**2 in units of mb-GeV**(-2),
C  SPIN(3) - the initial state electron spin (in the me rest frame),
C  LDIAG - logical flag to turn on diagnostic T-DR-M unpol xsection, 
C  LBF - logical flag to turn on diagnostic Brown+Feynman unpol xsection, 
C  NTRY - the number of event trials to generate with COMTN2,COMEGG,COMEEE
C         (COMEEE generates NTRY/20)
C
      COMMON /CONTRL/ EB,EPHOT,XME,XMG,KGMIN,ALPHA,PI,ROOT2,BARN,
     >                SPIN(3),LDIAG,LBF,NTRY   
      LOGICAL LDIAG,LBF
      REAL*8 KGMIN
C
      COMPLEX*16 EPP,EMPP,EMP,EPPP,SM,SP,LNK,LNT,G0,G0K,G0T
      REAL*8 P(4),S(4),PP(4),SPR(4),Q(4),QP(4),WGT(4)
      REAL*8 P1(4),P2(4),P1P(4,2),P2P(4,2),PLAB(4)
      REAL*8 KAP,K2,LNLAM,LNEGMN,ILAM,KT,JFAC
      COMPLEX*16 A1(6)
      COMPLEX*16 ELE(6,2,2,2),M0(6),M1(6),MAT0(2,2,2),MAT1(2,2,2)
C
C  RVEC contains the random numbers
C
        REAL*4 RVEC(2)
        integer, dimension(:), allocatable :: seed 
        integer, dimension(8) :: values  
        integer :: seed_size
        call date_and_time(VALUES=values) 
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        seed(:) = values(:) 
        call random_seed(put=seed) 
C
C  The center-of-mass energy**2 of the e-gamma system
C
      XME2=XME**2
      XMG2=XMG**2
      ARG=EB**2-XME2
      IF(ARG.LE.0.d0) ARG=0.d0
      PB=SQRT(ARG)
      S0=XME2+2.D0*(EB+PB)*EPHOT
C
C  The transformation to the lab frame
C
      PLAB(3)=(PB-EPHOT)
      PLAB(2)=0.d0
      PLAB(1)=0.d0
      PLAB(4)=SQRT(S0+PLAB(3)**2)
C
C  Define the initial state 4-vectors in the cm system
C
C  The electron momentum and spin
C
      EE=(S0+XME2)/(2.D0*SQRT(S0))
      PE=SQRT(EE**2-XME2)
      P(4)=EE
      P(3)=PE
      P(2)=0.D0
      P(1)=0.D0
      S(4)=PE/XME*SPIN(3)
      S(3)=EE/XME*SPIN(3)
      S(2)=SPIN(2)
      S(1)=SPIN(1)
C
C  The auxiliary momenta
C
      DO J=1,4
        P1(J)=(P(J)+XME*S(J))/2.d0
        P2(J)=(P(J)-XME*S(J))/2.d0
      ENDDO
C
C  The photon energy/momentum 
C
      EG=(S0-XME2)/(2.D0*SQRT(S0))
      PG=EG
      Q(4)=EG
      Q(3)=-PG
      Q(2)=0.D0
      Q(1)=0.D0
C
C  Useful variables (see T-D-M for definitions)
C
      KAP=-2.d0*(P(4)*Q(4)-P(3)*Q(3)-P(2)*Q(2)-P(1)*Q(1))/XME2
      LNK=DCMPLX(LOG(ABS(KAP)),-PI)
      G0K=G0(KAP)
      LNLAM=LOG(XMG/XME)
      LNEGMN=LOG(2.d0*KGMIN/XMG)-0.5d0
      C1=PI**2/6.d0
      ALPHA2=ALPHA**2
      ALPHPI=ALPHA/PI
C
C  Keep track of 2-2 cross section
C
      XS220=0.d0
      XS2202=0.d0
      XS221=0.d0
      XS2212=0.d0
      IF(LDIAG) THEN
        XSD0=0.D0
        XSD02=0.D0
        XSD1=0.D0
        XSD12=0.D0
      ENDIF
      IF(LBF) THEN
        XSF0=0.D0
        XSF02=0.D0
        XSF1=0.D0
        XSF12=0.D0
      ENDIF
C
C  Calculate the number of trials
C
      NTRIAL=NTRY
C
C  Calculate the phase space per trial
C
      PHASPC=2.D0*2.d0*PI/FLOAT(NTRIAL)
C
C  Calculate the cross section normalization
C
      CNORM=PHASPC*BARN/(64.d0*PI**2*(1.d0-KAP)*XME2)
C
C  Calculate the diagnostic cross section normalizations
C
      DNORM0=PHASPC*BARN*ALPHA2/(2.d0*(1.d0-KAP)*XME2)
      DNORM1=PHASPC*BARN*ALPHA**3/(2.d0*PI*(1.d0-KAP)*XME2)
C
C  Loop over many event trials
C
      DO IEPHOT=1,NTRIAL
C        IF(mod(IEPHOT,10000).EQ.0) print *, IEPHOT
C
C  Choose 2 random numbers
C
C        CALL RANMAR(RVEC,2)
        CALL RANDOM_NUMBER(RVEC(1))
        CALL RANDOM_NUMBER(RVEC(2))
C
C  Generate direction of the scattered electron
C
        COST=2.d0*RVEC(1)-1.d0
        SINT=SQRT(1.d0-COST**2)
        PHI=2.d0*PI*RVEC(2)
        COSY=SINT*SIN(PHI)
        COSX=SINT*COS(PHI)
C
C  Define the final state 4-vectors in the cm system
C
C  The electron momentum and spin
C
        PP(4)=EE
        PP(3)=PE*COST
        PP(2)=PE*COSY
        PP(1)=PE*COSX
C
C  There are two final state electron spins, SPR and -SPR
C
        SPR(4)=PE/XME
        SPR(3)=EE/XME*COST
        SPR(2)=EE/XME*COSY
        SPR(1)=EE/XME*COSX
C
C  The auxiliary momenta for the two spin states
C
        DO J=1,4
          P1P(J,1)=(PP(J)+XME*SPR(J))/2.d0
          P2P(J,1)=(PP(J)-XME*SPR(J))/2.d0
          P1P(J,2)=P2P(J,1)
          P2P(J,2)=P1P(J,1)
        ENDDO
C
C  The photon momentum 
C
        QP(4)=EG
        QP(3)=-PG*COST
        QP(2)=-PG*COSY
        QP(1)=-PG*COSX
C
C  Useful variables (all defined in Tsai,DeRaad,+Milton, eqs 1+12)
C
        TAU=2.d0*(P(4)*QP(4)-P(3)*QP(3)-P(2)*QP(2)-P(1)*QP(1))/XME2
        K2=-XME2*(KAP+TAU)/4.d0
        PK=-XME2*(KAP-TAU)/4.d0
        T=KAP+TAU
        KT=KAP*TAU
        D=T-KT
        LNT=DCMPLX(LOG(TAU),0.d0)
        G0T=G0(TAU)
        SINHY=SQRT(-T)/2.d0
        Y=LOG(SINHY+SQRT(SINHY**2+1.d0))
        TY=2.d0*Y
        COTHTY=1.d0/TANH(TY)
        COTHY=1.d0/TANH(Y)
        SINHTY=SINH(TY)
        YCSHTY=Y/SINHTY
        HY=LOG(2.d0*SINHY)-0.5d0*Y+(C1-DDILOG(EXP(-TY)))/TY
        H2Y=LOG(2.d0*SINHTY)-Y+(C1-DDILOG(EXP(-2.d0*TY)))/(2.d0*TY)
        ILAM=2.d0*(1.d0-TY*COTHTY)*LNLAM-4.d0*Y*COTHTY*(2.d0*HY-H2Y)
        JFAC=-ALPHPI*
     >       (2.d0*(1.d0-TY*COTHTY)*LNEGMN+4.d0*Y*COTHTY*(H2Y-1.d0))
C
C  Define the eps(l')*L_i*eps'(l) tensor for i=1,6 invariant matrix elements; 
C  l'=1,2 final photon helicity states;l=1,2 initial photon helicity states; 
C  and k=1,2 final electron helicity states (see Gongora+Stuart eqs 4.8-4.11)
C
        DO K=1,2
C
C i=1
C
          ELE(1,1,1,K)=K2*SM(Q,QP)/SP(Q,QP)*(SP(P2,P1P(1,K))
     >                 *SM(P1P(1,K),P2P(1,K))-SP(P1,P2)*SM(P1,P2P(1,K)))
          ELE(1,2,2,K)=(SP(Q,QP)/SM(Q,QP))**2*ELE(1,1,1,K)
          ELE(1,1,2,K)=DCMPLX(0.d0,0.d0)
          ELE(1,2,1,K)=DCMPLX(0.d0,0.d0)
C
C i=2
C
          ELE(2,1,1,K)=K2/2.d0*SM(Q,QP)/SP(Q,QP)*((SP(P1,P2)
     >                 *SM(P1P(1,K),P2P(1,K)))/XME2
     >                 *(SP(P1P(1,K),QP)*SM(P1,QP)
     >                   -SP(P1P(1,K),Q)*SM(P1,Q))
     >                 +(SP(P2,QP)*SM(P2P(1,K),QP)
     >                   -SP(P2,Q)*SM(P2P(1,K),Q)))
          ELE(2,2,2,K)=-(SP(Q,QP)/SM(Q,QP))**2*ELE(2,1,1,K)
          ELE(2,1,2,K)=DCMPLX(0.d0,0.d0)
          ELE(2,2,1,K)=DCMPLX(0.d0,0.d0)
C
C i=3
C
          ELE(3,1,1,K)=SM(Q,QP)/SP(Q,QP)*(SP(P1,P2)*SP(P1P(1,K),Q)
     >                 *SM(P1,Q)*SM(P1P(1,K),P2P(1,K))
     >                 -XME2*SP(P2,Q)*SM(P2P(1,K),Q))
     >                +2.d0*K2/SP(Q,QP)**2*(SP(Q,QP)*SP(P1,P2)*SM(P1,Q)*
     >                 SM(P2P(1,K),QP)+SP(P2,QP)*SP(P1P(1,K),Q)*SM(Q,QP)
     >                 *SM(P1P(1,K),P2P(1,K)))
     >                +(1.d0-PK/K2)*ELE(1,1,1,K)
          ELE(3,2,2,K)=SP(Q,QP)/SM(Q,QP)*(SP(P1,P2)*SP(P1P(1,K),Q)
     >                 *SM(P1,Q)*SM(P1P(1,K),P2P(1,K))
     >                 -XME2*SP(P2,Q)*SM(P2P(1,K),Q))
     >                -2.d0*K2/SM(Q,QP)**2*(SP(Q,QP)*SP(P1,P2)*SM(P1,QP)
     >                 *SM(P2P(1,K),Q)+SP(P2,Q)*SP(P1P(1,K),QP)*SM(Q,QP)
     >                 *SM(P1P(1,K),P2P(1,K)))
     >                +(1.d0-PK/K2)*ELE(1,2,2,K)
          ELE(3,1,2,K)=DCMPLX(0.d0,0.d0)
          ELE(3,2,1,K)=DCMPLX(0.d0,0.d0)
C
C i=4
C
          ELE(4,1,1,K)=DCMPLX(0.d0,0.d0)
          ELE(4,2,2,K)=DCMPLX(0.d0,0.d0)
          EPP=(SP(P1,QP)*SM(P1,Q)+SP(P2,QP)*SM(P2,Q))/SP(Q,QP)/ROOT2
          EMPP=SP(Q,QP)/SM(Q,QP)*EPP
          ELE(4,2,1,K)=(2.D0*ROOT2*K2*EPP)/(XME2*SM(Q,QP))*
     >                 (XME2*SP(P2,QP)*SM(P2P(1,K),Q)
     >                  -SP(P1,P2)*SP(P1P(1,K),QP)*SM(P1,Q)
     >                   *SM(P1P(1,K),P2P(1,K)))
          EMP=-(SP(P1,Q)*SM(P1,QP)+SP(P2,Q)*SM(P2,QP))/SM(Q,QP)/ROOT2
          EPPP=SM(Q,QP)/SP(Q,QP)*EMP
          ELE(4,1,2,K)=-(2.D0*ROOT2*K2*EMP)/(XME2*SP(Q,QP))*
     >                 (XME2*SP(P2,Q)*SM(P2P(1,K),QP)
     >                  -SP(P1,P2)*SP(P1P(1,K),Q)*SM(P1,QP)
     >                   *SM(P1P(1,K),P2P(1,K)))
C
C i=5
C
          ELE(5,1,1,K)=DCMPLX(0.d0,0.d0)
          ELE(5,2,2,K)=DCMPLX(0.d0,0.d0)
          ELE(5,2,1,K)=(K2/XME2)*EPP*EMPP*
     >                 (SP(P2,P1P(1,K))*SM(P1P(1,K),P2P(1,K))
     >                  -SP(P1,P2)*SM(P1,P2P(1,K)))
          ELE(5,1,2,K)=(EMP*EPPP)/(EPP*EMPP)*ELE(5,2,1,K)
C
C i=6
C
          ELE(6,1,1,K)=DCMPLX(0.d0,0.d0)
          ELE(6,2,2,K)=DCMPLX(0.d0,0.d0)
          ELE(6,2,1,K)=EPP*EMPP/XME2*
     >                 (SP(P1,P2)*SP(P1P(1,K),Q)*SM(P1,Q)*
     >                  SM(P1P(1,K),P2P(1,K))
     >                 -XME2*SP(P2,Q)*SM(P2P(1,K),Q))
     >                -PK/(2.d0*K2)*ELE(4,2,1,K)
          ELE(6,1,2,K)=EMP*EPPP/XME2*
     >                 (SP(P1,P2)*SP(P1P(1,K),Q)*SM(P1,Q)*
     >                  SM(P1P(1,K),P2P(1,K))
     >                 -XME2*SP(P2,Q)*SM(P2P(1,K),Q))
     >                -PK/(2.d0*K2)*ELE(4,1,2,K)
        ENDDO
C
C  Next, calculate the zeroth-order invariant amplitudes (G+S, eq 4.13)
C
        RM0=-32.d0*PI*ALPHA/(KAP*TAU)
        M0(1)=DCMPLX(RM0,0.d0)
        M0(2)=DCMPLX(RM0,0.d0)
        M0(3)=DCMPLX(0.d0,0.d0)
        M0(4)=DCMPLX(-RM0,0.d0)
        M0(5)=DCMPLX(0.d0,0.d0)
        M0(6)=DCMPLX(0.d0,0.d0)
C
C  Next, calculate the first-order helicity amplitudes
C  (Tsai, DeRaad, and Milton, eqs 6-11)
C
        A1(1)=ALPHA2*4.d0*SQRT(-T)/KAP*(-T*ILAM/KT-4.d0*Y**2/T-1.D0
     >        -2.D0*(1.D0/KAP+1.D0/TAU)+G0K-G0T+LNK*(2.D0*
     >        (2.D0*KAP+TAU-2.D0)/KAP*YCSHTY+(KAP-2.D0)/KAP/(KAP-1.D0))        
     >        +LNT*(2.D0*(KAP-2.D0)/TAU*YCSHTY
     >              +(3.D0*TAU-2.D0)/TAU/(TAU-1.D0)))       
        A1(2)=-ALPHA2*4.d0*SQRT(D)/KAP*(T*ILAM/KT+0.5D0+2.D0/KAP
     >        -(KAP-4.D0)/(2.D0*TAU)+2.D0/(T*D)*(KAP*(KAP+2.D0)
     >        -TAU*(KAP-2.D0))*Y**2+T*(KAP-1.D0)*G0K/(2.D0*D)
     >        +T*(KAP-TAU*(KAP-2.D0))*G0T/(2.D0*TAU*D)+LNT*(
     >        (KT*TAU+TAU*(KAP**2-6.D0*KAP+4.D0)-2.D0*KAP*(KAP-2.D0))
     >        *YCSHTY/(TAU*D)-(KAP+3.D0*TAU-2.D0)/(TAU*(TAU-1.D0)))
     >        +LNK*(-2.D0/KAP+(2.D0*TAU**2*(KAP-1.D0)+TAU*
     >        (3.D0*KAP**2-10.D0*KAP+4.D0)+KAP*(KAP**2-4.D0*KAP+4.D0))
     >        *YCSHTY/(KAP*D)))
C
C  A1(3) fixes overall sign error in eq 8 of T-DR-M
C
        A1(3)=ALPHA2*4.d0*SQRT(-T)/KAP*(T*(KAP-1.D0)*ILAM/KT
     >        +(3.D0*KAP-2.D0)/TAU-2.D0/KAP+2.D0-4.D0*Y**2/T-T*G0T/TAU
     >        +LNK*(-2.D0/KAP*(KAP**2-4.D0*KAP+2.D0+TAU*(KAP-1.D0))
     >        *YCSHTY-(3.D0*KAP-2.D0)/KAP)       
     >        -LNT*(2.D0/TAU*(KAP**2-3.D0*KAP+2.D0+KT)*YCSHTY
     >        +(3.D0*TAU*(KAP-1.D0)-4.D0*KAP+2.D0)/TAU/(TAU-1.D0)))               
        A1(4)=-ALPHA2*4.d0*SQRT(D)/KAP*(-(KAP**2+T)*ILAM/KT-TY/COTHY
     >        -2.D0/D*(KAP+TAU*(2.D0*KAP+1.D0)+TAU**2)*Y*HY*COTHY
     >        -0.5D0/(KAP-1.D0)-0.5D0*(KAP+1.D0)/(TAU-1.D0)-KAP/TAU
     >        -2.D0/KAP-2.D0/TAU-(TAU**2+TAU*(2.D0*KAP+1.D0)+2.D0*KAP**2        
     >        -KAP-8.D0)*Y**2/D+KAP/(2.D0*D)*(-KAP**2-2.D0*KAP+4.D0-TAU)
     >        *G0K+TAU/(2.D0*D)*(-TAU**2-2.D0*KT-KAP**2+3.D0*KAP+4.D0)
     >        *G0T+LNK*((TAU**2*(2.D0-KAP)
     >        +TAU*(-KAP**3-3.D0*KAP**2+8.D0*KAP-4.D0)
     >        +KAP*(-KAP**3+10.D0*KAP-4.D0))*YCSHTY/(KAP*D)
     >        +(5.D0*KAP**2-8.D0*KAP+4.D0)/(2.D0*KAP*(KAP-1.D0)**2))
     >       +LNT*((-KAP**3*(TAU-2.D0)-KAP**2*(TAU-1.D0)*(3.D0*TAU-2.D0)
     >        -KAP*(TAU-1.D0)*(3.D0*TAU**2-4.D0*TAU-4.D0)
     >        +TAU*(-TAU**3+2.D0*TAU**2+6.D0*TAU-4.D0))*YCSHTY/(TAU*D)
     >        -(KAP*(-3.D0*TAU**2+6.D0*TAU-4.D0)-2.D0*TAU**3+TAU**2
     >        +4.D0*TAU-4.D0)/(2.D0*TAU*(TAU-1.D0)**2)))
        A1(5)=-ALPHA2*4.d0*SQRT(-T)/KAP*(-D*ILAM/KT-D/(2.d0*TAU)
     >        *(4.D0/KAP+1.D0/(TAU-1.D0))+2.D0/D*(KAP*(2.D0*KAP-1.D0)
     >        -TAU*(TAU+1.D0))*Y*HY*COTHY-(KAP-TAU)*YCSHTY-1.D0/D*(
     >        TAU**2-TAU*(2.D0*KAP-1.D0)-2.D0*KAP**2+9.D0*KAP-8.D0)*Y**2
     >        -KAP/(2.D0*D)*(-TAU*(KAP-1.D0)-2.D0*KAP**2+7.D0*KAP-4.D0)
     >        *G0K+TAU/(2.D0*D)*(KAP*(TAU-2.D0)-TAU**2+4.D0)*G0T+LNK*(
     >        ((KAP**2-3.D0*KAP+2.D0)*TAU**2+TAU*(3.D0*KAP**3
     >        -10.D0*KAP**2+12.D0*KAP-4.D0)+KAP*(2.D0*KAP**3
     >        -11.D0*KAP**2+14.D0*KAP-4.D0))*YCSHTY/(KAP*D)
     >        -2.D0*(KAP-1.D0)/KAP)+LNT*((KAP**2*(TAU**2-2.D0*TAU+2.D0)
     >        -4.D0*KAP*(TAU**2-TAU+1.D0)-TAU*(TAU**3-2.D0*TAU**2
     >        -6.D0*TAU+4.D0))*YCSHTY/(TAU*D)-(TAU*(TAU-1.D0)*KAP
     >        -2.D0*TAU**3+TAU**2+4.D0*TAU-4.D0)
     >        /(2.D0*TAU*(TAU-1.D0)**2)))       
        A1(6)=-ALPHA2*4.d0/(KAP*SQRT(D))*(-D**2*ILAM/KT-D/2.D0*(-4.D0
     >        +4.D0/KAP+4.D0/TAU+TAU/(TAU-1.D0))+2.D0/D*(-KAP**2*(KAP**2
     >        -3.D0*KAP+1.D0)-KT*(2.D0*KAP**2-5.D0*KAP+2.D0)+TAU**2*(
     >        KAP-1.D0)-TAU**3)*Y*HY*COTHY-2.D0*(KAP*(KAP-3.D0)+TAU)
     >        *Y/COTHY+1.D0/D*(-KAP**4+KAP**3*(5.D0-2.D0*TAU)+KAP**2*(
     >        -2.D0*TAU**2+13.D0*TAU-11.D0)+KAP*(5.D0*TAU**2-20.D0*TAU
     >        +8.D0)+TAU*(-TAU**2-TAU+8.D0))*Y**2+KAP/(2.D0*D)*(-TAU**2
     >        *(KAP-1.D0)**2-TAU*(2.D0*KAP**3-10.D0*KAP**2+13.D0*KAP
     >        -4.D0)-KAP*(KAP**3-5.D0*KAP**2+8.D0*KAP-4.D0))*G0K
     >        +TAU/(2.D0*D)*(-KAP**2*(TAU**2-3.D0*TAU+3.D0)+KAP*(
     >        3.D0*TAU**2-7.D0*TAU+4.D0)-TAU*(TAU**2-4.D0))*G0T+LNK*(
     >        1.D0/(KAP*D)*(-TAU**3*(KAP-1.D0)**2*(KAP-2.D0)+TAU**2*(
     >        -3.D0*KAP**4+14.D0*KAP**3-24.D0*KAP**2+18.D0*KAP-4.D0)+KT*
     >        (-3.D0*KAP**4+19.D0*KAP**3-39.D0*KAP**2+30.D0*KAP-8.D0)
     >        +KAP**2*(-KAP**4+7.D0*KAP**3-16.D0*KAP**2+14.D0*KAP-4.D0))
     >        *YCSHTY+(KAP-2.D0)/(2.D0*KAP*(KAP-1.D0))*(TAU*(3.D0*KAP**2
     >        -5.D0*KAP+2.D0)+2.D0*KAP**3-3.D0*KAP**2+2.D0*KAP))+LNT*(
     >        1.D0/(TAU*D)*(KAP**3*(-TAU**3+3.D0*TAU**2-3.D0*TAU+2.D0)
     >        +KAP**2*(-TAU**4+8.D0*TAU**3-14.D0*TAU**2+10.D0*TAU-4.D0)
     >        +KT*(2.D0*TAU**3-13.D0*TAU**2+14.D0*TAU-8.D0)-TAU**2*(
     >        TAU**3-2.D0*TAU**2-6.D0*TAU+4.D0))*YCSHTY+(TAU*(
     >        2.D0*TAU**3-TAU**2-4.D0*TAU+4.D0)-KAP*(TAU-1.D0)*(
     >        3.D0*TAU**2-4.D0*TAU+4.D0))/(2.D0*TAU*(TAU-1.D0)**2)))
C
C  Now convert them into first-order invariant amplitudes by inverting
C  eq 2 of T-DR-M
C
         M1(1)=2.d0/(KAP*SQRT(-T))*(A1(1)+A1(3))
     >        -2.d0*(KAP**2+D+T)/(KAP*T*SQRT(D))*A1(2)
         M1(3)=-2.d0/(KAP*SQRT(-T))*(A1(1)+A1(3))
     >        +2.d0*(-KAP**2+D+T)/(KAP*T*SQRT(D))*A1(2)
         M1(2)=(2.d0-KAP)/KAP*M1(1)+2.d0/KAP*M1(3)
     >        -8.d0/SQRT(-T)**3*A1(3)       
         M1(5)=16.d0/D*(A1(5)/SQRT(-T)-A1(6)/SQRT(D))
         M1(4)=-(2.d0/(KAP*SQRT(D))*A1(4)
     >         +2.D0*(KAP**2-T)/(KAP*SQRT(D)**3)*A1(6)
     >         +(KAP**2-KT)/(4.D0*KAP**2)*M1(5))
         M1(6)=-(8.D0*KAP/SQRT(D)**3*A1(6)+2.D0*M1(4)+M1(5))
C
C  Form the zeroth-order and first-order matrix elements using G+S, eq 4.12
C
        DO J=1,2
          DO K=1,2
            DO L=1,2
              MAT0(L,K,J)=DCMPLX(0.d0,0.d0)
              MAT1(L,K,J)=DCMPLX(0.d0,0.d0)
              DO M=1,6
                MAT0(L,K,J)=MAT0(L,K,J)+ELE(M,L,K,J)*M0(M)
                MAT1(L,K,J)=MAT1(L,K,J)+ELE(M,L,K,J)*M1(M)
              ENDDO
              MAT0(L,K,J)=MAT0(L,K,J)/(2.d0*XME2**2)
              MAT1(L,K,J)=MAT1(L,K,J)/(2.d0*XME2**2)
            ENDDO
          ENDDO
        ENDDO
C
C  Zeroth order cross section
C
C  [WRGHT is the weight for a right-handed helicity (Left Circular Pol) photon
C  impinging on the electron and WLEFT is for a left-handed helicity (RCP)
C  photon]
C
        WLEFT0=0.d0
        WRGHT0=0.d0
        DO J=1,2
          DO L=1,2
            WRGHT0=WRGHT0+ABS(MAT0(L,1,J))**2
            WLEFT0=WLEFT0+ABS(MAT0(L,2,J))**2
          ENDDO
        ENDDO
        WRGHT0=WRGHT0*CNORM
        WLEFT0=WLEFT0*CNORM
C
C  Event weights 1+3 are the unpolarized cross sections (the averages
C  (of the left- and right-handed cross sections)
C
C  Event weights 2+4 are the polarized cross sections and are defined as the
C  the difference in cross sections for left-handed (negative helicity) 
C  and right-handed (positive helicity) photons impinging (traveling in the
C  -z direction) on the electron (traveling in the +z direction).  If the
C  electron spin direction is chosen to lie along the +z direction, these
C  weights are the difference in cross sections, sigma(Jz=3/2)-sigma(Jz=1/2).
C
        WGT(1)=(WLEFT0+WRGHT0)/2.d0
        WGT(2)=(WLEFT0-WRGHT0)/2.d0
        XS220=XS220+WGT(1)
        XS2202=XS2202+WGT(1)**2
C
C  First order cross section
C
        WLEFT1=0.d0
        WRGHT1=0.d0
        DO J=1,2
          DO L=1,2
            WRGHT1=WRGHT1+ABS(MAT0(L,1,J))**2*JFAC
     >           +2.d0*REAL(MAT0(L,1,J)*CONJG(MAT1(L,1,J)))
            WLEFT1=WLEFT1+ABS(MAT0(L,2,J))**2*JFAC
     >           +2.d0*REAL(MAT0(L,2,J)*CONJG(MAT1(L,2,J)))
          ENDDO
        ENDDO
        WRGHT1=WRGHT1*CNORM
        WLEFT1=WLEFT1*CNORM
        WGT(3)=(WLEFT1+WRGHT1)/2.d0
        WGT(4)=(WLEFT1-WRGHT1)/2.d0
        XS221=XS221+WGT(3)
        XS2212=XS2212+WGT(3)**2
C        print *, WGT(1),WGT(2),WGT(3),WGT(4)
C
C  Do diagnostic (unpolarized) xsections
C
        IF(LDIAG) THEN
          U=4.d0*(1.d0/KAP+1.d0/TAU)**2
     >      -4.d0*(1.d0/KAP+1.d0/TAU)-(KAP/TAU+TAU/KAP)
          DXS0=DNORM0*U
          DXS1=-DNORM1*(ILAM*U-4.d0*(1.d0-1.d0/KAP-1.D0/TAU)*Y/COTHY
     >   +2.D0*(4.D0+(TAU+2.D0)/KAP+(KAP+2.D0)/TAU)*Y*HY*COTHY
     >   +(3.D0*(KAP/TAU+TAU/KAP)-8.D0/KT+4.D0)*Y**2+
     >   LNK*((4.D0*KAP+3.D0*TAU+2.D0+(TAU**2+2.D0*TAU-24.D0)/KAP+
     >   2.D0*(KAP**2-12.D0)/TAU+16.D0/KT-8.D0*(TAU-2.D0)/KAP**2)*YCSHTY
     >   +1.D0+3.D0/TAU-7.D0/KT+(3.D0*TAU+16.D0)/(2.*KAP)
     >   +(3.D0*TAU-16.D0)/(2.D0*KAP**2)+(2.D0*KAP-TAU**2-KT*KAP)/(2.D0*
     >   KAP*KT*(KAP-1.D0))-(2.D0*KAP**2+TAU)/(2.D0*TAU*(KAP-1.D0)**2)) 
     >  +LNT*((4.D0*TAU+3.D0*KAP+2.D0+(KAP**2+2.D0*KAP-24.D0)/TAU+
     >   2.D0*(TAU**2-12.D0)/KAP+16.D0/KT-8.D0*(KAP-2.D0)/TAU**2)*YCSHTY
     >   +1.D0+3.D0/KAP-7.D0/KT+(3.D0*KAP+16.D0)/(2.*TAU)
     >   +(3.D0*KAP-16.D0)/(2.D0*TAU**2)+(2.D0*TAU-KAP**2-KT*TAU)/(2.D0*
     >   TAU*KT*(TAU-1.D0))-(2.D0*TAU**2+KAP)/(2.D0*KAP*(TAU-1.D0)**2))
     >   +(TAU+2.D0)/(2.D0*TAU*(KAP-1.D0))-2.D0*(TAU-4.D0)/KAP**2
     >   +(KAP+2.D0)/(2.D0*KAP*(TAU-1.D0))-2.D0*(KAP-4.D0)/TAU**2
     >   -(3.*TAU+22.D0)/(2.D0*KAP)-(3.D0*KAP+22.D0)/(2.D0*TAU)+16.D0/KT
     >+G0K*((KAP**2+KAP-3.D0)/TAU+2.D0/KAP+TAU/KAP**2+KAP+TAU/2.D0-1.D0)
     >+G0T*((TAU**2+TAU-3.D0)/KAP+2.D0/TAU+KAP/TAU**2+TAU+KAP/2.D0-1.D0)
     >     )
          WU0D=DXS0
          WU1D=DXS0*JFAC+DXS1
          XSD0=XSD0+WU0D
          XSD02=XSD02+WU0D**2
          XSD1=XSD1+WU1D
          XSD12=XSD12+WU1D**2
        ENDIF
        IF(LBF) THEN
C
C  Do the unpolarized cross section of Brown+Feynman, PR 85, 231 (1952)
C
C
          COSHY=SQRT(1.d0+SINHY**2)
          COSHTY=SQRT(1.d0+SINHTY**2)
          U=4.d0*(1.d0/KAP+1.d0/TAU)**2
     >      -4.d0*(1.d0/KAP+1.d0/TAU)-(KAP/TAU+TAU/KAP)
          FXS0=DNORM0*U
          FXS1=-ALPHA/PI*DNORM0*(
     >      2.d0*((1.d0-TY*COTHTY)*LNLAM-TY*COTHTY*(2.D0*HY-H2Y))*U
     >     +2.d0*(-4.d0*Y*SINHTY/KT*(2.d0-COSHTY)+TY*COTHY)*HY
     >     +REAL(LNK)*(4.d0*Y*COTHTY*(4.d0*COSHY**2/KT+(KAP-6.D0)/
     >      (2.D0*TAU*COSHTY)+4.D0/KAP**2-1.D0/KAP-TAU/(2.D0*KAP)
     >     -KAP/TAU-1.D0)+1.5D0*TAU/KAP**2+1.5D0*TAU/KAP+3.D0/TAU
     >     +1.D0-7.D0/KT+8.D0/KAP-8.D0/KAP**2+(2.D0*KAP-TAU**2-KAP*KT)/
     >      (2.D0*KAP*KT*(KAP-1.D0))-(2.D0*KAP**2+TAU)/
     >      (2.D0*TAU*(KAP-1.D0)**2))+(Y/SINHY)**2*(2.D0/KAP-7.D0*KAP/
     >      4.D0-0.75D0*TAU**2/KAP)-4.D0*Y/COTHY*(0.5D0-1.D0/KAP)
     >     +4.D0*(1.D0/KAP+1.D0/TAU)**2-12.D0/KAP-1.5*KAP/TAU
     >     -2.D0*KAP/TAU**2+(KAP/TAU+0.5D0)/(KAP-1.D0)+REAL(G0K)*(
     >      KAP**2/TAU+TAU/KAP**2+KAP/TAU+KAP+TAU/2.D0+2.D0/KAP
     >     -3.D0/TAU-1.D0)
     >     +REAL(LNT)*(4.d0*Y*COTHTY*(4.d0*COSHY**2/KT+(TAU-6.D0)/
     >      (2.D0*KAP*COSHTY)+4.D0/TAU**2-1.D0/TAU-KAP/(2.D0*TAU)
     >     -TAU/KAP-1.D0)+1.5D0*KAP/TAU**2+1.5D0*KAP/TAU+3.D0/KAP
     >     +1.D0-7.D0/KT+8.D0/TAU-8.D0/TAU**2+(2.D0*TAU-KAP**2-TAU*KT)/
     >      (2.D0*TAU*KT*(TAU-1.D0))-(2.D0*TAU**2+KAP)/
     >      (2.D0*KAP*(TAU-1.D0)**2))+(Y/SINHY)**2*(2.D0/TAU-7.D0*TAU/
     >      4.D0-0.75D0*KAP**2/TAU)-4.D0*Y/COTHY*(0.5D0-1.D0/TAU)
     >     +4.D0*(1.D0/TAU+1.D0/KAP)**2-12.D0/TAU-1.5*TAU/KAP
     >     -2.D0*TAU/KAP**2+(TAU/KAP+0.5D0)/(TAU-1.D0)+REAL(G0T)*(
     >      TAU**2/KAP+KAP/TAU**2+TAU/KAP+TAU+KAP/2.D0+2.D0/TAU
     >     -3.D0/KAP-1.D0))
          WU0F=FXS0
          WU1F=FXS0*JFAC+FXS1
          XSF0=XSF0+WU0F
          XSF02=XSF02+WU0F**2
          XSF1=XSF1+WU1F
          XSF12=XSF12+WU1F**2
        ENDIF
C
C  Boost back to the lab frame and accumulate event weights
C
        CALL LORENZ(1,PLAB,1,1,PP)
        CALL LORENZ(1,PLAB,1,1,QP)
        CALL WGTHST(1,1,PP,1,QP,0,PDUM,WGT)
C        print result of e-gamma
        print *, QP, WGT ,PP             !print photons
C        print *, PP, WGT     !check qiao  print electron
       
      ENDDO
      PRINT 1000, NTRIAL,XS220,SQRT(XS2202)
1000  FORMAT(/,2x,'***** COMTN2 Finished *****',/,
     >       1x,I8,' e+g=>e+g trials, unpolarized xs0 = ',
     >       E12.5,'+-',E12.5)
      PRINT 2000, XS221,SQRT(XS2212)
2000  FORMAT(1x,8x,'                  unpolarized xs1 = ',
     >       E12.5,'+-',E12.5)
      IF(LDIAG) THEN
        PRINT 3000, XSD0,SQRT(XSD02)
3000    FORMAT(1x,8x,'    diagnostic    unpolarized xs0 = ',
     >         E12.5,'+-',E12.5)
        PRINT 4000, XSD1,SQRT(XSD12)
4000    FORMAT(1x,8x,'    diagnostic    unpolarized xs1 = ',
     >         E12.5,'+-',E12.5)
       ENDIF
       IF(LBF) THEN
        PRINT 5000, XSF0,SQRT(XSF02)
5000    FORMAT(1x,8x,'    Brown+Feynman unpolarized xs0 = ',
     >         E12.5,'+-',E12.5)
        PRINT 6000, XSF1,SQRT(XSF12)
6000    FORMAT(1x,8x,'    Brown+Feynman unpolarized xs1 = ',
     >         E12.5,'+-',E12.5)
      ENDIF
      RETURN
      END
C
      SUBROUTINE COMEGG
C
C ***********************************************************************
C * This routine calculates the laboratory cross sections for E+G=>E+2G *
C * Compton scattering (order alpha3) using the expressions given in    *
C * Stuart and Gongora, Z.Phys. C42,617 (1989).                         *
C ***********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  The block /CONTRL/ contains the adjustable parameters needed for the
C  calculation and some useful physical constants: EB - e- beam energy in
C  GeV, EPHOT - photon energy in GeV, XME - electron mass (in GeV),
C  XMG - the fictitious photon mass (GeV), KGMIN - the minimum resolvable
C  photon energy (GeV), ALPHA - the fine structure constant, PI - 3.1416...
C  ROOT2 - SQRT(2), BARN - (hbar*c)**2 in units of mb-GeV**(-2),
C  SPIN(3) - the initial state electron spin (in the me rest frame),
C  LDIAG - logical flag to turn on diagnostic T-DR-M unpol xsection, 
C  LBF - logical flag to turn on diagnostic Brown+Feynman unpol xsection, 
C  NTRY - the number of event trials to generate with COMTN2,COMEGG,COMEEE
C         (COMEEE generates NTRY/20)
C
      COMMON /CONTRL/ EB,EPHOT,XME,XMG,KGMIN,ALPHA,PI,ROOT2,BARN,
     >                SPIN(3),LDIAG,LBF,NTRY   
      LOGICAL LDIAG,LBF
      REAL*8 KGMIN
C
C  The block /VECT/ contains four-vectors needed by the DXXX functions
C
      COMMON /VECT/ DCONST,XME2,P(4),PP(4),P1(4),P2(4),P1P(4,2),P2P(4,2)
     >             ,LSET1
      LOGICAL LSET1
      REAL*8 S(4),SPR(4),Q(4),QP(4),QPP(4),MQ(4),QOUT(4,2)
      EQUIVALENCE (QP(1),QOUT(1,1)), (QPP(1),QOUT(1,2))
      REAL*8 WLEFT(2),WRIGHT(2),PLAB(4),WGT(4)
      COMPLEX*16 MAT(2,2,2,2,2),DPPP,DPPM,DPMP,DMPP,DPMM,DMPM,DMMP,DMMM
C
C  RVEC contains the random numbers
C
      REAL*4 RVEC(5)
        integer, dimension(:), allocatable :: seed
        integer, dimension(8) :: values
        integer :: seed_size
        call date_and_time(VALUES=values)
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        seed(:) = values(:) 
        call random_seed(put=seed) 
C
C  The center-of-mass energy**2 of the e-gamma system
C
      XME2=XME**2
      ARG=EB**2-XME2
      IF(ARG.LE.0.d0) ARG=0.d0
      PB=SQRT(ARG)
      S0=XME2+2.D0*(EB+PB)*EPHOT
      ROOTS=SQRT(S0)
      SXME2=S0+XME2
      TROOTS=2.d0*ROOTS
C
C  The transformation to the electron rest frame can lower the min energy
C
      EGMN=XME/ROOTS*KGMIN
      GAME=SXME2/(2.d0*XME*ROOTS)
      GAMBE=(S0-XME2)/(2.d0*XME*ROOTS)
C
C  The transformation to the lab frame
C
      PLAB(3)=(PB-EPHOT)
      PLAB(2)=0.d0
      PLAB(1)=0.d0
      PLAB(4)=SQRT(S0+PLAB(3)**2)
C
C  Define the initial state 4-vectors in the cm system
C
C  The electron momentum and spin
C
      EE=(S0+XME2)/(2.D0*SQRT(S0))
      PE=SQRT(EE**2-XME2)
      P(4)=EE
      P(3)=PE
      P(2)=0.D0
      P(1)=0.D0
C
C  Construct the spin vector
C
      S(4)=PE/XME*SPIN(3)
      S(3)=EE/XME*SPIN(3)
      S(2)=SPIN(2)
      S(1)=SPIN(1)
C
C  The auxiliary momenta
C
      DO J=1,4
        P1(J)=(P(J)+XME*S(J))/2.d0
        P2(J)=(P(J)-XME*S(J))/2.d0
      ENDDO
C
C  The photon momentum and it's negative
C
      EG=(S0-XME2)/(2.D0*SQRT(S0))
      Q(4)=EG
      Q(3)=-EG
      Q(2)=0.D0
      Q(1)=0.D0
      DO J=1,4
        MQ(J)=-Q(J)
      ENDDO
C
C  Useful variables 
C
      ALPHA2=ALPHA**2
      TWOPI=2.d0*PI
      DCONST=2.d0*ROOT2*(4.d0*PI*ALPHA)**1.5d0
      PEMAX=(S0-XME**2)
     >     /(2.d0*SQRT(S0))
      PGMAX=PEMAX
      EEMAX=SQRT(PEMAX**2+XME2)
      EGMAX=PGMAX
      GG2=EGMAX+EGMN
      RG=EGMAX/EGMN
      CG=1.d0/(2.d0*LOG(RG))
      EE2=EEMAX+EGMN
      EE1=EE2-XME
      RE=EGMN/EE1
      CE=1.d0/LOG(EE1/EGMN)
C
C  Calculate the number of trials
C
      NTRIAL=NTRY
C
C  Calculate the phase space per trial
C
      PHASPC=(TWOPI)**2*2.d0/FLOAT(NTRIAL)
C
C  Calculate the cross section normalization (include ID particles factor)
C
      CNORM=PHASPC*BARN/(8.d0*TWOPI**5*4.d0*(EB+PB)*EPHOT)/2.d0
C
C  This program doesn't generate the tree-level cross section
C
      WGT(1)=0.d0
      WGT(2)=0.d0
C
C  Keep track of successful trials
C
      NGOOD=0
      XS23=0.d0
      XS232=0.d0
C
C  Loop over many event trials
C
      DO IEPHOT=1,NTRIAL
C
C  Choose 5 random numbers
C
C        CALL RANMAR(RVEC,5)
        CALL RANDOM_NUMBER(RVEC(1))
        CALL RANDOM_NUMBER(RVEC(2)) 
        CALL RANDOM_NUMBER(RVEC(3)) 
        CALL RANDOM_NUMBER(RVEC(4))
        CALL RANDOM_NUMBER(RVEC(5))
C
C  Generate energies for the electron and one of the photons
C
        R1=RVEC(1)
        R2=RVEC(2)
        EE=EE2-EE1*RE**R1
        XG=RG**(2.d0*R2-1.d0)
        EG=GG2*XG/(1.d0+XG)
C
C  Determine if these energies can balance energy/momentum
C
        PG=EG
        PE=SQRT(EE**2-XME2)
        IF(PE.EQ.0.d0.OR.PG.EQ.0.d0) GO TO 1000
        COSTEG=(SXME2+2.d0*EE*EG-TROOTS*(EE+EG))/(2.d0*PE*PG)
        IF(ABS(COSTEG).GT.1.d0) GO TO 1000
C
C  Generate direction of the scattered electron
C
        COST=2.d0*RVEC(3)-1.d0
        SINT=SQRT(1.d0-COST**2)
        PHI=TWOPI*RVEC(4)
        COSY=SINT*SIN(PHI)
        COSX=SINT*COS(PHI)
C
C  Define the final state 4-vectors in the cm system
C
C  The electron momentum and spin
C
        PP(4)=EE
        PP(3)=PE*COST
        PP(2)=PE*COSY
        PP(1)=PE*COSX
C
C  There are two final state electron spins, SPR and -SPR
C
        SPR(4)=PE/XME
        SPR(3)=EE/XME*COST
        SPR(2)=EE/XME*COSY
        SPR(1)=EE/XME*COSX
C
C  The auxiliary momenta for the two spin states
C
        DO J=1,4
          P1P(J,1)=(PP(J)+XME*SPR(J))/2.d0
          P2P(J,1)=(PP(J)-XME*SPR(J))/2.d0
          P1P(J,2)=P2P(J,1)
          P2P(J,2)=P1P(J,1)
        ENDDO
C
C  Generate azimuth of the 2nd photon in the frame with e- along the z-axis
C
        PHIG=TWOPI*RVEC(5)
        SINTEG=SQRT(1.d0-COSTEG**2)
C
C  The photon momentum 
C
        QPP(4)=EG
        QPP(3)=PG*COSTEG
        QPP(2)=PG*SINTEG*SIN(PHIG)
        QPP(1)=PG*SINTEG*COS(PHIG)
C
C  Rotate to the electron direction
C
        CALL NEWDIR(QPP,COST,PHI)
C
C  The last photon momentum follows from energy/mom conservation 
C
        QP(4)=ROOTS-PP(4)-QPP(4)
        QP(3)=-PP(3)-QPP(3)
        QP(2)=-PP(2)-QPP(2)
        QP(1)=-PP(1)-QPP(1)
C
C  Boost both photons back to the electron rest frame to check energies
C
        EEPP=GAME*QPP(4)-GAMBE*QPP(3)
        EEP=GAME*QP(4)-GAMBE*QP(3)
        IF(EEPP.LT.KGMIN) GO TO 1000
        IF(EEP.LT.KGMIN) GO TO 1000
C
C  Keep track of good event trials
C
        NGOOD=NGOOD+1
C
C  Useful variables 
C
C
C  Define the matrix element for I=1,2 initial gamma pol states, J,K=1,2
C  final gamma pol states and L=1,2 final electron pol states
C
C  Use only the second set of auxilliary momenta (which work for any inital
C  spin state), ISET=1 is a useful diagnostic
C
        DO ISET=2,2
          IF(ISET.EQ.1) THEN
            LSET1=.TRUE.
          ELSE
            LSET1=.FALSE.
          ENDIF
          DO L=1,2
C
            MAT(1,1,1,L,ISET)=DPPP(MQ,QP,QPP,L)+DPPP(MQ,QPP,QP,L)
     >                      +DPPP(QP,MQ,QPP,L)+DPPP(QPP,MQ,QP,L)
     >                      +DPPP(QP,QPP,MQ,L)+DPPP(QPP,QP,MQ,L)
C
            MAT(1,1,2,L,ISET)=DPPM(MQ,QP,QPP,L)+DPMP(MQ,QPP,QP,L)
     >                       +DPPM(QP,MQ,QPP,L)+DPMP(QP,QPP,MQ,L)
     >                       +DMPP(QPP,MQ,QP,L)+DMPP(QPP,QP,MQ,L)
C
            MAT(1,2,1,L,ISET)=DPMP(MQ,QP,QPP,L)+DPPM(MQ,QPP,QP,L)
     >                       +DMPP(QP,MQ,QPP,L)+DMPP(QP,QPP,MQ,L)
     >                       +DPPM(QPP,MQ,QP,L)+DPMP(QPP,QP,MQ,L)
C
            MAT(1,2,2,L,ISET)=DPMM(MQ,QP,QPP,L)+DPMM(MQ,QPP,QP,L)
     >                       +DMPM(QP,MQ,QPP,L)+DMMP(QP,QPP,MQ,L)
     >                       +DMPM(QPP,MQ,QP,L)+DMMP(QPP,QP,MQ,L)
C
            MAT(2,1,1,L,ISET)=DMPP(MQ,QP,QPP,L)+DMPP(MQ,QPP,QP,L)
     >                       +DPMP(QP,MQ,QPP,L)+DPPM(QP,QPP,MQ,L)
     >                       +DPMP(QPP,MQ,QP,L)+DPPM(QPP,QP,MQ,L)
C
            MAT(2,1,2,L,ISET)=DMPM(MQ,QP,QPP,L)+DMMP(MQ,QPP,QP,L)
     >                       +DPMM(QP,MQ,QPP,L)+DPMM(QP,QPP,MQ,L)
     >                       +DMMP(QPP,MQ,QP,L)+DMPM(QPP,QP,MQ,L)
C
            MAT(2,2,1,L,ISET)=DMMP(MQ,QP,QPP,L)+DMPM(MQ,QPP,QP,L)
     >                       +DMMP(QP,MQ,QPP,L)+DMPM(QP,QPP,MQ,L)
     >                       +DPMM(QPP,MQ,QP,L)+DPMM(QPP,QP,MQ,L)
C
            MAT(2,2,2,L,ISET)=DMMM(MQ,QP,QPP,L)+DMMM(MQ,QPP,QP,L)
     >                       +DMMM(QP,MQ,QPP,L)+DMMM(QP,QPP,MQ,L)
     >                       +DMMM(QPP,MQ,QP,L)+DMMM(QPP,QP,MQ,L)
C
          ENDDO
        ENDDO
C
C  Convert to cross sections 
C
C  [WRIGHT is the weight for a right-handed helicity (Left Circular Pol) photon
C  impinging on the electron and WLEFT is for a left-handed helicity (RCP)
C  photon]
C
        DO ISET=2,2
          WLEFT(ISET)=0.d0
          WRIGHT(ISET)=0.d0
          DO J=1,2
            DO K=1,2
              DO L=1,2
                WRIGHT(ISET)=WRIGHT(ISET)+ABS(MAT(1,J,K,L,ISET))**2
                WLEFT(ISET)=WLEFT(ISET)+ABS(MAT(2,J,K,L,ISET))**2
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        ISET=2
C
C  To convert to x-section, divide by generated energy densities
C  and multiply by constants 
C
        TNORM=CNORM*(EE2-EE)/CE/(CG*(1.d0/EG+1.d0/(GG2-EG)))
        WL=WLEFT(ISET)*TNORM
        WR=WRIGHT(ISET)*TNORM
C
C  Event weight 3 is the unpolarized cross section (the average
C  (of the left- and right-handed cross sections)
C
        WGT(3)=(WL+WR)/2.d0
C
C  Event weight 4 is the polarized cross section and is defined as the
C  the difference in cross section for left-handed (negative helicity) 
C  and right-handed (positive helicity) photons impinging (traveling in the
C  -z direction) on the electron (traveling in the +z direction).  If the
C  electron spin direction is chosen to lie along the +z direction, this
C  weight is the difference in cross sections, sigma(Jz=3/2)-sigma(Jz=1/2).
C
        WGT(4)=(WL-WR)/2.d0
        XS23=XS23+WGT(3)
        XS232=XS232+WGT(3)**2
C
C  Boost back to the lab frame and accumulate event weights
C
        CALL LORENZ(1,PLAB,1,1,PP)
        CALL LORENZ(1,PLAB,1,2,QOUT)
        CALL WGTHST(1,1,PP,2,QOUT,0,PDUM,WGT)
        print *, QOUT, WGT,PP
C        print *,PP,WGT      !check qiao
1000    CONTINUE
      ENDDO
      PRINT 2000, NGOOD,NTRIAL,XS23,SQRT(XS232)
2000  FORMAT(/,2x,'***** COMEEG Finished *****',/,1x,I7,'/',I7,
     >' good e+g=>e+2g trials, total xs = ',E12.5,'+-',E12.5)
      RETURN
      END
C
      SUBROUTINE COMEEE
C
C **********************************************************************
C * This routine calculates the laboratory cross sections for E+G=>3E  *
C * Compton scattering (order alpha3) using the techniques given in    *
C * Stuart and Gongora, Z.Phys. C42,617 (1989).                        *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  The block /CONTRL/ contains the adjustable parameters needed for the
C  calculation and some useful physical constants: EB - e- beam energy in
C  GeV, EPHOT - photon energy in GeV, XME - electron mass (in GeV),
C  XMG - the fictitious photon mass (GeV), KGMIN - the minimum resolvable
C  photon energy (GeV), ALPHA - the fine structure constant, PI - 3.1416...
C  ROOT2 - SQRT(2), BARN - (hbar*c)**2 in units of mb-GeV**(-2),
C  SPIN(3) - the initial state electron spin (in the me rest frame),
C  LDIAG - logical flag to turn on diagnostic T-DR-M unpol xsection, 
C  LBF - logical flag to turn on diagnostic Brown+Feynman unpol xsection, 
C  NTRY - the number of event trials to generate with COMTN2,COMEGG,COMEEE
C         (COMEEE generates NTRY/20)
C
      COMMON /CONTRL/ EB,EPHOT,XME,XMG,KGMIN,ALPHA,PI,ROOT2,BARN,
     >                SPIN(3),LDIAG,LBF,NTRY   
      LOGICAL LDIAG,LBF
      REAL*8 KGMIN
      REAL*8 P(4),S(4,2),Q(4),PP(4),SP(4,2),PPP(4),SPP(4,2)
     >      ,PB(4),SB(4,2),MP(4),POUT(4,2),QH(4)
      EQUIVALENCE (PP(1),POUT(1,1)), (PPP(1),POUT(1,2))
      REAL*8 WLEFT(2),WRIGHT(2),PLAB(4),WGT(4)
      COMPLEX*16 MAT(2,2,2,2,2)
      COMPLEX*16 D1P,D1M,D2P,D2M
C
C  RVEC contains the random numbers
C
      REAL*4 RVEC(5)
C
C  The center-of-mass energy**2 of the e-gamma system
C
      XME2=XME**2
      ARG=EB**2-XME2
      IF(ARG.LE.0.d0) ARG=0.d0
      PBEAM=SQRT(ARG)
      S0=XME2+2.D0*(EB+PBEAM)*EPHOT
      ROOTS=SQRT(S0)
C
C  Return if there isn't enough energy to create the e-e+e- final state
C
        print *, ROOTS, 3.d0* XME
      IF(ROOTS.LT.3.d0*XME) RETURN
      SXME2=S0+XME2
      TROOTS=2.d0*ROOTS
C
C  The transformation to the lab frame
C
      PLAB(3)=(PBEAM-EPHOT)
      PLAB(2)=0.d0
      PLAB(1)=0.d0
      PLAB(4)=SQRT(S0+PLAB(3)**2)
C
C  Define the initial state 4-vectors in the cm system
C
C  The electron momentum and spin
C
      EE=(S0+XME2)/(2.D0*SQRT(S0))
      PE=SQRT(EE**2-XME2)
      P(4)=EE
      P(3)=PE
      P(2)=0.D0
      P(1)=0.D0
C
C  MP is the negative of P
C
      DO J=1,4
        MP(J)=-P(J)
      ENDDO
C
C  Construct the spin vector
C
      S(4,1)=PE/XME*SPIN(3)
      S(3,1)=EE/XME*SPIN(3)
      S(2,1)=SPIN(2)
      S(1,1)=SPIN(1)
C
C  Generate the spin-flipped vector
C
      DO J=1,4
        S(J,2)=-S(J,1)
      ENDDO
C
C  The photon momentum
C
      EG=(S0-XME2)/(2.D0*SQRT(S0))
      Q(4)=EG
      Q(3)=-EG
      Q(2)=0.D0
      Q(1)=0.D0
C
C  Useful variables 
C
      ALPHA2=ALPHA**2
      TWOPI=2.d0*PI
      DCONST=2.d0*ROOT2*(4.d0*PI*ALPHA)**1.5d0
      PEMAX=SQRT((S0-9.d0*XME2)*(S0-XME2))/(2.d0*SQRT(S0))
      EEMAX=SQRT(PEMAX**2+XME2)
      EEMIN=XME
      DE=EEMAX-EEMIN
C
C  Calculate the number of trials
C
      NTRIAL=NTRY/20
C
C  Calculate the phase space per trial
C
      PHASPC=(TWOPI)**2*2.d0*DE**2/FLOAT(NTRIAL)
C
C  Calculate the cross section normalization (include ID particles factor)
C
      CNORM=PHASPC*BARN/(8.d0*TWOPI**5*4.d0*(EB+PBEAM)*EPHOT)/2.d0
C
C  This program doesn't generate the tree-level cross section
C
      WGT(1)=0.d0
      WGT(2)=0.d0
C
C  Keep track of successful trials
C
      NGOOD=0
      XS23=0.d0
      XS232=0.d0
C
C  Loop over many event trials
C
      DO IEPHOT=1,NTRIAL
C
C  Choose 5 random numbers
C
C        CALL RANMAR(RVEC,5)
        CALL RANDOM_NUMBER(RVEC(1))
        CALL RANDOM_NUMBER(RVEC(2)) 
        CALL RANDOM_NUMBER(RVEC(3))
        CALL RANDOM_NUMBER(RVEC(4))
        CALL RANDOM_NUMBER(RVEC(5))
C
C  Generate energies for the identical electrons
C
        EE=EEMIN+DE*RVEC(1)
        EEP=EEMIN+DE*RVEC(2)
C
C  Determine if these energies can balance energy/momentum
C
        PE=SQRT(EE**2-XME2)
        PEP=SQRT(EEP**2-XME2)
        IF(PE.EQ.0.d0.OR.PEP.EQ.0.d0) GO TO 1000
        COSTEE=(SXME2+2.d0*EE*EEP-TROOTS*(EE+EEP))/(2.d0*PE*PEP)
        IF(ABS(COSTEE).GT.1.d0) GO TO 1000
C
C  Generate direction of the scattered electron
C
        COST=2.d0*RVEC(3)-1.d0
        SINT=SQRT(1.d0-COST**2)
        PHI=TWOPI*RVEC(4)
        COSY=SINT*SIN(PHI)
        COSX=SINT*COS(PHI)
C
C  Define the final state 4-vectors in the cm system
C
C  The electron momentum and spin
C
        PP(4)=EE
        PP(3)=PE*COST
        PP(2)=PE*COSY
        PP(1)=PE*COSX
C
C  There are two final state electron spins, SPR and -SPR
C
        SP(4,1)=PE/XME
        SP(3,1)=EE/XME*COST
        SP(2,1)=EE/XME*COSY
        SP(1,1)=EE/XME*COSX
C
C  Generate the spin-flipped vector
C
        DO J=1,4
          SP(J,2)=-SP(J,1)
C
C  Choose a random definition of the photon auxilliary vector (P1P)
C
          QH(J)=(PP(J)+XME*SP(J,1))/2.d0
        ENDDO
C
C  Generate azimuth of the 2nd e- in the frame with 1st e- along the z-axis
C
        PHIP=TWOPI*RVEC(5)
        SINTEE=SQRT(1.d0-COSTEE**2)
C
C  The momentum four-vector 
C
        PPP(4)=EEP
        PPP(3)=PEP*COSTEE
        PPP(2)=PEP*SINTEE*SIN(PHIP)
        PPP(1)=PEP*SINTEE*COS(PHIP)
C
C  Rotate to the electron direction
C
        CALL NEWDIR(PPP,COST,PHI)
C
C  The positron momentum follows from energy/mom conservation 
C
        PB(4)=ROOTS-PP(4)-PPP(4)
        PB(3)=-PP(3)-PPP(3)
        PB(2)=-PP(2)-PPP(2)
        PB(1)=-PP(1)-PPP(1)
C
C  Construct spin vectors for both electron and positron
C
C  first the electron
C
        SPP(4,1)=PEP/XME
        SPP(3,1)=EEP/XME*PPP(3)/PEP
        SPP(2,1)=EEP/XME*PPP(2)/PEP
        SPP(1,1)=EEP/XME*PPP(1)/PEP
C
C  and finally the positron
C
        EPOS=PB(4)
        PPOS=SQRT(PB(1)**2+PB(2)**2+PB(3)**2)
        SB(4,1)=PPOS/XME
        SB(3,1)=EPOS/XME*PB(3)/PPOS
        SB(2,1)=EPOS/XME*PB(2)/PPOS
        SB(1,1)=EPOS/XME*PB(1)/PPOS
C
C  Generate the spin-flipped vectors
C
        DO J=1,4
          SPP(J,2)=-SPP(J,1)
          SB(J,2)=-SB(J,1)
        ENDDO
C
C  Keep track of good event trials
C
        NGOOD=NGOOD+1

        print *, NGOOD
C
C  Useful variables 
C
C
C  Define the matrix element MAT(I,J,K,L,M) for I=1,2 initial gamma pol
C  states, J=1,2 initial electron spin states, K,L=1,2 final electron pol
C  states and M=1,2 final positron pol states
C
        DO J=1,1
          DO K=1,2
            DO L=1,2
              DO M=1,2
C
                MAT(1,J,K,L,M)=D1P(MP,S(1,J),PB,SB(1,M),
     >                             PP,SP(1,K),PPP,SPP(1,L),Q,QH)
     >                        -D1P(MP,S(1,J),PB,SB(1,M),
     >                             PPP,SPP(1,L),PP,SP(1,K),Q,QH)
     >                        -D1P(PB,SB(1,M),MP,S(1,J),
     >                             PP,SP(1,K),PPP,SPP(1,L),Q,QH)
     >                        +D1P(PB,SB(1,M),MP,S(1,J),
     >                             PPP,SPP(1,L),PP,SP(1,K),Q,QH)
     >                        +D2P(MP,S(1,J),PB,SB(1,M),
     >                             PP,SP(1,K),PPP,SPP(1,L),Q,QH)
     >                        -D2P(MP,S(1,J),PB,SB(1,M),
     >                             PPP,SPP(1,L),PP,SP(1,K),Q,QH)
     >                        -D2P(PB,SB(1,M),MP,S(1,J),
     >                             PP,SP(1,K),PPP,SPP(1,L),Q,QH)
     >                        +D2P(PB,SB(1,M),MP,S(1,J),
     >                             PPP,SPP(1,L),PP,SP(1,K),Q,QH)
C
                MAT(2,J,K,L,M)=D1M(MP,S(1,J),PB,SB(1,M),
     >                             PP,SP(1,K),PPP,SPP(1,L),Q,QH)
     >                        -D1M(MP,S(1,J),PB,SB(1,M),
     >                             PPP,SPP(1,L),PP,SP(1,K),Q,QH)
     >                        -D1M(PB,SB(1,M),MP,S(1,J),
     >                             PP,SP(1,K),PPP,SPP(1,L),Q,QH)
     >                        +D1M(PB,SB(1,M),MP,S(1,J),
     >                             PPP,SPP(1,L),PP,SP(1,K),Q,QH)
     >                        +D2M(MP,S(1,J),PB,SB(1,M),
     >                             PP,SP(1,K),PPP,SPP(1,L),Q,QH)
     >                        -D2M(MP,S(1,J),PB,SB(1,M),
     >                             PPP,SPP(1,L),PP,SP(1,K),Q,QH)
     >                        -D2M(PB,SB(1,M),MP,S(1,J),
     >                             PP,SP(1,K),PPP,SPP(1,L),Q,QH)
     >                        +D2M(PB,SB(1,M),MP,S(1,J),
     >                             PPP,SPP(1,L),PP,SP(1,K),Q,QH)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
C
C  Convert to cross sections
C
C  [WRIGHT is the weight for a right-handed helicity (Left Circular Pol) photon
C  impinging on the electron and WLEFT is for a left-handed helicity (RCP)
C  photon]
C
        DO I=1,1
          WLEFT(I)=0.d0
          WRIGHT(I)=0.d0
          DO K=1,2
            DO L=1,2
              DO M=1,2
                WRIGHT(I)=WRIGHT(I)+ABS(MAT(I,I,K,L,M))**2
                WLEFT(I)=WLEFT(I)+ABS(MAT(3-I,I,K,L,M))**2
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        ISET=1
C
C  To convert to x-section, divide by generated energy densities
C  and multiply by constants 
C
        WL=WLEFT(ISET)*CNORM
        WR=WRIGHT(ISET)*CNORM
C
C  Event weight 3 is the unpolarized cross section (the average
C  (of the left- and right-handed cross sections)
C
        WGT(3)=(WL+WR)/2.d0
C
C  Event weight 4 is the polarized cross section and is defined as the
C  the difference in cross section for left-handed (negative helicity) 
C  and right-handed (positive helicity) photons impinging (traveling in the
C  -z direction) on the electron (traveling in the +z direction).  If the
C  electron spin direction is chosen to lie along the +z direction, this
C  weight is the difference in cross sections, sigma(Jz=3/2)-sigma(Jz=1/2).
C
        WGT(4)=(WL-WR)/2.d0
        XS23=XS23+WGT(3)
        XS232=XS232+WGT(3)**2
C
C  Boost back to the lab frame and accumulate event weights
C
        CALL LORENZ(1,PLAB,1,2,POUT)
        CALL LORENZ(1,PLAB,1,1,PB)
        CALL WGTHST(1,2,POUT,0,QOUT,1,PB,WGT)
1000    CONTINUE
      ENDDO
      PRINT 2000, NGOOD,NTRIAL,XS23,SQRT(XS232)
2000  FORMAT(/,2x,'***** COMEEE Finished *****',/,1x,I7,'/',I7,
     >' good e+g=>3e trials, total xs = ',E12.5,'+-',E12.5)
      RETURN
      END
C
      FUNCTION DPPP(Q,QP,QPP,L)
C
C **********************************************************************
C * This function implements the amplitude D+++ as defined in Stuart   *
C * and Gongora, Z.Phys. C42,617 (1989).                               *
C *   Parameters:   Q(4) - REAL*8 4-vector                             *
C *                QP(4) - REAL*8 4-vector                             *
C *               QPP(4) - REAL*8 4-vector                             *
C *                    L - Integer index for e- final state spin       *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DPPP,SP,SM
      COMPLEX*16 ANS
      REAL*8 Q(4),QP(4),QPP(4)
      INTEGER L
C
C  The block /VECT/ contains four-vectors needed by the DXXX functions
C
      COMMON /VECT/ DCONST,XME2,P(4),PP(4),P1(4),P2(4),P1P(4,2),P2P(4,2)
     >             ,LSET1
      LOGICAL LSET1
C
C  Calculate amplitudes for either set of auxiliary momenta
C
      IF(LSET1) THEN
        ANS=-SP(P1,P2)*SP(P2,P1P(1,L))*SM(P1,Q)*SM(P1P(1,L),P2P(1,L))        
     >    /(XME2*SP(P2,Q)*SP(P2,QP)*SP(P2,QPP))
     >    *(SM(QPP,P1P(1,L))*SP(P1P(1,L),P2)
     >      +SM(QPP,P2P(1,L))*SP(P2P(1,L),P2))
     >    *(SM(QP,P1)*SP(P1,P2)-SM(QP,Q)*SP(Q,P2))
      ELSE
        ANS=-SP(P2,P1P(1,L))*SM(P2P(1,L),QPP)
     >    /(SP(P1P(1,L),Q)*SP(P1P(1,L),QP)*SP(P1P(1,L),QPP))
     >    *(SP(P1P(1,L),P2P(1,L))*SM(P2P(1,L),QP)
     >      +SP(P1P(1,L),QPP)*SM(QPP,QP))
     >    *(SP(P1P(1,L),P1)*SM(P1,Q)+SP(P1P(1,L),P2)*SM(P2,Q))
      ENDIF
C
C  Calculate the normalization
C
      PQ=P(4)*Q(4)-P(3)*Q(3)-P(2)*Q(2)-P(1)*Q(1)
      PPQPP=PP(4)*QPP(4)-PP(3)*QPP(3)-PP(2)*QPP(2)-PP(1)*QPP(1)
      ANS=ANS*DCONST/(-4.d0*PQ*PPQPP)
C
C  Return the answer
C
      DPPP=ANS
      RETURN
      END
C
      FUNCTION DMMM(Q,QP,QPP,L)
C
C **********************************************************************
C * This function implements the amplitude D--- as defined in Stuart   *
C * and Gongora, Z.Phys. C42,617 (1989).                               *
C *   Parameters:   Q(4) - REAL*8 4-vector                             *
C *                QP(4) - REAL*8 4-vector                             *
C *               QPP(4) - REAL*8 4-vector                             *
C *                    L - Integer index for e- final state spin       *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DMMM,SP,SM
      COMPLEX*16 ANS
      REAL*8 Q(4),QP(4),QPP(4)
      INTEGER L
C
C  The block /VECT/ contains four-vectors needed by the DXXX functions
C
      COMMON /VECT/ DCONST,XME2,P(4),PP(4),P1(4),P2(4),P1P(4,2),P2P(4,2)
     >             ,LSET1
      LOGICAL LSET1
C
C  Calculate amplitudes for either set of auxiliary momenta
C
      IF(LSET1) THEN
        ANS=-SP(P2,Q)*SM(P1,P2P(1,L))
     >    /(SM(P1,Q)*SM(P1,QP)*SM(P1,QPP))
     >    *(SP(QPP,P1P(1,L))*SM(P1P(1,L),P1)
     >      +SP(QPP,P2P(1,L))*SM(P2P(1,L),P1))
     >    *(SP(QP,P2)*SM(P2,P1)-SP(QP,Q)*SM(Q,P1))
      ELSE
        ANS=-SP(P1,P2)*SP(P1P(1,L),QPP)*SM(P1,P2P(1,L))
     >    *SM(P1P(1,L),P2P(1,L))
     >    /(XME2*SM(P2P(1,L),Q)*SM(P2P(1,L),QP)*SM(P2P(1,L),QPP))
     >    *(SM(P2P(1,L),P1P(1,L))*SP(P1P(1,L),QP)
     >      +SM(P2P(1,L),QPP)*SP(QPP,QP))
     >    *(SM(P2P(1,L),P1)*SP(P1,Q)+SM(P2P(1,L),P2)*SP(P2,Q))
      ENDIF
C
C  Calculate the normalization
C
      PQ=P(4)*Q(4)-P(3)*Q(3)-P(2)*Q(2)-P(1)*Q(1)
      PPQPP=PP(4)*QPP(4)-PP(3)*QPP(3)-PP(2)*QPP(2)-PP(1)*QPP(1)
      ANS=ANS*DCONST/(-4.d0*PQ*PPQPP)
C
C  Return the answer
C
      DMMM=ANS
      RETURN
      END
C
      FUNCTION DPPM(Q,QP,QPP,L)
C
C **********************************************************************
C * This function implements the amplitude D++- as defined in Stuart   *
C * and Gongora, Z.Phys. C42,617 (1989).                               *
C *   Parameters:   Q(4) - REAL*8 4-vector                             *
C *                QP(4) - REAL*8 4-vector                             *
C *               QPP(4) - REAL*8 4-vector                             *
C *                    L - Integer index for e- final state spin       *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DPPM,SP,SM
      COMPLEX*16 ANS,AUX1,AUX2,AUX4,AUX5,AUX6
      REAL*8 Q(4),QP(4),QPP(4)
      INTEGER L
C
C  The block /VECT/ contains four-vectors needed by the DXXX functions
C
      COMMON /VECT/ DCONST,XME2,P(4),PP(4),P1(4),P2(4),P1P(4,2),P2P(4,2)
     >             ,LSET1
      LOGICAL LSET1
C
C  Calculate amplitudes for either set of auxiliary momenta
C
      IF(LSET1) THEN
        AUX1=-(SM(P1,Q)*SP(Q,P2)+SM(P1,QP)*SP(QP,P2))
        AUX2=SM(QP,P1)*SP(P1,P2)-SM(QP,Q)*SP(Q,P2)
        ANS=-SP(P1,P2)*SM(P1,Q)
     >    /(SP(P2,Q)*SP(P2,QP)*SM(P1,QPP))
     >    *(SP(P1P(1,L),QPP)*SM(P1P(1,L),P2P(1,L))*AUX1*AUX2/XME2
     >      -SP(P2,QPP)*SM(P1,P2P(1,L))*AUX2)
      ELSE
        AUX4=SM(P2P(1,L),QPP)*SP(QPP,P1P(1,L))
        AUX5=SM(QP,P2P(1,L))*SP(P2P(1,l),P1P(1,L))
     >    +SM(QP,QPP)*SP(QPP,P1P(1,L))
        AUX6=SP(P1P(1,L),P1)*SM(P1,Q)+SP(P1P(1,L),P2)*SM(P2,Q)
        ANS=-SP(P1P(1,L),QPP)*SM(P1P(1,L),P2P(1,L))
     >    /(SP(P1P(1,L),Q)*SP(P1P(1,L),QP)*SM(P2P(1,L),QPP))
     >    *(SP(P1,P2)*SM(P1,Q)*AUX4*AUX5/XME2
     >      -SP(P2,P1P(1,L))*SM(Q,QP)*AUX4
     >      +SP(P2,P1P(1,L))*SM(P2P(1,L),QP)*AUX6)
      ENDIF
C
C  Calculate the normalization
C
      PQ=P(4)*Q(4)-P(3)*Q(3)-P(2)*Q(2)-P(1)*Q(1)
      PPQPP=PP(4)*QPP(4)-PP(3)*QPP(3)-PP(2)*QPP(2)-PP(1)*QPP(1)
      ANS=ANS*DCONST/(-4.d0*PQ*PPQPP)
C
C  Return the answer
C
      DPPM=ANS
      RETURN
      END
C
      FUNCTION DPMP(Q,QP,QPP,L)
C
C **********************************************************************
C * This function implements the amplitude D+-+ as defined in Stuart   *
C * and Gongora, Z.Phys. C42,617 (1989).                               *
C *   Parameters:   Q(4) - REAL*8 4-vector                             *
C *                QP(4) - REAL*8 4-vector                             *
C *               QPP(4) - REAL*8 4-vector                             *
C *                    L - Integer index for e- final state spin       *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DPMP,SP,SM
      COMPLEX*16 ANS,AUX1,AUX2,AUX3,AUX4,AUX5,AUX6
      REAL*8 Q(4),QP(4),QPP(4)
      INTEGER L
C
C  The block /VECT/ contains four-vectors needed by the DXXX functions
C
      COMMON /VECT/ DCONST,XME2,P(4),PP(4),P1(4),P2(4),P1P(4,2),P2P(4,2)
     >             ,LSET1
      LOGICAL LSET1
C
C  Calculate amplitudes for either set of auxiliary momenta
C
      IF(LSET1) THEN
        AUX1=SM(QPP,P1P(1,L))*SP(P1P(1,L),QP)
     >      +SM(QPP,P2P(1,L))*SP(P2P(1,L),QP)
        AUX2=-(SP(P2,Q)*SM(Q,P1)+SP(P2,QP)*SM(QP,P1))
        AUX3=-SM(P1,Q)*SP(Q,P2)
        ANS=SP(P1,P2)*SM(P1,Q)
     >    /(SP(P2,Q)*SM(P1,QP)*SP(P2,QPP))
     >    *(SP(P2,P1P(1,L))*SM(P1P(1,L),P2P(1,L))*AUX1*AUX3/XME2
     >     +SP(P2,QP)*SM(P2P(1,L),QPP)*(AUX3-AUX2)
     >     +SP(P2,QP)*SP(P2,P1P(1,L))*SM(P1,QPP)*SM(P1P(1,L),P2P(1,L)))
      ELSE
        AUX4=SP(P1P(1,L),QPP)*SM(QPP,P2P(1,L))
        AUX5=SP(QP,P1)*SM(P1,Q)+SP(QP,P2)*SM(P2,Q)
        AUX6=SM(P2P(1,L),QP)*SP(QP,P1P(1,L))
     >    +SM(P2P(1,L),QPP)*SP(QPP,P1P(1,L))
        ANS=SM(P2P(1,L),QPP)
     >    /(SP(P1P(1,L),Q)*SM(P2P(1,L),QP)*SP(P1P(1,L),QPP))
     >    *(SP(P2,P1P(1,L))*AUX4*AUX5
     >      +SP(P1,P2)*SP(P1P(1,L),QP)*SM(P1,Q)*(AUX6-AUX4)
     >      +XME2*SP(P2,P1P(1,L))*SP(P1P(1,L),QP)*SM(P2P(1,L),Q))
      ENDIF
C
C  Calculate the normalization
C
      PQ=P(4)*Q(4)-P(3)*Q(3)-P(2)*Q(2)-P(1)*Q(1)
      PPQPP=PP(4)*QPP(4)-PP(3)*QPP(3)-PP(2)*QPP(2)-PP(1)*QPP(1)
      ANS=ANS*DCONST/(-4.d0*PQ*PPQPP)
C
C  Return the answer
C
      DPMP=ANS
      RETURN
      END
C
      FUNCTION DMPP(Q,QP,QPP,L)
C
C **********************************************************************
C * This function implements the amplitude D-++ as defined in Stuart   *
C * and Gongora, Z.Phys. C42,617 (1989).                               *
C *   Parameters:   Q(4) - REAL*8 4-vector                             *
C *                QP(4) - REAL*8 4-vector                             *
C *               QPP(4) - REAL*8 4-vector                             *
C *                    L - Integer index for e- final state spin       *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DMPP,SP,SM
      COMPLEX*16 ANS,AUX1,AUX2,AUX3,AUX4,AUX5
      REAL*8 Q(4),QP(4),QPP(4)
      INTEGER L
C
C  The block /VECT/ contains four-vectors needed by the DXXX functions
C
      COMMON /VECT/ DCONST,XME2,P(4),PP(4),P1(4),P2(4),P1P(4,2),P2P(4,2)
     >             ,LSET1
      LOGICAL LSET1
C
C  Calculate amplitudes for either set of auxiliary momenta
C
      IF(LSET1) THEN
        AUX1=SP(P2,P1)*SM(P1,QP)-SP(P2,Q)*SM(Q,QP)
        AUX2=SM(QPP,P1P(1,L))*SP(P1P(1,L),P2)
     >     +SM(QPP,P2P(1,L))*SP(P2P(1,L),P2)
        AUX3=-SP(P2,Q)*SM(Q,P1)
        ANS=SP(P2,Q)
     >    /(SM(P1,Q)*SP(P2,QP)*SP(P2,QPP))
     >    *(SM(P2P(1,L),QPP)*AUX1*AUX3
     >     -SP(P2,P1P(1,L))*SM(P1,QP)*SM(P1P(1,L),P2P(1,L))*AUX2
     >     -SP(P2,P1P(1,L))*SM(QP,QPP)*SM(P1P(1,L),P2P(1,L))*AUX3)
      ELSE
        AUX4=SP(P1P(1,L),P2P(1,L))*SM(P2P(1,L),QP)
     >    +SP(P1P(1,L),QPP)*SM(QPP,QP)
        AUX5=SP(P1P(1,L),QP)*SM(QP,P2P(1,L))
     >    +SP(P1P(1,l),QPP)*SM(QPP,P2P(1,L))
        ANS=SM(P2P(1,L),QPP)
     >    /(SM(P2P(1,L),Q)*SP(P1P(1,L),QP)*SP(P1P(1,L),QPP))
     >    *(SP(P2,Q)*AUX4*AUX5
     >      +SP(P1,P2)*SP(P1P(1,L),Q)*SM(P1,P2P(1,L))*AUX4)
      ENDIF
C
C  Calculate the normalization
C
      PQ=P(4)*Q(4)-P(3)*Q(3)-P(2)*Q(2)-P(1)*Q(1)
      PPQPP=PP(4)*QPP(4)-PP(3)*QPP(3)-PP(2)*QPP(2)-PP(1)*QPP(1)
      ANS=ANS*DCONST/(-4.d0*PQ*PPQPP)
C
C  Return the answer
C
      DMPP=ANS
      RETURN
      END
C
      FUNCTION DMMP(Q,QP,QPP,L)
C
C **********************************************************************
C * This function implements the amplitude D--+ as defined in Stuart   *
C * and Gongora, Z.Phys. C42,617 (1989).                               *
C *   Parameters:   Q(4) - REAL*8 4-vector                             *
C *                QP(4) - REAL*8 4-vector                             *
C *               QPP(4) - REAL*8 4-vector                             *
C *                    L - Integer index for e- final state spin       *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DMMP,SP,SM
      COMPLEX*16 ANS,AUX1,AUX2,AUX4,AUX5,AUX6
      REAL*8 Q(4),QP(4),QPP(4)
      INTEGER L
C
C  The block /VECT/ contains four-vectors needed by the DXXX functions
C
      COMMON /VECT/ DCONST,XME2,P(4),PP(4),P1(4),P2(4),P1P(4,2),P2P(4,2)
     >             ,LSET1
      LOGICAL LSET1
C
C  Calculate amplitudes for either set of auxiliary momenta
C
      IF(LSET1) THEN
        AUX1=-(SP(P2,Q)*SM(Q,P1)+SP(P2,QP)*SM(QP,P1))
        AUX2=SP(QP,P2)*SM(P2,P1)-SP(QP,Q)*SM(Q,P1)
        ANS=-SP(P2,Q)
     >    /(SM(P1,Q)*SM(P1,QP)*SP(P2,QPP))
     >   *(SM(P2P(1,L),QPP)*AUX1*AUX2
     >   -SP(P2,P1P(1,L))*SM(P1,QPP)*SM(P1P(1,L),P2P(1,L))*AUX2) 
      ELSE  
        AUX4=SP(P1P(1,L),QPP)*SM(QPP,P2P(1,L))
        AUX5=SP(QP,P1P(1,L))*SM(P1P(1,L),P2P(1,L))
     >      +SP(QP,QPP)*SM(QPP,P2P(1,L))
        AUX6=SM(P2P(1,L),P1)*SP(P1,Q)+SM(P2P(1,L),P2)*SP(P2,Q)
        ANS=-SM(P2P(1,L),QPP)
     >    /(SM(P2P(1,L),Q)*SM(P2P(1,L),QP)*SP(P1P(1,L),QPP))
     >    *(SP(P2,Q)*AUX4*AUX5
     >      -SP(P1,P2)*SP(Q,QP)*SM(P1,P2P(1,L))*AUX4
C     >      +SP(P1,P2)*SP(Q,QP)*SM(P1,P2P(1,L))*AUX4
     >      +SP(P1,P2)*SP(P1P(1,L),QP)*SM(P1,P2P(1,L))*AUX6)
      ENDIF      
C
C  Calculate the normalization
C
      PQ=P(4)*Q(4)-P(3)*Q(3)-P(2)*Q(2)-P(1)*Q(1)
      PPQPP=PP(4)*QPP(4)-PP(3)*QPP(3)-PP(2)*QPP(2)-PP(1)*QPP(1)
      ANS=ANS*DCONST/(-4.d0*PQ*PPQPP)
C
C  Return the answer
C
      DMMP=ANS
      RETURN
      END
C
      FUNCTION DMPM(Q,QP,QPP,L)
C
C **********************************************************************
C * This function implements the amplitude D-+- as defined in Stuart   *
C * and Gongora, Z.Phys. C42,617 (1989).                               *
C *   Parameters:   Q(4) - REAL*8 4-vector                             *
C *                QP(4) - REAL*8 4-vector                             *
C *               QPP(4) - REAL*8 4-vector                             *
C *                    L - Integer index for e- final state spin       *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DMPM,SP,SM
      COMPLEX*16 ANS,AUX1,AUX2,AUX3,AUX4,AUX5,AUX6
      REAL*8 Q(4),QP(4),QPP(4)
      INTEGER L
C
C  The block /VECT/ contains four-vectors needed by the DXXX functions
C
      COMMON /VECT/ DCONST,XME2,P(4),PP(4),P1(4),P2(4),P1P(4,2),P2P(4,2)
     >             ,LSET1
      LOGICAL LSET1
C
C  Calculate amplitudes for either set of auxiliary momenta
C
      IF(LSET1) THEN
        AUX1=SP(QPP,P1P(1,L))*SM(P1P(1,L),QP)
     >     +SP(QPP,P2P(1,L))*SM(P2P(1,L),QP)
        AUX2=-(SM(P1,Q)*SP(Q,P2)+SM(P1,QP)*SP(QP,P2))
        AUX3=-SP(P2,Q)*SM(Q,P1)
        ANS=SP(P2,Q)
     >    /(SM(P1,Q)*SP(P2,QP)*SM(P1,QPP))
     >    *(SM(P1,P2P(1,L))*AUX1*AUX3
     >    +SP(P1P(1,L),QPP)*SM(P1,QP)*SM(P1P(1,L),P2P(1,L))*(AUX3-AUX2)         
     >    +XME2*SP(P2,QPP)*SM(P1,QP)*SM(P1,P2P(1,L)))
      ELSE
        AUX4=SM(P2P(1,L),QPP)*SP(QPP,P1P(1,L))
        AUX5=SM(QP,P1)*SP(P1,Q)+SM(QP,P2)*SP(P2,Q)
        AUX6=SP(P1P(1,L),QP)*SM(QP,P2P(1,L))
     >      +SP(P1P(1,L),QPP)*SM(QPP,P2P(1,L))
        ANS=SP(P1P(1,L),QPP)*SM(P1P(1,L),P2P(1,L))
     >    /(SM(P2P(1,L),Q)*SP(P1P(1,L),QP)*SM(P2P(1,L),QPP))
     >    *(SP(P1,P2)*SM(P1,P2P(1,L))*AUX4*AUX5/XME2
     >     +SP(P2,Q)*SM(P2P(1,L),QP)*(AUX6-AUX4)
     >     +SP(P1,P2)*SP(P1P(1,L),Q)*SM(P1,P2P(1,L))*SM(P2P(1,L),QP)) 
      ENDIF         
C
C  Calculate the normalization
C
      PQ=P(4)*Q(4)-P(3)*Q(3)-P(2)*Q(2)-P(1)*Q(1)
      PPQPP=PP(4)*QPP(4)-PP(3)*QPP(3)-PP(2)*QPP(2)-PP(1)*QPP(1)
      ANS=ANS*DCONST/(-4.d0*PQ*PPQPP)
C
C  Return the answer
C
      DMPM=ANS
      RETURN
      END
C
      FUNCTION DPMM(Q,QP,QPP,L)
C
C **********************************************************************
C * This function implements the amplitude D+-- as defined in Stuart   *
C * and Gongora, Z.Phys. C42,617 (1989).                               *
C *   Parameters:   Q(4) - REAL*8 4-vector                             *
C *                QP(4) - REAL*8 4-vector                             *
C *               QPP(4) - REAL*8 4-vector                             *
C *                    L - Integer index for e- final state spin       *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DPMM,SP,SM
      COMPLEX*16 ANS,AUX1,AUX2,AUX3,AUX4,AUX5
      REAL*8 Q(4),QP(4),QPP(4)
      INTEGER L
C
C  The block /VECT/ contains four-vectors needed by the DXXX functions
C
      COMMON /VECT/ DCONST,XME2,P(4),PP(4),P1(4),P2(4),P1P(4,2),P2P(4,2)
     >             ,LSET1
      LOGICAL LSET1
C
C  Calculate amplitudes for either set of auxiliary momenta
C
      IF(LSET1) THEN
        AUX1=SM(P1,P2)*SP(P2,QP)-SM(P1,Q)*SP(Q,QP)
        AUX2=SP(QPP,P1P(1,L))*SM(P1P(1,L),P1)
     >      +SP(QPP,P2P(1,L))*SM(P2P(1,L),P1)
        AUX3=-SM(P1,Q)*SP(Q,P2)
        ANS=SP(P1,P2)*SM(P1,Q)
     >    /(SP(P2,Q)*SM(P1,QP)*SM(P1,QPP))
     >    *(SP(P1P(1,L),QPP)*SM(P1P(1,L),P2P(1,L))*AUX1*AUX3/XME2
     >    -SP(P2,QP)*SM(P1,P2P(1,L))*AUX2
     >    -SP(QP,QPP)*SM(P1,P2P(1,L))*AUX3)
      ELSE
        AUX4=SM(P2P(1,L),P1P(1,L))*SP(P1P(1,L),QP)
     >      +SM(P2P(1,L),QPP)*SP(QPP,QP)
        AUX5=SM(P2P(1,L),QP)*SP(QP,P1P(1,L))
     >      +SM(P2P(1,L),QPP)*SP(QPP,P1P(1,L))
        ANS=SP(P1P(1,L),QPP)*SM(P1P(1,L),P2P(1,L))
     >    /(SP(P1P(1,L),Q)*SM(P2P(1,L),QP)*SM(P2P(1,L),QPP))
     >    *(SP(P1,P2)*SM(P1,Q)*AUX4*AUX5/XME2
     >      +SP(P2,P1P(1,L))*SM(P2P(1,L),Q)*AUX4)
      ENDIF
C
C  Calculate the normalization
C
      PQ=P(4)*Q(4)-P(3)*Q(3)-P(2)*Q(2)-P(1)*Q(1)
      PPQPP=PP(4)*QPP(4)-PP(3)*QPP(3)-PP(2)*QPP(2)-PP(1)*QPP(1)
      ANS=ANS*DCONST/(-4.d0*PQ*PPQPP)
C
C  Return the answer
C
      DPMM=ANS
      RETURN
      END
C
      FUNCTION D1P(P,S,PB,SB,PP,SPR,PPP,SPP,K,QH)
C
C **********************************************************************
C * This function implements the amplitude D1+ for the production of   *
C * the 3e final state.                                                *
C *   Parameters:   P(4) - REAL*8 mom 4-vector (initial e-/final e+)   *
C *                 S(4) - REAL*8 spin 4-vector for above              *
C *                PB(4) - REAL*8 mom 4-vector (initial e-/final e+)   *
C *                SB(4) - REAL*8 spin 4-vector for above              *
C *                PP(4) - REAL*8 mom 4-vector (final e-)              *
C *               SPR(4) - REAL*8 spin 4-vector for above              *
C *               PPP(4) - REAL*8 mom 4-vector (final e-)              *
C *               SPP(4) - REAL*8 spin 4-vector for above              *
C *                 K(4) - REAL*8 mom 4-vector (initial photon)        *
C *                QH(4) - REAL*8 auxilliary momentum 4-vector         *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 D1P,SP,SM
      COMPLEX*16 ANS
      REAL*8 P1(4),P2(4),P1P(4),P2P(4),P1PP(4),P2PP(4),P1B(4),P2B(4),
     >       Q(4)
      REAL*8 P(4),S(4),PB(4),SB(4),PP(4),SPR(4),PPP(4),SPP(4),K(4),
     >       QH(4)
C
C  The block /CONTRL/ contains the adjustable parameters needed for the
C  calculation and some useful physical constants: EB - e- beam energy in
C  GeV, EPHOT - photon energy in GeV, XME - electron mass (in GeV),
C  XMG - the fictitious photon mass (GeV), KGMIN - the minimum resolvable
C  photon energy (GeV), ALPHA - the fine structure constant, PI - 3.1416...
C  ROOT2 - SQRT(2), BARN - (hbar*c)**2 in units of mb-GeV**(-2),
C  SPIN(3) - the initial state electron spin (in the me rest frame),
C  LDIAG - logical flag to turn on diagnostic T-DR-M unpol xsection, 
C  LBF - logical flag to turn on diagnostic Brown+Feynman unpol xsection, 
C  NTRY - the number of event trials to generate with COMTN2,COMEGG,COMEEE
C         (COMEEE generates NTRY/20)
C
      COMMON /CONTRL/ EB,EPHOT,XME,XMG,KGMIN,ALPHA,PI,ROOT2,BARN,
     >                SPIN(3),LDIAG,LBF,NTRY   
      LOGICAL LDIAG,LBF
      REAL*8 KGMIN
      LOGICAL FCALL
      DATA FCALL/.TRUE./
C
C  Initialize things
C
      IF(FCALL) THEN
        DCONST=2.d0*ROOT2*(4.d0*PI*ALPHA)**1.5d0
        XME2=XME**2
        FCALL=.FALSE.
      ENDIF
C
C  Calculate the light-like momenta
C
        DO J=1,4
          P1(J)=(P(J)+XME*S(J))/2.d0
          P2(J)=(P(J)-XME*S(J))/2.d0
          P1B(J)=(PB(J)+XME*SB(J))/2.d0
          P2B(J)=(PB(J)-XME*SB(J))/2.d0
          P1P(J)=(PP(J)+XME*SPR(J))/2.d0
          P2P(J)=(PP(J)-XME*SPR(J))/2.d0
          P1PP(J)=(PPP(J)+XME*SPP(J))/2.d0
          P2PP(J)=(PPP(J)-XME*SPP(J))/2.d0
        ENDDO
C
C  Calculate the amplitude
C
      ANS=(SM(P2P,K)*SP(P1,P2)*SP(QH,P1B)*SM(P2PP,P2)
     >     +SP(QH,P1P)*SM(P1P,P2P)*SM(K,P2PP)*SP(P1B,P1)
     >     +SM(P2P,K)*SP(P1B,P1)*(SP(QH,P1P)*SM(P1P,P2PP)
     >      +SP(QH,P2P)*SM(P2P,P2PP)-SP(QH,K)*SM(K,P2PP))
     >     +SP(QH,P1P)*SM(P1P,P2P)*SP(P1,P2)*SM(P2PP,P2)/XME2
     >      *(SM(K,P1P)*SP(P1P,P1B)+SM(K,P2P)*SP(P2P,P1B))
     >    +SM(P1PP,P2PP)*SP(P2B,P1B)/XME2*(
     >     SM(P2P,K)*SP(P1,P2)*SP(QH,P1PP)*SM(P2B,P2)
     >     +SP(QH,P1P)*SM(P1P,P2P)*SM(K,P2B)*SP(P1PP,P1)
     >     +SM(P2P,K)*SP(P1PP,P1)*(SP(QH,P1P)*SM(P1P,P2B)
     >      +SP(QH,P2P)*SM(P2P,P2B)-SP(QH,K)*SM(K,P2B))
     >     +SP(QH,P1P)*SM(P1P,P2P)*SP(P1,P2)*SM(P2B,P2)/XME2
     >      *(SM(K,P1P)*SP(P1P,P1PP)+SM(K,P2P)*SP(P2P,P1PP))))/SP(QH,K)
C
C  Calculate the normalization
C
      DO I=1,4
        Q(I)=PPP(I)+PB(I)
      ENDDO
      PPK=PP(4)*K(4)-PP(3)*K(3)-PP(2)*K(2)-PP(1)*K(1)
      Q2=Q(4)**2-Q(3)**2-Q(2)**2-Q(1)**2
      ANS=ANS*DCONST/(-2.d0*PPK*Q2)
C
C  Return the answer
C
      D1P=ANS
      RETURN
      END
C
      FUNCTION D1M(P,S,PB,SB,PP,SPR,PPP,SPP,K,QH)
C
C **********************************************************************
C * This function implements the amplitude D1- for the production of   *
C * the 3e final state.                                                *
C *   Parameters:   P(4) - REAL*8 mom 4-vector (initial e-/final e+)   *
C *                 S(4) - REAL*8 spin 4-vector for above              *
C *                PB(4) - REAL*8 mom 4-vector (initial e-/final e+)   *
C *                SB(4) - REAL*8 spin 4-vector for above              *
C *                PP(4) - REAL*8 mom 4-vector (final e-)              *
C *               SPR(4) - REAL*8 spin 4-vector for above              *
C *               PPP(4) - REAL*8 mom 4-vector (final e-)              *
C *               SPP(4) - REAL*8 spin 4-vector for above              *
C *                 K(4) - REAL*8 mom 4-vector (initial photon)        *
C *                QH(4) - REAL*8 auxilliary momentum 4-vector         *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 D1M,SP,SM
      COMPLEX*16 ANS
      REAL*8 P1(4),P2(4),P1P(4),P2P(4),P1PP(4),P2PP(4),P1B(4),P2B(4),
     >       Q(4)
      REAL*8 P(4),S(4),PB(4),SB(4),PP(4),SPR(4),PPP(4),SPP(4),K(4),
     >       QH(4)
C
C  The block /CONTRL/ contains the adjustable parameters needed for the
C  calculation and some useful physical constants: EB - e- beam energy in
C  GeV, EPHOT - photon energy in GeV, XME - electron mass (in GeV),
C  XMG - the fictitious photon mass (GeV), KGMIN - the minimum resolvable
C  photon energy (GeV), ALPHA - the fine structure constant, PI - 3.1416...
C  ROOT2 - SQRT(2), BARN - (hbar*c)**2 in units of mb-GeV**(-2),
C  SPIN(3) - the initial state electron spin (in the me rest frame),
C  LDIAG - logical flag to turn on diagnostic T-DR-M unpol xsection, 
C  LBF - logical flag to turn on diagnostic Brown+Feynman unpol xsection, 
C  NTRY - the number of event trials to generate with COMTN2,COMEGG,COMEEE
C         (COMEEE generates NTRY/20)
C
      COMMON /CONTRL/ EB,EPHOT,XME,XMG,KGMIN,ALPHA,PI,ROOT2,BARN,
     >                SPIN(3),LDIAG,LBF,NTRY   
      LOGICAL LDIAG,LBF
      REAL*8 KGMIN
      LOGICAL FCALL
      DATA FCALL/.TRUE./
C
C  Initialize things
C
      IF(FCALL) THEN
        DCONST=2.d0*ROOT2*(4.d0*PI*ALPHA)**1.5d0
        XME2=XME**2
        FCALL=.FALSE.
      ENDIF
C
C  Calculate the light-like momenta
C
        DO J=1,4
          P1(J)=(P(J)+XME*S(J))/2.d0
          P2(J)=(P(J)-XME*S(J))/2.d0
          P1B(J)=(PB(J)+XME*SB(J))/2.d0
          P2B(J)=(PB(J)-XME*SB(J))/2.d0
          P1P(J)=(PP(J)+XME*SPR(J))/2.d0
          P2P(J)=(PP(J)-XME*SPR(J))/2.d0
          P1PP(J)=(PPP(J)+XME*SPP(J))/2.d0
          P2PP(J)=(PPP(J)-XME*SPP(J))/2.d0
        ENDDO
C
C  Calculate the amplitude
C
      ANS=(SM(P2P,QH)*SP(P1,P2)*SP(K,P1B)*SM(P2PP,P2)
     >     +SP(K,P1P)*SM(P1P,P2P)*SM(QH,P2PP)*SP(P1B,P1)
     >     +SM(P2P,QH)*SP(P1B,P1)*(SP(K,P1P)*SM(P1P,P2PP)
     >      +SP(K,P2P)*SM(P2P,P2PP))
     >     +SP(K,P1P)*SM(P1P,P2P)*SP(P1,P2)*SM(P2PP,P2)/XME2
     >      *(SM(QH,P1P)*SP(P1P,P1B)+SM(QH,P2P)*SP(P2P,P1B)
     >     -SM(QH,K)*SP(K,P1B))+SM(P1PP,P2PP)*SP(P2B,P1B)/XME2*(
     >      SM(P2P,QH)*SP(P1,P2)*SP(K,P1PP)*SM(P2B,P2)
     >      +SP(K,P1P)*SM(P1P,P2P)*SM(QH,P2B)*SP(P1PP,P1)
     >      +SM(P2P,QH)*SP(P1PP,P1)*(SP(K,P1P)*SM(P1P,P2B)
     >       +SP(K,P2P)*SM(P2P,P2B))
     >      +SP(K,P1P)*SM(P1P,P2P)*SP(P1,P2)*SM(P2B,P2)/XME2
     >       *(SM(QH,P1P)*SP(P1P,P1PP)+SM(QH,P2P)*SP(P2P,P1PP)
     >      -SM(QH,K)*SP(K,P1PP))))/SM(K,QH)
C
C  Calculate the normalization
C
      DO I=1,4
        Q(I)=PPP(I)+PB(I)
      ENDDO
      PPK=PP(4)*K(4)-PP(3)*K(3)-PP(2)*K(2)-PP(1)*K(1)
      Q2=Q(4)**2-Q(3)**2-Q(2)**2-Q(1)**2
      ANS=ANS*DCONST/(-2.d0*PPK*Q2)
C
C  Return the answer
C
      D1M=ANS
      RETURN
      END
C
      FUNCTION D2P(P,S,PB,SB,PP,SPR,PPP,SPP,K,QH)
C
C **********************************************************************
C * This function implements the amplitude D2+ for the production of   *
C * the 3e final state.                                                *
C *   Parameters:   P(4) - REAL*8 mom 4-vector (initial e-/final e+)   *
C *                 S(4) - REAL*8 spin 4-vector for above              *
C *                PB(4) - REAL*8 mom 4-vector (initial e-/final e+)   *
C *                SB(4) - REAL*8 spin 4-vector for above              *
C *                PP(4) - REAL*8 mom 4-vector (final e-)              *
C *               SPR(4) - REAL*8 spin 4-vector for above              *
C *               PPP(4) - REAL*8 mom 4-vector (final e-)              *
C *               SPP(4) - REAL*8 spin 4-vector for above              *
C *                 K(4) - REAL*8 mom 4-vector (initial photon)        *
C *                QH(4) - REAL*8 auxilliary momentum 4-vector         *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 D2P,SP,SM
      COMPLEX*16 ANS
      REAL*8 P1(4),P2(4),P1P(4),P2P(4),P1PP(4),P2PP(4),P1B(4),P2B(4),
     >       Q(4)
      REAL*8 P(4),S(4),PB(4),SB(4),PP(4),SPR(4),PPP(4),SPP(4),K(4),
     >       QH(4)
C
C  The block /CONTRL/ contains the adjustable parameters needed for the
C  calculation and some useful physical constants: EB - e- beam energy in
C  GeV, EPHOT - photon energy in GeV, XME - electron mass (in GeV),
C  XMG - the fictitious photon mass (GeV), KGMIN - the minimum resolvable
C  photon energy (GeV), ALPHA - the fine structure constant, PI - 3.1416...
C  ROOT2 - SQRT(2), BARN - (hbar*c)**2 in units of mb-GeV**(-2),
C  SPIN(3) - the initial state electron spin (in the me rest frame),
C  LDIAG - logical flag to turn on diagnostic T-DR-M unpol xsection, 
C  LBF - logical flag to turn on diagnostic Brown+Feynman unpol xsection, 
C  NTRY - the number of event trials to generate with COMTN2,COMEGG,COMEEE
C         (COMEEE generates NTRY/20)
C
      COMMON /CONTRL/ EB,EPHOT,XME,XMG,KGMIN,ALPHA,PI,ROOT2,BARN,
     >                SPIN(3),LDIAG,LBF,NTRY   
      LOGICAL LDIAG,LBF
      REAL*8 KGMIN
      LOGICAL FCALL
      DATA FCALL/.TRUE./
C
C  Initialize things
C
      IF(FCALL) THEN
        DCONST=2.d0*ROOT2*(4.d0*PI*ALPHA)**1.5d0
        XME2=XME**2
        FCALL=.FALSE.
      ENDIF
C
C  Calculate the light-like momenta
C
        DO J=1,4
          P1(J)=(P(J)+XME*S(J))/2.d0
          P2(J)=(P(J)-XME*S(J))/2.d0
          P1B(J)=(PB(J)+XME*SB(J))/2.d0
          P2B(J)=(PB(J)-XME*SB(J))/2.d0
          P1P(J)=(PP(J)+XME*SPR(J))/2.d0
          P2P(J)=(PP(J)-XME*SPR(J))/2.d0
          P1PP(J)=(PPP(J)+XME*SPP(J))/2.d0
          P2PP(J)=(PPP(J)-XME*SPP(J))/2.d0
        ENDDO
C
C  Calculate the amplitude
C
      ANS=(-SM(P1P,P2P)*SP(QH,P1)*SP(P1P,P1B)*SM(P2PP,K)
     >     +SP(P1,P2)*SM(K,P2)*SM(P2P,P2PP)*SP(P1B,QH)
     >     -SP(QH,P1)*SM(P2P,P2PP)*(SP(P1B,P1)*SM(P1,K)
     >      +SP(P1B,P2)*SM(P2,K))
     >     -SM(P1P,P2P)*SP(P1,P2)*SM(K,P2)*SP(P1P,P1B)/XME2
     >      *(SM(P2PP,K)*SP(K,QH)-SM(P2PP,P1)*SP(P1,QH)
     >      -SM(P2PP,P2)*SP(P2,QH))+SM(P1PP,P2PP)*SP(P2B,P1B)/XME2*(
     >     -SM(P1P,P2P)*SP(QH,P1)*SP(P1P,P1PP)*SM(P2B,K)
     >     +SP(P1,P2)*SM(K,P2)*SM(P2P,P2B)*SP(P1PP,QH)
     >     -SP(QH,P1)*SM(P2P,P2B)*(SP(P1PP,P1)*SM(P1,K)
     >      +SP(P1PP,P2)*SM(P2,K))
     >     -SM(P1P,P2P)*SP(P1,P2)*SM(K,P2)*SP(P1P,P1PP)/XME2
     >      *(SM(P2B,K)*SP(K,QH)-SM(P2B,P1)*SP(P1,QH)
     >      -SM(P2B,P2)*SP(P2,QH))))/SP(QH,K)
C
C  Calculate the normalization
C
      DO I=1,4
        Q(I)=PPP(I)+PB(I)
      ENDDO
      PK=P(4)*K(4)-P(3)*K(3)-P(2)*K(2)-P(1)*K(1)
      Q2=Q(4)**2-Q(3)**2-Q(2)**2-Q(1)**2
      ANS=ANS*DCONST/(-2.d0*PK*Q2)
C
C  Return the answer
C
      D2P=ANS
      RETURN
      END
C
      FUNCTION D2M(P,S,PB,SB,PP,SPR,PPP,SPP,K,QH)
C
C **********************************************************************
C * This function implements the amplitude D2- for the production of   *
C * the 3e final state.                                                *
C *   Parameters:   P(4) - REAL*8 mom 4-vector (initial e-/final e+)   *
C *                 S(4) - REAL*8 spin 4-vector for above              *
C *                PB(4) - REAL*8 mom 4-vector (initial e-/final e+)   *
C *                SB(4) - REAL*8 spin 4-vector for above              *
C *                PP(4) - REAL*8 mom 4-vector (final e-)              *
C *               SPR(4) - REAL*8 spin 4-vector for above              *
C *               PPP(4) - REAL*8 mom 4-vector (final e-)              *
C *               SPP(4) - REAL*8 spin 4-vector for above              *
C *                 K(4) - REAL*8 mom 4-vector (initial photon)        *
C *                QH(4) - REAL*8 auxilliary momentum 4-vector         *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 D2M,SP,SM
      COMPLEX*16 ANS
      REAL*8 P1(4),P2(4),P1P(4),P2P(4),P1PP(4),P2PP(4),P1B(4),P2B(4),
     >       Q(4)
      REAL*8 P(4),S(4),PB(4),SB(4),PP(4),SPR(4),PPP(4),SPP(4),K(4),
     >       QH(4)
C
C  The block /CONTRL/ contains the adjustable parameters needed for the
C  calculation and some useful physical constants: EB - e- beam energy in
C  GeV, EPHOT - photon energy in GeV, XME - electron mass (in GeV),
C  XMG - the fictitious photon mass (GeV), KGMIN - the minimum resolvable
C  photon energy (GeV), ALPHA - the fine structure constant, PI - 3.1416...
C  ROOT2 - SQRT(2), BARN - (hbar*c)**2 in units of mb-GeV**(-2),
C  SPIN(3) - the initial state electron spin (in the me rest frame),
C  LDIAG - logical flag to turn on diagnostic T-DR-M unpol xsection, 
C  LBF - logical flag to turn on diagnostic Brown+Feynman unpol xsection, 
C  NTRY - the number of event trials to generate with COMTN2,COMEGG,COMEEE
C         (COMEEE generates NTRY/20)
C
      COMMON /CONTRL/ EB,EPHOT,XME,XMG,KGMIN,ALPHA,PI,ROOT2,BARN,
     >                SPIN(3),LDIAG,LBF,NTRY   
      LOGICAL LDIAG,LBF
      REAL*8 KGMIN
      LOGICAL FCALL
      DATA FCALL/.TRUE./
C
C  Initialize things
C
      IF(FCALL) THEN
        DCONST=2.d0*ROOT2*(4.d0*PI*ALPHA)**1.5d0
        XME2=XME**2
        FCALL=.FALSE.
      ENDIF
C
C  Calculate the light-like momenta
C
        DO J=1,4
          P1(J)=(P(J)+XME*S(J))/2.d0
          P2(J)=(P(J)-XME*S(J))/2.d0
          P1B(J)=(PB(J)+XME*SB(J))/2.d0
          P2B(J)=(PB(J)-XME*SB(J))/2.d0
          P1P(J)=(PP(J)+XME*SPR(J))/2.d0
          P2P(J)=(PP(J)-XME*SPR(J))/2.d0
          P1PP(J)=(PPP(J)+XME*SPP(J))/2.d0
          P2PP(J)=(PPP(J)-XME*SPP(J))/2.d0
        ENDDO
C
C  Calculate the amplitude
C
      ANS=(-SM(P1P,P2P)*SP(K,P1)*SP(P1P,P1B)*SM(P2PP,QH)
     >     +SP(P1,P2)*SM(QH,P2)*SM(P2P,P2PP)*SP(P1B,K)
     >     +SP(K,P1)*SM(P2P,P2PP)*(SP(P1B,K)*SM(K,QH)
     >      -SP(P1B,P1)*SM(P1,QH)-SP(P1B,P2)*SM(P2,QH))
     >     +SM(P1P,P2P)*SP(P1,P2)*SM(QH,P2)*SP(P1P,P1B)/XME2
     >      *(SM(P2PP,P1)*SP(P1,K)+SM(P2PP,P2)*SP(P2,K))
     >      +SM(P1PP,P2PP)*SP(P2B,P1B)/XME2*(
     >      -SM(P1P,P2P)*SP(K,P1)*SP(P1P,P1PP)*SM(P2B,QH)
     >      +SP(P1,P2)*SM(QH,P2)*SM(P2P,P2B)*SP(P1PP,K)
     >      +SP(K,P1)*SM(P2P,P2B)*(SP(P1PP,K)*SM(K,QH)
     >       -SP(P1PP,P1)*SM(P1,QH)-SP(P1PP,P2)*SM(P2,QH))
     >      +SM(P1P,P2P)*SP(P1,P2)*SM(QH,P2)*SP(P1P,P1PP)/XME2
     >       *(SM(P2B,P1)*SP(P1,K)+SM(P2B,P2)*SP(P2,K))))/SM(K,QH)
C
C  Calculate the normalization
C
      DO I=1,4
        Q(I)=PPP(I)+PB(I)
      ENDDO
      PK=P(4)*K(4)-P(3)*K(3)-P(2)*K(2)-P(1)*K(1)
      Q2=Q(4)**2-Q(3)**2-Q(2)**2-Q(1)**2
      ANS=ANS*DCONST/(-2.d0*PK*Q2)
C
C  Return the answer
C
      D2M=ANS
      RETURN
      END
C
      FUNCTION SP(P1,P2)
C
C **********************************************************************
C * This function implements the spinor product S+ of Kleiss and       *
C * Stirling defined in Stuart and Gongora, Z.Phys. C42,617 (1989).    *
C *   Parameters: P1(4) - REAL*8 4-vector                              *
C *               P2(4) - REAL*8 4-vector                              *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 SP
      COMPLEX*16 ANS,I
      REAL*8 P1(4), P2(4)
      DATA I/(0.d0,1.d0)/
C
C  Make certain that both 4-vectors have positive energies
C
      IF(P1(4).GT.0.d0.AND.P2(4).GT.0.d0) THEN
        ARG1=P2(4)-P2(1)
        ARG2=P1(4)-P1(1)
        ANS=DCMPLX(P1(2),P1(3))*SQRT(ARG1/ARG2)
     >     -DCMPLX(P2(2),P2(3))*SQRT(ARG2/ARG1)
      ELSEIF(P1(4).LT.0.d0.AND.P2(4).GT.0.d0) THEN
        ARG1=P2(4)-P2(1)
        ARG2=-(P1(4)-P1(1))
        ANS=-I*(DCMPLX(P1(2),P1(3))*SQRT(ARG1/ARG2)
     >         +DCMPLX(P2(2),P2(3))*SQRT(ARG2/ARG1))
      ELSEIF(P1(4).GT.0.d0.AND.P2(4).LT.0.d0) THEN
        ARG1=-(P2(4)-P2(1))
        ARG2=P1(4)-P1(1)
        ANS=I*(DCMPLX(P1(2),P1(3))*SQRT(ARG1/ARG2)
     >         +DCMPLX(P2(2),P2(3))*SQRT(ARG2/ARG1))
      ELSE
        ARG1=-(P2(4)-P2(1))
        ARG2=-(P1(4)-P1(1))
        ANS=DCMPLX(P1(2),P1(3))*SQRT(ARG1/ARG2)
     >         -DCMPLX(P2(2),P2(3))*SQRT(ARG2/ARG1)
      ENDIF
      SP=ANS
      RETURN
      END
C
      FUNCTION SM(P1,P2)
C
C **********************************************************************
C * This function implements the spinor product S- of Kleiss and       *
C * Stirling defined in Stuart and Gongora, Z.Phys. C42,617 (1989).    *
C *   Parameters: P1(4) - REAL*8 4-vector                              *
C *               P2(4) - REAL*8 4-vector                              *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL SP
      COMPLEX*16 SM,SP
      COMPLEX*16 ANS,I
      REAL*8 P1(4), P2(4), R1(4), R2(4)
      DATA I/(0.d0,1.d0)/
C
C  Make certain that both 4-vectors have positive energies
C
      IF(P1(4).GT.0.d0.AND.P2(4).GT.0.d0) THEN
C
C  Answer is defined in terms of SP
C
        ANS=-CONJG(SP(P1,P2))
      ELSEIF(P1(4).LT.0.d0.AND.P2(4).GT.0.d0) THEN
        DO J=1,4
          R1(J)=-P1(J)
        ENDDO
        ANS=-I*CONJG(SP(R1,P2))
      ELSEIF(P1(4).GT.0.d0.AND.P2(4).LT.0.d0) THEN
        DO J=1,4
          R2(J)=-P2(J)
        ENDDO
        ANS=-I*CONJG(SP(P1,R2))
      ELSE
        DO J=1,4
          R1(J)=-P1(J)
          R2(J)=-P2(J)
        ENDDO
        ANS=CONJG(SP(R1,R2))
      ENDIF
      SM=ANS
      RETURN
      END
C
      FUNCTION G0(X)
C
C **********************************************************************
C * This function implements the G0(X) as defined in Milton,Tsai and  *
C * DeRaad, PR D6, 1428 (1972) and Brown+Feynman, PR 85, 231 (1952).   *
C *   Parameter: X - REAL*8                                            *
C **********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 G0
      REAL*8 X,U,PI,C1,REG0,IMG0
C      LOGICAL FCALL
C      DATA FCALL/.TRUE./
C      IF(FCALL) THEN
        PI=ACOS(-1.d0)
        C1=PI**2/6.d0
C        print *, PI, C1
C        FCALL=.FALSE.
C      ENDIF
      U=1.d0-X
      REG0=2.d0/X*(C1-DDILOG(U))
      IF(X.LT.0.d0) THEN
        IMG0=-2.d0/X*PI*LOG(U)
        G0=DCMPLX(REG0,IMG0)
      ELSE
        G0=DCMPLX(REG0,0.d0)
      ENDIF
      RETURN
      END
C
      SUBROUTINE NEWDIR(V,COST,PHI)
C
C **************************************************
C  THIS ROUTINE ROTATES THE Z AXIS TO THE DIRECTION
C  ACOS(COST),PHI.
C **************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(3),ROTMAT(3,3),TMP(3)
      DATA COSLST,PHILST,CPHI,SPHI/0.,0.,1.,0./
      IF((COST.EQ.COSLST).AND.(PHILST.EQ.PHI)) GO TO 100
      COSLST=COST
      SINT=DSQRT(1.-COST**2)
      PHILST=PHI
      CPHI=DCOS(PHI)
      SPHI=DSIN(PHI)
      ROTMAT(1,1)=COST*CPHI
      ROTMAT(2,1)=COST*SPHI
      ROTMAT(3,1)=-SINT
      ROTMAT(1,2)=-SPHI
      ROTMAT(2,2)=CPHI
      ROTMAT(3,2)=0.
      ROTMAT(1,3)=SINT*CPHI
      ROTMAT(2,3)=SINT*SPHI
      ROTMAT(3,3)=COST
100   CONTINUE
      DO 200 I=1,3
      TMP(I)=V(1)*ROTMAT(I,1)+V(2)*ROTMAT(I,2)+V(3)*ROTMAT(I,3)
200   CONTINUE
      DO 300 I=1,3
      V(I)=TMP(I)
300   CONTINUE
      RETURN
      END
C
      SUBROUTINE LORENZ(IQ,PCMQ,I1,I2,PCMTR)
C ***********************************************************************
C * LORENTZ TRANSFORMATION FROM THE REST SYSTEM  (DSQRT(Q2),0) TO THE   *
C * BOOSTED SYSTEM PCMQ(*,IQ) FOR THE VECTORS PCMTR(*,I1)...PCMTR(*,I2) *
C ***********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PCMQ(4,18),PCMTR(4,18)
      Q0=PCMQ(4,IQ)
      Q2=Q0**2-PCMQ(1,IQ)**2-PCMQ(2,IQ)**2-PCMQ(3,IQ)**2
      IF(Q2.LE.0.)GOTO 500
      Q=DSQRT(Q2)
      DO I=I1,I2
        E=PCMTR(4,I)
        PCMTR(4,I)=(Q0*E+PCMQ(1,IQ)*PCMTR(1,I)+PCMQ(2,IQ)*PCMTR(2,I)
     >            +PCMQ(3,IQ)*PCMTR(3,I))/Q
        F=(E+PCMTR(4,I))/(Q0+Q)
        DO J=1,3
          PCMTR(J,I)=PCMTR(J,I)+PCMQ(J,IQ)*F
        ENDDO
      ENDDO
      RETURN
500   WRITE(6,501)Q2
501   FORMAT(18H LORENZ ERROR, Q2=,E14.4)
      STOP
      END
C
      SUBROUTINE WGTHST(IFLAG,NEM,PP,NGAM,QP,NEP,PB,WGT)
C
C ***********************************************************************
C * This routine is an interface to the 2=>2 and 2=>3 MC generators     *
C * COMTN2, COMEEG, COMEEE.  The parameters are:                       *
C *   IFLAG - determines the function of this routine:                  *
C *              = 0 - initialize everything                            *
C *              = 1 - accumulate event weights                         *
C *              = 2 - print the results                                *
C *     NEM - the number of electrons in the final state (1 or 2)       *
C * PP(4,2) - the 4-vectors of the scattered electrons in the lab frame *
C *    NGAM - the number of photons in the final state (0-2)            *
C * QP(4,2) - the lab-frame 4-vectors of the scattered photons          *
C *     NEP - the number of positrons in the final state (0 or 1)       *
C *   PB(4) - the 4-vector of the scattered positron in the lab frame   *
C *  WGT(4) - the event weights for the event as follows:               *
C *              = 1 tree-level unpolarized cross section               *
C *              = 2 tree-level polarized xsection (sig(-)-sig(+))/2    *
C *              = 3 order-alpha unpolarized cross section              *
C *              = 4 order-alpha polarized cross section                *
C ***********************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  Define histograms for unpolarized and polarized cross sections in
C  zeroth order (i=1,2) and for the order-alpha corrections (i=3,4).
C  (Keep full error matrix to accurately assess uncertainties.)
C  Histogram the longitudinal asymmetry in momentum space (PSUM),
C  keep track of the full 4x4 error matrix for each bin (PSUM2).
C  Histogram the projected vertical angle in photon energy space (TYSUM),
C  keep track of the full 4x4 error matrix for each bin (TYSUM2).
C  Do energy-weighted sum of photon energy (ESUM) with the lab frame cut
C  on maximum photon angle TMAX.  The 4x4 error matrix is kept in ESUM2.
C  Do energy distribution of the backscattered gammas (EGSUM) for a single
C  photon device with lab frame cut on maximum photon angle TMAX.  The 4x4 
C  error matrix is kept in EGSUM2
C
      PARAMETER (NB=100)
      DIMENSION PSUM(4,NB),PSUM2(4,4,NB),ESUM(4),ESUM2(4,4)
     >   ,TYSUM(4,NB),TYSUM2(4,4,NB),EGSUM(4,NB),EGSUM2(4,4,NB)
      LOGICAL LONG,TRANS
C
C  The block /CONTRL/ contains the adjustable parameters needed for the
C  calculation and some useful physical constants: EB - e- beam energy in
C  GeV, EPHOT - photon energy in GeV, XME - electron mass (in GeV),
C  XMG - the fictitious photon mass (GeV), KGMIN - the minimum resolvable
C  photon energy (GeV), ALPHA - the fine structure constant, PI - 3.1416...
C  ROOT2 - SQRT(2), BARN - (hbar*c)**2 in units of mb-GeV**(-2),
C  SPIN(3) - the initial state electron spin (in the me rest frame),
C  LDIAG - logical flag to turn on diagnostic T-DR-M unpol xsection, 
C  LBF - logical flag to turn on diagnostic Brown+Feynman unpol xsection, 
C  NTRY - the number of event trials to generate with COMTN2,COMEGG,COMEEE
C         (COMEEE generates NTRY/20)
C
      COMMON /CONTRL/ EB,EPHOT,XME,XMG,KGMIN,ALPHA,PI,ROOT2,BARN,
     >                SPIN(3),LDIAG,LBF,NTRY   
      LOGICAL LDIAG,LBF
      REAL*8 KGMIN
      DIMENSION DDADX(4),DDAADX(4)
      DIMENSION PP(4,2),QP(4,2),PB(4),WGT(4)
C
      IF(IFLAG.EQ.0) THEN
C
C *************************
C * Initialize everything *
C *************************
C
C  TMAX is the maximum accepted angle of the photons in the laboratory
C  frame (for the photon energy asymmetry)
C  
        TMAX=1.15d0/10.d2
C
C  Set up histogramming flags
C
        IF(SPIN(3).NE.0.d0) THEN
          LONG=.TRUE.
        ELSE
          LONG=.FALSE.
        ENDIF
        SPINT=SQRT(SPIN(1)**2+SPIN(2)**2)
        IF(SPINT.GT.0.d0) THEN
          TRANS=.TRUE.
C
C  Calculate the components of the transverse polarization vector
C
          PTX=SPIN(1)/SPINT
          PTY=SPIN(2)/SPINT
        ELSE
          TRANS=.FALSE.
        ENDIF
C
C  Zero the histograms
C
        DO K=1,NB
          DO I=1,4
            PSUM(I,K)=0.d0
            TYSUM(I,K)=0.d0
            EGSUM(I,K)=0.d0
            DO J=1,4
              PSUM2(I,J,K)=0.d0
              TYSUM2(I,J,K)=0.d0
              EGSUM2(I,J,K)=0.d0
            ENDDO
          ENDDO
        ENDDO
        DO I=1,4
          ESUM(I)=0.d0
          DO J=1,4
            ESUM2(I,J)=0.d0
          ENDDO
        ENDDO
C
C  Work out histogram binning
C
        PBIN=EB/FLOAT(NB)
C
      ELSEIF(IFLAG.EQ.1) THEN
C
C ****************************
C * Accumulate event weights *
C ****************************
C
C  Normalize the event weights for the appropriate analysis
C
        WU0=WGT(1)
        WU1=WGT(3)
        IF(LONG) THEN
C
C  divide out partial polarization to get asymmetries for 100% pol
C
          WZ0=WGT(2)/SPIN(3)
          WZ1=WGT(4)/SPIN(3)
        ELSE
          WZ0=0.d0
          WZ1=0.d0
        ENDIF
C
C  Loop over all e- in the final state
C 
        DO IEM=1,NEM
C
C  Compute the energy bin
C
          ELAB=PP(4,IEM)
          IELAB=INT(ELAB/PBIN)+1
          IF(IELAB.LT.1) IELAB=1
          IF(IELAB.GT.NB) IELAB=NB
C
C  Sum the total event weights by electron energy bin
C
          PSUM(1,IELAB)=PSUM(1,IELAB)+WU0
          PSUM(2,IELAB)=PSUM(2,IELAB)+WZ0
          PSUM(3,IELAB)=PSUM(3,IELAB)+WU1
          PSUM(4,IELAB)=PSUM(4,IELAB)+WZ1
C
C  Do the full 4x4 error matrix for each energy bin
C
          PSUM2(1,1,IELAB)=PSUM2(1,1,IELAB)+WU0*WU0
          PSUM2(1,2,IELAB)=PSUM2(1,2,IELAB)+WU0*WZ0
          PSUM2(1,3,IELAB)=PSUM2(1,3,IELAB)+WU0*WU1
          PSUM2(1,4,IELAB)=PSUM2(1,4,IELAB)+WU0*WZ1
          PSUM2(2,2,IELAB)=PSUM2(2,2,IELAB)+WZ0*WZ0
          PSUM2(2,3,IELAB)=PSUM2(2,3,IELAB)+WZ0*WU1
          PSUM2(2,4,IELAB)=PSUM2(2,4,IELAB)+WZ0*WZ1
          PSUM2(3,3,IELAB)=PSUM2(3,3,IELAB)+WU1*WU1
          PSUM2(3,4,IELAB)=PSUM2(3,4,IELAB)+WU1*WZ1
          PSUM2(4,4,IELAB)=PSUM2(4,4,IELAB)+WZ1*WZ1
        ENDDO
C
C  Check for photons
C
        IF(NGAM.GT.0) THEN
C
C  Loop over photons
C
          DO K=1,NGAM
C
C  Compute the gamma angle wrt the beam
C
            EQLAB=QP(4,K)
            PTQ=SQRT(QP(1,K)**2+QP(2,K)**2)
            TQLAB=PTQ/EQLAB
C
C  Compute the energy bin
C
            IQLAB=INT(EQLAB/PBIN)+1
            IF(IQLAB.LT.1) IQLAB=1
            IF(IQLAB.GT.NB) IQLAB=NB
C
C  Check to see if it's accepted by the detector
C
            IF(TQLAB.LT.TMAX) THEN
              IF(LONG) THEN
C
C  Do the energy-weighted cross sections
C
                WUE0=WU0*EQLAB
                WZE0=WZ0*EQLAB
                WUE1=WU1*EQLAB
                WZE1=WZ1*EQLAB
C
C  Accumulate each type
C
                ESUM(1)=ESUM(1)+WUE0
                ESUM(2)=ESUM(2)+WZE0
                ESUM(3)=ESUM(3)+WUE1
                ESUM(4)=ESUM(4)+WZE1
C
C  Do the full 4x4 error matrix for each energy bin
C
                ESUM2(1,1)=ESUM2(1,1)+WUE0*WUE0
                ESUM2(1,2)=ESUM2(1,2)+WUE0*WZE0
                ESUM2(1,3)=ESUM2(1,3)+WUE0*WUE1
                ESUM2(1,4)=ESUM2(1,4)+WUE0*WZE1
                ESUM2(2,2)=ESUM2(2,2)+WZE0*WZE0
                ESUM2(2,3)=ESUM2(2,3)+WZE0*WUE1
                ESUM2(2,4)=ESUM2(2,4)+WZE0*WZE1
                ESUM2(3,3)=ESUM2(3,3)+WUE1*WUE1
                ESUM2(3,4)=ESUM2(3,4)+WUE1*WZE1
                ESUM2(4,4)=ESUM2(4,4)+WZE1*WZE1

C
C  Sum the total event weights by photon energy bin
C
                EGSUM(1,IQLAB)=EGSUM(1,IQLAB)+WU0
                EGSUM(2,IQLAB)=EGSUM(2,IQLAB)+WZ0
                EGSUM(3,IQLAB)=EGSUM(3,IQLAB)+WU1
                EGSUM(4,IQLAB)=EGSUM(4,IQLAB)+WZ1
C
C  Do the full 4x4 error matrix for each energy bin
C
                EGSUM2(1,1,IQLAB)=EGSUM2(1,1,IQLAB)+WU0*WU0
                EGSUM2(1,2,IQLAB)=EGSUM2(1,2,IQLAB)+WU0*WZ0
                EGSUM2(1,3,IQLAB)=EGSUM2(1,3,IQLAB)+WU0*WU1
                EGSUM2(1,4,IQLAB)=EGSUM2(1,4,IQLAB)+WU0*WZ1
                EGSUM2(2,2,IQLAB)=EGSUM2(2,2,IQLAB)+WZ0*WZ0
                EGSUM2(2,3,IQLAB)=EGSUM2(2,3,IQLAB)+WZ0*WU1
                EGSUM2(2,4,IQLAB)=EGSUM2(2,4,IQLAB)+WZ0*WZ1
                EGSUM2(3,3,IQLAB)=EGSUM2(3,3,IQLAB)+WU1*WU1
                EGSUM2(3,4,IQLAB)=EGSUM2(3,4,IQLAB)+WU1*WZ1
                EGSUM2(4,4,IQLAB)=EGSUM2(4,4,IQLAB)+WZ1*WZ1
              ENDIF
            ENDIF
          ENDDO
          IF(TRANS) THEN
            EQLAB=0.d0
            TYXE=0.d0
C
C  Loop over the photons in each event
C
            DO K=1,NGAM
C
C  Consider only those photons which hit the detector
C
              TQLAB=SQRT(QP(1,K)**2+QP(2,K)**2)/QP(4,K)
              IF(TQLAB.LT.TMAX) THEN
C
C  Compute the total gamma energy hitting the detector
C
                EQLAB=EQLAB+QP(4,K)
C
C  Compute the transverse angle (along the pol dir) times the gamma energy
C
                TYXE=TYXE+PTX*QP(1,K)+PTY*QP(2,K)
              ENDIF
            ENDDO
C
C  Calculate the energy-weighted average projected vertical angle
C  (we assume that the HERA tranverse polarimeter does this)
C
            TYXE=TYXE/EQLAB
C
C  Compute photon energy bin
C
            X=EQLAB
            IX=INT(X/PBIN)+1
            IF(IX.LT.1) IX=1
            IF(IX.GT.NB) IX=NB
            WUE0=WU0
            WUE1=WU1
C
C  Compute mean projected angle difference for the helicity states
C  (divide out partial polarization to get asymmetries for 100% pol)
C
            WTY0=2.d0*WGT(2)*TYXE/SPINT
            WTY1=2.d0*WGT(4)*TYXE/SPINT
            TYSUM(1,IX)=TYSUM(1,IX)+WUE0
            TYSUM(2,IX)=TYSUM(2,IX)+WTY0
            TYSUM(3,IX)=TYSUM(3,IX)+WUE1
            TYSUM(4,IX)=TYSUM(4,IX)+WTY1
C
C  Accumulate the full 4x4 error matrix for each energy bin
C
            TYSUM2(1,1,IX)=TYSUM2(1,1,IX)+WUE0*WUE0
            TYSUM2(1,2,IX)=TYSUM2(1,2,IX)+WUE0*WTY0
            TYSUM2(1,3,IX)=TYSUM2(1,3,IX)+WUE0*WUE1
            TYSUM2(1,4,IX)=TYSUM2(1,4,IX)+WUE0*WTY1
            TYSUM2(2,2,IX)=TYSUM2(2,2,IX)+WTY0*WTY0
            TYSUM2(2,3,IX)=TYSUM2(2,3,IX)+WTY0*WUE1
            TYSUM2(2,4,IX)=TYSUM2(2,4,IX)+WTY0*WTY1
            TYSUM2(3,3,IX)=TYSUM2(3,3,IX)+WUE1*WUE1
            TYSUM2(3,4,IX)=TYSUM2(3,4,IX)+WUE1*WTY1
            TYSUM2(4,4,IX)=TYSUM2(4,4,IX)+WTY1*WTY1
          ENDIF
        ENDIF
C
      ELSEIF(IFLAG.EQ.2) THEN
C
C *******************************
C * Calculate and print results *
C *******************************
C
C
C  First, symmetrize the error matrices
C
        DO I=1,3
          L=I+1
          DO J=L,4
            ESUM2(J,I)=ESUM2(I,J)
            DO K=1,NB
              PSUM2(J,I,K)=PSUM2(I,J,K)
              TYSUM2(J,I,K)=TYSUM2(I,J,K)
              EGSUM2(J,I,K)=EGSUM2(I,J,K)
            ENDDO
          ENDDO
        ENDDO
C
C  Do the unpolarized cross section first
C
        PRINT 3500
3500    FORMAT(/,10x,'***** Unpolarized Cross Section for e- *****',/)
        PRINT 4000
4000    FORMAT(1x,'  Energy Bin   cross section0   (1-0)/0 (%)')
        DO I=1,NB
          EE=(I-0.5)*PBIN
          IF(PSUM(1,I).GT.0.D0) THEN
            CHANGE=PSUM(3,I)/PSUM(1,I)
C
C  The correct errors involve correlations between the various x-sections
C  (COMTN2 generates 4 weights/trial....correlations improve the precision
C   on many quantities)
C
            SCHANG=CHANGE*SQRT(PSUM2(1,1,I)/PSUM(1,I)**2
     >            +PSUM2(3,3,I)/PSUM(3,I)**2
     >            -2.d0*PSUM2(1,3,I)/(PSUM(1,I)*PSUM(3,I)))
          ELSE
            CHANGE=0.d0
            SCHANG=0.d0
          ENDIF
          PRINT 5000, EE,PSUM(1,I),SQRT(PSUM2(1,1,I)),
     >                change*100.,schang*100. 
5000      FORMAT(1x,5x,f7.2,2x,f6.1,
     >           '+-',f6.1,2x,f6.3,'+-',f6.4)
        ENDDO
C
C  Next the longitudinal asymmetry of the electrons
C
        IF(LONG) THEN
          PRINT 5500
5500      FORMAT(/,10x,
     >    '***** Longitudinal Polarization Asymmetry for e- *****',/)        
          PRINT 6000
6000     FORMAT(1x,'  Energy Bin     Asymmetry0      (1-0)/0 (%)     ',
     >              '(1-0)*100')       
          DO I=1,NB
            EE=(I-0.5)*PBIN
            IF(PSUM(1,I).GT.0.D0) THEN
C
C  The order-alpha**2 longitudinal asymmetry A(0) and its uncertainty
C
              ASY0=PSUM(2,I)/PSUM(1,I)
              SASY0=ASY0*SQRT(PSUM2(1,1,I)/PSUM(1,I)**2
     >              +PSUM2(2,2,I)/PSUM(2,I)**2
     >              -2.d0*PSUM2(1,2,I)/(PSUM(1,I)*PSUM(2,I)))
              TOTXS=PSUM(1,I)+PSUM(3,I)
              DEN=PSUM(1,I)*TOTXS
C
C  DIFF is the difference between PSUM(2,I)/PSUM(1,I) (ie A(0)) the
C  fully corrected asymmetry (PSUM(2,I)+PSUM(4,I))/(PSUM(1,I)+PSUM(3,I)).
C
              DIFF=(PSUM(1,I)*PSUM(4,I)-PSUM(2,I)*PSUM(3,I))/DEN
C
C  Calculate the derivatives of DIFF given in equation 40 
C
              DDADX(1)=(-PSUM(1,I)**2*PSUM(4,I)
     >                +2.d0*PSUM(1,I)*PSUM(2,I)*PSUM(3,I)
     >                +PSUM(3,I)**2*PSUM(2,I))/DEN**2
              DDADX(2)=-PSUM(3,I)/DEN
              DDADX(3)=-(PSUM(2,I)+PSUM(4,I))/TOTXS**2
              DDADX(4)=1.d0/TOTXS
C
C  Calculate the uncertainty on DIFF from equation 40 
C
              SDIFF=0.d0
              DO K=1,4
                DO J=1,4
                  SDIFF=SDIFF+DDADX(K)*DDADX(J)*PSUM2(K,J,I)
                ENDDO
              ENDDO
              SDIFF=SQRT(SDIFF)
C
C  CHANGE is the fractional difference, DIFF/A(0).
C
              CHANGE=DIFF/ASY0
C
C  Calculate the derivatives of CHANGE given in equation 40 
C
              DDAADX(1)=PSUM(3,I)*(PSUM(2,I)+PSUM(4,I))
     >                 /(PSUM(2,I)*TOTXS**2)        
              DDAADX(2)=-PSUM(1,I)*PSUM(4,I)/(PSUM(2,I)**2*TOTXS)
              DDAADX(3)=-PSUM(1,I)*(PSUM(2,I)+PSUM(4,I))
     >                 /(PSUM(2,I)*TOTXS**2)
              DDAADX(4)=PSUM(1,I)/(PSUM(2,I)*TOTXS)
C
C  Calculate the uncertainty on CHANGE from equation 40 
C
              SCHANG=0.d0
              DO K=1,4
                DO J=1,4
                  SCHANG=SCHANG+DDAADX(K)*DDAADX(J)*PSUM2(K,J,I)
                ENDDO
              ENDDO
              SCHANG=SQRT(SCHANG)
            ELSE
              ASY0=0.d0
              SASY0=0.d0
              DIFF=0.d0
              SDIFF=0.d0
              CHANGE=0.d0
              SCHANG=0.d0
            ENDIF
            PRINT 7000, EE,asy0,sasy0,change*100.,schang*100.,
     >                  diff*100.,sdiff*100
7000        FORMAT(1x,5x,f7.2,2x,f6.3,'+-',f6.3,2x,f6.3,'+-',f6.4,2x,
     >             f6.3,'+-',f6.4)
          ENDDO
        ENDIF
C
C  Next, the projected vertical angle distribution vs photon energy
C
        IF(TRANS) THEN
          PRINT 11050
11050     FORMAT(/,10x,
     >           '***** Transverse Polarization Asymmetry *****',/)
          PRINT 11120
11120     FORMAT(1x,' Energy Bin    Vertical Angle0 (rad)  (1-0)/0 (%)')
          DO I=1,NB
            X=(I-0.5)*PBIN
            IF(TYSUM(1,I).GT.0.D0) THEN
C
C  The order-alpha**2 mean projected angle TY0 and its uncertainty
C
              TY0=TYSUM(2,I)/TYSUM(1,I)
              STY0=TY0*SQRT(TYSUM2(1,1,I)/TYSUM(1,I)**2
     >              +TYSUM2(2,2,I)/TYSUM(2,I)**2
     >              -2.d0*TYSUM2(1,2,I)/(TYSUM(1,I)*TYSUM(2,I)))
              TOTXS=TYSUM(1,I)+TYSUM(3,I)
              DEN=TYSUM(1,I)*TOTXS
C
C  DIFF is the difference between TYSUM(2,I)/TYSUM(1,I) (ie Ty(0)) the
C  fully corrected shift (TYSUM(2,I)+TYSUM(4,I))/(TYSUM(1,I)+TYSUM(3,I)).
C
              DIFF=(TYSUM(1,I)*TYSUM(4,I)-TYSUM(2,I)*TYSUM(3,I))/DEN
C
C  Calculate the derivatives of DIFF given in equation 40 
C
              DDADX(1)=(-TYSUM(1,I)**2*TYSUM(4,I)
     >                +2.d0*TYSUM(1,I)*TYSUM(2,I)*TYSUM(3,I)
     >                +TYSUM(3,I)**2*TYSUM(2,I))/DEN**2
              DDADX(2)=-TYSUM(3,I)/DEN
              DDADX(3)=-(TYSUM(2,I)+TYSUM(4,I))/TOTXS**2
              DDADX(4)=1.d0/TOTXS
C
C  Calculate the uncertainty on DIFF from equation 40 
C
              SDIFF=0.d0
              DO K=1,4
                DO J=1,4
                  SDIFF=SDIFF+DDADX(K)*DDADX(J)*TYSUM2(K,J,I)
                ENDDO
              ENDDO
              SDIFF=SQRT(SDIFF)
C
C  CHANGE is the fractional difference, DIFF/Ty(0).
C
              CHANGE=DIFF/TY0
C
C  Calculate the derivatives of CHANGE given in equation 40 
C
              DDAADX(1)=TYSUM(3,I)*(TYSUM(2,I)+TYSUM(4,I))
     >                 /(TYSUM(2,I)*TOTXS**2)        
              DDAADX(2)=-TYSUM(1,I)*TYSUM(4,I)/(TYSUM(2,I)**2*TOTXS)
              DDAADX(3)=-TYSUM(1,I)*(TYSUM(2,I)+TYSUM(4,I))
     >                 /(TYSUM(2,I)*TOTXS**2)
              DDAADX(4)=TYSUM(1,I)/(TYSUM(2,I)*TOTXS)
C
C  Calculate the uncertainty on CHANGE from equation 40 
C
              SCHANG=0.d0
              DO K=1,4
                DO J=1,4
                  SCHANG=SCHANG+DDAADX(K)*DDAADX(J)*TYSUM2(K,J,I)
                ENDDO
              ENDDO
              SCHANG=SQRT(SCHANG)
            ELSE
              TY0=0.d0
              STY0=0.d0
              DIFF=0.d0
              SDIFF=0.d0
              CHANGE=0.d0
              SCHANG=0.d0
            ENDIF
            PRINT 11140, X,TY0,STY0,change*100.,schang*100.
11140       FORMAT(5x,f7.2,2x,d10.3,'+-',d10.3,2x,f6.3,'+-',f6.3)
          ENDDO
        ENDIF
C
C  Finally, do the total photon energy asymmetry
C
        IF(LONG) THEN
          IF(ESUM(1).GT.0.D0) THEN
C
C  The order-alpha**2 energy asymmetry ASY0 and its uncertainty
C
            ASY0=ESUM(2)/ESUM(1)
            SASY0=ASY0*SQRT(ESUM2(1,1)/ESUM(1)**2
     >            +ESUM2(2,2)/ESUM(2)**2
     >            -2.d0*ESUM2(1,2)/(ESUM(1)*ESUM(2)))
            TOTXS=ESUM(1)+ESUM(3)
            DEN=ESUM(1)*TOTXS
C
C  DIFF is the difference between ESUM(2)/ESUM(1) (ie ASY(0)) the
C  fully corrected asymmetry (ESUM(2)+ESUM(4))/(ESUM(1)+ESUM(I)).
C
            DIFF=(ESUM(1)*ESUM(4)-ESUM(2)*ESUM(3))/DEN
C
C  Calculate the derivatives of DIFF given in equation 40 
C
            DDADX(1)=(-ESUM(1)**2*ESUM(4)
     >              +2.d0*ESUM(1)*ESUM(2)*ESUM(3)
     >              +ESUM(3)**2*ESUM(2))/DEN**2
            DDADX(2)=-ESUM(3)/DEN
            DDADX(3)=-(ESUM(2)+ESUM(4))/TOTXS**2
            DDADX(4)=1.d0/TOTXS
C
C  Calculate the uncertainty on DIFF from equation 40 
C
            SDIFF=0.d0
            DO K=1,4
              DO J=1,4
                SDIFF=SDIFF+DDADX(K)*DDADX(J)*ESUM2(K,J)
              ENDDO
            ENDDO
            SDIFF=SQRT(SDIFF)
C
C  CHANGE is the fractional difference, DIFF/ASY0.
C
            CHANGE=DIFF/ASY0
C
C  Calculate the derivatives of CHANGE given in equation 40 
C
            DDAADX(1)=ESUM(3)*(ESUM(2)+ESUM(4))/(ESUM(2)*TOTXS**2)        
            DDAADX(2)=-ESUM(1)*ESUM(4)/(ESUM(2)**2*TOTXS)
            DDAADX(3)=-ESUM(1)*(ESUM(2)+ESUM(4))
     >               /(ESUM(2)*TOTXS**2)
            DDAADX(4)=ESUM(1)/(ESUM(2)*TOTXS)
C
C  Calculate the uncertainty on CHANGE from equation 40 
C
            SCHANG=0.d0
            DO K=1,4
              DO J=1,4
                SCHANG=SCHANG+DDAADX(K)*DDAADX(J)*ESUM2(K,J)
              ENDDO
            ENDDO
            SCHANG=SQRT(SCHANG)
          ELSE
            ASY0=0.d0
            SASY0=0.d0
            DIFF=0.d0
            SDIFF=0.d0
            CHANGE=0.d0
            SCHANG=0.d0
          ENDIF
          PRINT 11500, tmax,asy0,sasy0,change*100.,schang*100.,
     >                diff*100.,sdiff*100
11500     FORMAT(/,5x,
     >    '***** Energy-weighted asymmetry, theta < ',f6.4,' rad *****',        
     >         /,/,1x,'asy0 = ',f6.4,'+-',f6.4,', (1-0)/0 = ',f6.4,'+-'        
     >        ,f6.4,' (%), (1-0)x100 = ',f6.4,'+-',f6.4)
        ENDIF
C
C  Do the unpolarized cross section first
C
        PRINT 12000
12000   FORMAT(/,10x,'***** Unpolarized Photon Cross Section *****',/)
        PRINT 4000
        DO I=1,NB
          EQ=(I-0.5)*PBIN
          IF(EGSUM(1,I).GT.0.D0) THEN
            CHANGE=EGSUM(3,I)/EGSUM(1,I)
C
C  The correct errors involve correlations between the various x-sections
C  (COMTN2 generates 4 weights/trial....correlations improve the precision
C   on many quantities)
C
            SCHANG=CHANGE*SQRT(EGSUM2(1,1,I)/EGSUM(1,I)**2
     >            +EGSUM2(3,3,I)/EGSUM(3,I)**2
     >            -2.d0*EGSUM2(1,3,I)/(EGSUM(1,I)*EGSUM(3,I)))
          ELSE
            CHANGE=0.d0
            SCHANG=0.d0
          ENDIF
          PRINT 5000, EQ,EGSUM(1,I),SQRT(EGSUM2(1,1,I)),
     >                change*100.,schang*100. 
        ENDDO
C
C  Next the longitudinal asymmetry of the photons
C
        IF(LONG) THEN
          PRINT 12500
12500     FORMAT(/,10x,
     >    '***** Longitudinal Photon Cross Section Asymmetry *****',/)        
          PRINT 6000      
          DO I=1,NB
            EQ=(I-0.5)*PBIN
            IF(EGSUM(1,I).GT.0.D0) THEN
C
C  The order-alpha**2 longitudinal asymmetry A(0) and its uncertainty
C
              ASY0=EGSUM(2,I)/EGSUM(1,I)
              SASY0=ASY0*SQRT(EGSUM2(1,1,I)/EGSUM(1,I)**2
     >              +EGSUM2(2,2,I)/EGSUM(2,I)**2
     >              -2.d0*EGSUM2(1,2,I)/(EGSUM(1,I)*EGSUM(2,I)))
              TOTXS=EGSUM(1,I)+EGSUM(3,I)
              DEN=EGSUM(1,I)*TOTXS
C
C  DIFF is the difference between EGSUM(2,I)/EGSUM(1,I) (ie A(0)) and the
C  fully corrected asymmetry (EGSUM(2,I)+EGSUM(4,I))/(EGSUM(1,I)+EGSUM(3,I)).
C
              DIFF=(EGSUM(1,I)*EGSUM(4,I)-EGSUM(2,I)*EGSUM(3,I))/DEN
C
C  Calculate the derivatives of DIFF given in equation 40 
C
              DDADX(1)=(-EGSUM(1,I)**2*EGSUM(4,I)
     >                +2.d0*EGSUM(1,I)*EGSUM(2,I)*EGSUM(3,I)
     >                +EGSUM(3,I)**2*EGSUM(2,I))/DEN**2
              DDADX(2)=-EGSUM(3,I)/DEN
              DDADX(3)=-(EGSUM(2,I)+EGSUM(4,I))/TOTXS**2
              DDADX(4)=1.d0/TOTXS
C
C  Calculate the uncertainty on DIFF from equation 40 
C
              SDIFF=0.d0
              DO K=1,4
                DO J=1,4
                  SDIFF=SDIFF+DDADX(K)*DDADX(J)*EGSUM2(K,J,I)
                ENDDO
              ENDDO
              SDIFF=SQRT(SDIFF)
C
C  CHANGE is the fractional difference, DIFF/A(0).
C
              CHANGE=DIFF/ASY0
C
C  Calculate the derivatives of CHANGE given in equation 40 
C
              DDAADX(1)=EGSUM(3,I)*(EGSUM(2,I)+EGSUM(4,I))
     >                 /(EGSUM(2,I)*TOTXS**2)        
              DDAADX(2)=-EGSUM(1,I)*EGSUM(4,I)/(EGSUM(2,I)**2*TOTXS)
              DDAADX(3)=-EGSUM(1,I)*(EGSUM(2,I)+EGSUM(4,I))
     >                 /(EGSUM(2,I)*TOTXS**2)
              DDAADX(4)=EGSUM(1,I)/(EGSUM(2,I)*TOTXS)
C
C  Calculate the uncertainty on CHANGE from equation 40 
C
              SCHANG=0.d0
              DO K=1,4
                DO J=1,4
                  SCHANG=SCHANG+DDAADX(K)*DDAADX(J)*EGSUM2(K,J,I)
                ENDDO
              ENDDO
              SCHANG=SQRT(SCHANG)
            ELSE
              ASY0=0.d0
              SASY0=0.d0
              DIFF=0.d0
              SDIFF=0.d0
              CHANGE=0.d0
              SCHANG=0.d0
            ENDIF
            PRINT 7000, EQ,asy0,sasy0,change*100.,schang*100.,
     >                  diff*100.,sdiff*100
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END
