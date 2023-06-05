c Remove the comments on the following lines to enable smearing of the
c two-body kinematics due to Fermi motion in the case of a nuclear target.
#define TARGET_FERMI_MOTION on
c added by RTJ 12-6-99
 
c     version 1.0    c. 1970 for FNAL E8
c
c             2.0
c
c          through   Descriptions lost in the mists of time.
c
c             5.0  
c
C             6.0    FROM YU ZHO (FNAL E621)
C
C             7.0    ADDED MULTIPLE MODE CAPABILITY,
C                    PRIMITIVE NEUTRAL HIT GENERATION,
C                    CONVERTED TO HBOOK (FNAL E756)
C
C             8.0    ALL ROUTINES TO DOUBLE PRECISION,
C                    PRIMITIVE RADIATIVE MODE RECONSTRUCTION
C      
C                    8.1   RADIATIVE CASCADE RECONSTRUCTION ADDED,
C                          "THREE PHOTON" RECONSTRUCTION ADDED,
C                          STRUCTURE (BUT NOT FUNCTION) FOR TRIGGER
C                          SIMULATION ADDED,
C                          LOSS TRACING ADDED.
c
c             9.0    CEBAF-RadPhi-phase1 (Design)
c                    Removed all charged particle capability.
c                    Installed hybrid-method photon shower generation.
c                    Installed cluster-finding pattern recognition.
c                    Installed multi-photon reconstruction.
c                    Installed photon beam (removed other options)
c                    Removed unstable recoil capability (require proton recoil)
c                    Converted to centimeters (From fathoms)
c                    reverted to single precision (for compatability)
c                 
c                    9.1   Structural modifications to eliminate
c                          useless statement labels and modularize
c                          code functions.
c                          Multiple shower generation options:
c                          MCTYPE=1, exact photon energy deposited
c                                    in exactly 1 block.
c                          MCTYPE=2, photon energy "smeared" as
c                                    a/sqrt(e) + b and deposited in one block.
c                          MCTYPE=3, hybrid monte-carlo, clusters looked up
c                                    from tables obtained from e852 data
c                          MCTYPE=4, clusters looked up from tables generated
c                                    by GEANT.
C                    9.2   MCTYPE=0, clusters as MCTYPE=4, writes events
c                                    for MCTYPE=5,no analysis
c                          MCTYPE=5, reads MCTYPE=4 events, does analysis.
c
c             10.0   TJNAF Hall-D version.
c                    Undid Version 9.x changes except structural changes
c                    Wrote c-driver 
c
c     integer getone
c     do j=1,1000000000
c        i=getone(x)
c        if(i.eq.1) stop
c     enddo
c     stop
c     end
      integer function getone(flags,emin,emax,flname) 
      implicit none
      integer  flags(*)   !control variables from driver
      real     emin,emax  !range of beam photon energies
      character*100 flname
      COMMON/PAWC/IA(100000)
      real    ia       
      include 'mcommon.inc'
      integer madc        !number of blocks with data
      integer kadc(2,625) !address of block with data and block data
      real    pbeam(3)    !beam 3-vector momentum
      integer MODE/54/       
      integer i           !loop index
      integer jtrig       !number of decays to generate
      real    accp        !accetance of a decay mode
      external brehm
      real     brehm
      integer  isave
      real     x
      logical  enter /.false./
c
c     return codes:
c
c     getone=0   normal event
c     getone=1   finished with all events
c     getone=2   new type of event 
c
      getone=0
      if(.not.enter) then
         enter=.true.
         if(flags(2).eq.1) then
            open(unit=20,file=flname,form='unformatted')
         endif
c
c        signal for first event of this type
c
         getone=2
         ntrig=-1
         CALL HLIMIT(100000)
c
c        make histograms used to pick values from distributions
c
         call hbook1(9998,'cth$',100,-1.,1.,0.)
         call hbfun1(9999,'Brem Spectrum$',100,emin,emax,brehm)
c
c        setup entire program run.
c     
         CALL START
         do i=1,999
            call mcread(jtrig)
            if(jtrig.ne.0) then
               if(flags(1).eq.0) PRINT 10,jtrig 
10             FORMAT('  START GENERATING EVENTS ASKED FOR',I10)
               CALL BEGMODE(I)
               CALL MCLIST (flags)
               isave=i+1
               go to 5
            endif
         enddo
      endif
 5    continue
c
c     have we finished this mode?
c
      if(ntrig.eq.jtrig) then
         CALL ENDMODE(I,ntrig)      
         if(imade(i).ne.0) then
            accp=100.*float(igot(i))/float(imade(i))
            print 99,accp
 99         format(' acceptance = ',f6.2,'%')
         endif
         NTRIG=-1  
         do i=isave,999
c
c           loop over modes until we find one asking for events
c
            call mcread(jtrig)
c
c           test for completion, do output.
c
            if(jtrig.eq.-99) then
               print*,' normal completion'
               call hrput(0,'mcany.hst','N') 
               getone=1
               return
            endif
            if(jtrig.ne.0) then
               if(flags(1).eq.0) PRINT 10,jtrig 
               CALL BEGMODE(I)
               CALL MCLIST(flags)
               getone=2
               isave=i+1
               go to 7
            endif
         enddo
      endif
 7    continue
c
c     start loop over production/decay modes
c
c
c     know how many decays in this mode to generate
c
C         
C     start loop over events
C         
      if (ntrig.lt.jtrig) then
         NTRIG=NTRIG+1
c
c        fill the event.
c
         CALL GETUM(I,madc,kadc,pbeam)
         if(mctype.ne.0) then !v9.2
c
c           compare monte-carlo event to reconstructed event
c 
            CALL MCUSER(I)
         endif !v9.2
         if(flags(4).ne.0) call mcdump
         if(flags(2).ne.0) call mcwrite(20)
      endif
      return
      END 
      SUBROUTINE BEGMODE(MODE)
      implicit none
      integer mode
C
C     INITIALIZE FOR THE CURRENT DECAY MODE.
C     PRINT OUT DESCRIPTION OF DECAY MODE
C
C
C     SET UP HISTOGRAMS.
C
      call hbook1(100*mode+71,'Smallest Photon Angle$',
     $              100,0.,0.3,0.)
      call hbook1(100*mode+73,'Largest Photon Angle$',
     $              100,0.,1.5,0.)
      call hbook1(100*mode+75,'Smallest Photon Energy$',
     $              100,0.,1.0,0.)
      call hbook1(100*mode+77,'Largest Photon Energy$',
     $              100,0.,5.0,0.)
      call hbook1(100*mode+79,'Charged Particle Angles$',
     $              100,0.,0.5,0.)
      call hbook1(100*mode+81,'Charged Particle Momenta$',
     $              100,0.,5.0,0.)
      call hbook1(100*mode+83,'Minimum Photon Opening angle$',
     $              100,0.,0.10,0.)
      call hbook1(100*mode+85,'beta,pions$',
     $              100,0.95,1.0,0.)
      call hbook1(100*mode+87,'beta,kaons$',
     $              100,0.95,1.0,0.) 
      call hbook1(100*mode+89,'d-beta,pions$',
     $              100,0.,0.5,0.)
      call hbook1(100*mode+91,'d-beta,kaons$',
     $              100,0.,0.5,0.)
      RETURN        
      END 
      SUBROUTINE ENDMODE(MODE,ntrig)
      implicit none
      integer mode
      integer ntrig
      integer idum
      real    cont(100),sum
      integer i,j,k
      real    ec(6)
      real    ev(6) /95.,90.,80.,75.,50.,25./
      real    ex(6)
C
C     SET UP FOR NEXT DECAY MODE, OUTPUT FOR PREVIOUS.
C
      call hbook1(100*mode+72,'Smallest Photon Angle$',
     $              100,0.,0.3,0.)
      call hbook1(100*mode+74,'Largest Photon Angle$',
     $              100,0.,1.5,0.)
      call hbook1(100*mode+76,'Smallest Photon Energy$',
     $              100,0.,1.,0.)
      call hbook1(100*mode+78,'Largest Photon Energy$',
     $              100,0.,5.,0.)
      call hbook1(100*mode+80,'Charged Particle Angles$',
     $              100,0.,0.5,0.)
      call hbook1(100*mode+82,'Charged Particle Momenta$',
     $              100,0.,5.,0.)
      call hbook1(100*mode+84,'Minimum Photon Opening angle$',
     $              100,0.,0.10,0.)
      call hbook1(100*mode+86,'beta,pions$',
     $              100,0.95,1.0,0.)
      call hbook1(100*mode+88,'beta,kaons$',
     $              100,0.95,1.0,0.)
      call hbook1(100*mode+90,'d-beta,pions$',
     $              100,0.,0.5,0.)
      call hbook1(100*mode+92,'d-beta,kaons$',
     $              100,0.,0.5,0.)
      do j=71,91,2
         call hunpak(100*mode+j,cont,' ',0)
         sum=0.
         do i=1,100
            sum=sum+cont(i)
         enddo
         if(sum.ge.0.) then
            cont(1)=cont(1)/sum
            do i=2,100
               cont(i)=cont(i-1)+cont(i)/sum
            enddo
            if((j.eq.71).or.(j.eq.75).or.
     $         (j.eq.91)) then
               do i=1,100
                  cont(i)=1.-cont(i)
               enddo
            endif
            call hpak(100*mode+j+1,cont)
            do k=1,6
               ec(k)=99999.
            enddo
            do i=1,100
               do k=1,6
                  if(abs(cont(i)-ev(k)).le.ec(k)) then
                     ec(k)=abs(cont(i)-ev(k))
                     call hix(100*mode+j+1,i,ex(k))
                     ex(k)=ex(k)*1000.
                  endif
               enddo
            enddo
c            write(10,1) 100*mode+j+1
c 1          format(1x,i5)
c            write(10,2) ex
c 2          format(6(1x,f5.0))
         endif
      enddo
      RETURN        
      END 
      SUBROUTINE MCUSER(MODE)
      implicit none
      include 'mcommon.inc'
      integer mode,i,j
      integer kgam,ipgam(50)
      real    anmin,ang
      real    anmax
      real    enmax,enmin
      real    beta
      real    dbeta
c
c     loop over particles, find the type you want
c
      kgam=0
      DO I=1,50  
         IF(ID(I).NE.0) THEN
            if((i.lt.jrecoil).and.(i.ne.1)) then
               if(XYZV(ITRKEND(I),3).eq.0.) then
c
c                 forward system final state particle
c
                  if(id(i).eq.21) then
c
c                    forward system final state photon
c
                     kgam=kgam+1
                     ipgam(kgam)=i
                  endif
               endif
            endif
         ENDIF      
      enddo
      if(kgam.ge.1) then
         anmax=-99999.
         anmin=99999.
         enmin=99999.
         enmax=-99999.
         do i=1,kgam
            ang=acos(pl(iph(i),3)/pl(iph(i),4))
            if(ang.le.anmin) anmin=ang
            if(ang.ge.anmax) anmax=ang
            if(pl(iph(i),4).le.enmin) enmin=pl(iph(i),4)
            if(pl(iph(i),4).ge.enmin) enmax=pl(iph(i),4)
         enddo
         call hf1(100*mode+75,enmin,1.)
         if(enmin.le.0.100) return
         call hf1(100*mode+71,anmin,1.)
         call hf1(100*mode+73,anmax,1.)
         call hf1(100*mode+77,enmax,1.)
      endif
      if(kgam.ge.2) then
         anmin=99999.
         do i=1,kgam-1
            do j=i+1,kgam
               ang=(pl(ipgam(i),1)*pl(ipgam(j),1)
     $             +pl(ipgam(i),2)*pl(ipgam(j),2)
     $             +pl(ipgam(i),3)*pl(ipgam(j),3))/
     $             (pl(ipgam(i),4)*pl(ipgam(j),4))
               ang=acos(ang)
               if(ang.le.anmin) anmin=ang
            enddo
         enddo
         call hf1(100*mode+83,anmin,1.)
      endif
      do i=1,50
         if(id(i).ne.0) then
            ang=(180./3.1415928)*acos(pl(i,3)/pl(i,5))
            if(ang.ge.30.) then
               if((id(i).eq.15).or.(id(i).eq.16)) then
c
c                 charged pions
c
                  beta=pl(i,5)/pl(i,4)
                  call hf1(100*mode+85,beta,1.)
                  dbeta=0.5*(rms(13)**2-rms(15)**2)/(pl(i,5)**2)
                  call hf1(100*mode+89,dbeta,1.)
               endif
               if((id(i).eq.13).or.(id(i).eq.14)) then
c
c                 charged kaons
c
                  beta=pl(i,5)/pl(i,4)
                  call hf1(100*mode+87,beta,1.)
                  dbeta=0.5*(rms(13)**2-rms(15)**2)/(pl(i,5)**2)
                  call hf1(100*mode+91,dbeta,1.)
               endif
            endif
         endif
      enddo
      RETURN 
      END 

      SUBROUTINE DECAY2(I1,I2,I3)       
C
C   DECAY OF PARTICLE I1 IN ITS C. M. S. INTO  PARTICLES I2 
C    AND I3 
c (written by the dreaded Bob Handler, henct the creative variable names)
C
      implicit none
      include 'mcommon.inc'
      real    A(5) /5*0./
      complex ROOT(4)
      real    C(4)
      real    D(4)  
      real    TWOPI/6.2831856/
      real    SQR5/2.236068/
      real    TOT/0./
      integer nid,nid1,nid2,nid3
      integer i1,i2,i3
      real    rm1sq,rm2sq,rm3sq,e2sq,e3sq
      real    dum,dumsq
      real    term1,psq,pcenm
      real    theta,costh,sinth
      real    phi,cosph,sinph
      real    y,z
      real    ranf
      external ranf
      NID1=ID(I1)   
      NID2=ID(I2)   
      NID3=ID(I3)
      RM1SQ=RMS(NID1)*RMS(NID1)         
      RM2SQ=RMS(NID2)*RMS(NID2)         
      RM3SQ=RMS(NID3)*RMS(NID3)         
      DUM=RM1SQ+RM2SQ-RM3SQ   
      DUMSQ=DUM*DUM 
      TERM1=DUMSQ/(4.*RM1SQ)  
      PSQ=TERM1-RM2SQ         
      PCENM=SQRT(PSQ)         
      PHI=TWOPI*RANf()        
      COSPH=COS(PHI)
      SINPH=SIN(PHI)
c
choose here (choose what Bob?)
c
      if((nid1.eq.31).or.(nid1.eq.27).or.(nid.eq.26)) then
c
c        have requested a phi,omega or rho decay.
c
c 76      continue
         costh=-1.+2.*ranf()
c         y=ranf()
c         z=0.5*(1+costh**2)
c         if(y.gt.z) go to 76
c         if(ranf().gt.0.5) costh=-costh
c         costh=exp(alog(ranf())/3.)
         call hf1(9998,costh,1.)
      else
 77      continue
         COSTH=-1.+2.*RANf()
         if(abs(costh).ge.1.) go to 77   
      endif  
      SINTH=SQRT(1.-COSTH*COSTH)        
      PCM(I2,3)=PCENM*COSTH   
      PCM(I2,1)=PCENM*SINTH*COSPH       
      PCM(I2,2)=PCENM*SINTH*SINPH       
      PCM(I3,1)=-PCM(I2,1)    
      PCM(I3,2)=-PCM(I2,2)    
      PCM(I3,3)=-PCM(I2,3)    
      E2SQ=PSQ+RM2SQ
      E3SQ=PSQ+RM3SQ
      PCM(I2,4)=SQRT(E2SQ)    
      PCM(I3,4)=SQRT(E3SQ)    
      PCM(I2,5)=PCENM         
      PCM(I3,5)=PCENM         
      RETURN        
      END 

      SUBROUTINE DECAY3(I1,I2,I3,I4)  
      implicit none
      include 'mcommon.inc'  
      REAL    M,M1,M2,M3,M12,M13MAX,M13MIN,M13,M23,TMP(3)      
      real    TWOPI/6.28318/ 
      real    qf,elow12,ehig12,delta12,e13min,elow,ehig
      integer i1,i2,i3,i4
      real    r11,r13,r21,r23,r31,r33 !a rotation matrix
      real    alf,sina,cosa,sa,ca     !angles for the rotation and sin,cosine
      real    phi,sp,cp               !angles for the rotation and sin,cosine
      real    st,ct                   !angles for the rotation and sin,cosine
      real    sinb,cosb               !angles for the rotation and sin,cosine
      real    p1,p2,p3                !cm momenta
      real    e1,e2,e3                !cm energies
      real    q1,q3
      real    delta
      real    ranf
      external ranf
      real    pip(5),pim(5),pi0(5)
      integer i
c
c     check for special cases
c
      if(id(i1).eq.25) then
c
c        eta decay.
c
         if((id(i2).eq.17)
     $     .and.(id(i2).eq.17)
     $     .and.(id(i3).eq.17)) then
c
c           eta -> 3pi0 decay
c
            call eta0_decay(pip,pim,pi0)
            do i=1,5
               pcm(i2,i)=pi0(i)
               pcm(i3,i)=pip(i)
               pcm(i4,i)=pim(i)
            enddo
            return
         else
c
c           eta -> pi+ pi- pi0 decay
c
            call eta_decay(pip,pim,pi0)
            if(id(i2).eq.17) then
               do i=1,5
                  pcm(i2,i)=pi0(i)
                  pcm(i3,i)=pip(i)
                  pcm(i4,i)=pim(i)
               enddo
            endif
            if(id(i3).eq.17) then
               do i=1,5
                  pcm(i3,i)=pi0(i)
                  pcm(i2,i)=pip(i)
                  pcm(i4,i)=pim(i)
               enddo
            endif
            if(id(i4).eq.17) then
               do i=1,5
                  pcm(i4,i)=pi0(i)
                  pcm(i3,i)=pip(i)
                  pcm(i2,i)=pim(i)
               enddo
            endif
            return
         endif
      endif
c
      if(id(i1).eq.27) then
c
c        omega -> pi+ pi- pi0 decay
c
         call omega_decay(pip,pim,pi0)
         if(id(i2).eq.17) then
            do i=1,5
               pcm(i2,i)=pi0(i)
               pcm(i3,i)=pip(i)
               pcm(i4,i)=pim(i)
            enddo
         endif
         if(id(i3).eq.17) then
            do i=1,5
               pcm(i3,i)=pi0(i)
               pcm(i2,i)=pip(i)
               pcm(i4,i)=pim(i)
            enddo
         endif
         if(id(i4).eq.17) then
            do i=1,5
               pcm(i4,i)=pi0(i)
               pcm(i3,i)=pip(i)
               pcm(i2,i)=pim(i)
            enddo
         endif
         return
      endif
      M =RMS(ID(I1))
      M1=RMS(ID(I2))
      M2=RMS(ID(I3))
      M3=RMS(ID(I4))

      QF=M1**2+M2**2+M3**2+M**2         
      ELOW12=(M1+M2)**2       
      EHIG12=(M-M3)**2        
      DELTA12=EHIG12-ELOW12   
      ELOW=(M2+M3)**2         
      EHIG=(M-M1)**2
      DELTA=EHIG-ELOW         
      E13MIN=(M1+M3)**2 
C         
C GET TWO BODY EFFECTIVE MASSES APPROXIMATLY INSIDE DALITZ BOUNDRY.   
1     CONTINUE 
      M12=ELOW12+RANf()*DELTA12         
      M23=ELOW+RANf()*DELTA
      if(m12*m23.le.0.) go to 1
      M13=QF-M12-M23
      IF(M13.LE.E13MIN) GO TO 1
C         
C DO DETAILED CHECK TO SEE IF THIS PAIR OF EFFECTIVE MASSES IS ALLOWED.         
      E3=(M**2-M12-M3**2)/(2.*SQRT(M12))
      E1=(M12+M1**2-M2**2)/(2.*SQRT(M12))         
      Q1=SQRT(E1**2-M1**2)    
      Q3=SQRT(E3**2-M3**2)    
      M13MAX=(E1+E3)**2-(Q1-Q3)**2      
      M13MIN=(E1+E3)**2-(Q1+Q3)**2      
c      PRINT 101,M13MIN,M13,M13MAX       
101   FORMAT(' MIN,CURRENT,MAX ',3(1X,F12.5))     
      IF((M13.GT.M13MAX).OR.(M13.LT.M13MIN)) GO TO 1
C         
C HAVE THE TWO BODY EFFECTIVE MASSES, GET THE ONE BODY ENERGIES.      
      E3=(M**2+M3**2-M12)/(2.*M)        
      E2=(M**2+M2**2-M13)/(2.*M)        
      E1=M-E2-E3    
C         
C COMPUTE THE MAGNITIUDE OF THE MOMENTUM OF EACH PARTICLE.  
      P3=SQRT(E3**2-M3**2)    
      P2=SQRT(E2**2-M2**2)    
      P1=SQRT(E1**2-M1**2)    
C         
C GET THE DECAY ANGLES.       
      COSB=(P1**2+P3**2-P2**2)/(2.*P1*P3)         
      if(abs(cosb).ge.1.) go to 1
      SINB=SQRT(1.-COSB**2)   
      COSA=(P1-P3*COSB)/P2    
      if(abs(cosa).ge.1.) go to 1
      SINA=SQRT(1.-COSA**2)   
C         
C RANDOMIZE THE ORIENTATION OF THE EVENT.         
      PHI=TWOPI*RANf()        
      CP=COS(PHI)   
      SP=SIN(PHI)   
      CT=2.*RANf()-1.         
      ST=SQRT(1.-CT**2)       
      ALF=TWOPI*RANf()        
      CA=COS(ALF)   
      SA=SIN(ALF)   
C         
      R11=CA*CT*CP-SA*SP      
      R13=CA*ST     
C         
      R21=-SA*CT*CP-CA*SP     
      R23=-SA*ST    
C         
      R31=-ST*CP    
      R33=CT        
C         
C ROTATE THE EVENT. 
      TMP(1)=P1     
      PCM(I2,1)=R11*TMP(1)    
      PCM(I2,2)=R21*TMP(1)    
      PCM(I2,3)=R31*TMP(1)    
      PCM(I2,4)=E1  
      PCM(I2,5)=P1  
C         
      TMP(1)=-P2*COSA         
      TMP(3)=-P2*SINA         
      PCM(I3,1)=R11*TMP(1)+R13*TMP(3)   
      PCM(I3,2)=R21*TMP(1)+R23*TMP(3)   
      PCM(I3,3)=R31*TMP(1)+R33*TMP(3)   
      PCM(I3,4)=E2  
      PCM(I3,5)=P2  
C         
      TMP(1)=-P3*COSB         
      TMP(3)=P3*SINB
      PCM(I4,1)=R11*TMP(1)+R13*TMP(3)   
      PCM(I4,2)=R21*TMP(1)+R23*TMP(3)   
      PCM(I4,3)=R31*TMP(1)+R33*TMP(3)   
      PCM(I4,4)=E3  
      PCM(I4,5)=P3  
      RETURN        
      END 

      SUBROUTINE EVENT(MODE,pbeam,ebeam)  
      implicit none
      include 'mcommon.inc'
      integer mode
      real    pbeam(*),ebeam
      real    pi /3.1415928/
      logical enter /.false./  
      real    PLK(5)
      real    rmss(60)        
      logical mflag(15)
      real    ptmp1(3),ptmp2(3),ptmp3(3)
      integer i,j,k,l,m,nv,nid
      real    eta,gbct,rdn,dum,dpos,dzpos
      real    xtt
      real    min_xtt_cut
      parameter (min_xtt_cut=0.5e-3)
      real    ranf
      external ranf
      logical  cata
      if(.not.enter) then
         enter=.true.
         do i=1,60
            rmss(i)=rms(i)
         enddo
      endif    
      CALL UZERO(PL,1,750)    
      CALL UZERO(XYZV,1,150)   
      CALL UZERO(NODEC,1,50)  
1     CONTINUE
 11   continue
c
c     choose the masses of the particles in the decay chains.
c
      do i=1,ntot
         if(width(id(i)).ne.0.) then
 12         continue
            k=id(i)
            rms(k)=rmss(k)
     $       +0.5*width(k)*tan(pi*(ranf()-0.5))
            if(abs(rms(k)-rmss(k)).ge.5.*width(k)) go to 12
         endif
      enddo
c
c     check for energy conservation in this set of mass choices
c
      do i=1,nvee
         if(idv(i).eq.2) then
            xtt=rms(id(lb(i,1)))-rms(id(lb(i,2)))-rms(id(lb(i,3)))
         else
            xtt=rms(id(lb(i,1)))-rms(id(lb(i,2)))-rms(id(lb(i,3)))
     $         -rms(id(lb(i,4)))
         endif
         if(xtt.lt.min_xtt_cut) go to 11
      enddo
C         
C     GENERATE PRIMARY PARTICLE AT THE TARGET     
c
      CALL PKAON(MODE,pbeam,ebeam,cata)
      if(cata) go to 11
      DO NV=1,NVEE
         I=LB(NV,1)
         J=LB(NV,2)
         K=LB(NV,3)
         L=LB(NV,4)
C
C        IF NODEC=1,NO DECAY IS GENERATED         
C
         IF(NODEC(I).EQ.1) THEN
            NODEC(J)=1        
            NODEC(K)=1        
            IF(IDV(NV).EQ.3) NODEC(L)=1 
         ELSE  
C
C           GENERATE THE DECAY DISTANCE FOR PARTICLE I      
C
            NID=ID(I)         
            ETA=PL(I,5)/RMS(NID)        
            GBCT=ETA*PTAU(NID)
            RDN=RANf()        
            DUM=-aLOG(RDN)    
            DPOS=GBCT*DUM     
            DZPOS=DPOS*PL(I,3)/PL(I,5)  
            XYZV(J,3)=XYZV(I,3)+DZPOS   
            XYZV(K,3)=XYZV(J,3)         
C         
C         
            XYZV(J,1)=XYZV(I,1)+DZPOS*PL(I,1)/PL(I,3)       
            XYZV(K,1)=XYZV(J,1)         
            XYZV(J,2)=XYZV(I,2)+DZPOS*PL(I,2)/PL(I,3)       
            XYZV(K,2)=XYZV(J,2)         
C           DECAY OF PARTICLE I
            IF(IDV(NV).EQ.2) then
               do m=1,3
                  ptmp1(m)=pl(i,m) 
                  ptmp2(m)=ptmp1(m)
               enddo
               call rotptoz(ptmp1,ptmp1,ptmp2)
               do m=1,3
                  pl(i,m)=ptmp2(m)
               enddo
               CALL DECAY2(I,J,K) 
               do m=1,3
                  ptmp2(m)=pl(j,m)
               enddo
               call rotztop(ptmp1,ptmp2,ptmp3)
               do m=1,3
                  pl(j,m)=ptmp3(m)
               enddo
               do m=1,3
                  ptmp2(m)=pl(k,m)
               enddo
               call rotztop(ptmp1,ptmp2,ptmp3)
               do m=1,3
                  pl(k,m)=ptmp3(m)
                  pl(i,m)=ptmp1(m)
               enddo
            endif
            IF(IDV(NV).EQ.3) CALL DECAY3(I,J,K,L) 
            CALL TRANS(I,J)
            CALL TRANS(I,K)   
            IF(L.NE.0) THEN   
               CALL TRANS(I,L)
               DO M=1,3    
                  XYZV(L,M)=XYZV(J,M)   
               enddo
            ENDIF   
        ENDIF      
      enddo   
      RETURN        
      END  

      SUBROUTINE GETUM(MODE,madc,kadc,pbeam)
c
c     generate final state, showers in the glass and reconstruct.
c
      implicit none
      include 'mcommon.inc'
      integer mode
c
c     returned variables
c
      integer madc           !number of ADC channels with data
      integer kadc(2,*)      !ADC data
      real    pbeam(*)
c
c     locals
c
      common/tracked/xfrnt,yfrnt,xback,yback,nfront,nback
      real    XFRNT(10)
      real    YFRNT(10)
      real    XBACK(10)
      real    YBACK(10)
      integer nfront
      integer nback
      real    ebeam
      real    vertex(3) /0.,0.,0./
      integer inoise
      integer ipois
      external ipois
c
      real    cc(625)     /625*0.001/
      real    ccgod(625)  /625*0.001/
      integer ipeds(625)  /625*16/
c
c     return variables
c
c
c     cluster finder return variables
c
      integer nclus
      integer nbclus(8,25)
      real    ebclus(8,25)
      integer mbclus(8)
      real    eclus(8)
      real    xclus(8)
      real    yclus(8)
      real    sigmac(8)
c
c     for event analysis
c
      real    etotc         !total energy in the clusters
      logical fidcut        !true if satisfied the fiducial cut
      real    efmass(28)    !effective masses of photon pairs
      real    rgamma(3,8)   !photon cluster positions
      real    pgamma(3,8)   !photon vector momenta
      real    eprime(8)     !fitted photon energies
      integer icase         !topology index
      integer nmp           !number of parents to 2 photons
      integer itype(4)      !type of the parents
      integer ipoint(4,2)   !pointers to the parents
      integer nomg          !number of omega to pi-gamma decays
      integer iomg(2,2)     !pointes for above
      real    chi6
      integer ichi(2)
      real    pi1(3)
      real    epi1
      real    pi2(3)
      real    epi2
      real    pi3(3)
      real    epi3
c
c     cluster finder input variables
c
      integer nblk
      real    eblk(500)
      integer iblk(500)
c
c     for kaon fitting
c
      real    pk(3),pf(3)
      real    r0(3) /0.,0.,0./
      real*8  sigma(4)
      real*8  chi8
      real    sf
      real    s0
      real    sigmaz /0./
      real    chis(5),st(5)
c
c     for loading showers
c
      real    eout(28,28)
c
      real    eedge,ehole,etotal,rmax,rtemp
c
      LOGICAL ENTER /.FALSE./
      logical nofit
      integer i,j           !loop indicies
      real    xgls,ygls     !scratch projected positions
      integer irow,icol,ich !block addresses
      real    sig           !the width of a distribution
      real    gausran       !gaussian random number generator
      real    xx,yy         !scratch hf2 arguments
      external gausran
      real    dummy         !name says it all.
      real    wt
      real    ang
      real    anmin,anmax,enmin,enmax
      IF(.NOT.ENTER) THEN
         enter=.true.
         call hbook1(30,'edge energy$',
     $               100,0.,2.,0.)
         call hbook1(40,'hole energy$',
     $               100,0.,2.,0.)
         call hbook1(50,'total energy$',
     $               100,0.,10.,0.)
         call hbook1(51,'Energy of missing photons MeV$',
     $               100,0.,500.,0.)
         call hbook1(52,'Angle of missing photons degrees$',
     $               100,25.,100.,0.)
         call hbook2(53,'Angle vs Energy$',
     $               100,0.,500.,100,25.,100.,0.)
         call hbook2(250,'Energy Deposition$',
     $        100,-70.,70.,100,-70.,70.,0.)
         if(mctype.eq.0) then !start v9.2
c
c           open output file
c
            open(unit=55,file='mctype0.dat',STATUS='new', 
     X       FORM='unformatted')      
         endif 
         if(mctype.eq.5) then
c
c           open input file
c
#ifndef LINUX
            open(unit=55,file='mctype0.dat',READONLY,STATUS='old', 
     X       FORM='unformatted')      
#else /* g77 does not support READONLY */
            open(unit=55,file='mctype0.dat',STATUS='old', 
     X       FORM='unformatted')      
#endif
         endif !end v9.2
      ENDIF
c
c     fill the decay tree.
c
      CALL EVENT(MODE,pbeam,ebeam)
      RETURN
      END

      SUBROUTINE LOREN4(DIR,P4IN,P4OUT)
c (another Bob Handler Special)
      implicit none
      real    DIR(4),P4IN(4),P4OUT(4)
      real    pcm2,epbeat,prod,onmcm,epbeta
      integer i
      PCM2=DIR(1)*DIR(1)+DIR(2)*DIR(2)+DIR(3)*DIR(3)       
      ONMCM=1./SQRT(DIR(4)*DIR(4)-PCM2)
      EPBETA=P4IN(1)*DIR(1)+P4IN(2)*DIR(2)+P4IN(3)*DIR(3)  
      PROD=EPBETA*(DIR(4)*ONMCM-1.)/PCM2-P4IN(4)*ONMCM     
      P4OUT(4)=ONMCM*(P4IN(4)*DIR(4)-EPBETA)     
      DO I=1,3  
         P4OUT(I)=P4IN(I)+DIR(I)*PROD   
      enddo  
      RETURN       
      END

      SUBROUTINE MCDUMP
      implicit none
      include 'mcommon.inc'      
      real    ZEND(50) 
      integer i,k,j,ii,l
      real    mass
      DO I=1,50
         ZEND(I)=999999.
      enddo
      ZEND(LB(1,1))=XYZV(LB(1,2),3)
      DO I=1,NVEE
         K=IDV(I)+1
         DO J=2,K
            II=LB(I,J)
            DO L=1,NVEE
               IF(L.NE.I) THEN
                  IF(NODEC(LB(L,2)).EQ.0) THEN
                     IF(II.EQ.LB(L,1)) ZEND(II)=XYZV(LB(L,2),3)      
                  ELSE
                     ZEND(II)=999999.
                  ENDIF
               ENDIF
            enddo
         enddo
      enddo
C
C     PRINTS OUT MONTE-CARLO AND DATA MOMENTA
C
      PRINT 2
2     FORMAT('  NAME',13x,'PX',11X,'PY',11X,'PZ',10X,' E')
      k=10000+ibeam
      print 99,(beamdump(j),j=1,4),k
 99   FORMAT('  Beam     ',4(1X,F12.5),1x,i6)
      DO I=1,50
         IF(ID(I).NE.0) THEN
            if(i.eq.jrecoil) then
               k=20000+id(i)
               mass=sqrt(pl(i,4)**2
     $                  -pl(i,1)**2
     $                  -pl(i,2)**2
     $                  -pl(i,3)**2 )
               if(XYZV(ITRKEND(I),3).eq.0) then
c
c                recoil particle is final state particle,
c                override the tag.
c
                  k=70000+id(i)
               endif
               print*,' recoil mass ',mass
            endif
            if(i.eq.1) then
               k=30000+id(i)
               if(XYZV(ITRKEND(I),3).eq.0) then
c
c                forward particle is final state particle,
c                override the tag.
c
                 k=50000+id(i)
               endif
            endif
            if((i.lt.jrecoil).and.(i.ne.1)) then
               if(XYZV(ITRKEND(I),3).ne.0.) then
                  k=40000+id(i)
               else
                  k=50000+id(i)
               endif
            endif
            if(i.gt.jrecoil) then
               if(XYZV(ITRKEND(I),3).ne.0.) then
                  k=60000+id(i)
               else
                  k=70000+id(i)
               endif
            endif
            PRINT 100,DNAME(ID(I)),(PL(I,J),J=1,4),k
         ENDIF      
      enddo
100   FORMAT(1x,a10,4(1X,F12.5),1x,i6)  
      DO I=1,NTOT
         PRINT 201,DNAME(ID(I)),XYZV(I,3),XYZV(ITRKEND(I),3)
      enddo
201   FORMAT(1X,A10,' TRACK FROM ',F9.5,' TO ',F9.5)
      RETURN        
      END

      SUBROUTINE MCwrite(lun)
      implicit none
      include 'mcommon.inc'      
      real    ZEND(50) 
      integer i,k,j,ii,l,lun   
C         
C     writes OUT MONTE-CARLO AND DATA MOMENTA     
C     
      write(lun) ntot+1
      k=10000+ibeam
      write(lun) (beamdump(j),j=1,4),k
c
c     flags as follows:   10,000 series:   Beam particle.
c                         20,000           Primary recoil particle
c                         30,000           Forward primary
c                         40,000           Forward intermediate
c                         50,000           Forward final
c                         60,000           Recoil intermediate
c                         70,000           Recoil final
c
      DO I=1,50  
         IF(ID(I).NE.0) THEN
            if(i.eq.jrecoil) then
               k=20000+id(i)
               if(XYZV(ITRKEND(I),3).eq.0.) then
                  k=70000+id(i)
               endif
            endif
            if(i.eq.1) then
               k=30000+id(i)
               if(XYZV(ITRKEND(I),3).eq.0.) then
                  k=50000+id(i)
               endif
            endif
            if((i.lt.jrecoil).and.(i.ne.1)) then
               if(XYZV(ITRKEND(I),3).eq.0.) then
                  k=50000+id(i)
               else
                  k=40000+id(i)
               endif
            endif
            if(i.gt.jrecoil) then
               if(XYZV(ITRKEND(I),3).eq.0.) then
                  k=70000+id(i)
               else
                  k=60000+id(i)
               endif
            endif
            write(lun) (PL(I,J),J=1,4),(XYZV(I,j),j=1,3),k
         ENDIF      
      enddo
      RETURN        
      END 

      SUBROUTINE MCgive(bret,pret,vret,kd,jflag,itot)
      implicit none
      include 'mcommon.inc'      
      real    ZEND(50) 
      integer i,k,j,ii,l,lun   
      integer jflag(50),kd(50)
      real    bret(4),pret(50,5),vret(50,3)
      integer itot,kk
C         
C     returns MONTE-CARLO MOMENTA etc.    
C     
c
c     flags as follows:   10,000 series:   Beam particle.
c                         20,000           Primary recoil particle
c                         30,000           Forward primary
c                         40,000           Forward intermediate
c                         50,000           Forward final
c                         60,000           Recoil intermediate
c                         70,000           Recoil final
c
      do i=1,4
         bret(i)=beamdump(i)
      enddo
      itot=0
      DO I=1,50 
         kd(i)=id(i) 
         do j=1,5
            pret(i,j)=pl(i,j)
         enddo
         do j=1,3
            vret(i,j)=xyzv(i,j)
         enddo
         IF(ID(I).NE.0) THEN
            itot=itot+1
            if(i.eq.jrecoil) then
               k=20000
               if(XYZV(ITRKEND(I),3).eq.0.) then
                  k=70000
               endif
            endif
            if(i.eq.1) then
               k=30000
               if(XYZV(ITRKEND(I),3).eq.0.) then
                  k=50000
               endif
            endif
            if((i.lt.jrecoil).and.(i.ne.1)) then
               if(XYZV(ITRKEND(I),3).eq.0.) then
                  k=50000
               else
                  k=40000
               endif
            endif
            if(i.gt.jrecoil) then
               if(XYZV(ITRKEND(I),3).eq.0.) then
                  k=70000
               else
                  k=60000
               endif
            endif
            jflag(i)=k
         else
            jflag(i)=0
         ENDIF      
      enddo
      RETURN        
      END 

      SUBROUTINE MCLIST(flags)
      implicit none
      include 'mcommon.inc'
      integer flags(*)
      integer i,j
      DATA DNAME/' Omega    ',' Cascade- ',' Cascade0 ',' Sigma-   ', 
     $           ' Sigma0   ',' Sigma+   ',' Lambda   ',' Neutron  ', 
     $           ' Proton   ',' P-BAR    ',' K-long   ',' K-short  ', 
     $           ' K-minus  ',' K-plus   ',' Pi-plus  ',' Pi-minus ', 
     $           ' Pi-zero  ',' Positron ',' Electron ',' Neutrino ',
     $           ' Gamma    ',' Mu-plus  ',' Mu-minus ',' x(2200)  ',
     $           ' Eta      ',' Rho(770) ',' Omega783 ',' Etap(958)',
     $           ' f0(975)  ',' a0(980)  ',' phi(1020)',' h1(1170) ',
     $           ' b1(1235) ',' a1(1260) ',' f2(1270) ',' f1(1285) ',
     $           ' eta(1295)',' a2(1320) ',' omg(1390)',' f1(1420) ',
     $           ' eta(1440)',' K*0(892) ',' K*0(1350)',' Delta+   ',
     $           ' Sigma1385',' Rho1700  ',' N(1440)  ',' N(1535)  ',
     $           ' Delta0   ',' x(2000)  ',' f0(500)  ', 
     $         9*'          '/
C         
C     PRINTS OUT DESCRIPTION OF CURRENT DECAY MODE.         
C     SETS UP POINTERS TO FIND END OF TRACKS.
C
      if(flags(1).eq.0) then
         print*,dname(ibeam),'+',dname(itgt),'==>',
     $     dname(id(1)),'+',dname(irecoil)
         print*,' Produced with t-dist as ',tcoef
         PRINT 1
      endif
1     FORMAT(' ***************************************')    
      DO I=1,50
         ITRKEND(I)=0
      enddo
      DO I=1,NVEE
         IF(LB(I,1).NE.0) THEN
            ITRKEND(LB(I,1))=LB(I,2)
            if(id(lb(i,1)).eq.irecoil) then
               if(flags(1).eq.0) print 1
            endif
            IF(IDV(I).EQ.2) THEN        
               if(flags(1).eq.0) PRINT 100,(DNAME(ID(LB(I,J))),J=1,3)         
            ELSE    
               if(flags(1).eq.0) PRINT 101,(DNAME(ID(LB(I,J))),J=1,4)         
            ENDIF   
         ENDIF      
      enddo    
100   FORMAT(1X,A10,' ==> ',A10,' + ',A10)        
101   FORMAT(1X,A10,' ==> ',A10,' + ',A10,' + ',A10)        
      if(flags(1).eq.0) PRINT 1       
      RETURN        
      END 

      SUBROUTINE MCREAD(jtrig)
      implicit none
      include 'mcommon.inc'
      integer jtrig
      integer mtot
C         
c
c   read the file describing the decays to be generated.
c
C   DESCRIPTION OF EVENT IS READ
c   then the number of events to be generated
c   then the coefficient in dN/dt=ext(alfa*t)
c   then idv, the number of particles at each vertex in the decay tree
c   then a list of the particle (by type) in the event
c   
c
      real    A(10)  
      integer i,j,m
C         
C  BELOW ARE REST MASSES (RMS), LIFETIMES (PTAU), AND CHARGES (PCH)    
C  OF THE PARTICLES CONTAINED IN THESE TABLES
C  ALL UNITS ARE IN GEVS AND centiMETERS 
c  
C    IDENTIFICATION IS MADE AS FOLLOWS  
C       1=OMEGA               9=PROTON             17=PIZERO
C       2=CASCADE MINUS      10=ANTI-PROTON        18=EPLUS 
C       3=CASCADE ZERO       11=KLONG              19=EMINUS
C       4=SIGMA MINUS        12=KSHORT             20=NEUTRINO        
C       5=SIGMA ZERO         13=KMINUS             21=GAMMA 
C       6=SIGMA PLUS         14=KPLUS              22=MUPLUS
C       7=LAMBDA             15=PIPLUS             23=MUMINUS         
C       8=NEUTRON            16=PIMINUS            24="spare"
c
c      25=eta                33=b1(1235)           41=eta(1440)
c      26=rho                34=a1(1260)           42=K*0(892)
c      27=omega              35=f2(1270)           43=K*0(1350)
c      28=eta'               36=f1(1285)           44=DELTA PLUS
c      29=f0(975)            37=eta(1295)          45=SIGMA(1385)
c      30=a0(980)            38=a2(1320)           46=rho(1700)
c      31=phi                39=omega(1390)        47=NSTAR(1440)
c      32=h1(1170)           40=f1(1420)           48=NSTAR(1535)
c
c      49=DELTA ZERO
c      50=X(2000)
c      51=f0(500)
c
      DATA RMS/1.67242  ,1.32132,1.3149  ,1.19734  ,1.19248  ,
     $         1.18936  ,1.1156 ,0.93957 ,0.9382796,0.9382796,
     $         0.4977   ,0.4977 ,0.49367 ,0.493668 ,0.1395669,  
     $         0.1395669,0.13496,0.000511,0.000511 ,0.       ,
     $         0.       ,0.10566,0.10566 ,2.200    ,0.54745  ,
     $         0.7681   ,0.78195,0.95775 ,0.9741   ,0.9827   ,
     $         1.019413 ,1.17   ,1.232   ,1.26     ,1.275    ,
     $         1.282    ,1.295  ,1.3184  ,1.391    ,1.4253   ,
     $         1.440    ,0.89645,1.425   ,1.2349   ,1.3872   ,
     $         1.691    ,1.440  ,1.535   ,1.2336   ,2.000    ,
     $         0.700    ,9*0./  
c   
      DATA PTAU/  2.5, 4.92, 8.69, 4.44,   0.,
     $            2.4, 7.89,  -1.,  -1.,  -1.,
     $          1554.,2.675,370.9,370.9,780.4,
     $          780.4,0.   , -1. , -1. , -1. ,
     $           -1. ,-1.  , -1. ,  3. ,36*0./   
c   
      DATA PCH/-1.,-1., 0.,-1., 0.,
     $          1., 0., 0., 1.,-1.,
     $          0., 0.,-1., 1., 1.,
     $         -1., 0., 1.,-1., 0.,
     $          0., 1.,-1., 0.,36*0./
c      data width/25*0.,.1515,0.,0.,.047,.057,0.,.36,.155,.4,
c     &.185,.024,0.035,0.110,.224,0.0553,0.06,0.05024,0.250,
c     & 0.1311, 0.0358, 0.215, 0.350, 0.150,12*0./
      data width/   0.,    0.,    0.,   0.,   0.,   0.,   0.,   0.,
     &              0.,    0.,    0.,   0.,   0.,   0.,   0.,   0.,
     &              0.,    0.,    0.,   0.,   0.,   0.,   0.,   0.,
c    &              0.,  .151,  .008,   0., .047, .057, .004,  .36,
     &              0.,  .151,  .008,   0., .047, .100, .004,  .36,
     &           0.155,   0.4,  0.85,0.024,0.035,0.110,0.224,.0553,
     &            0.06,.05024, 0.250,.1311,.0358,0.215,0.350,0.150,
     &            0.113,100.0, 0.500, 9*0.  /
C         
C  READ IN DESCRIPTION OF EVENTS TO BE GENERATED  
C         
      READ (18,10,END=9999)(A(I),I=1,10)         
   10 FORMAT(BZ,10A6)      
      read(18,*,end=9999) jtrig
      read(18,*,end=9999) tcoef
      read(18,*,end=9999) ibeam,itgt
 11   format(bz,5i4)
C  INITIALIZE BY READING IN PARAMETERS  
      READ (18,30,END=9999)(IDV(I),I=1,20)        
   30 FORMAT(BZ,20I4)
      READ (18,50,END=9999)(ID(I),I=1,50)        
   50 FORMAT(BZ,50I4)         
      NVEE=0        
      DO I=1,20   
         M=IDV(I)      
         IF(M.EQ.0) GO TO 72    
         READ (18,60,END=9999)(LB(I,J),J=1,4)       
   60    FORMAT(BZ,4I4)
         NVEE=NVEE+1   
      enddo 
   72 CONTINUE   
      NTOT=0        
      NCH=0         
      DO I=1,50 
         IF(ID(I).NE.0) THEN  
            NTOT=NTOT+1       
         ENDIF  
      enddo
      mtot=ntot
      irecvee=nvee+1
      READ (18,50,END=9999)(jD(I),I=1,50)  
      irecoil=jd(1)
      if(jd(2).ne.0) then
c
c        the recoil particle has a decay tree specified
c        Append its decay tree to the event tree.
c
         READ (18,30,END=9999)(jDV(I),I=1,20)     
         DO I=1,20
            M=jDV(I)
            IF(M.EQ.0) GO TO 82
            READ (18,60,END=9999)(mB(I,J),J=1,4)
         enddo 
   82    CONTINUE
         do i=1,50
            if(jd(i).ne.0) then
               ntot=ntot+1
               id(ntot)=jd(i)
            endif
         enddo
         do i=1,20
            if(jdv(i).ne.0) then
               nvee=nvee+1
               idv(nvee)=jdv(i)
               do j=1,idv(nvee)+1
                  lb(nvee,j)=mb(i,j)+mtot
               enddo
             endif
         enddo
      else
c
c        recoil particle is stable
c
         ntot=ntot+1
         id(ntot)=jd(1)
      endif
      jrecoil=mtot+1
c
c     initialization for later stuff.
c
      NTRIG=0
c
c     count photons and setup pointers
c
      MPHOTON=0     
      NEUTRON=0
      DO I=1,NTOT         
         IF(ID(I).EQ.21) THEN 
            MPHOTON=MPHOTON+1 
            IPH(MPHOTON)=I    
         ENDIF      
         IF(ID(I).EQ.8) THEN
            NEUTRON=I
         ENDIF
      enddo
      RETURN 
 9999 jtrig=-99
      return
      END 


      SUBROUTINE PKAON(MODE,pb,eb,cata) 
      implicit none
      include 'mcommon.inc'
      integer mode
      real    pbeam(4),psave(5),pin(4),pout(4),pdir(4)
      real    pb(3),eb
      real    rs
      real    PI /3.14159273/ 
      real    xxx
      real    rho
      real    phi
      real    t
      real    eneut
      real    xtt
      integer i,j
      real    ranf
      external ranf
      integer itried
      logical cata
      itried=0
      cata=.false.
c
c     get the momentum of the beam 'photon'.
c
      call tgtmom(pbeam)
c
c     setup initial point in space of the 'phi'.
c
 1    continue
c
c     gaussian distribution in (x,y)
c
      xxx=ranf()
      rho=sqrt(-2.*alog(xxx))
      phi=2.*ranf()*pi 
      XYZV(1,1)=xt+sigxt*rho*cos(phi)
      XYZV(1,2)=yt+sigyt*rho*sin(phi)
      if((xyzv(1,1).lt.xlow).or.(xyzv(1,1).gt.xhig).or.
     $   (xyzv(1,2).lt.ylow).or.(xyzv(1,2).gt.yhig)) go to 1
c
c     flat distribution in z
c
      XYZV(1,3)=zlow+ranf()*(zhig-zlow)
c
c     fill the vector momentum of the 'phi'.
c
 2    continue
      itried=itried+1
      if(itried.ge.1000) then
c
c        calling routine has setup a recoil or forward
c        particle that is below threashold, dump the event
c        and start over.
c
         cata=.true.
         return
      endif
      xxx=ranf()
      if(xxx.eq.0.) xxx=ranf()
      t=alog(xxx)/tcoef
      pl(1,4)=(t+rms(itgt)**2-rms(irecoil)**2+2.*pbeam(4)*rms(itgt))/
     $                 (2.*rms(itgt))
      eneut=pbeam(4)-pl(1,4)+rms(itgt)
c
c     demand physical t (physical recoil particle momentum)
c
      if(eneut.le.rms(irecoil)) go to 2
      pl(1,3)=(2.*pbeam(4)*pl(1,4)+t-rms(id(1))**2
     $     -rms(ibeam)**2)/(2.*pbeam(3))
      xtt=(pl(1,4)**2-rms(id(1))**2-pl(1,3)**2)
c
c     demand physical t (physical 'phi' momentum)
c
      if(xtt.le.0.) go to 2
      xtt=sqrt(xtt)
      phi=2.*pi*ranf()
      pl(1,2)=xtt*cos(phi)
      if(pl(1,2).ne.pl(1,2)) then
         print*,tcoef
         print*,'t,eta(4)',t,pl(1,4)
         print*,'eneut,eta(3)',eneut,pl(1,3)
         print*,'phi,etaperp',phi,xtt
      endif
      pl(1,1)=xtt*sin(phi)
      pl(1,5)=sqrt(pl(1,1)**2+pl(1,2)**2+pl(1,3)**2)
      prec(1)=pbeam(1)-pl(1,1)
      prec(2)=pbeam(2)-pl(1,2)
      prec(3)=pbeam(3)-pl(1,3)
      prec(4)=eneut
      prec(5)=sqrt(prec(1)**2+prec(2)**2+prec(3)**2)

#if defined TARGET_FERMI_MOTION
      call genFermi(pbeam,tcoef,rms(id(1)),rms(irecoil),psave,prec)
      do i=1,5
        pl(1,i)=psave(i) 
      enddo
#endif

      do i=1,5
         pl(jrecoil,i)=prec(i)
         if(i.le.3) then
            xyzv(jrecoil,i)=xyzv(1,i)
         endif
      enddo
      rdec(1,1)=xyzv(1,1)
      rdec(1,2)=xyzv(1,2)
      rdec(1,3)=xyzv(1,3)
      pb(1)=pbeam(1)
      pb(2)=pbeam(2)
      pb(3)=pbeam(3)
      eb=pbeam(4)
      return
      END 

      SUBROUTINE PROJEC(XI,YI,ZI,U,V,W,XO,YO,ZO)  
      implicit none
      real XI,YI,ZI,U,V,W,XO,YO,ZO
      XO = XI + (ZO-ZI)*U/W   
      YO = YI + (ZO-ZI)*V/W   
      RETURN        
      END 

      SUBROUTINE START
      implicit none
      include 'mcommon.inc'
      character*17 geover
C         
C     READ DATA CARDS  (told you this baby was old!)  
C
c     open input files.
c
c     detector,beam geometry and monte-carlo configuration file.
c
#ifndef LINUX
      open(unit=17,file='cebaf.geom',READONLY,STATUS='old',
     X FORM='formatted')      
#else /* g77 does not support READONLY */
      open(unit=17,file='cebaf.geom',STATUS='old',
     X FORM='formatted')      
#endif
c
c     list of requested/available decay trees.
c   
#ifndef LINUX
      open(unit=18,file='decayt.dat',READONLY,STATUS='old',
     X FORM='formatted')  
#else /* g77 does not support READONLY */
      open(unit=18,file='decayt.dat',STATUS='old',
     X FORM='formatted')  
#endif
c       
      read(17,3) geover
      read(17,*) xt,yt,zt
      read(17,*) sigxt,sigyt
      read(17,*) xlow,xhig
      read(17,*) ylow,yhig
      read(17,*) zlow,zhig
      read(17,*) jspike,pfix,sigmap
      read(17,*) zglass
      read(17,*) holex,holey
      read(17,*) xgmin,xgmax
      read(17,*) ygmin,ygmax
      read(17,*) pmin,pmax
      read(17,*) mctype,statterm,floorterm
      read(17,*) donoise
      read(17,*) barnoise
3     format(a17)
      return
      end

      SUBROUTINE TGTMOM(pbeam)
      implicit none
      include 'mcommon.inc'
      real    pbeam(4)
      real    p
      real    dummy
      real    ranf
      external ranf   
      real    gausran
      external gausran
      real    hrndm1
      external hrndm1
      IF(JSPIKE.EQ.0) THEN   
c
c        beam photon momentum flat between [pmin,pmax]
c 
         P=PMIN+RANf()*(PMAX-PMIN)
         pbeam(3)=p      
         pbeam(4)=sqrt(p**2+rms(ibeam)**2)
      ELSE IF(JSPIKE.EQ.1) THEN  
c
c        beam photon momentum gaussianly distributed about pfix
c       
         P=PFIX+sigmap*gausran(dummy)     
         pbeam(3)=p      
         pbeam(4)=sqrt(p**2+rms(ibeam)**2) 
      ELSE IF(JSPIKE.EQ.2) THEN   
c
c        beam photon momentum distributed by contents of histogram 9999
c      
         P=HRNDM1(9999,DUMMY) 
         pbeam(3)=p      
         pbeam(4)=sqrt(p**2+rms(ibeam)**2)
      ENDIF
      pbeam(1)=0.
      pbeam(2)=0.
      beamdump(1)=pbeam(1)
      beamdump(2)=pbeam(2)
      beamdump(3)=pbeam(3)
      beamdump(4)=pbeam(4)
      RETURN
      END 

      SUBROUTINE TRANS(I1,J1)
      implicit none
      include 'mcommon.inc'
      integer i1,j1
      integer kk,m,m1,m2,k1,k2,k3,k4,n
      real PDIR(4),PIN(4),POUT(4)     
C     USE LOREN4 TO TRANSFORM THE MOMENTUM OF PARTICLE J1 FROM        
C     THE REST FRAME OF PARTICLE I1 TO LAB. VIA VARIOUS REST FRAMES   
      DO KK=1,4   
         POUT(KK)=PCM(J1,KK)  
         PCM(1,KK)=PL(1,KK)
         pcm(lb(irecvee,1),kk)=pl(lb(irecvee,1),kk) 
      enddo
      DO M=1,NVEE
         IF(I1.EQ.LB(M,1)) GO TO 30     
      enddo
25    M=M2
30    CONTINUE 
      M1=LB(M,1)    
      DO K1=1,4  
         PDIR(K1)=-PCM(M1,K1) 
         PIN(K1)=POUT(K1)     
      enddo
      PDIR(4)=-PDIR(4)        
      CALL LOREN4(PDIR,PIN,POUT)        
      IF((M.EQ.1).or.(m.eq.irecvee)) then
         DO K4=1,4  
            PL(J1,K4)=POUT(K4)   
         enddo
         PL(J1,5)=SQRT(PL(J1,1)**2+PL(J1,2)**2+PL(J1,3)**2) 
         return
      endif   
      DO K2=1,M  
         N=IDV(K2)+1
         DO K3=2,N         
            IF(M1.EQ.LB(K2,K3)) GO TO 55
         enddo
      enddo
55    M2=K2         
      GO TO 25      
      END 

      real FUNCTION DOT(A,B)
      implicit none
      REAL A(3),B(3)     
      DOT=A(1)*B(1)+A(2)*B(2)+A(3)*B(3) 
      RETURN        
      END 

      real FUNCTION EFMASQ(I,J)    
      implicit none
      include 'mcommon.inc'
      integer i,j
      EFMASQ=(PL(I,4)+PL(J,4))**2       
     $      -(PL(I,3)+PL(J,3))**2       
     $      -(PL(I,2)+PL(J,2))**2       
     $      -(PL(I,1)+PL(J,1))**2       
      RETURN        
      END 

      real FUNCTION EMISSQ(I,J)    
      implicit none
      include 'mcommon.inc'
      integer i,j
      EMISSQ=(PL(I,4)-PL(J,4))**2       
     $      -(PL(I,3)-PL(J,3))**2       
     $      -(PL(I,2)-PL(J,2))**2       
     $      -(PL(I,1)-PL(J,1))**2       
      RETURN        
      END 

      real FUNCTION GAUSRAN(DUM) 
      implicit none
      integer igauss
      real    dum
      real    ranf
      external ranf
C         
C   THIS FUNCTION RETURNS A RANDOM NUMBER DISTRIBUTED       
C   WITH A GAUSSIAN WITH SIGMA = 1      
C   APOCRYPHAL.     
C         
 10   GAUSRAN=0.    
      DO IGAUSS=1,12       
         GAUSRAN=GAUSRAN+RANf()         
      enddo
      GAUSRAN=GAUSRAN-6.      
      RETURN        
      END 

      subroutine rotptoz(pp,pd,pout)
c
c     rotate 3-vector (pd) such that 3-vector pp is the new z-axis
c
      implicit none
c
c     input variables
c
      real              pp(3)   !specify new z-axis on old coordinate syatem
      real              pd(3)   !3-vector to be rotated
c
c     output variables
c
      real              pout(3) !rotated 3-vector
c
      real              ltf(3,3)
c
c     scratch 
c
      real              ppm,st,sp,ct,cp
c
c rotate p to z
c
      ppm=sqrt(pp(1)**2+pp(2)**2+pp(3)**2)
      ct=pp(3)/ppm
      st=sqrt(pp(1)**2+pp(2)**2)/ppm
      cp=pp(1)/(st*ppm)
      sp=pp(2)/(st*ppm)
c
      ltf(1,1)=ct*cp
      ltf(1,2)=ct*sp
      ltf(1,3)=-st
c
      ltf(2,1)=-sp
      ltf(2,2)=cp
      ltf(2,3)=0.
c
      ltf(3,1)=st*cp
      ltf(3,2)=st*sp
      ltf(3,3)=ct
c
      pout(1)=ltf(1,1)*pd(1)+ltf(1,2)*pd(2)+ltf(1,3)*pd(3)
      pout(2)=ltf(2,1)*pd(1)+ltf(2,2)*pd(2)+ltf(2,3)*pd(3)
      pout(3)=ltf(3,1)*pd(1)+ltf(3,2)*pd(2)+ltf(3,3)*pd(3)
      return
      end
      subroutine rotztop(pp,pd,pout)
c
c     inverse of rotptoz
c
      implicit none
c
c     input variables
c
      real              pp(3)
      real              pd(3)
c
c     output variables
c
c
      real              pout(3)
c
      real              ltf(3,3)
c
c     scratch 
c
      real              ppm,st,sp,ct,cp

c
c rotate z to p
c
      ppm=sqrt(pp(1)**2+pp(2)**2+pp(3)**2)
      ct=pp(3)/ppm
      st=sqrt(pp(1)**2+pp(2)**2)/ppm
      cp=pp(1)/(st*ppm)
      sp=pp(2)/(st*ppm)
c
      ltf(1,1)=ct*cp
      ltf(1,2)=-sp
      ltf(1,3)=st*cp
c
      ltf(2,1)=ct*sp
      ltf(2,2)=cp
      ltf(2,3)=st*sp
c
      ltf(3,1)=-st
      ltf(3,2)=0.
      ltf(3,3)=ct
c
      pout(1)=ltf(1,1)*pd(1)+ltf(1,2)*pd(2)+ltf(1,3)*pd(3)
      pout(2)=ltf(2,1)*pd(1)+ltf(2,2)*pd(2)+ltf(2,3)*pd(3)
      pout(3)=ltf(3,1)*pd(1)+ltf(3,2)*pd(2)+ltf(3,3)*pd(3)
      return
      end
c

      real function brehm(x)
      brehm=1./x
      return
      end
      SUBROUTINE omega_decay(PiPlus,PiMinus,PiZero)
      implicit none
c
c     omega -> pi+ pi- pi0, with proper dalitz plot distribution
c
c     output variables are:  PiPlus,PiMinus,PiZero
c     the momenta of the pions in the CM of the omega
c
      real    PiPlus(5),PiMinus(5),PiZero(5)
c
c     components 1-3 momentum, 4, energy, 5 magnitude of p
c
      real    z,zmax
      real    omax /4.2332921e-02/
      REAL    M,M1,M2,M3,M12,M13MAX,M13MIN,M13,M23,TMP(3)      
      real    TWOPI/6.28318/ 
      real    qf,elow12,ehig12,delta12,e13min,elow,ehig
      real    r11,r13,r21,r23,r31,r33 !a rotation matrix
      real    alf,sina,cosa,sa,ca     !angles for the rotation and sin,cosine
      real    phi,sp,cp               !angles for the rotation and sin,cosine
      real    st,ct                   !angles for the rotation and sin,cosine
      real    sinb,cosb               !angles for the rotation and sin,cosine
      real    p1,p2,p3                !cm momenta
      real    e1,e2,e3                !cm energies
      real    q1,q3
      real    delta
      real    ranf
      external ranf   
      logical enter /.false./
c
      if(.not.enter) then
         enter=.true.
         M =0.78195
         M1=0.1395679
         M2=0.1395679
         M3=0.1349743
         QF=M1**2+M2**2+M3**2+M**2       
         ELOW12=(M1+M2)**2
         EHIG12=(M-M3)**2
         DELTA12=EHIG12-ELOW12
         ELOW=(M2+M3)**2
         EHIG=(M-M1)**2
         DELTA=EHIG-ELOW
         E13MIN=(M1+M3)**2
      endif
C
C GET TWO BODY EFFECTIVE MASSES APPROXIMATLY INSIDE DALITZ BOUNDRY.
C
1     CONTINUE
      M12=ELOW12+RANf()*DELTA12
      M23=ELOW+RANf()*DELTA
      if(m12*m23.le.0.) go to 1
      M13=QF-M12-M23
      IF(M13.LE.E13MIN) GO TO 1
C
C DO DETAILED CHECK TO SEE IF THIS PAIR OF EFFECTIVE MASSES IS ALLOWED.
C
      E3=(M**2-M12-M3**2)/(2.*SQRT(M12))
      E1=(M12+M1**2-M2**2)/(2.*SQRT(M12))  
      Q1=SQRT(E1**2-M1**2)
      Q3=SQRT(E3**2-M3**2)
      M13MAX=(E1+E3)**2-(Q1-Q3)**2
      M13MIN=(E1+E3)**2-(Q1+Q3)**2
101   FORMAT(' MIN,CURRENT,MAX ',3(1X,F12.5))
      IF((M13.GT.M13MAX).OR.(M13.LT.M13MIN)) GO TO 1
C
C HAVE THE TWO BODY EFFECTIVE MASSES, GET THE ONE BODY ENERGIES.
C 
      E3=(M**2+M3**2-M12)/(2.*M)
      E2=(M**2+M2**2-M13)/(2.*M)
      E1=M-E2-E3
C
C COMPUTE THE MAGNITIUDE OF THE MOMENTUM OF EACH PARTICLE.
C
      P3=SQRT(E3**2-M3**2)
      P2=SQRT(E2**2-M2**2)
      P1=SQRT(E1**2-M1**2)
C
C GET THE DECAY ANGLES.
C
      COSB=(P1**2+P3**2-P2**2)/(2.*P1*P3)   
      if(abs(cosb).ge.1.) go to 1
      SINB=SQRT(1.-COSB**2)
      COSA=(P1-P3*COSB)/P2
      if(abs(cosa).ge.1.) go to 1
      SINA=SQRT(1.-COSA**2)
c
c     omega decay, take vector nature into account.
c
      z=(sina*p1*p2/omax)**2
      if(z.le.ranf()) go to 1
C
C RANDOMIZE THE ORIENTATION OF THE EVENT.
C
      PHI=TWOPI*RANf()
      CP=COS(PHI)
      SP=SIN(PHI)
      CT=2.*RANf()-1.   
      ST=SQRT(1.-CT**2)
      ALF=TWOPI*RANf()
      CA=COS(ALF)
      SA=SIN(ALF)
C
      R11=CA*CT*CP-SA*SP
      R13=CA*ST
C
      R21=-SA*CT*CP-CA*SP
      R23=-SA*ST
C
      R31=-ST*CP
      R33=CT
C
C ROTATE THE EVENT.
      TMP(1)=P1
c
      PiPlus(1)=R11*TMP(1)  
      PiPlus(2)=R21*TMP(1)
      PiPlus(3)=R31*TMP(1)
      PiPlus(4)=E1
      PiPlus(5)=P1
C
      TMP(1)=-P2*COSA
      TMP(3)=-P2*SINA
c
      PiMinus(1)=R11*TMP(1)+R13*TMP(3)
      PiMinus(2)=R21*TMP(1)+R23*TMP(3)
      PiMinus(3)=R31*TMP(1)+R33*TMP(3)
      PiMinus(4)=E2
      PiMinus(5)=P2
C
      TMP(1)=-P3*COSB
      TMP(3)=P3*SINB
      PiZero(1)=R11*TMP(1)+R13*TMP(3)
      PiZero(2)=R21*TMP(1)+R23*TMP(3)
      PiZero(3)=R31*TMP(1)+R33*TMP(3)
      PiZero(4)=E3
      PiZero(5)=P3
c
      RETURN
      END
      SUBROUTINE eta_decay(PiPlus,PiMinus,PiZero)
      implicit none
c
c     eta -> pi+ pi- pi0, with proper dalitz plot distribution
c
c     output variables are:  PiPlus,PiMinus,PiZero
c     the momenta of the pions in the CM of the eta
c
      real    PiPlus(5),PiMinus(5),PiZero(5)
c
c     components 1-3 momentum, 4, energy, 5 magnitude of p
c
      real    z,zmax,t,x,y
      REAL    M,M1,M2,M3,M12,M13MAX,M13MIN,M13,M23,TMP(3)      
      real    TWOPI/6.28318/ 
      real    qf,elow12,ehig12,delta12,e13min,elow,ehig
      real    r11,r13,r21,r23,r31,r33 !a rotation matrix
      real    alf,sina,cosa,sa,ca     !angles for the rotation and sin,cosine
      real    phi,sp,cp               !angles for the rotation and sin,cosine
      real    st,ct                   !angles for the rotation and sin,cosine
      real    sinb,cosb               !angles for the rotation and sin,cosine
      real    p1,p2,p3                !cm momenta
      real    e1,e2,e3                !cm energies
      real    q1,q3
      real    delta
      real    ranf
      external ranf   
      real    alpha33
      real    q33
      logical enter /.false./
c
      if(.not.enter) then
         enter=.true.
         M =0.54745
         M3=0.1395679
         M2=0.1395679
         M1=0.1349743
         QF=M1**2+M2**2+M3**2+M**2       
         ELOW12=(M1+M2)**2
         EHIG12=(M-M3)**2
         DELTA12=EHIG12-ELOW12
         ELOW=(M2+M3)**2
         EHIG=(M-M1)**2
         DELTA=EHIG-ELOW
         E13MIN=(M1+M3)**2
         q33=(2.+M1/M3)/(M-M1-2.*M3)
         alpha33=-0.40
      endif
C
C GET TWO BODY EFFECTIVE MASSES APPROXIMATLY INSIDE DALITZ BOUNDRY.
C
1     CONTINUE
      M12=ELOW12+RANf()*DELTA12
      M23=ELOW+RANf()*DELTA
      if(m12*m23.le.0.) go to 1
      M13=QF-M12-M23
      IF(M13.LE.E13MIN) GO TO 1
C
C DO DETAILED CHECK TO SEE IF THIS PAIR OF EFFECTIVE MASSES IS ALLOWED.
C
      E3=(M**2-M12-M3**2)/(2.*SQRT(M12))
      E1=(M12+M1**2-M2**2)/(2.*SQRT(M12))  
      Q1=SQRT(E1**2-M1**2)
      Q3=SQRT(E3**2-M3**2)
      M13MAX=(E1+E3)**2-(Q1-Q3)**2
      M13MIN=(E1+E3)**2-(Q1+Q3)**2
101   FORMAT(' MIN,CURRENT,MAX ',3(1X,F12.5))
      IF((M13.GT.M13MAX).OR.(M13.LT.M13MIN)) GO TO 1
C
C HAVE THE TWO BODY EFFECTIVE MASSES, GET THE ONE BODY ENERGIES.
C 
      E3=(M**2+M3**2-M12)/(2.*M)
      E2=(M**2+M2**2-M13)/(2.*M)
      E1=M-E2-E3
C
C COMPUTE THE MAGNITIUDE OF THE MOMENTUM OF EACH PARTICLE.
C
      P3=SQRT(E3**2-M3**2)
      P2=SQRT(E2**2-M2**2)
      P1=SQRT(E1**2-M1**2)
C
C GET THE DECAY ANGLES.
C
      COSB=(P1**2+P3**2-P2**2)/(2.*P1*P3)   
      if(abs(cosb).ge.1.) go to 1
      SINB=SQRT(1.-COSB**2)
      COSA=(P1-P3*COSB)/P2
      if(abs(cosa).ge.1.) go to 1
      SINA=SQRT(1.-COSA**2)
c
c     eta decay, take alpha parameter into account.
c
      t=e1-m1
      x=q33*t - 1.
      y=alpha33*x +1.
      if(y.le.ranf()*(abs(alpha33)+1.)) go to 1
C
C RANDOMIZE THE ORIENTATION OF THE EVENT.
C
      PHI=TWOPI*RANf()
      CP=COS(PHI)
      SP=SIN(PHI)
      CT=2.*RANf()-1.   
      ST=SQRT(1.-CT**2)
      ALF=TWOPI*RANf()
      CA=COS(ALF)
      SA=SIN(ALF)
C
      R11=CA*CT*CP-SA*SP
      R13=CA*ST
C
      R21=-SA*CT*CP-CA*SP
      R23=-SA*ST
C
      R31=-ST*CP
      R33=CT
C
C ROTATE THE EVENT.
      TMP(1)=P1
c
      PiZero(1)=R11*TMP(1)  
      PiZero(2)=R21*TMP(1)
      PiZero(3)=R31*TMP(1)
      PiZero(4)=E1
      PiZero(5)=P1
C
      TMP(1)=-P2*COSA
      TMP(3)=-P2*SINA
c
      PiMinus(1)=R11*TMP(1)+R13*TMP(3)
      PiMinus(2)=R21*TMP(1)+R23*TMP(3)
      PiMinus(3)=R31*TMP(1)+R33*TMP(3)
      PiMinus(4)=E2
      PiMinus(5)=P2
C
      TMP(1)=-P3*COSB
      TMP(3)=P3*SINB
      PiPlus(1)=R11*TMP(1)+R13*TMP(3)
      PiPlus(2)=R21*TMP(1)+R23*TMP(3)
      PiPlus(3)=R31*TMP(1)+R33*TMP(3)
      PiPlus(4)=E3
      PiPlus(5)=P3
c
      RETURN
      END

      SUBROUTINE eta0_decay(Pi1,Pi2,Pi3)
      implicit none
c
c     eta -> 3 pi0, with proper dalitz plot distribution
c
c     output variables are:  Pi1,Pi2,Pi3
c     the momenta of the pions in the CM of the eta
c
      real    Pi1(5),Pi2(5),Pi3(5)
c
c     components 1-3 momentum, 4, energy, 5 magnitude of p
c
      real    z,zmax,t,x,y
      REAL    M,M1,M2,M3,M12,M13MAX,M13MIN,M13,M23,TMP(3)      
      real    TWOPI/6.28318/ 
      real    qf,elow12,ehig12,delta12,e13min,elow,ehig
      real    r11,r13,r21,r23,r31,r33 !a rotation matrix
      real    alf,sina,cosa,sa,ca     !angles for the rotation and sin,cosine
      real    phi,sp,cp               !angles for the rotation and sin,cosine
      real    st,ct                   !angles for the rotation and sin,cosine
      real    sinb,cosb               !angles for the rotation and sin,cosine
      real    p1,p2,p3                !cm momenta
      real    e1,e2,e3                !cm energies
      real    q1,q3
      real    delta
      real    ranf
      external ranf   
      real    g
      real    q,qx,qy,del
      logical enter /.false./
c
      if(.not.enter) then
         enter=.true.
         M =0.54745
         M3=0.1349743
         M2=0.1349743
         M1=0.1349743
         QF=M1**2+M2**2+M3**2+M**2       
         ELOW12=(M1+M2)**2
         EHIG12=(M-M3)**2
         DELTA12=EHIG12-ELOW12
         ELOW=(M2+M3)**2
         EHIG=(M-M1)**2
         DELTA=EHIG-ELOW
         E13MIN=(M1+M3)**2
         q=(M-M1-2.*M3)
         qx=sqrt(3.)/(2.*0.54745*q)
         qy=     3. /(2.*0.54745*q)
         del=(0.54745 - 0.1349743)**2
c
c        the parameter describing the deviation from phase space
c
         g=-0.15
      endif
C
C GET TWO BODY EFFECTIVE MASSES APPROXIMATLY INSIDE DALITZ BOUNDRY.
C
1     CONTINUE
      M12=ELOW12+RANf()*DELTA12
      M23=ELOW+RANf()*DELTA
      if(m12*m23.le.0.) go to 1
      M13=QF-M12-M23
      IF(M13.LE.E13MIN) GO TO 1
C
C DO DETAILED CHECK TO SEE IF THIS PAIR OF EFFECTIVE MASSES IS ALLOWED.
C
      E3=(M**2-M12-M3**2)/(2.*SQRT(M12))
      E1=(M12+M1**2-M2**2)/(2.*SQRT(M12))  
      Q1=SQRT(E1**2-M1**2)
      Q3=SQRT(E3**2-M3**2)
      M13MAX=(E1+E3)**2-(Q1-Q3)**2
      M13MIN=(E1+E3)**2-(Q1+Q3)**2
101   FORMAT(' MIN,CURRENT,MAX ',3(1X,F12.5))
      IF((M13.GT.M13MAX).OR.(M13.LT.M13MIN)) GO TO 1
C
C HAVE THE TWO BODY EFFECTIVE MASSES, GET THE ONE BODY ENERGIES.
C 
      E3=(M**2+M3**2-M12)/(2.*M)
      E2=(M**2+M2**2-M13)/(2.*M)
      E1=M-E2-E3
C
C COMPUTE THE MAGNITIUDE OF THE MOMENTUM OF EACH PARTICLE.
C
      P3=SQRT(E3**2-M3**2)
      P2=SQRT(E2**2-M2**2)
      P1=SQRT(E1**2-M1**2)
C
C GET THE DECAY ANGLES.
C
      COSB=(P1**2+P3**2-P2**2)/(2.*P1*P3)   
      if(abs(cosb).ge.1.) go to 1
      SINB=SQRT(1.-COSB**2)
      COSA=(P1-P3*COSB)/P2
      if(abs(cosa).ge.1.) go to 1
      SINA=SQRT(1.-COSA**2)
c
c     eta decay, take the g-parameter into account.
c
      x=qx*(m12-m13)
      y=qy*(del-m23) -1.
      t=1.+g*(x*x+y*y)
      if(ranf().gt.t) go to 1
C
C RANDOMIZE THE ORIENTATION OF THE EVENT.
C
      PHI=TWOPI*RANf()
      CP=COS(PHI)
      SP=SIN(PHI)
      CT=2.*RANf()-1.   
      ST=SQRT(1.-CT**2)
      ALF=TWOPI*RANf()
      CA=COS(ALF)
      SA=SIN(ALF)
C
      R11=CA*CT*CP-SA*SP
      R13=CA*ST
C
      R21=-SA*CT*CP-CA*SP
      R23=-SA*ST
C
      R31=-ST*CP
      R33=CT
C
C ROTATE THE EVENT.
      TMP(1)=P1
c
      Pi3(1)=R11*TMP(1)  
      Pi3(2)=R21*TMP(1)
      Pi3(3)=R31*TMP(1)
      Pi3(4)=E1
      Pi3(5)=P1
C
      TMP(1)=-P2*COSA
      TMP(3)=-P2*SINA
c
      Pi2(1)=R11*TMP(1)+R13*TMP(3)
      Pi2(2)=R21*TMP(1)+R23*TMP(3)
      Pi2(3)=R31*TMP(1)+R33*TMP(3)
      Pi2(4)=E2
      Pi2(5)=P2
C
      TMP(1)=-P3*COSB
      TMP(3)=P3*SINB
      Pi1(1)=R11*TMP(1)+R13*TMP(3)
      Pi1(2)=R21*TMP(1)+R23*TMP(3)
      Pi1(3)=R31*TMP(1)+R33*TMP(3)
      Pi1(4)=E3
      Pi1(5)=P3
c
      RETURN
      END
