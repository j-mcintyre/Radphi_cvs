      common/accepta/imade(60),igot(60)
      integer imade
      integer igot
      COMMON/BMDEF/XT,YT,ZT,sigxt,sigyt,xlow,xhig,ylow,yhig,zlow,zhig
      real    XT
      real    YT
      real    ZT
      real    sigxt
      real    sigyt
      real    xlow
      real    xhig
      real    ylow
      real    yhig
      real    zlow
      real    zhig
      COMMON/bmtype/JSPIKE,PMIN,PMAX,PFIX,sigmap 
      integer jspike
      real    pmin
      real    pmax 
      real    pfix
      real    sigmap      
      common/glass/zglass,holex,holey,xgmax,xgmin,ygmax,ygmin
      real    zglass
      real    holex
      real    holey
      real    xgmax
      real    xgmin
      real    ygmax
      real    ygmin
      COMMON/LOCATE/XYZV(50,3)
      real    xyzv
      COMMON/MOMS/PL(50,5),PCM(50,5),PLD(50,5) 
      real    pl
      real    pcm
      real    pld   
      COMMON/NEUTRAL/MPHOTON,IPH(10),NEUTRON 
      integer mphoton
      integer iph
      integer neutron
      COMMON/NOISE/donoise,barnoise
      logical donoise
      real    barnoise
      COMMON/PARDAT/RMS(60),PTAU(60),PCH(60),width(60)      
      real    rms
      real    ptau
      real    pch
      real    width
      COMMON/PARNAM/DNAME(60) 
      CHARACTER*10 DNAME      
      COMMON/PARTID/IDV(20),IDC(50),ID(50),LB(20,4)
      integer idv
      integer idc
      integer id
      integer lb 
      COMMON/recID/jDV(20),jDC(50),jD(50),mB(20,4)
      integer jdv
      integer jdc
      integer jd
      integer mb 
      common/recoil/prec(5),rdec(20,3),irecoil,jrecoil,irecvee
      real    prec
      real    rdec
      integer irecoil
      integer jrecoil
      integer irecvee
      common/shtype/mctype,statterm,floorterm
      integer mctype
      real    statterm
      real    floorterm
      COMMON/TALLY/NTRIG,IFLAG,NVEE,NCH,NTOT,NODEC(50)
      integer ntrig
      integer iflag
      integer nvee
      integer nch
      integer ntot
      integer nodec      
      common/tco/tcoef
      real    tcoef
      COMMON/TRKEND/ITRKEND(50)
      integer itrkend

      common/initial/ibeam,itgt,beamdump(4)
      integer ibeam,itgt
      real    beamdump

