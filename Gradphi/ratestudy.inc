      real Xupv,Yupv,Zupv,Nnupv,Ncupv,Enupv,Ecupv,Xnupv,Xcupv
      real Xcpv,Ycpv,Zcpv,Nncpv,Nccpv,Encpv,Eccpv,Xncpv,Xccpv
      real Xbsd,Ybsd,Zbsd,Nnbsd,Ncbsd,Enbsd,Ecbsd,Xnbsd,Xcbsd
      common /ntuple1/Nnupv,Ncupv,Enupv,Ecupv,Xnupv,Xcupv,
     +                Xupv,Yupv,Zupv,Tupv,
     +                Nncpv,Nccpv,Encpv,Eccpv,Xncpv,Xccpv,
     +                Xcpv,Ycpv,Zcpv,Tcpv,
     +                Nnbsd(3),Ncbsd(3),Enbsd(3),Ecbsd(3),
     +                Xnbsd(3),Xcbsd(3),Xbsd(3),Ybsd(3),Zbsd(3),Tbsd(3)
      integer ntvars(50)
      equivalence (ntvars(1),Nnupv)

      integer ntID
      parameter (ntID=10)
      character ntnames*350
      parameter (ntnames='Nnupv:R,Ncupv:R,'
     +                // 'Enupv:R,Ecupv:R,'
     +                // 'Xnupv:R,Xcupv:R,'
     +                // 'Xupv:R,Yupv:R,Zupv:R,Tupv:R,'
     +                // 'Nncpv:R,Nccpv:R,'
     +                // 'Encpv:R,Eccpv:R,'
     +                // 'Xncpv:R,Xccpv:R,'
     +                // 'Xcpv:R,Ycpv:R,Zcpv:R,Tcpv:R,'
     +                // 'Nnbsd(3):R,Ncbsd(3):R,'
     +                // 'Enbsd(3):R,Ecbsd(3):R,'
     +                // 'Xnbsd(3):R,Xcbsd(3):R,'
     +                // 'Xbsd(3):R,Ybsd(3):R,Zbsd(3):R,Tbsd(3):R')
