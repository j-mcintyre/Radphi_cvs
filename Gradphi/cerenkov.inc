      Integer nbins
      Parameter (nbins=51)
      Real Ephot(nbins)
      Real eff0(nbins),eff1(nbins),effPMT(nbins)
      Real rindexLGD(nbins),alenLGD(nbins)
      Real rindexWrap(nbins),alenWrap(nbins)
      Real rindexGap(nbins),alenGap(nbins)
      Real rindexAIR(nbins),alenAIR(nbins)
      Real rindexPMT(nbins),alenPMT(nbins)
      Common /Cerenkov/Ephot,eff0,eff1,effPMT,rindexLGD,alenLGD,
     +                 rindexWrap,alenWrap,rindexGap,alenGap,
     +                 rindexAIR,alenAIR,rindexPMT,alenPMT
      save /Cerenkov/
