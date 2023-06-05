c
c  These parameters describe the background particle flux that is
c  superimposed on top of the standard generation.
c
c  author: Richard Jones
c  date: March 22, 2005
c
c  Notes:
c  -------
c  ebeam:	energy of the electron beam (GeV)
c  bgrate: 	the rate (1/ns) of beam photons above the threshold BCUTE
c          	for discrete bremsstrahlung to generate during the ADC gate;
c          	default value of 0 represents a background-free simulation.
c  beamorigin:	starting coordinates x,y,z (cm) for beam electrons which 
c		are incident on the bremsstrahlung radiator.
c  bggate:	time (ns) of the start (1) and end (2) of the ADC gate,
c		during which random background is simulated.

      real ebeam
      real bgrate
      real beamorigin(3)
      real bggate(2)
      common /backgrounds/ebeam,bgrate,beamorigin,bggate
