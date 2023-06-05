mcwrap: Monte Carlo event generator for peripheral production reactions
author: Scott Teige (teige@dustbunny.physics.indiana.edu)
version: 1.0

This program was adapted for the Radphi experiment at Jefferson Lab
(e94-016) which detects all-neutral decays of mesons produced in the reaction

    gamma p  -->  p X  -->  p N(gamma)

The generator produces 4-vectors for the final-state particles under the
assumption of a chain of 2-body decays.  Intermediate states in the
decay chain are described by a Breit-Wigner lineshape.  The three-body
decays of the omega, eta and etaprime are also supported in the code.

The code was adapted for Radphi from a previous version used by E-852 at
the BNL-MPS.  The major addition was the incorporation of Fermi motion in
the initial state that is required when using a nuclear target.  The
new routine genFermi.F was introduced to the package, and it is called
from within fffnmc.F when the conditional TARGET_FERMI_MOTION is defined
at the head of the code.  The Fermi motion of the initial nucleon is
treated in the harmonic oscillator model with the rms momentum of bound
nucleons adjusted by the parameter pFermiRMS.  This parameter should be
adjusted to describe the target nucleus.

The following table lists the set of particles supported in this version
of mcwrap.  Baryons are capitalized and mesons are in lower case.

1. Omega-
2. Cascade-
3. Cascade0
4. Sigma-
5. Sigma0
6. Sigma+
7. Lamba
8. Neutron
9. Proton
10. antiProton
11. klong
12. kshort
13. k-
14. k+
15. pi+
16. pi-
17. pi0
18. positron
19. electron
20. neutrino
21. gamma
22. mu+
23. mu-
24. x(2200)
25. eta
26. rho0(770)
27. omega(783)
28. etaprime(958)
29. f0(975)
30. a0(980)
31. phi(1020)
32. h1(1170)
33. b1(1285)
34. a1(1260)
35. f2(1270)
36. f1(1285)
37. eta(1295)
38. a2(1320)
39. omega(1390)
40. f1(1420)
41. eta(1440)
42. k0*(892)
43. k0*(1350)
44. Delta+
45. Sigma(1385)
46. rho0(1700)
47. N*(1440)
48. N*(1535)
49. Delta0
50. x(2000)
51. f0(500)

The reaction to be simulated is described in an ascii file called decayt.dat
which consists of a sequence of sections of the following form.

===start of new reaction section in decayt.dat===
 Text title for this reaction - note the leading space
N                             - number of events of this type to generate
b                             - t-distribution is generated as exp(b*t)dt
   ibeam itarget              - particle id for beam and target
    k                         - number of unstable nodes in meson tree
   m1   m2   ...  mk          - multiplicities for unstable nodes 1..k
   f1   f2   f3   f4  ... fn  - particle id of all meson vertex tree nodes
   u1  d11  d12  ...          - node list descending from parent node u1
   u2  d21  d22  ...          - node list descending from parent node u2
   ...
   uk  dk1  dk2  ...          - node list descending from parent node uk
   f1   f2   f3   f4  ... fn  - particle id of all baryon vertex tree nodes 
   m1   m2  ...   mk          - number of unstable nodes in baryon tree
   u1  d11  d12  ...          - node list descending from parent node u1
   u2  d21  d22  ...          - node list descending from parent node u2
   ...
   uk  dk1  dk2  ...          - node list descending from parent node uk
   
The reaction is viewed as a tree-structured heirarchy of 2-body and 3-body
decays under two vertices: the meson vertex and the baryon vertex.  The
two trees are described in parallel structures, with one subtlety that the
multiplicity line follows the nodelist line for the baryons and preceeds
it for the mesons.  For leaf nodes no multiplicity entry is present.
