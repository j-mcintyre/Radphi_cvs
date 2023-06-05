geo. version 1.1
0. 0. 0.
0.26 0.26
-2.0 2.0
-2.0 2.0
0. 2.62
2  8. 0.01
100.
6. 6.
-38. 38.
-38. 38.
3.0 3.8
2  0.05 0.02
.false.
50.
------------Everything below this line is a comment--------RTJ

The format of the data records above is as follows, line-by-line.

geover      (char[]) : label indicating the version of this file
xt,yt,zt    (float)  : coordinates of center of target upstream face
sigxt,sigyt (float)  : rms gaussian beam spot dimensions on target
xlow,xhig   (float)  : cutoff limits for primary vertex position
ylow,yhig   (float)  : ditto 
zlow,zhig   (float)  : coordinates of upstream,downstream target faces
jspike,pfix,sigmap   : three cases
   jspike=0 (int)    : photon beam energy flat on [pmin,pmax] (see below)
   jspike=1 (int)    : photon beam energy distributed as a Gaussian with
                       pfix,sigmap (float) as mean and rms
   jspike=2 (int)    : photon beam energy distributed according to a
                       bremsstrahlung spectrum (histogram 9999)
zglass      (float)  : unused
holex,holey (float)  : unused
xgmin,xgmax (float)  : unused
ygmin,ygmax (float)  : unused
pmin,pmax   (float)  : limits of photon beam spectrum if jspike=0 (see above)
mctype,statterm,floorterm : obsolete (leave these values as they are)
donoise     (boolean): unused
barnoise    (float)  : unused
