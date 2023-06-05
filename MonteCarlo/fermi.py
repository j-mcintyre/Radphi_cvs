import ROOT
import math
import re

hbarc = 0.197 # GeV

hmag = ROOT.TH1D("hmag", "model nucleon momentum distribution for 9Be",
                 200, 0, 2)
hcom = ROOT.TH1D("hcom", "model nucleon momentum distribution for 9Be",
                 200, -1, 1)

for line in open("fermi.txt"):
   p = line.rstrip().split()
   hcom.Fill(float(p[0]))
   hmag.Fill(float(p[4]))
hmag.GetXaxis().SetTitle("p magnitude (GeV/c)")
hmag.GetYaxis().SetTitle("counts")
hcom.GetXaxis().SetTitle("p component (GeV/c)")
hcom.GetYaxis().SetTitle("counts")
hmag.Draw()

table = []
for line in open("9Be_nucleon_momentum_dist.txt"):
   """
    9.6    3.939E-6  3.317E-7    3.815E-6  3.323E-7
   """
   if "spin" in line:
      break
   fline = line.rstrip().split()
   try:
      vals = [float(fline[i]) for i in range(0,5)]
   except:
      continue
   table.append([vals[i] for i in range(0,5)])

nx = len(table)
print "nx=",nx
x0 = hbarc * (table[0][0] - (table[1][0] - table[0][0])/2)
x1 = hbarc * (table[nx-1][0] + (table[nx-1][0] - table[nx-2][0])/2)
hmagp = ROOT.TH1D("hmagp", "theory proton momentum distribution for 9Be",
                  nx, x0, x1)
hmagn = ROOT.TH1D("hmagn", "theory neutron momentum distribution for 9Be",
                  nx, x0, x1)
for i in range(0, nx):
   p = table[i][0] * hbarc
   hmagp.SetBinContent(i+1, table[i][3] * p**2)
   hmagp.SetBinError(i+1, table[i][4] * p**2)
   hmagn.SetBinContent(i+1, table[i][1] * p**2)
   hmagn.SetBinError(i+1, table[i][2] * p**2)

hmagp.GetXaxis().SetTitle("p magnitude (GeV/c)")
hmagp.GetYaxis().SetTitle("counts")
hmagn.GetXaxis().SetTitle("p magnitude (GeV/c)")
hmagn.GetYaxis().SetTitle("counts")
hmagp.SetLineColor(ROOT.kRed)
hmagn.SetLineColor(ROOT.kBlue)
hmagp.Scale(hmag.GetMaximum()/hmagp.GetMaximum())
hmagn.Scale(hmag.GetMaximum()/hmagn.GetMaximum())

def gaus2(var, par):
   return (par[0] * math.exp(-0.5 * (var[0] / par[1])**2) * var[0]**2 +
           par[2] * math.exp(-0.5 * (var[0] / par[3])**2) * var[0]**2)

def fit_gaus2(h):
   fgaus2 = ROOT.TF1("fgaus2", gaus2, 0, 2, 4)
   fgaus2.SetParameter(0, h.GetMaximum())
   fgaus2.SetParameter(1, 0.10)
   fgaus2.SetParameter(2, h.GetMaximum() / 20)
   fgaus2.SetParameter(3, 0.20)
   fitres = h.Fit(fgaus2, "s")
   return fitres
