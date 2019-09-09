import sys
import ROOT
from ROOT import TFile
import math
ROOT.gStyle.SetOptStat(1111111)
ROOT.gStyle.SetOptFit(1111)

#ROOT.gStyle.SetStatX(0.9)
#ROOT.gStyle.SetStatY(0.7)


def drawText(text, x, y):        
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.06)
    latex.DrawLatex(x,y, "#font[42]{%s}"%text)

# Efficiency
def setupFit():
   fitFunc = ROOT.TF1("f1", "[0]*TMath::Erf((x-[1])/[2]) + [3]", 0, 1000)     #Error function scaled to [0,1]
   fitFunc.SetParNames("Norm.", "Edge", "Resolution", "Y-Offset")
   fitFunc.SetParameters(0.5, 5, 20, 0.5)
   fitParams = {}
   for x in xrange(0, 4):
     fitParams['var%s'%x] = fitFunc.GetParameter(x)
     fitParams['var%s'%x] = fitFunc.GetParError(x)
     print 'before fit', x, fitFunc.GetParameter(x), fitFunc.GetParError(x)
   return fitFunc
     
def makeEffPlot(passed,total):
   a = passed.Clone()
   b = total.Clone()
   eff = ROOT.TEfficiency(a,b)
   eff.SetTitle("Efficiency Plot")
   eff.SetMarkerStyle(10)
   eff.SetMarkerSize(1)
   eff.SetLineWidth(1)
   return eff
   
def outputFit():
   fitParams = {}
   for x in xrange(0, 4):
      fitParams['var%s'%x] = fitFunc.GetParameter(x)
      fitParams['var%s'%x] = fitFunc.GetParError(x)
      print x, fitFunc.GetParameter(x), fitFunc.GetParError(x)   
   return

#   hfile = TFile('jeioutfile.root', 'RECREATE', 'Output file')
#   hfile = TFile('jeioutfile.root', 'READ')
hfile = TFile(sys.argv[1], 'READ')

hist_dR_mu_jet = hfile.Get("hist_dR_mu_jet")
hist_dR_eg_jet = hfile.Get("hist_dR_eg_jet")
hist_dR2_mu_jet = hfile.Get("hist_dR2_mu_jet")
hist_dR2_eg_jet = hfile.Get("hist_dR2_eg_jet")
hist_dR2dR_mu_jet = hfile.Get("hist_dR2dR_mu_jet")
hist_dR2dR_eg_jet = hfile.Get("hist_dR2dR_eg_jet")

hist_dPhi_mu_jet = hfile.Get("hist_dPhi_mu_jet")
hist_dPhi_eg_jet = hfile.Get("hist_dPhi_eg_jet")
hist_dEta_mu_jet = hfile.Get("hist_dEta_mu_jet")
hist_dEta_eg_jet = hfile.Get("hist_dEta_eg_jet")

hist_dPhi_eg_jet_fine = hfile.Get("hist_dPhi_eg_jet_fine")
hist_dEta_eg_jet_fine = hfile.Get("hist_dEta_eg_jet_fine")
hist_dR2_eg_jet_fine = hfile.Get("hist_dR2_eg_jet_fine")

hist_mu_pt_all = hfile.Get("mu_pt_all")
hist_mu_pt_iso = hfile.Get("mu_pt_iso")
hist_eg_pt_all = hfile.Get("eg_pt_all")
hist_eg_pt_iso = hfile.Get("eg_pt_iso")
hist_jet_et_all = hfile.Get("jet_et_all")
hist_nJet0     = hfile.Get("nJet0")
hist_nJet0_mu  = hfile.Get("nJet0_mu")
hist_nJet0_eg  = hfile.Get("nJet0_eg")
hist_nMu0      = hfile.Get("nMu0")
hist_nEG0      = hfile.Get("nEG0")



#canv = ROOT.TCanvas("dR", "dR",  700, 500)
#canv = ROOT.TCanvas("dR", "dR",  1300, 900)
canv = ROOT.TCanvas("dR", "dR",  1300, 950)
pad = ROOT.TPad("pad", "", 0, 0, 1, 1)
pad.SetFillStyle(0)
pad.Divide(2,4)
pad.Draw()
pad.cd(1)
ROOT.gPad.SetLogy()
hist_dR_mu_jet.Draw()
pad.cd(3)
ROOT.gPad.SetLogy(0)
hist_dR_mu_jet.Draw()
ROOT.gPad.Update()

pad.cd(2)
ROOT.gPad.SetLogy()
hist_dR_eg_jet.Draw()
pad.cd(4)
ROOT.gPad.SetLogy(0)
hist_dR_eg_jet.Draw()

pad.cd(5)
hist_dPhi_mu_jet.Draw()
pad.cd(7)
hist_dEta_mu_jet.Draw()
pad.cd(6)
hist_dPhi_eg_jet.Draw()
pad.cd(8)
hist_dEta_eg_jet.Draw()



ROOT.gPad.Update()

canv2 = ROOT.TCanvas("iso", "iso", 1300, 950)
pad2 = ROOT.TPad("pad2", "", 0, 0, 1, 1)
pad2.SetFillStyle(0)
pad2.Divide(2,3)
pad2.Draw()
pad2.cd(1)
ROOT.gPad.SetLogy()
hist_mu_pt_all.Draw()
hist_mu_pt_iso.Draw("SAME")
hist_mu_pt_iso.SetLineColor(2)
pad2.cd(3)
ROOT.gPad.SetLogy(0)
hist_mu_pt_all.Draw()
hist_mu_pt_iso.Draw("SAME")
hist_mu_pt_iso.SetLineColor(2)
pad2.cd(2)
ROOT.gPad.SetLogy()
hist_eg_pt_all.Draw()
hist_eg_pt_iso.Draw("SAME")
hist_eg_pt_iso.SetLineColor(2)
pad2.cd(4)
ROOT.gPad.SetLogy(0)
hist_eg_pt_all.Draw()
hist_eg_pt_iso.Draw("SAME")
hist_eg_pt_iso.SetLineColor(2)
pad2.cd(5)
ROOT.gPad.SetLogy(0)
hist_jet_et_all.Draw()


canv3 = ROOT.TCanvas("nobj", "nobj", 1300, 950)
pad3 = ROOT.TPad("pad3", "", 0, 0, 1, 1)
pad3.SetFillStyle(0)
pad3.Divide(2,3)
pad3.Draw()
pad3.cd(1)
ROOT.gPad.SetLogy(0)
hist_nMu0.Draw()
pad3.cd(2)
hist_nEG0.Draw()
pad3.cd(3)
hist_nJet0.Draw()
pad3.cd(4)
hist_nJet0_mu.Draw()
pad3.cd(5)
hist_nJet0_eg.Draw()

ROOT.gPad.Update()


canv4 = ROOT.TCanvas("dR2", "dR2", 1300, 950)
pad4 = ROOT.TPad("pad4", "", 0, 0, 1, 1)
pad4.SetFillStyle(0)
pad4.Divide(2,2)
pad4.Draw()
pad4.cd(1)
ROOT.gPad.SetLogy(0)
hist_dR2_mu_jet.Draw()
pad4.cd(3)
hist_dR2dR_mu_jet.Draw('box')
pad4.cd(2)
hist_dR2_eg_jet.Draw()
pad4.cd(4)
hist_dR2dR_eg_jet.Draw('box')

canv5 = ROOT.TCanvas("eg delta fine", "eg delta fine", 1300, 950)
pad5 = ROOT.TPad("pad4", "", 0, 0, 1, 1)
pad5.SetFillStyle(0)
pad5.Divide(1,3)
pad5.Draw()
ROOT.gPad.SetLogy(0)
pad5.cd(1)
hist_dPhi_eg_jet_fine.Draw()
pad5.cd(2)
hist_dEta_eg_jet_fine.Draw()
pad5.cd(3)
hist_dR2_eg_jet_fine.Draw()

ROOT.gPad.Update()
