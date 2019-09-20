import ROOT
from ROOT import TFile
import math
ROOT.gStyle.SetOptStat(1111111)
ROOT.gStyle.SetOptFit(1111)

#ROOT.gStyle.SetStatX(0.9)
#ROOT.gStyle.SetStatY(0.7)

#path1 = "/afs/hephy.at/data/mzarucki02/nanoAOD/DegenerateStopAnalysis/postProcessing/processing_RunII_v5_1/nanoAOD_v5_1-0/Summer16_05Feb2018/met200_ht200_isr90/incLep/"
#path2 = "T2tt_mStop_500_mLSP_325/"
#path1 = "/afs/hephy.at/data/mzarucki02/nanoAOD/DegenerateStopAnalysis/postProcessing/processing_RunII_v6_6/nanoAOD_v6_6-1/Run2018_14Dec2018/inc/oneLepTight/"
#path2 = "MET_Run2018D_14Dec2018/"
#print "first file is ", path1+path2+"T2tt_mStop_500_mLSP_325.root"
#tree = ROOT.TChain("Events")
#tree.Add(path1+path2+"T2tt_mStop_500_mLSP_325.root")/afs/cern.ch/work/m/mjeitler/public/WorkingArea/CMSSW_9_4_13/src/Workspace/DegenerateStopAnalysis/plotsManfred//afs/cern.ch/work/m/mjeitler/public/WorkingArea/CMSSW_9_4_13/src/Workspace/DegenerateStopAnalysis/plotsManfred
#tree.Add(path1+path2+"MET_Run2018D_14Dec2018_9.root")
#filename = '/afs/cern.ch/work/m/mjeitler/public/WorkingArea/CMSSW_9_4_13/src/Workspace/DegenerateStopAnalysis/plotsManfred/T2tt_mStop_500_mLSP_325.root'
#f = ROOT.TFile("root://eoscms//eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/TEAshiftNtuples/ZeroBias2018E-week45-l1t-integration-v101p0-CMSSW-10_2_1/ZeroBias/crab_ZeroBias2018E-week45-l1t-integration-v101p0-CMSSW-10_2_1__325308_ZeroBias_Run2018E-v1/181105_165655/0000/L1Ntuple_8.root")
#filename = '/afs/cern.ch/work/m/mjeitler/public/WorkingArea/CMSSW_9_4_13/src/Workspace/DegenerateStopAnalysis/plotsManfred/L1Ntuple_8.root'
#f = ROOT.TFile(filename)

#tree = f.Get("l1UpgradeTree").Get("L1UpgradeTree")
#uGTTree = f.Get("l1uGTTree").Get("L1uGTTree")
#caloTree = f.Get("l1CaloTowerTree").Get("L1CaloTowerTree")
#eventTree = f.Get("l1EventTree").Get("L1EventTree")



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
   
def deltaPhi(phi1, phi2):
    pi = 3.141592654
    dphi = phi2-phi1
    if  dphi > pi:
        dphi -= 2.0*pi
    if dphi <= -pi:
        dphi += 2.0*pi
    return abs(dphi)   

hfile = TFile('outfile.root', 'RECREATE', 'Output file')

hist_dR_mu_jet = ROOT.TH1D("hist_dR_mu_jet", "hist_dR_mu_jet", 100, 0., 10.)
hist_dR_eg_jet = ROOT.TH1D("hist_dR_eg_jet", "hist_dR_eg_jet", 100, 0., 10.)
hist_dR2_mu_jet = ROOT.TH1D("hist_dR2_mu_jet", "hist_dR2_mu_jet", 100, 0., 100.)
hist_dR2_eg_jet = ROOT.TH1D("hist_dR2_eg_jet", "hist_dR2_eg_jet", 100, 0., 100.)
hist_dR2_eg_jet_fine = ROOT.TH1D("hist_dR2_eg_jet_fine", "hist_dR2_eg_jet_fine", 10000, 0., 1.)
hist_dR2dR_mu_jet = ROOT.TH2D("hist_dR2dR_mu_jet", "hist_dR2dR_mu_jet", 20, 0., 10., 20, 0., 100.)
hist_dR2dR_eg_jet = ROOT.TH2D("hist_dR2dR_eg_jet", "hist_dR2dR_eg_jet", 20, 0., 10., 20, 0., 100.)

hist_dPhi_mu_jet = ROOT.TH1D("hist_dPhi_mu_jet", "hist_dPhi_mu_jet", 70, 0., 3.15)
hist_dPhi_eg_jet = ROOT.TH1D("hist_dPhi_eg_jet", "hist_dPhi_eg_jet", 70, 0., 3.15)
hist_dPhi_eg_jet_fine = ROOT.TH1D("hist_dPhi_eg_jet_fine", "hist_dPhi_eg_jet_fine", 10000, 0., 1.)
hist_dEta_mu_jet = ROOT.TH1D("hist_dEta_mu_jet", "hist_dEta_mu_jet", 100, -7., 7.)
hist_dEta_eg_jet = ROOT.TH1D("hist_dEta_eg_jet", "hist_dEta_eg_jet", 100, -7., 7.)
hist_dEta_eg_jet_fine = ROOT.TH1D("hist_dEta_eg_jet_fine", "hist_dEta_eg_jet_fine", 10000, -1., 1.)
hist_dPhidEta_mu_jet = ROOT.TH2D("hist_dPhidEta_mu_jet", "hist_dPhidEta_mu_jet", 60, -4.,4., 20, 0., 3.15)
hist_dPhidEta_eg_jet = ROOT.TH2D("hist_dPhidEta_eg_jet", "hist_dPhidEta_eg_jet", 60, -4.,4., 20, 0., 3.15)

hist_mu_pt_all = ROOT.TH1D("mu_pt_all", "mu_pt_all", 200, 0, 100)
hist_mu_pt_iso = ROOT.TH1D("mu_pt_iso", "mu_pt_iso", 200, 0, 100)
hist_eg_pt_all = ROOT.TH1D("eg_pt_all", "eg_pt_all", 200, 0, 100)
hist_eg_pt_iso = ROOT.TH1D("eg_pt_iso", "eg_pt_iso", 200, 0, 100)
hist_jet_et_all = ROOT.TH1D("jet_et_all", "jet_et_all", 202, -10, 1000)
hist_nJet0     = ROOT.TH1D("nJet0", "nJet0", 20, 0., 20.)
hist_nJet0_mu  = ROOT.TH1D("nJet0_mu", "nJet0_mu", 20, 0., 20.)
hist_nJet0_eg  = ROOT.TH1D("nJet0_eg", "nJet0_eg", 20, 0., 20.)
hist_nMu0      = ROOT.TH1D("nMu0", "nMu0", 20, 0., 20.)
hist_nEG0      = ROOT.TH1D("nEG0", "nEG0", 20, 0., 20.)


tree1 = ROOT.TChain("l1UpgradeTree/L1UpgradeTree")
#for i in range(25,31):
#   print 'i=', i
#  324747_ZeroBias:
#   tree1.Add("root://eoscms//eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/TEAshiftNtuples/ZeroBias2018D-week43-l1t-integration-v101p0-CMSSW-10_2_1/ZeroBias/crab_ZeroBias2018D-week43-l1t-integration-v101p0-CMSSW-10_2_1__324747_ZeroBias_Run2018D-v1/181024_154920/0000/L1Ntuple_{}.root".format(i))
#   tree1.Add("root://eoscms//eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/TEAshiftNtuples/NanoDST2018D-week36-l1t-integration-v100p0-CMSSW-10_2_1/L1Accept/crab_NanoDST2018D-week36-l1t-integration-v100p0-CMSSW-10_2_1__322079_L1Accept_Run2018D-v1/180908_191159/0000/L1Ntuple_{}.root".format(i))
tree1.Add("/afs/cern.ch/work/m/mjeitler/public/WorkingArea/CMSSW_10_2_1/src/L1Ntuple18000.root")

print ' entries in tree1: ', tree1.GetEntries()
tree2 = ROOT.TChain("l1uGTTree/L1uGTTree")
#for i in range(25,31):
#   print 'i=', i
#  324747_ZeroBias:
#   tree2.Add("root://eoscms//eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/TEAshiftNtuples/ZeroBias2018D-week43-l1t-integration-v101p0-CMSSW-10_2_1/ZeroBias/crab_ZeroBias2018D-week43-l1t-integration-v101p0-CMSSW-10_2_1__324747_ZeroBias_Run2018D-v1/181024_154920/0000/L1Ntuple_{}.root".format(i))
#   tree2.Add("root://eoscms//eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/TEAshiftNtuples/NanoDST2018D-week36-l1t-integration-v100p0-CMSSW-10_2_1/L1Accept/crab_NanoDST2018D-week36-l1t-integration-v100p0-CMSSW-10_2_1__322079_L1Accept_Run2018D-v1/180908_191159/0000/L1Ntuple_{}.root".format(i))
tree2.Add("/afs/cern.ch/work/m/mjeitler/public/WorkingArea/CMSSW_10_2_1/src/L1Ntuple18000.root")

print ' entries in tree2: ', tree2.GetEntries()
tree3 = ROOT.TChain("l1EventTree/L1EventTree")
#for i in range(25,31):
#   print 'i=', i
#  324747_ZeroBias:
#   tree3.Add("root://eoscms//eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/TEAshiftNtuples/ZeroBias2018D-week43-l1t-integration-v101p0-CMSSW-10_2_1/ZeroBias/crab_ZeroBias2018D-week43-l1t-integration-v101p0-CMSSW-10_2_1__324747_ZeroBias_Run2018D-v1/181024_154920/0000/L1Ntuple_{}.root".format(i))
#   tree3.Add("root://eoscms//eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/TEAshiftNtuples/NanoDST2018D-week36-l1t-integration-v100p0-CMSSW-10_2_1/L1Accept/crab_NanoDST2018D-week36-l1t-integration-v100p0-CMSSW-10_2_1__322079_L1Accept_Run2018D-v1/180908_191159/0000/L1Ntuple_{}.root".format(i))
tree3.Add("/afs/cern.ch/work/m/mjeitler/public/WorkingArea/CMSSW_10_2_1/src/L1Ntuple18000.root")

print ' entries in tree3: ', tree3.GetEntries()


jcnt=0
trgcnt=0
#print 'test 1'
#for idx, (l1UpgradeTree/L1UpgradeTree), (l1uGTTree/L1uGTTree) in zip(range(tree1.GetEntries()), tree1, tree2):
#for idx, L1UpgradeTree, L1uGTTree in zip(range(tree1.GetEntries()), tree1, tree2):
#for idx, l1UpgradeTree.L1UpgradeTree, l1uGTTree.L1uGTTree in zip(range(tree1.GetEntries()), tree1, tree2):
#for idx, (l1UpgradeTree/L1UpgradeTree), (l1uGTTree/L1uGTTree) in zip(range(100), tree1, tree2):
for idx, tree1, tree2, tree3 in zip(range(tree1.GetEntries()), tree1, tree2, tree3):
  jcnt+=1
#  print 'test 2'
#  if jcnt > 1000: break   
#  (l1UpgradeTree/L1UpgradeTree).GetEntry(idx)
#  (l1uGTTree/L1uGTTree).GetEntry(idx)
#  L1UpgradeTree.GetEntry(idx)
#  L1uGTTree.GetEntry(idx)
  tree1.GetEntry(idx)
  tree2.GetEntry(idx)
  tree3.GetEntry(idx)
#  print idx, 'nMuons', l1UpgradeTree.L1UpgradeTree.L1Upgrade.nMuons
#  print idx, 'nMuons', L1Upgrade.nMuons
  if (jcnt < 10): 
    print idx, 'nEGs', tree1.L1Upgrade.nEGs, ' event=', tree3.Event.event, \
      bool(tree2.L1uGT.getAlgoDecisionInitial()[9])

#    for ii in range (0, 512):
#      if bool(uGTTree.L1uGT.getAlgoDecisionInitial()[ii]):
#        print 'Algo=', ii, ' is true'
 
# Mateusz degstop trigger must be fulfilled:
#  if not (bool(tree2.L1uGT.getAlgoDecisionInitial()[128])): 
#    continue
# Mu3 must be fulfilled:
#  if not (bool(tree2.L1uGT.getAlgoDecisionInitial()[9])): 
#    continue
# SingleJet90 must be fulfilled (there is no SingleJet100):
#  if not (bool(tree2.L1uGT.getAlgoDecisionInitial()[311])): 
#    continue
  trgcnt+=1
#  print 'Algo[9] = ', bool(uGTTree.L1uGT.getAlgoDecisionInitial()[9])


#jcnt=0
#for uGTiev in uGTTree:
#  jcnt+=1
#  if jcnt > 1: break
#  if (jcnt == 1):   print( dir(uGTiev.L1uGT) )
#  print jcnt, uGTiev.L1uGT.getL1MenuUUID()  #, uGTiev.L1uGT.m_orbitNr
#  for ii in range (0, 512):
#   print 'Algo=', ii, bool(uGTiev.L1uGT.getAlgoDecisionInitial()[ii])

# print 'jcnt==', jcnt

#icnt=0
#for iev in tree:
#  icnt=icnt+1
#  if icnt > 10: break
#  if icnt < 1:
#    print icnt, 'nEGs=', iev.L1Upgrade.nEGs, ' nJets=', iev.L1Upgrade.nJets, \
#    ' nMuons=', iev.L1Upgrade.nMuons
#  for j in range (0, iev.L1Upgrade.nEGs):
#    print '   j=', j, ' Et=', iev.L1Upgrade.egEt[j], ' Eta=', iev.L1Upgrade.egEta[j], \
#     ' Phi=', iev.L1Upgrade.egPhi[j]
#    print '          IEt=', iev.L1Upgrade.egIEt[j], ' IEta=', iev.L1Upgrade.egIEta[j], \
#      ' IPhi=', iev.L1Upgrade.egIPhi[j]
  imucnt=0
  for ii in range (0, tree1.L1Upgrade.nMuons): 
    mujetfound=0  
    if (tree1 .L1Upgrade.muonBx[ii] != 0):
      continue
    imucnt +=1
    ijetcnt=0        
#    for jj in range (0, min(1, tree1.L1Upgrade.nJets)):
    for jj in range (0, tree1.L1Upgrade.nJets):
      if (tree1.L1Upgrade.jetBx[jj] != 0):
        continue
      if (tree1.L1Upgrade.jetEt[jj] < 20.):
        continue
#      print  'jet no:', jj, ' jet Et =', tree1.L1Upgrade.jetEt[jj]
      if (ijetcnt > 0):
      	break	
      ijetcnt +=1
      dPhi = deltaPhi(tree1.L1Upgrade.muonPhi[ii], tree1.L1Upgrade.jetPhi[jj])
      dEta = tree1.L1Upgrade.muonEta[ii] - tree1.L1Upgrade.jetEta[jj]
      hist_dPhi_mu_jet.Fill(dPhi)
      hist_dEta_mu_jet.Fill(dEta)
      hist_dPhidEta_mu_jet.Fill(dEta, dPhi)
#      dR2 =  (tree1.L1Upgrade.muonPhi[ii] - tree1.L1Upgrade.jetPhi[jj])**2 + \
#             (tree1.L1Upgrade.muonEta[ii] - tree1.L1Upgrade.jetEta[jj])**2 
      dR2 =  dPhi**2 + dEta**2
      dR = math.sqrt(dR2)
#      dR = math.sqrt( (tree1.L1Upgrade.muonPhi[ii] - tree1.L1Upgrade.jetPhi[jj])**2 + \
#                      (tree1.L1Upgrade.muonEta[ii] - tree1.L1Upgrade.jetEta[jj])**2 )
      hist_dR_mu_jet.Fill(dR)
      hist_dR2_mu_jet.Fill(dR2)
      hist_dR2dR_mu_jet.Fill(dR, dR2)
#      dPhi = tree1.L1Upgrade.muonPhi[ii] - tree1.L1Upgrade.jetPhi[jj]

      if (dR < 0.3): 
        mujetfound +=1
#    print ' ijetcnt =', ijetcnt
    hist_nJet0_mu.Fill(ijetcnt)
#    if (iev.L1Upgrade.nMuons > 0):
    hist_mu_pt_all.Fill(tree1.L1Upgrade.muonEt[ii])
    if (mujetfound == 0):
      hist_mu_pt_iso.Fill(tree1.L1Upgrade.muonEt[ii])
#  print ' imucnt =', imucnt
  hist_nMu0.Fill(imucnt)
    		
  iegcnt=0		      
  for ii in range (0, tree1.L1Upgrade.nEGs): 
    egjetfound=0
    if (tree1.L1Upgrade.egBx[ii] != 0):
      continue
    iegcnt +=1
    ijetcnt=0        
#    for jj in range (0, min(2, tree1.L1Upgrade.nJets)):
    for jj in range (0, tree1.L1Upgrade.nJets):
      if (tree1.L1Upgrade.jetBx[jj] != 0):
        continue
      if (tree1.L1Upgrade.jetEt[jj] < 20.):
        continue
      if (ijetcnt > 0):
      	break		

      dPhi = deltaPhi(tree1.L1Upgrade.egPhi[ii], tree1.L1Upgrade.jetPhi[jj])
      dEta = tree1.L1Upgrade.egEta[ii] - tree1.L1Upgrade.jetEta[jj]
      dR2 =  dPhi**2 + dEta**2

      hist_dPhi_eg_jet_fine.Fill(dPhi)     
      hist_dEta_eg_jet_fine.Fill(dEta) 
      hist_dR2_eg_jet_fine.Fill(dR2)	
     
#      if (dR2 < 0.0001):
      if ( (dPhi < 0.1) and (abs(dEta) < 0.1) ):	
        continue
	
      ijetcnt +=1      
      dR = math.sqrt(dR2)
#      print 'dR2=', dR2, ' dR=', dR, ' dPhi=', dPhi, tree1.L1Upgrade.egPhi[ii], tree1.L1Upgrade.jetPhi[jj],\
#       ' dEta=', dEta, tree1.L1Upgrade.egEta[ii], tree1.L1Upgrade.jetEta[jj]
            
      hist_dPhi_eg_jet.Fill(dPhi)     
      hist_dEta_eg_jet.Fill(dEta) 
      hist_dPhidEta_eg_jet.Fill(dEta, dPhi)     
      hist_dR_eg_jet.Fill(dR)	
      hist_dR2_eg_jet.Fill(dR2)	
      hist_dR2dR_eg_jet.Fill(dR, dR2)	
	      
      if (dR < 0.3): 
        egjetfound +=1
#      break
    hist_nJet0_eg.Fill(ijetcnt)
    hist_eg_pt_all.Fill(tree1.L1Upgrade.egEt[ii])
    if (egjetfound == 0):
      hist_eg_pt_iso.Fill(tree1.L1Upgrade.egEt[ii])
  hist_nEG0.Fill(iegcnt)
 
 
  ijetcnt=0  
  maxjetet=-10
#  for jj in range (0, min(1, tree1.L1Upgrade.nJets)):
  for jj in range (0, tree1.L1Upgrade.nJets):
    if (tree1.L1Upgrade.jetBx[jj] != 0):
      continue
    ijetcnt +=1
    if (tree1.L1Upgrade.jetEt[jj] > maxjetet):
      maxjetet = tree1.L1Upgrade.jetEt[jj]
#    print 'jet: ', jj, ijetcnt, maxjetet

  hist_nJet0.Fill(ijetcnt)
  hist_jet_et_all.Fill(maxjetet)
  


print 'used events up to ', jcnt-1, ' and found ', trgcnt, ' triggers'

hfile.Write()

