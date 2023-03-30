import numpy as np
import ROOT

histname = "h_pT"
#histname = "h_E"
#histname = "h_eta"
hist = "fevt/%s"%histname

t_Hgg = ROOT.TFile('~/eos/IMGs/HighLumiPileUp_ROOT/H125GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv2_FEVTDEBUG_nXXX_IMG.root', 'READ')
#t_Hgg = ROOT.TFile('~/eos/IMGs/HighLumiPileUp_ROOT/SingleHiggsPt10to80_Eta14_pythia8_HighLumiPileUpv2_FEVTDEBUG_nXXX_IMG.root', 'READ')
h_Hgg = ROOT.gDirectory.Get(hist)
print " >> Hgg entries:",h_Hgg.GetEntries()

#t_Ggg = ROOT.TFile('~/eos/IMGs/TEST_eta23/PromptDiPhoton_PtHat50_MGG80toInf_Pt25_Eta23_13TeV_TuneCUETP8M1_FEVTDEBUG_TEST_IMG.root', 'READ')
#t_Ggg = ROOT.TFile('~/eos/IMGs/TEST_eta23/PromptDiPhoton_PtHat10_MGG80toInf_Pt25_Eta23_13TeV_TuneCUETP8M1_FEVTDEBUG_TEST_IMG.root', 'READ')
#t_Ggg = ROOT.TFile('~/eos/IMGs/TEST_eta23/PromptDiPhoton_PtHat10_MGG80toInf_Pt25_Eta23_13TeV_TuneCUETP8M1_FEVTDEBUG_TEST_IMG.root', 'READ')
t_Ggg = ROOT.TFile('~/eos/IMGs/HighLumiPileUp_ROOT/PromptDiPhoton_PtHat45_MGG80toInf_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUp_FEVTDEBUG_nXXX_IMG.root', 'READ')
#t_Ggg = ROOT.TFile('~/eos/IMGs/HighLumiPileUp_ROOT/PromptDiPhoton_MGG80toInf_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUp_FEVTDEBUG_nXXX_IMG.root', 'READ')
#t_Ggg = ROOT.TFile('~/eos/IMGs/HighLumiPileUp_ROOT/DoublePhotonGaussPt55_StdDev20_HighLumiPileUpv3_FEVTDEBUG_nXXX_IMG.root', 'READ')
h_Ggg = ROOT.gDirectory.Get(hist)
print " >> gg entries:",h_Ggg.GetEntries()

c = ROOT.TCanvas("c","c",600,600)
c.SetBorderSize(0);
c.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
#ROOT.gPad.SetLogy()
c.cd()

ROOT.gPad.SetLeftMargin(0.15)

ymax = 0.18
h_Hgg.SetLineColor(2)
h_Hgg.GetYaxis().SetTitleOffset(2.)
h_Hgg.GetYaxis().SetTitle("N_{#gamma} / N_{#gamma,tot}")
h_Hgg.GetXaxis().SetTitleOffset(1.2)
h_Hgg.GetXaxis().SetTitle("p_{T} [GeV]")
h_Hgg.Scale(1./h_Hgg.Integral())
h_Hgg.Draw("")
h_Hgg.GetXaxis().SetRangeUser(0.,350.)
ROOT.gPad.Update()
h_Ggg.Scale(1./h_Ggg.Integral())

if histname == 'h_eta':
    ymax = 0.06
    h_Hgg.GetXaxis().SetTitle("#eta")
    h_Hgg.GetXaxis().SetRangeUser(-1.5,1.5)
h_Hgg.GetYaxis().SetRangeUser(0.,ymax)
h_Ggg.Draw("SAME")

leg = ROOT.TLegend(0.65,0.7,0.85,0.85)
leg.AddEntry(h_Hgg,"H#rightarrow#gamma#gamma","LP")
leg.AddEntry(h_Ggg,"#gamma#gamma","LP")
leg.SetBorderSize(0)
leg.Draw()

c.Print('%s.png'%histname)
