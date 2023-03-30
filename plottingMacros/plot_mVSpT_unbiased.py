#!/usr/bin/env python
import os
from glob import glob

import numpy as np

import ROOT
from ROOT import TCanvas, TPad, TFile, TPaveText, TLegend
from ROOT import gBenchmark, gStyle, gROOT, TStyle
from ROOT import TH1D, TF1, TGraphErrors, TMultiGraph

from math import sqrt

from array import array

import tdrstyle
tdrstyle.setTDRStyle()

import CMS_lumi

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_13TeV = '13 TeV'
CMS_lumi.writeExtraText = 1
#CMS_lumi.extraText = 'Preliminary'
CMS_lumi.extraText = 'Simulation'

iPos    = 0
iPeriod = 0

gStyle.SetOptFit(0)

def loadcanvas(name):
  canvas = TCanvas(name,name,400,20,1400,1000)
  canvas.SetFillColor(0)
  canvas.SetBorderMode(0)
  canvas.SetFrameFillStyle(0)
  canvas.SetFrameBorderMode(0)
  canvas.SetTickx(0)
  canvas.SetTicky(0)
  return canvas

def loadlegend(top, bottom, left, right):
  relPosX    = 0.001
  relPosY    = 0.005
  posX = 1 - right - relPosX*(1-left-right)
  posY = 1 - top - relPosY*(1-top-bottom)
  legendOffsetX = 0.0
  legendOffsetY = - 0.05
  textSize   = 0.05
  textFont   = 60
  legendSizeX = 0.4
  legendSizeY = 0.2
  legend = TLegend(posX-legendSizeX+legendOffsetX,posY-legendSizeY+legendOffsetY,posX+legendOffsetX,posY+legendOffsetY)
  legend.SetTextSize(textSize)
  legend.SetLineStyle(0)
  legend.SetBorderSize(0)
  return legend

histos={}
eosDir='root://cmsxrootd.fnal.gov//store/user/ddicroce/test'
#eosDir='/eos/uscms/store/user/ddicroce/test'

files_ = []
firstfile = True
filelist = '/uscms/home/ddicroce/nobackup/TauClassifier/CMSSW_10_2_20_UL/src/MLAnalyzer/list_HTauTau_unbiased.txt'
with open(filelist) as list_:
    content = list_.readlines()
paths = [x.strip() for x in content] 
print(paths)

for path in paths:
  #print(path)
  files_.append( TFile.Open(path) )
  #print(file)
  tmp_2d = files_[-1].Get('fevt/h_a_m_pT')
  tmp_m  = files_[-1].Get('fevt/h_jet_ma')
  tmp_pt = files_[-1].Get('fevt/h_jet_pta')
  if (firstfile):
    histos['mVSpT'] = tmp_2d.Clone('mVSpT')
    histos['mass']  = tmp_m.Clone('m')
    histos['pt']    = tmp_pt.Clone('pt')
    firstfile = False
  if not (firstfile):
    histos['mVSpT'].Add(tmp_2d)
    histos['mass'].Add(tmp_m)
    histos['pt'].Add(tmp_pt)

print (histos['mVSpT'].GetNbinsX())
print (histos['mVSpT'].GetNbinsY())
binx = []
biny = []
binz = []

histos['mVSpT_inverted'] = histos['mVSpT'].Clone('mVSpT_inverted')
binint = histos['mVSpT'].Integral()
binmax = 0
for iBinX in range(histos['mVSpT'].GetNbinsX()):
  #binx.append(histos['mVSpT'].GetXaxis().GetBinUpEdge(iBinX+1))
  for iBinY in range(histos['mVSpT'].GetNbinsY()):
    binz.append(histos['mVSpT'].GetBinContent(iBinX+1,iBinY+1))
    if (histos['mVSpT'].GetBinContent(iBinX+1,iBinY+1) > binmax):
      binmax = histos['mVSpT'].GetBinContent(iBinX+1,iBinY+1)
print(binmax)

histos['mVSpT_ratio'] = histos['mVSpT'].Clone('mVSpT_ratio')
for iBinX in range(histos['mVSpT_ratio'].GetNbinsX()):
  for iBinY in range(histos['mVSpT_ratio'].GetNbinsY()):
    if (histos['mVSpT'].GetBinContent(iBinX+1,iBinY+1) == 0): continue   #Avoids division by 0, in the case that the bin content was 0
    histos['mVSpT_ratio'].SetBinContent(iBinX+1, iBinY+1, ((1/binmax)*(histos['mVSpT'].GetBinContent(iBinX+1,iBinY+1))) )
    histos['mVSpT_inverted'].SetBinContent(iBinX+1, iBinY+1, (1/binmax)*(binint/histos['mVSpT'].GetBinContent(iBinX+1,iBinY+1)))

histos['mass_inverted'] = histos['mass'].Clone('mass_inverted')
massint = histos['mass'].Integral()
massmax = 0
for iBinX in range(histos['mass'].GetNbinsX()):
  binx.append(histos['mass'].GetXaxis().GetBinUpEdge(iBinX+1))
  if (massmax < massint/histos['mass'].GetBinContent(iBinX+1)): 
    massmax = massint/histos['mass'].GetBinContent(iBinX+1)
for iBinX in range(histos['mass_inverted'].GetNbinsX()):
  if (histos['mass'].GetBinContent(iBinX+1) == 0): continue  #Avoids division by 0, in the case that the bin content was 0
  histos['mass_inverted'].SetBinContent(iBinX+1, (1/massmax)*(massint/histos['mass'].GetBinContent(iBinX+1)))

histos['pt_inverted'] = histos['pt'].Clone('pt_inverted')
ptint = histos['pt'].Integral()
ptmax = 0
for iBinX in range(histos['pt'].GetNbinsX()):
  biny.append(histos['pt'].GetXaxis().GetBinUpEdge(iBinX+1))
  if (ptmax < ptint/histos['pt'].GetBinContent(iBinX+1)):
    ptmax = ptint/histos['pt'].GetBinContent(iBinX+1)
for iBinX in range(histos['pt_inverted'].GetNbinsX()):
  if (histos['pt'].GetBinContent(iBinX+1) == 0): continue  #Avoids division by 0, in the case that the bin content was 0
  histos['pt_inverted'].SetBinContent(iBinX+1, (1/ptmax)*(ptint/histos['pt'].GetBinContent(iBinX+1)))
  
print(binx)
print(biny)
print(binz)

canvas = loadcanvas("c1")
canvas.cd()
histos['mVSpT'].GetXaxis().SetTitle("m^{a} (GeV)")
histos['mVSpT'].GetYaxis().SetTitle("p_{T}^{a} (GeV)")
histos['mVSpT'].SetMinimum(0)
histos['mVSpT'].Draw('COLZ TEXT')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
canvas.SaveAs('a_massVspT_unbiased.root')
canvas.SaveAs('a_massVspT_unbiased.png')

canvas = loadcanvas("c2")
canvas.cd()
legend = loadlegend(canvas.GetTopMargin(), canvas.GetBottomMargin(), canvas.GetLeftMargin(), canvas.GetRightMargin())
histos['mass'].SetLineColor(2)
histos['mass'].SetLineWidth(3)
histos['mass'].SetXTitle("m^{a} (GeV)")
histos['mass'].SetYTitle("Jets")
histos['mass'].SetMinimum(0)
histos['mass'].Draw('COLZ TEXT')
legend.AddEntry(histos['mass'], 'Biased','lf')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
legend.Draw()
canvas.SaveAs('a_m_unbiased.root')
canvas.SaveAs('a_m_unbiased.png')

canvas = loadcanvas("c3")
canvas.cd()
legend = loadlegend(canvas.GetTopMargin(), canvas.GetBottomMargin(), canvas.GetLeftMargin(), canvas.GetRightMargin())
histos['pt'].SetLineColor(2)
histos['pt'].SetLineWidth(3)
histos['pt'].SetXTitle("p_{T}^{a} (GeV)")
histos['pt'].SetYTitle("Jets")
histos['pt'].SetMinimum(0)
histos['pt'].Draw('COLZ TEXT')
legend.AddEntry(histos['pt'], 'Biased','lf')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
legend.Draw()
canvas.SaveAs('a_pt_unbiased.root')
canvas.SaveAs('a_pt_unbiased.png')

canvas = loadcanvas("c4")
canvas.cd()
histos['mVSpT_ratio'].GetXaxis().SetTitle("m^{a} (GeV)")
histos['mVSpT_ratio'].GetYaxis().SetTitle("p_{T}^{a} (GeV)")
histos['mVSpT_ratio'].SetMinimum(0)
histos['mVSpT_ratio'].Draw('COLZ TEXT')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
canvas.SaveAs('a_massVspT_ratio_unbiased.root')
canvas.SaveAs('a_massVspT_ratio_unbiased.png')

canvas = loadcanvas("c5")
canvas.cd()
legend = loadlegend(canvas.GetTopMargin(), canvas.GetBottomMargin(), canvas.GetLeftMargin(), canvas.GetRightMargin())
histos['mass_inverted'].SetLineColor(2)
histos['mass_inverted'].SetLineWidth(3)
histos['mass_inverted'].SetXTitle("m^{a} (GeV)")
histos['mass_inverted'].SetYTitle("Jets")
histos['mass_inverted'].Draw('COLZ TEXT')
legend.AddEntry(histos['mass_inverted'], 'Biased','lf')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
legend.Draw()
canvas.SaveAs('a_m_inverted_unbiased.root')
canvas.SaveAs('a_m_inverted_unbiased.png')

canvas = loadcanvas("c6")
canvas.cd()
legend = loadlegend(canvas.GetTopMargin(), canvas.GetBottomMargin(), canvas.GetLeftMargin(), canvas.GetRightMargin())
histos['pt_inverted'].SetLineColor(2)
histos['pt_inverted'].SetLineWidth(3)
histos['pt_inverted'].SetXTitle("p_{T}^{a} (GeV)")
histos['pt_inverted'].SetYTitle("Jets")
histos['pt_inverted'].Draw('COLZ TEXT')
legend.AddEntry(histos['pt_inverted'], 'Biased','lf')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
legend.Draw()
canvas.SaveAs('a_pt_inverted_unbiased.root')
canvas.SaveAs('a_pt_inverted_unbiased.png')

canvas = loadcanvas("c7")
canvas.cd()
histos['mVSpT_inverted'].GetXaxis().SetTitle("m^{a} (GeV)")
histos['mVSpT_inverted'].GetYaxis().SetTitle("p_{T}^{a} (GeV)")
histos['mVSpT_inverted'].SetMinimum(0)
histos['mVSpT_inverted'].Draw('COLZ TEXT')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
canvas.SaveAs('a_massVspT_inverted_unbiased.root')
canvas.SaveAs('a_massVspT_inverted_unbiased.png')
