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
eosDir='/eos/user/d/ddicroce/ML/TauClassifier/HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_unbiased/HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_unbiased_v5/201026_213424/0000'

firstfile = True

files_ = []
paths = ['%s'%(path) for path in glob('%s/*root'%(eosDir))]
for path in paths:
  #print(path)
  files_.append( TFile(path) )
  #print(infile)
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

binmax = 0
for iBinX in range(histos['mVSpT'].GetNbinsX()):
  for iBinY in range(histos['mVSpT'].GetNbinsY()):
    if (histos['mVSpT'].GetBinContent(iBinX+1,iBinY+1) > binmax):
      binmax = histos['mVSpT'].GetBinContent(iBinX+1,iBinY+1)
histos['mVSpT_ratio'] = histos['mVSpT'].Clone('mVSpT_ratio')
for iBinX in range(histos['mVSpT_ratio'].GetNbinsX()):
  for iBinY in range(histos['mVSpT_ratio'].GetNbinsY()):
    histos['mVSpT_ratio'].SetBinContent(iBinX+1, iBinY+1, ((1/binmax)*(histos['mVSpT'].GetBinContent(iBinX+1,iBinY+1))) )
  
canvas = loadcanvas("c1")
canvas.cd()
histos['mVSpT'].GetXaxis().SetTitle("m^{a} (GeV)")
histos['mVSpT'].GetYaxis().SetTitle("p_{T}^{a} (GeV)")
histos['mVSpT'].Draw('COLZ TEXT')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
canvas.SaveAs('a_massVspT_unbiased.png')

canvas = loadcanvas("c2")
canvas.cd()
legend = loadlegend(canvas.GetTopMargin(), canvas.GetBottomMargin(), canvas.GetLeftMargin(), canvas.GetRightMargin())
histos['mass'].SetLineColor(2)
histos['mass'].SetLineWidth(3)
histos['mass'].SetXTitle("m^{a} GeV")
histos['mass'].SetYTitle("Jets")
histos['mass'].Draw('COLZ TEXT')
legend.AddEntry(histos['mass'], 'Unbiased','lf')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
legend.Draw()
canvas.SaveAs('a_m_unbiased.png')

canvas = loadcanvas("c3")
canvas.cd()
legend = loadlegend(canvas.GetTopMargin(), canvas.GetBottomMargin(), canvas.GetLeftMargin(), canvas.GetRightMargin())
histos['pt'].SetLineColor(2)
histos['pt'].SetLineWidth(3)
histos['pt'].SetXTitle("p_{T}^{a} GeV")
histos['pt'].SetYTitle("Jets")
histos['pt'].Draw('COLZ TEXT')
legend.AddEntry(histos['pt'], 'Unbiased','lf')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
legend.Draw()
canvas.SaveAs('a_pt_unbiased.png')

canvas = loadcanvas("c4")
canvas.cd()
histos['mVSpT_ratio'].GetXaxis().SetTitle("m^{a} (GeV)")
histos['mVSpT_ratio'].GetYaxis().SetTitle("p_{T}^{a} (GeV)")
histos['mVSpT_ratio'].Draw('COLZ TEXT')
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.Update()
canvas.SaveAs('a_massVspT_ratio_unbiased.png')
