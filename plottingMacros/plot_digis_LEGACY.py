import ROOT

sw = ROOT.TStopwatch()
sw.Start()

nDigis = 10

# Init histos
h_adc_pho = []
h_adc_ele = []
#xMax = [312,312,312,312,1896,2292,2160,1764,1236,972]
xMax = [312,312,312,312,2292,2292,2292,2292,2292,2292]
for i in range(nDigis):
  h_label = "h_adc_pho%d"%i
  h_adc_pho.append(ROOT.TH1F(h_label,h_label,132,180.,xMax[i]))
  h_adc_pho[-1].Sumw2()
  h_label = "h_adc_ele%d"%i
  h_adc_ele.append(ROOT.TH1F(h_label,h_label,132,180.,xMax[i]))
  h_adc_ele[-1].Sumw2()

# Read into histograms
print " >> Doing photons..."
chain_in_pho = ROOT.TChain("fevt/RHTree")
chain_in_pho.Add("histo_pho_.root")
for i in range(nDigis):
  branch_label = "EB_adc%d >> h_adc_pho%d"%(i,i)
  cut_label = "EB_adc%d > 0."%(i)
  chain_in_pho.Draw(branch_label,cut_label)
print " >> Doing electrons..."
chain_in_ele = ROOT.TChain("fevt/RHTree")
chain_in_ele.Add("histo_ele_.root")
for i in range(nDigis):
  branch_label = "EB_adc%d >> h_adc_ele%d"%(i,i)
  cut_label = "EB_adc%d > 0."%(i)
  chain_in_ele.Draw(branch_label,cut_label)

# Draw
c_adc = []
pho_color = 4
ele_color = 2
for i in range(nDigis):
  c_label = "c_adc%d"%(i)
  c_adc.append(ROOT.TCanvas(c_label,c_label,600,400))
  ROOT.gPad.SetLogy(1)
  h_adc_pho[i].SetLineColor(pho_color)
  h_adc_pho[i].Draw("HIST")
  h_adc_ele[i].SetLineColor(ele_color)
  h_adc_ele[i].Draw("HIST SAME")
  c_adc[-1].Print("HISTS/%s.png"%(c_label))
  c_adc[-1].Print("HISTS/%s.eps"%(c_label))

sw.Stop()
print " >> Real time:",sw.RealTime()/60.
