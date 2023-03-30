import ROOT

sw = ROOT.TStopwatch()
sw.Start()

##### photons #####


h_E_pho = ROOT.TH1F("h_E_pho","h_E_pho",100,0.,100.)
h_E_pho.Sumw2()
h_t_pho = ROOT.TH1F("h_t_pho","h_t_pho",150,-150.,150.)
h_t_pho.Sumw2()
h_Er_pho = ROOT.TH1F("h_Er_pho","h_Er_pho",100,0.,100.)
h_Er_pho.Sumw2()
h_tr_pho = ROOT.TH1F("h_tr_pho","h_tr_pho",150,-150.,150.)
h_tr_pho.Sumw2()

h_E_ele = ROOT.TH1F("h_E_ele","h_E_ele",100,0.,100.)
h_E_ele.Sumw2()
h_t_ele = ROOT.TH1F("h_t_ele","h_t_ele",150,-150.,150.)
h_t_ele.Sumw2()
h_Er_ele = ROOT.TH1F("h_Er_ele","h_Er_ele",100,0.,100.)
h_Er_ele.Sumw2()
h_tr_ele = ROOT.TH1F("h_tr_ele","h_tr_ele",150,-150.,150.)
h_tr_ele.Sumw2()

print " >> Doing photons..."
chain_in_pho = ROOT.TChain("fevt/RHTree")
chain_in_pho.Add("histo_Pho.root")
chain_in_pho.Draw("EBenergy >> h_E_pho","EBenergy > 0.")
chain_in_pho.Draw("EBenergyRed >> h_Er_pho","EBenergyRed > 0.")
chain_in_pho.Draw("EBtime >> h_t_pho","EBtime != 0.")
chain_in_pho.Draw("EBtimeRed >> h_tr_pho","EBtimeRed != 0.")

print " >> Doing electrons..."
chain_in_ele = ROOT.TChain("fevt/RHTree")
chain_in_ele.Add("histo_Ele.root")
chain_in_ele.Draw("EBenergy >> h_E_ele","EBenergy > 0.")
chain_in_ele.Draw("EBenergyRed >> h_Er_ele","EBenergyRed > 0.")
chain_in_ele.Draw("EBtime >> h_t_ele","EBtime != 0.")
chain_in_ele.Draw("EBtimeRed >> h_tr_ele","EBtimeRed != 0.")

pho_color = 4
ele_color = 2

c_E = ROOT.TCanvas("c_E","c_E",600,400)
ROOT.gPad.SetLogy(1)
h_E_pho.SetLineColor(pho_color)
h_E_pho.Draw("HIST")
h_E_ele.SetLineColor(ele_color)
h_E_ele.Draw("HIST SAME")
c_E.Print("HISTS/c_E.png")
c_E.Print("HISTS/c_E.eps")

c_t = ROOT.TCanvas("c_t","c_t",600,400)
ROOT.gPad.SetLogy(1)
h_t_pho.SetLineColor(pho_color)
h_t_pho.Draw("HIST")
h_t_ele.SetLineColor(ele_color)
h_t_ele.Draw("SAME HIST")
c_t.Print("HISTS/c_t.png")
c_t.Print("HISTS/c_t.eps")

c_Er = ROOT.TCanvas("c_Er","c_Er",600,400)
ROOT.gPad.SetLogy(1)
h_Er_pho.SetLineColor(pho_color)
h_Er_pho.Draw("")
h_Er_ele.SetLineColor(ele_color)
h_Er_ele.Draw("SAME")
c_Er.Print("HISTS/c_Er.png")
c_Er.Print("HISTS/c_Er.eps")

c_tr = ROOT.TCanvas("c_tr","c_tr",600,400)
ROOT.gPad.SetLogy(1)
h_tr_pho.SetLineColor(pho_color)
h_tr_pho.Draw("")
h_tr_ele.SetLineColor(ele_color)
h_tr_ele.Draw("SAME")
c_tr.Print("HISTS/c_tr.png")
c_tr.Print("HISTS/c_tr.eps")

sw.Stop()
print " >> Real time:",sw.RealTime()/60.
