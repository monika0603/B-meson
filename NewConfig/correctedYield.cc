#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include "myloop.h"
#include "plotDressing.h"
#include <iostream>

using namespace RooFit;

// General fitting options
#define NUMBER_OF_CPU       1
#define DO_MINOS            kTRUE
// 0 - w/o DISPLAY
// 1 - w/  DISPLAY
#define DISPLAY             1

#define SOURCE              "myloop_reference.root"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda

void correctedYield(int channel = 1, int beamSpotErrEstimate = 1, int comparison = 1)
{
    double ptBins[] = {10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100.};
    int ptBinsSize = (sizeof(ptBins)/sizeof(*ptBins))-1;
    cout<<"pT bin size = "<<ptBinsSize<<'\t'<<ptBins[0]<<endl;
    
    TH1D *bPtHist = new TH1D("bPtHist", "B-meson p_{T} distribution", ptBinsSize, &ptBins[0]);
    
    double mass_min, mass_max, mass_peak;
    int nbins;
    TString ntuple_name = "", xaxis_title = "", yaxis_title = "";
    
    switch (channel) {
        case 1:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = BP_MASS;
            nbins = 50;
            ntuple_name = "ntkp";
            xaxis_title = "M_{J/#psi K^{#pm}} [GeV]";
            break;
        case 2:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = B0_MASS;
            nbins = 50;
            ntuple_name = "ntkstar";
            xaxis_title = "M_{J/#psi K^{#pm}#pi^{#mp}} [GeV]";
            break;
        case 3:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = B0_MASS;
            nbins = 50;
            ntuple_name = "ntks";
            xaxis_title = "M_{J/#psi K^{0}_{S}} [GeV]";
            break;
        case 4:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = BS_MASS;
            nbins = 50;
            ntuple_name = "ntphi";
            xaxis_title = "M_{J/#psi K^{#pm}K^{#mp}} [GeV]";
            break;
        case 5:
            mass_min = 3.6; mass_max = 4.0;
            mass_peak = PSI2S_MASS;
            nbins = 80;
            ntuple_name = "ntmix";
            xaxis_title = "M_{J/#psi #pi^{#pm}#pi^{#mp}} [GeV]";
            break;
        case 6:
            mass_min = 5.3; mass_max = 6.3;
            mass_peak = LAMBDAB_MASS;
            nbins = 50;
            ntuple_name = "ntlambda";
            xaxis_title = "M_{J/#psi #Lambda} [GeV]";
            break;
    }
    
    TFile *fin = new TFile(SOURCE);
    TTree *tin = (TTree*)fin->Get(ntuple_name);
    
    ReducedBranches br;
    br.setbranchadd(tin);
    
    int n_br_queued = 0;
    ReducedBranches br_queue[32];
    
    for (int evt=0;evt<tin->GetEntries();evt++) {
        tin->GetEntry(evt);
        
        if (channel==1) { // cuts for B+ -> J/psi K+
            switch (beamSpotErrEstimate) {
                case 1:
                    if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
                    break;
                case 2:
                    if (br.hltbook[HLT_Dimuon16_Jpsi_v1]!=1) continue;
                    break;
                case 3:
                    if (br.hltbook[HLT_Dimuon0_Jpsi_Muon_v1]!=1) continue;
                    break;
                case 4:
                    if (br.hltbook[HLT_Dimuon10_Jpsi_Barrel_v1]!=1) continue;
                    break;
            }
            if (br.vtxprob<=0.1) continue;
            if (br.tk1pt<=1.0) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            
            if (br.mass >= 5.16 && br.mass <= 5.365)
                bPtHist->Fill(br.pt);
            
        }else
        if (channel==2) { // cuts for B0 -> J/psi K*
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-KSTAR_MASS)>=0.05) continue;
            
            TLorentzVector v4_tk1, v4_tk2;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-PHI_MASS)<=0.01) continue;
            
            if (n_br_queued==0) {
                memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
                n_br_queued++;
            }else
            if (br.run == br_queue[n_br_queued-1].run && br.event == br_queue[n_br_queued-1].event) { // same event
                memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
                n_br_queued++;
                if (n_br_queued>=32) printf("Warning: maximum queued branches reached.\n");
            }
            
            if (br.run != br_queue[n_br_queued-1].run || br.event != br_queue[n_br_queued-1].event || evt==tin->GetEntries()-1) {
                for (int i=0; i<n_br_queued; i++) {
                    
                    bool isBestKstarMass = true;
                    for (int j=0; j<n_br_queued; j++) {
                        if (j==i) continue;
                        if (br_queue[i].mu1idx==br_queue[j].mu1idx &&
                            br_queue[i].mu2idx==br_queue[j].mu2idx &&
                            br_queue[i].tk1idx==br_queue[j].tk1idx &&
                            br_queue[i].tk2idx==br_queue[j].tk2idx) {
                        
                            if (fabs(br_queue[j].tktkmass-KSTAR_MASS)<fabs(br_queue[i].tktkmass-KSTAR_MASS)) {
                                isBestKstarMass = false;
                                continue;
                            }
                        }
                    }
                                 
                  //  if (isBestKstarMass) _nt->Fill(&br_queue[i].mass);
                }
                
                n_br_queued = 0;
                memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
                n_br_queued++;
            }
        }else
        if (channel==3) { // cuts for B0 -> J/psi Ks
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-KSHORT_MASS)>=0.015) continue;
                
           // _nt->Fill(&br.mass);
        }else
        if (channel==4) { // cuts for Bs -> J/psi phi
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-PHI_MASS)>=0.010) continue;
            
            TLorentzVector v4_tk1, v4_tk2;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
                
          //  _nt->Fill(&br.mass);
        }else
        if (channel==5) { // cuts for psi(2S)/X(3872) -> J/psi pipi
            if (br.vtxprob<=0.2) continue;
            if (fabs(br.tk1eta)>=1.6) continue;
            if (fabs(br.tk2eta)>=1.6) continue;
            
           // _nt->Fill(&br.mass);
        }else
        if (channel==6) {
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-LAMBDA_MASS)>=0.015) continue;
            
            TLorentzVector v4_tk1, v4_tk2;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSHORT_MASS)<=0.015) continue;
            
            
          //  _nt->Fill(&br.mass);
        }
    }
    fin->Close();
    
    //------- Reading acceptance * efficiency numbers -----------------
    
    TFile *fileEffAcc = new TFile("./EffAcc.root","r");
    TH1D *h_ptEffAcc = (TH1D*)fileEffAcc->Get("h_bp_pt_AccEff");
    cout<<"Number of bins in acceptance*efficiency = "<<h_ptEffAcc->GetNbinsX()<<endl;
    
    //------- Published 7 TeV numbers----------------------------------
    double cross_section7TeV[] = {4.07, 1.47, 0.412, 0.181, 0.042};
    double cross_section7TeVstatE[] = {0.47, 0.13, 0.041, 0.015, 0.007}; //Statistical errors
    double cross_section7TeVE[] = {0.31, 0.09, 0.026, 0.012, 0.004}; //Systematic errors
    double ptBins_7TeV[] = {5., 10., 13., 17., 24., 30.};
    int ptBinsSize_7TeV = (sizeof(ptBins_7TeV)/sizeof(*ptBins_7TeV))-1;
    
    TH1D *bPtHist_7TeV = new TH1D("bPtHist_7TeV", "7 TeV B-meson p_{T} distribution", ptBinsSize_7TeV, &ptBins_7TeV[0]);
    
    for(int j = 1; j <= bPtHist_7TeV->GetNbinsX(); j++)
    {
        cout<<j<<'\t'<<cross_section7TeV[j-1]<<endl;
        bPtHist_7TeV->SetBinContent(j,cross_section7TeV[j-1]);
        bPtHist_7TeV->SetBinError(j,cross_section7TeVE[j-1]);
    }
    
    //------- FONLL 13 TeV numbers----------------------------------
    
    double FONLL_13TeV[] = {3.3063e+06, 3.3497e+05, 6.4092e+04, 1.7701e+04, 6.1757e+03, 2.5320e+03, 1.1675e+03,
    5.8851e+02, 3.1817e+02, 1.8195e+02};
    double FONLL_13TeVE[] = {2.2612e+06, 1.6605e+05, 2.5186e+04, 5.872e+03, 1.8076e+03, 6.738e+02, 2.8859e+02,
    1.3737e+02, 7.0989e+01, 3.917e+01};
    double ptBins_13TeV[] = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
    int ptBinsSize_13TeV = (sizeof(ptBins_13TeV)/sizeof(*ptBins_13TeV))-1;
    
    TH1D *hFONLL_13TeV = new TH1D("hFONLL_13TeV", "13 TeV B-meson FONLL distribution", ptBinsSize_13TeV, &ptBins_13TeV[0]);
    
    for(int j = 1; j <= hFONLL_13TeV->GetNbinsX(); j++)
    {
        hFONLL_13TeV->SetBinContent(j,FONLL_13TeV[j-1]);
    }
    
    //------- ATLAS 13 TeV numbers----------------------------------
    
    double ATLAS_7TeV[] = {103.4, 36.03, 15.33, 6.056, 1.814, 0.3477, 0.06244,
    0.006099};
    double ATLAS_7TeVE[] = {7.6, 2.32, 0.98, 0.376, 0.115, 0.028, 0.00526, 6.66e-14};
    double ptBins_ATLAS7TeV[] = {9., 13., 16., 20., 25., 35., 50., 70., 120.};
    int ptBinsSize_ATLAS7TeV = (sizeof(ptBins_ATLAS7TeV)/sizeof(*ptBins_ATLAS7TeV))-1;
    
    TH1D *hATLAS7TeV = new TH1D("hATLAS7TeV", "7 TeV B+ meson ATLAS cross-section", ptBinsSize_ATLAS7TeV, &ptBins_ATLAS7TeV[0]);
    
    for(int j = 1; j <= hATLAS7TeV->GetNbinsX(); j++)
    {
       // cout<<j<<'\t'<<ATLAS_7TeV[j-1]<<endl;
        hATLAS7TeV->SetBinContent(j,ATLAS_7TeV[j-1]);
        hATLAS7TeV->SetBinError(j,ATLAS_7TeVE[j-1]);
    }
    //-----------------------------------------------------------------

    double etaMax = 2.4;
    double etaMin = -2.4;
    double etaWid = etaMax - etaMin;
    double nEvents = 739847.;
    double BR = (1.027e-03) * (5.961e-02);
    double Lint = 0.; //47 (pb)-1
    double FF = 0.402;
    
    switch (comparison) {
        case 1: //7 TeV data comparison
            Lint = 47e06;
            yaxis_title = "d#sigma/dp_{T} (pp #rightarrow B^{+}X; |y|<2.4) [#mub/GeV]";
            break;
            
        case 2:
            Lint = 47; //FONLL comparison
            yaxis_title = "d#sigma/dp_{T} (pp #rightarrow B^{+}X; |y|<2.4) [pb/GeV]";
            break;
    }
    
    bPtHist->Sumw2();
    double y[6];
    
    cout<<"--- CMS 13 TeV results ---- "<<endl;
    for(int i = 1; i <= bPtHist->GetNbinsX(); i++)
    {
        double binWid = bPtHist->GetBinWidth(i);
        double accEff = h_ptEffAcc->GetBinContent(i);
        
        switch (comparison) {
            case 1: //7 TeV data comparison
                bPtHist->SetBinContent(i, bPtHist->GetBinContent(i)/binWid/BR/accEff/Lint/2.0);
                bPtHist->SetBinError(i, bPtHist->GetBinError(i)/binWid/BR/accEff/Lint/2.0);
                break;
                
            case 2:
                bPtHist->SetBinContent(i, bPtHist->GetBinContent(i)/binWid/BR/accEff/Lint/FF/2.0);
                bPtHist->SetBinError(i, bPtHist->GetBinError(i)/binWid/BR/accEff/Lint/FF/2.0);
                break;
        }
        cout<<i<<'\t'<<bPtHist->GetBinContent(i)<<endl;
        if(i>4){
            y[i-5] = bPtHist->GetBinContent(i)/FONLL_13TeV[i-3];
             cout<<bPtHist->GetBinContent(i)<<'\t'<<FONLL_13TeV[i-3]<<'\t'<<y[i-5]<<'\t'<<i-5<<endl;
        }

    }

#if DISPLAY
    TCanvas *c1 = canvasDressing("c1");
    c1->SetLogy();
    c1->cd();
    
 /*   TPad *pp1=0, *pp1_1=0;
    pp1 = new TPad("p1","p1",0,0.34,1,1,0,0,0);
    pp1->SetBottomMargin(0.14);
    pp1->SetTopMargin(0.05*(1/0.72));
    pp1->SetLeftMargin(0.15);
    pp1->SetRightMargin(0.035);
    
    pp1->Draw();
    pp1->cd();
    pp1->SetNumber(1);
    pp1->SetLogy();*/
    bPtHist->SetTitle("");
    bPtHist->GetXaxis()->SetTitle("p_{T}[GeV]");
    bPtHist->GetXaxis()->SetLabelFont(42);
    bPtHist->GetXaxis()->SetLabelOffset(0.01);
    bPtHist->GetXaxis()->SetTitleSize(0.06);
    bPtHist->GetXaxis()->SetTitleOffset(1.09);
    bPtHist->GetXaxis()->SetLabelFont(42);
    bPtHist->GetXaxis()->SetLabelSize(0.055);
    bPtHist->GetXaxis()->SetTitleFont(42);
    bPtHist->GetYaxis()->SetTitle(yaxis_title);
    bPtHist->GetYaxis()->SetLabelFont(42);
    bPtHist->GetYaxis()->SetLabelOffset(0.01);
    bPtHist->GetYaxis()->SetTitleOffset(1.14);
    bPtHist->GetYaxis()->SetTitleSize(0.06);
    bPtHist->GetYaxis()->SetTitleFont(42);
    bPtHist->GetYaxis()->SetLabelFont(42);
    bPtHist->GetYaxis()->SetLabelSize(0.055);
    bPtHist->SetMarkerStyle(20);
    bPtHist->SetMarkerColor(kBlack);
    bPtHist->SetLineColor(kRed-7);
    bPtHist->SetLineWidth(4);
  //  bPtHist->SetMinimum(100.);
  //  bPtHist->SetMaximum(10e6);
    bPtHist->Draw();

    switch (comparison) {
        case 1: //7 TeV data comparison
            bPtHist_7TeV->SetMarkerStyle(21);
            bPtHist_7TeV->SetMarkerColor(kBlack);
            bPtHist_7TeV->SetLineWidth(4);
            bPtHist_7TeV->SetLineColor(kYellow);
            bPtHist_7TeV->Draw("same");
            break;
            
        case 2: //FONLL comparison
            TBox *box1 = new TBox(10,2.3867e+06,20,4.6479e+06);
            box1->SetFillColor(kOrange);
            box1->Draw();
            
            TBox *box2 = new TBox(20,2.6535e+05,30,4.3140e+05);
            box2->SetFillColor(kOrange);
            box2->Draw();
            
            TBox *box3 = new TBox(30,5.2986e+04,40,7.8172e+04);
            box3->SetFillColor(kOrange);
            box3->Draw();
            
            TBox *box4 = new TBox(40,1.5014e+04,50,2.0886e+04);
            box4->SetFillColor(kOrange);
            box4->Draw();
            
            TBox *box5 = new TBox(50,5.3260e+03,60,7.1336e+03);
            box5->SetFillColor(kOrange);
            box5->Draw();
            
            TBox *box6 = new TBox(60,2.2088e+03,70,2.8827e+03);
            box6->SetFillColor(kOrange);
            box6->Draw();
            
            TBox *box7 = new TBox(70,1.0270e+03,80,1.3156e+03);
            box7->SetFillColor(kOrange);
            box7->Draw();
            
            TBox *box8 = new TBox(80,5.2087e+02,90,6.5824e+02);
            box8->SetFillColor(kOrange);
            box8->Draw();
            
            TBox *box9 = new TBox(90,2.8291e+02,100,3.5390e+02);
            box9->SetFillColor(kOrange);
            box9->Draw();
            
            hFONLL_13TeV->SetMarkerStyle(21);
            hFONLL_13TeV->SetMarkerColor(kBlack);
            hFONLL_13TeV->Draw("P same");
            break;
    }
  
   /* c1->cd();
    
    pp1_1 = new TPad("p1_1","p1_1",0,0.0,1,0.34,0,0,0);
    pp1_1->SetBottomMargin(0.2);
    pp1_1->SetTopMargin(0.04);
    pp1_1->SetLeftMargin(0.15);
    pp1_1->SetRightMargin(0.04);

    pp1_1->Draw();
    pp1_1->cd();
    pp1_1->SetNumber(2);

    int n = 6;
    double x[6] = {35.,45.,55.,65.,75.,90.};
    TGraph *grDataFONLL = new TGraph(n,x,y);
    grDataFONLL->SetTitle("");
    grDataFONLL->GetXaxis()->SetTitle("p_{T}[GeV]");
    grDataFONLL->GetXaxis()->SetLabelFont(42);
    grDataFONLL->GetXaxis()->SetLabelOffset(0.01);
    grDataFONLL->GetXaxis()->SetTitleSize(0.09);
    grDataFONLL->GetXaxis()->SetTitleOffset(0.99);
    grDataFONLL->GetXaxis()->SetLabelFont(42);
    grDataFONLL->GetXaxis()->SetLabelSize(0.1);
    grDataFONLL->GetXaxis()->SetTitleFont(42);
    grDataFONLL->GetYaxis()->SetTitle("Ratio(Data/FONLL)");
    grDataFONLL->GetYaxis()->SetLabelFont(42);
    grDataFONLL->GetYaxis()->SetLabelOffset(0.01);
    grDataFONLL->GetYaxis()->SetTitleOffset(0.59);
    grDataFONLL->GetYaxis()->SetTitleSize(0.11);
    grDataFONLL->GetYaxis()->SetTitleFont(42);
    grDataFONLL->GetYaxis()->SetLabelFont(42);
    grDataFONLL->GetYaxis()->SetLabelSize(0.1);
    grDataFONLL->Draw("AC*");*/
    
  //hATLAS7TeV->Draw("same");
  //  LegendpTSpectrum();
    
#endif
    
}
