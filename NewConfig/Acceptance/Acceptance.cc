#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TLatex.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "../Bfinder/interface/format.h"
#include "plotDressing.h"
#include <iostream>

void Acceptance()
{
    TChain *root = new TChain("demo/root");
    
    root->Add("DumpGenInfo_all.root");
    
    int n_entries = root->GetEntries();
    printf("Going to process %d entries.\n",n_entries);
    
    //-----------------------------------------------------------------
    // setting memory addresses to the branches
    // create new output trees
    
    GenInfoBranches *GenInfo = new GenInfoBranches;
    GenInfo->setbranchadd(root);
    
    //------------------------------------------------------------------
    // Defining histograms
    double ptBins[] = {10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100.};
    int ptBinsSize = (sizeof(ptBins)/sizeof(*ptBins))-1;
    cout<<"pT bin size = "<<ptBinsSize<<'\t'<<ptBins[0]<<endl;
    
    TFile *fout = new TFile("acceptance.root","recreate");
    
    TH1D* h_bp_pt      = new TH1D("h_bp_pt","Acceptance numerator B+ p_{T}",40,10.,100.);
    TH1D* h_bp_y       = new TH1D("h_bp_y","Acceptance numerator B+ rapidity",40,-2.4,2.4);
    
    TH1D* h_bpfilter_pt      = new TH1D("h_bpfilter_pt","Filter numerator for correction factor",40,10.,100.);
    TH1D* h_bpfilter_y       = new TH1D("h_bpfilter_y","Filter numerator for correction factor",40,-2.4,2.4);
    
    TH2D* h_bp_pt_y_num    = new TH2D("h_bp_pt_y_num","Numerator B+ pt vs rapidity",80,-2.4,2.4,80,10.,100.);
    TH2D* h_bp_pt_y_den    = new TH2D("h_bp_pt_y_den","Denominator B+ pt vs rapidity",80,-2.4,2.4,80,10.,100.);

    TH1D* h_bp_pt_norm = new TH1D("h_bp_pt_norm","generator B+ p_{T}",40,10.,100.);
    TH1D* h_bp_y_norm  = new TH1D("h_bp_y_norm","generator B+ rapidity",40,-2.4,2.4);
    
    h_bp_pt->Sumw2();
    h_bp_y->Sumw2();
    h_bpfilter_pt->Sumw2();
    h_bpfilter_y->Sumw2();
    
    char histoName[200];
    char histoTitle[200];
    std::map<std::string,TH2D*> h_mu1pt_mu2pt;
    
    for(int kPt=0; kPt<ptBinsSize; kPt++) {
        sprintf(histoName, "hSignalPtBin%d", kPt);
        sprintf(histoTitle, "Muon pT correlation for %5.2f < p_{T} < %5.2f ", ptBins[kPt], ptBins[kPt+1]);
        h_mu1pt_mu2pt[histoName] = new TH2D(histoName, histoTitle, 200,0,50,200,0,50);
    
        h_mu1pt_mu2pt[histoName]->GetXaxis()->SetTitle("p_{T}^{#mu_{1}} [GeV]");
        h_mu1pt_mu2pt[histoName]->GetXaxis()->SetLabelFont(42);
        h_mu1pt_mu2pt[histoName]->GetXaxis()->SetLabelOffset(0.01);
        h_mu1pt_mu2pt[histoName]->GetXaxis()->SetTitleSize(0.06);
        h_mu1pt_mu2pt[histoName]->GetXaxis()->SetTitleOffset(1.09);
        h_mu1pt_mu2pt[histoName]->GetXaxis()->SetLabelFont(42);
        h_mu1pt_mu2pt[histoName]->GetXaxis()->SetLabelSize(0.055);
        h_mu1pt_mu2pt[histoName]->GetXaxis()->SetTitleFont(42);
        h_mu1pt_mu2pt[histoName]->GetYaxis()->SetTitle("p_{T}^{#mu_{2}} [GeV]");
        h_mu1pt_mu2pt[histoName]->GetYaxis()->SetLabelFont(42);
        h_mu1pt_mu2pt[histoName]->GetYaxis()->SetLabelOffset(0.01);
        h_mu1pt_mu2pt[histoName]->GetYaxis()->SetTitleOffset(1.14);
        h_mu1pt_mu2pt[histoName]->GetYaxis()->SetTitleSize(0.06);
        h_mu1pt_mu2pt[histoName]->GetYaxis()->SetTitleFont(42);
        h_mu1pt_mu2pt[histoName]->GetYaxis()->SetLabelFont(42);
        h_mu1pt_mu2pt[histoName]->GetYaxis()->SetLabelSize(0.055);
    }
    
    for (int evt=0; evt<n_entries; evt++) {
        if (evt%1000==0 || evt==n_entries-1) printf("processing %d/%d (%.2f%%).\n",evt,n_entries-1,(double)evt/(double)(n_entries-1)*100.);
        
        root->GetEntry(evt);
        
        // Look for indices of the whole decay tree
        for (int idx = 0; idx < GenInfo->size; idx++) {
            
            if (abs(GenInfo->pdgId[idx])==521) { // B+ find
                int idx_bp   = idx;
                int idx_jpsi = GenInfo->da1[idx_bp];
                int idx_kp   = GenInfo->da2[idx_bp];
                int idx_mu1  = GenInfo->da1[idx_jpsi];
                int idx_mu2  = GenInfo->da2[idx_jpsi];
                
                if (GenInfo->pdgId[idx_jpsi]!=443) continue; // not J/psi
                if (abs(GenInfo->pdgId[idx_kp])!=321) continue; // not K+-
                if (abs(GenInfo->pdgId[idx_mu1])!=13) continue; // not mu+-
                if (abs(GenInfo->pdgId[idx_mu2])!=13) continue; // not mu+-

                TLorentzVector v4_bp, v4_mu1, v4_mu2;
                v4_bp.SetPtEtaPhiM(GenInfo->pt[idx_bp],GenInfo->eta[idx_bp],GenInfo->phi[idx_bp],GenInfo->mass[idx_bp]);
                v4_mu1.SetPtEtaPhiM(GenInfo->pt[idx_mu1],GenInfo->eta[idx_mu1],GenInfo->phi[idx_mu1],GenInfo->mass[idx_mu1]);
                v4_mu2.SetPtEtaPhiM(GenInfo->pt[idx_mu2],GenInfo->eta[idx_mu2],GenInfo->phi[idx_mu2],GenInfo->mass[idx_mu2]);

                h_bp_pt_norm->Fill(GenInfo->pt[idx_bp]);
                h_bp_y_norm->Fill(v4_bp.Rapidity());
                h_bp_pt_y_den->Fill(v4_bp.Rapidity(),GenInfo->pt[idx_bp]);
                
                bool muon1Acc = ((fabs(GenInfo->eta[idx_mu1])<1.3 && GenInfo->pt[idx_mu1]>3.3) ||
                                 (fabs(GenInfo->eta[idx_mu1])>1.3 && fabs(GenInfo->eta[idx_mu1])<2.2 && v4_mu1.P()>2.9) ||
                                 (fabs(GenInfo->eta[idx_mu1])>2.2 && fabs(GenInfo->eta[idx_mu1])<2.4 && GenInfo->pt[idx_mu1]>0.8));
                bool muon2Acc = ((fabs(GenInfo->eta[idx_mu2])<1.3 && GenInfo->pt[idx_mu2]>3.3) ||
                                 (fabs(GenInfo->eta[idx_mu2])>1.3 && fabs(GenInfo->eta[idx_mu2])<2.2 && v4_mu2.P()>2.9) ||
                                 (fabs(GenInfo->eta[idx_mu2])>2.2 && fabs(GenInfo->eta[idx_mu2])<2.4 && GenInfo->pt[idx_mu2]>0.8));

                bool kplusAcc = (GenInfo->pt[idx_kp]>0.8 && fabs(GenInfo->eta[idx_kp])<2.4);
                
                bool muon1Filter = fabs(GenInfo->eta[idx_mu1])<2.4 && GenInfo->pt[idx_mu1]>2.8;
                bool muon2Filter = fabs(GenInfo->eta[idx_mu2])<2.4 && GenInfo->pt[idx_mu2]>2.8;
                
                if (muon1Filter && muon2Filter) {
                    h_bpfilter_pt->Fill(GenInfo->pt[idx_bp]);
                    h_bpfilter_y->Fill(v4_bp.Rapidity());
                }
                
                if (muon1Acc && muon2Acc && kplusAcc){
                    h_bp_pt->Fill(GenInfo->pt[idx_bp]);
                    h_bp_y->Fill(v4_bp.Rapidity());
                    h_bp_pt_y_num->Fill(v4_bp.Rapidity(), GenInfo->pt[idx_bp]);
                    
                    char histogramName[200];
                    for(int kPt=0; kPt<ptBinsSize; kPt++) {
                        sprintf(histogramName, "hSignalPtBin%d", kPt);
                        if(GenInfo->pt[idx_bp] >= ptBins[kPt] && GenInfo->pt[idx_bp] < ptBins[kPt+1])
                        {
                           h_mu1pt_mu2pt[histogramName]->Fill(GenInfo->pt[idx_mu1], GenInfo->pt[idx_mu2]);
                        }
                    }
                }
            }
        }
        
    } // end of evt loop
    
    TH1D *h_bp_pt_acc = (TH1D*)h_bp_pt->Clone("h_bp_pt_acc");
    h_bp_pt_acc->Divide(h_bp_pt_norm);
    
    for (int i=1; i<=h_bp_pt_acc->GetNbinsX(); i++) {
        double p = h_bp_pt_acc->GetBinContent(i);
        double N = h_bp_pt_norm->GetBinContent(i);
        h_bp_pt_acc->SetBinError(i,sqrt(p*(1.-p)/N)); // binominal error
        if((sqrt(p*(1.-p)/N)) == 0) cout<<i<<'\t'<<p<<'\t'<<N<<'\t'<<sqrt(p*(1.-p)/N)<<endl;
    }
    
  //  TCanvas *c1 = canvasDressing("c1");
    h_bp_pt_acc->SetTitle("");
    h_bp_pt_acc->GetXaxis()->SetTitle("p_{T} [GeV]");
    h_bp_pt_acc->GetXaxis()->SetLabelFont(42);
    h_bp_pt_acc->GetXaxis()->SetLabelOffset(0.01);
    h_bp_pt_acc->GetXaxis()->SetTitleSize(0.06);
    h_bp_pt_acc->GetXaxis()->SetTitleOffset(1.09);
    h_bp_pt_acc->GetXaxis()->SetLabelFont(42);
    h_bp_pt_acc->GetXaxis()->SetLabelSize(0.055);
    h_bp_pt_acc->GetXaxis()->SetTitleFont(42);
    h_bp_pt_acc->GetYaxis()->SetTitle("Acceptance");
    h_bp_pt_acc->GetYaxis()->SetLabelFont(42);
    h_bp_pt_acc->GetYaxis()->SetLabelOffset(0.01);
    h_bp_pt_acc->GetYaxis()->SetTitleOffset(1.14);
    h_bp_pt_acc->GetYaxis()->SetTitleSize(0.06);
    h_bp_pt_acc->GetYaxis()->SetTitleFont(42);
    h_bp_pt_acc->GetYaxis()->SetLabelFont(42);
    h_bp_pt_acc->GetYaxis()->SetLabelSize(0.055);
    
    TH1D *h_bp_y_acc = (TH1D*)h_bp_y->Clone("h_bp_y_acc");
    h_bp_y_acc->Divide(h_bp_y_norm);
    
    for (int i=1; i<=h_bp_y_acc->GetNbinsX(); i++) {
        double p = h_bp_y_acc->GetBinContent(i);
        double N = h_bp_y_norm->GetBinContent(i);
        h_bp_y_acc->SetBinError(i,sqrt(p*(1.-p)/N)); // binominal error
    }
    
  //  TCanvas *c2 = canvasDressing("c2");
    h_bp_y_acc->SetTitle("");
    h_bp_y_acc->GetXaxis()->SetTitle("Rapidity");
    h_bp_y_acc->GetXaxis()->SetLabelFont(42);
    h_bp_y_acc->GetXaxis()->SetLabelOffset(0.01);
    h_bp_y_acc->GetXaxis()->SetTitleSize(0.06);
    h_bp_y_acc->GetXaxis()->SetTitleOffset(1.09);
    h_bp_y_acc->GetXaxis()->SetLabelFont(42);
    h_bp_y_acc->GetXaxis()->SetLabelSize(0.055);
    h_bp_y_acc->GetXaxis()->SetTitleFont(42);
    h_bp_y_acc->GetYaxis()->SetTitle("Acceptance");
    h_bp_y_acc->GetYaxis()->SetLabelFont(42);
    h_bp_y_acc->GetYaxis()->SetLabelOffset(0.01);
    h_bp_y_acc->GetYaxis()->SetTitleOffset(1.14);
    h_bp_y_acc->GetYaxis()->SetTitleSize(0.06);
    h_bp_y_acc->GetYaxis()->SetTitleFont(42);
    h_bp_y_acc->GetYaxis()->SetLabelFont(42);
    h_bp_y_acc->GetYaxis()->SetLabelSize(0.055);
    
    TH1D *h2d_bp_y_pt_acc = (TH1D*)h_bp_pt_y_num->Clone("h2d_bp_y_pt_acc");
    h2d_bp_y_pt_acc->Divide(h_bp_pt_y_den);
    h2d_bp_y_pt_acc->Sumw2();
    
    h2d_bp_y_pt_acc->SetTitle("");
    h2d_bp_y_pt_acc->GetXaxis()->SetTitle("Rapidity");
    h2d_bp_y_pt_acc->GetXaxis()->SetLabelFont(42);
    h2d_bp_y_pt_acc->GetXaxis()->SetLabelOffset(0.01);
    h2d_bp_y_pt_acc->GetXaxis()->SetTitleSize(0.06);
    h2d_bp_y_pt_acc->GetXaxis()->SetTitleOffset(1.09);
    h2d_bp_y_pt_acc->GetXaxis()->SetLabelFont(42);
    h2d_bp_y_pt_acc->GetXaxis()->SetLabelSize(0.055);
    h2d_bp_y_pt_acc->GetXaxis()->SetTitleFont(42);
    h2d_bp_y_pt_acc->GetYaxis()->SetTitle("Transverse momentum");
    h2d_bp_y_pt_acc->GetYaxis()->SetLabelFont(42);
    h2d_bp_y_pt_acc->GetYaxis()->SetLabelOffset(0.01);
    h2d_bp_y_pt_acc->GetYaxis()->SetTitleOffset(1.14);
    h2d_bp_y_pt_acc->GetYaxis()->SetTitleSize(0.06);
    h2d_bp_y_pt_acc->GetYaxis()->SetTitleFont(42);
    h2d_bp_y_pt_acc->GetYaxis()->SetLabelFont(42);
    h2d_bp_y_pt_acc->GetYaxis()->SetLabelSize(0.055);
    
    // For correction factors //
    TH1D *h_bp_pt_corr = (TH1D*)h_bpfilter_pt->Clone("h_bp_pt_corr");
    h_bp_pt_corr->Divide(h_bp_pt);
    h_bp_pt_corr->Sumw2();
    
   /* for (int i=1; i<=h_bp_pt_corr->GetNbinsX(); i++) {
        double p = h_bp_pt_corr->GetBinContent(i);
        double N = h_bp_pt->GetBinContent(i);
        h_bp_pt_corr->SetBinError(i,sqrt(p*(1.-p)/N)); // binominal error
    }*/
    
    TH1D *h_bp_y_corr = (TH1D*)h_bpfilter_y->Clone("h_bp_y_corr");
    h_bp_y_corr->Divide(h_bp_y);
    h_bp_y_corr->Sumw2();
    
    // Acceptance factors to match official MC sample //
    TH1D *h_bp_pt_AccPrime = (TH1D*)h_bpfilter_pt->Clone("h_bp_pt_AccPrime");
    h_bp_pt_AccPrime->Divide(h_bp_pt_norm);
    h_bp_pt_AccPrime->Sumw2();
    
    fout->Write();
    fout->Close();
    
    delete GenInfo;
}
