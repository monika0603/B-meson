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
#include "myloop_mc.h"
#include "plotDressing.h"

using namespace RooFit;

// General fitting options
#define NUMBER_OF_CPU       1
#define DO_MINOS            kTRUE
// 0 - w/o DISPLAY
// 1 - w/  DISPLAY
#define DISPLAY             1

#define SOURCE              "myloop_mc.root"

void eval_eff(double pt_min = 10., double pt_max = 100., double *res = NULL)
{
    double mass_min = 5.18, mass_max = 6.0;
    double mass_peak = BP_MASS;
    int nbins = 50;
    
    RooRealVar mass("mass","mass",mass_min,mass_max);
    
    TFile *fout = new TFile("myfitter.root","recreate");
    TNtupleD *_nt = new TNtupleD("_nt","_nt","mass");
    
    TFile *fin = new TFile(SOURCE);
    TTree *tin = (TTree*)fin->Get("ntkp_mc");
    
    ReducedBranches br;
    br.setbranchadd(tin);
    
    int norm_count = 0;
    int current_run = -1, current_event = -1;
    
    for (int evt=0;evt<tin->GetEntries();evt++) {
        tin->GetEntry(evt);
        
        if (fabs(br.geny)>2.4) continue;
        if (br.genpt<pt_min || br.genpt>=pt_max) continue; // within the pt binning

        if (br.event != current_event || br.run != current_run) {
            current_run = br.run;
            current_event = br.event;
            norm_count++;
        }
        
        //if (br.hltbook[HLT_Dimuon16_Jpsi_v1]!=1) continue;
        if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1) continue;
        if (br.vtxprob<=0.1) continue;
        if (br.tk1pt<=1.0) continue;
        if (br.lxy/br.errxy<=3.0) continue;
        if (br.cosalpha2d<=0.99) continue;
        
        _nt->Fill(&br.mass);
    }
    fin->Close();
    
    RooDataSet *data = new RooDataSet("data","data",_nt,RooArgSet(mass));
    
    double n_signal_initial = data->sumEntries(TString::Format("abs(mass-%g)<0.015",mass_peak))
                            - data->sumEntries(TString::Format("abs(mass-%g)<0.030&&abs(mass-%g)>0.015",mass_peak,mass_peak));

    double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
    
    //-----------------------------------------------------------------
    // signal PDF 
    
    RooRealVar m_mean("m_mean","m_mean",mass_peak,mass_min,mass_max);
    RooRealVar m_sigma1("m_sigma1","m_sigma1",0.015,0.001,0.050);
    RooRealVar m_sigma2("m_sigma2","m_sigma2",0.030,0.001,0.100);
    RooRealVar m_fraction("m_fraction","m_fraction",0.5);
    
    if (pt_min>=60.) {
        m_fraction.setVal(1.);
        m_sigma2.setConstant(kTRUE);
    }
    
    RooGaussian m_gaussian1("m_gaussian1","m_gaussian1",mass,m_mean,m_sigma1);
    RooGaussian m_gaussian2("m_gaussian2","m_gaussian2",mass,m_mean,m_sigma2);
    
    RooAddPdf pdf_m_signal("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));
    
    //-----------------------------------------------------------------
    // combinatorial background PDF (exponential)
    
    RooRealVar m_exp("m_exp","m_exp",-0.3,-10.,+10.);
    RooExponential pdf_m_combinatorial_exp("pdf_m_combinatorial_exp","pdf_m_combinatorial_exp",mass,m_exp);
    
    //-----------------------------------------------------------------
    // full model
    
    RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries()*2.);
    RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
    
    RooAddPdf *model = new RooAddPdf("model","model",
        RooArgList(pdf_m_signal, pdf_m_combinatorial_exp),
        RooArgList(n_signal, n_combinatorial));
    
    model->fitTo(*data,Minos(DO_MINOS),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));
    
    if (res!=NULL) {
        res[0] = n_signal.getVal();
        res[1] = n_signal.getError();
        res[2] = n_signal.getErrorLo();
        res[3] = n_signal.getErrorHi();
        res[4] = norm_count;
    }

#if DISPLAY
    TCanvas *c1 = canvasDressing("c1");
    
    // Display
    RooPlot* frame_m = mass.frame();
    
    TH1D* histo_data = (TH1D*)data->createHistogram("histo_data", mass, Binning(nbins,mass_min,mass_max));
    histo_data->Sumw2(false);
    histo_data->SetBinErrorOption(TH1::kPoisson);
    histo_data->SetMarkerStyle(20);
    histo_data->SetMarkerSize(0.8);
    histo_data->SetLineColor(kBlack);
    for (int i=1; i<=nbins; i++)
        if (histo_data->GetBinContent(i)==0) histo_data->SetBinError(i,0.);
    
    data->plotOn(frame_m,Binning(nbins),Invisible());
    
    model->plotOn(frame_m,Precision(2E-4));
    model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_signal),LineColor(kRed),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kRed), VLines(), DrawOption("F"));
    model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_combinatorial_exp),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    
    frame_m->SetTitle("");
    frame_m->GetXaxis()->SetTitle("M_{J/#psi K^{#pm}} [GeV]");
    frame_m->GetXaxis()->SetLabelFont(42);
    frame_m->GetXaxis()->SetLabelOffset(0.01);
    frame_m->GetXaxis()->SetTitleSize(0.06);
    frame_m->GetXaxis()->SetTitleOffset(1.09);
    frame_m->GetXaxis()->SetLabelFont(42);
    frame_m->GetXaxis()->SetLabelSize(0.055);
    frame_m->GetXaxis()->SetTitleFont(42);
    frame_m->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass_max-mass_min)*1000./nbins));
    frame_m->GetYaxis()->SetLabelFont(42);
    frame_m->GetYaxis()->SetLabelOffset(0.01);
    frame_m->GetYaxis()->SetTitleOffset(1.14);
    frame_m->GetYaxis()->SetTitleSize(0.06);
    frame_m->GetYaxis()->SetTitleFont(42);
    frame_m->GetYaxis()->SetLabelFont(42);
    frame_m->GetYaxis()->SetLabelSize(0.055);
    
    frame_m->Draw();
    histo_data->Draw("Esame");
    
    c1->SaveAs(TString::Format("fig/bp_efficiency_bin_%g_%g.pdf",pt_min,pt_max));
    
#endif
    
    fout->Write();
    fout->Close();
}

void myfitter_eff()
{
    double res[5];
    double pt_bins[11] = {10.,15.,20.,25.,30.,40.,50.,60.,70.,80.,100.};
    TH1D* h_bp_efficiency = new TH1D("h_bp_efficiency","Efficiencies of B+",10,pt_bins);
    
    for (int i=0; i<10; i++) {
        eval_eff(pt_bins[i],pt_bins[i+1],res);
        
        double width = pt_bins[i+1]-pt_bins[i];
        
        h_bp_efficiency->SetBinContent(i+1,res[0]/res[4]);
        h_bp_efficiency->SetBinError(i+1,res[1]/res[4]);
    }
    
    TFile *fout = new TFile("myfitter_eff.root","recreate");
    h_bp_efficiency->Write();
    fout->Close();
}