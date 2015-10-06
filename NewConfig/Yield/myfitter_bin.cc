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
#include <TGraphErrors.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include "../myloop.h"
#include "../plotDressing.h"

using namespace RooFit;

// General fitting options
#define NUMBER_OF_CPU       1
#define DO_MINOS            kTRUE
// 0 - w/o DISPLAY
// 1 - w/  DISPLAY
#define invDISPLAY             0
#define yieldDISPLAY           1

#define SOURCE              "myloop_reference.root"

void eval_yield(double pt_min = 10., double pt_max = 100., double *res = NULL, int bkg_choice = 1)
{
    double mass_min = 5.18, mass_max = 6.0;
    double mass_peak = BP_MASS;
    int nbins = 41;
    
    RooRealVar mass("mass","mass",mass_min,mass_max);
    
    TFile *fout = new TFile("myfitter.root","recreate");
    TNtupleD *_nt = new TNtupleD("_nt","_nt","mass");
    
    TFile *fin = new TFile(SOURCE);
    TTree *tin = (TTree*)fin->Get("ntkp");
    
    ReducedBranches br;
    br.setbranchadd(tin);
    
    for (int evt=0;evt<tin->GetEntries();evt++) {
        tin->GetEntry(evt);
        
        //if (br.hltbook[HLT_Dimuon16_Jpsi_v1]!=1) continue;
        if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
        if (br.vtxprob<=0.1) continue;
        if (br.tk1pt<=1.0) continue;
        if (br.lxy/br.errxy<=3.0) continue;
        if (br.cosalpha2d<=0.99) continue;
        
        if (br.pt<pt_min || br.pt>=pt_max) continue; // within the pt binning
        
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
    
    RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries()*2.);
    RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
    
    RooAddPdf *model = 0;
    
    //-----------------------------------------------------------------
    // combinatorial background PDF (Bernstein Polynomial)

    RooRealVar m1_poly("m1_poly", "coeficient of mass term", 0.5, 0., 10.);
    RooRealVar m2_poly("m2_poly", "coeficient of mass^2 term", 0.5, 0., 10.);
    RooRealVar m3_poly("m3_poly", "coeficient of mass^3 term", 0.5, 0., 10.);
    
    // f(x) = 1+c_1x+c_2x^2
  //  RooPolynomial pdf_m_combinatorial_poly("pdf_m_combinatorial_poly", "cubic polynomial", mass, RooArgList(m1_poly) );
    RooBernstein pdf_m_combinatorial_bern("pdf_m_combinatorial_bern", "Bernstein polynomial", mass, RooArgList(RooConst(1.),m1_poly, m2_poly) );
    
    //-----------------------------------------------------------------
    // combinatorial background PDF (exponential)
    
    RooRealVar m_exp("m_exp","m_exp",-0.3,-10.,+10.);
    RooExponential pdf_m_combinatorial_exp("pdf_m_combinatorial_exp","pdf_m_combinatorial_exp",mass,m_exp);

    // Full model
    switch (bkg_choice){
        
        case 1:
            model = new RooAddPdf("model","model",
                                             RooArgList(pdf_m_signal, pdf_m_combinatorial_bern),
                                             RooArgList(n_signal, n_combinatorial));
            break;
            
        case 2:
            model = new RooAddPdf("model","model",
                                  RooArgList(pdf_m_signal, pdf_m_combinatorial_exp),
                                  RooArgList(n_signal, n_combinatorial));
            break;
    }
    
    model->fitTo(*data,Minos(DO_MINOS),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));
    
    if (res!=NULL) {
        res[0] = n_signal.getVal();
        res[1] = n_signal.getError();
        res[2] = n_signal.getErrorLo();
        res[3] = n_signal.getErrorHi();
    }
    
    cout<<"****************************************"<<endl;
    cout<<"*******"<< n_signal.getVal()<<endl;
    cout<<"****************************************"<<endl;
    
#if invDISPLAY
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
    c1->SaveAs(TString::Format("fig/bp_yield_bin_%g_%g.pdf",pt_min,pt_max));
#endif
    
    fout->Write();
    fout->Close();
}

void myfitter_bin()
{
    double res[4];
    double ptBins[] = {10.,15.,20.,25.,30.,40.,50.,60.,70.,80.,100.};
    int ptBinsSize = (sizeof(ptBins)/sizeof(*ptBins))-1;
    cout<<"pT bin size = "<<ptBinsSize<<'\t'<<ptBins[0]<<endl;
    double y[9];
    double ey[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    
    TH1D *h_dummy = new TH1D("h_dummy", "Dummy histo", ptBinsSize, &ptBins[0]);
    TH1D *h_bp_yields_expBkg = new TH1D("h_bp_yields_expBkg", "B-meson p_{T} distribution", ptBinsSize, &ptBins[0]);
    
    TH1D *h_bp_yields_polyBkg = new TH1D("h_bp_yields_polyBkg", "B-meson p_{T} distribution", ptBinsSize, &ptBins[0]);
    
    TString yaxis_title = "d#sigma/dp_{T} (pp #rightarrow B^{+}X; |y|<2.4) [pb/GeV]";
    
    //------ Reading efficiency times acceptance -----------------
    TFile *fileEffAcc = new TFile("./EffAcc.root","r");
    TH1D *h_ptEffAcc = (TH1D*)fileEffAcc->Get("h_bp_pt_AccEff");
    cout<<"Number of bins in acceptance*efficiency = "<<h_ptEffAcc->GetNbinsX()<<endl;
    
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
    
    //--------------------------------------------------------------
    double BR = (1.027e-03) * (5.961e-02);
    double Lint = 47; //47 (pb)-1
    double FF = 0.402;
    
    for (int i=1; i<=h_bp_yields_expBkg->GetNbinsX(); i++) {
        eval_yield(ptBins[i-1],ptBins[i],res, 2);
        
        double binWid = h_bp_yields_expBkg->GetBinWidth(i);
        double accEff = h_ptEffAcc->GetBinContent(i);
    
        h_bp_yields_expBkg->SetBinContent(i,res[0]/binWid/BR/accEff/Lint/2.0/FF);
        h_bp_yields_expBkg->SetBinError(i,res[1]/binWid/BR/accEff/Lint/2.0/FF);
        
        ey[i-1] = sqrt( h_bp_yields_expBkg->GetBinError(i)/h_bp_yields_expBkg->GetBinContent(i) + FONLL_13TeVE[i-1]/FONLL_13TeV[i-1]);
        y[i-1] = h_bp_yields_expBkg->GetBinContent(i)/FONLL_13TeV[i-1];
    }
    
    for (int i=1; i<=h_bp_yields_polyBkg->GetNbinsX(); i++) {
        eval_yield(ptBins[i-1],ptBins[i],res, 1);
        
        double binWid = h_bp_yields_polyBkg->GetBinWidth(i);
        double accEff = h_ptEffAcc->GetBinContent(i);
        
        h_bp_yields_polyBkg->SetBinContent(i,res[0]/binWid/BR/accEff/Lint/2.0/FF);
        h_bp_yields_polyBkg->SetBinError(i,res[1]/binWid/BR/accEff/Lint/2.0/FF);
    }
    
    TFile *fout = new TFile("myfitter_bin.root","recreate");
    h_bp_yields_expBkg->Write();
    fout->Close();
    
#if yieldDISPLAY
    TCanvas *c2 = canvasDressing("c2");
    c2->SetLogy();
    c2->cd();
    
    TPad *pp1=0, *pp1_1=0;
    pp1 = new TPad("p1","p1",0,0.34,1,1,0,0,0);
    pp1->SetBottomMargin(0.14);
    pp1->SetTopMargin(0.05*(1/0.72));
    pp1->SetLeftMargin(0.15);
    pp1->SetRightMargin(0.035);
     
    pp1->Draw();
    pp1->cd();
    pp1->SetNumber(1);
    pp1->SetLogy();
    
    TH1D *hDummy = new TH1D("hDummy", "", 1, 10, 100);
    
    hDummy->SetTitle("");
    hDummy->GetXaxis()->SetTitle("p_{T}[GeV]");
    hDummy->GetXaxis()->SetLabelFont(42);
    hDummy->GetXaxis()->SetLabelOffset(0.01);
    hDummy->GetXaxis()->SetTitleSize(0.06);
    hDummy->GetXaxis()->SetTitleOffset(1.09);
    hDummy->GetXaxis()->SetLabelFont(42);
    hDummy->GetXaxis()->SetLabelSize(0.055);
    hDummy->GetXaxis()->SetTitleFont(42);
    hDummy->GetYaxis()->SetTitle(yaxis_title);
    hDummy->GetYaxis()->SetLabelFont(42);
    hDummy->GetYaxis()->SetLabelOffset(0.01);
    hDummy->GetYaxis()->SetTitleOffset(1.14);
    hDummy->GetYaxis()->SetTitleSize(0.06);
    hDummy->GetYaxis()->SetTitleFont(42);
    hDummy->GetYaxis()->SetLabelFont(42);
    hDummy->GetYaxis()->SetLabelSize(0.055);
    hDummy->SetMinimum(100.);
    hDummy->SetMaximum(10e6);
    hDummy->Draw();
    
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
    hFONLL_13TeV->Draw("p same");
    
    h_bp_yields_expBkg->SetMarkerStyle(20);
    h_bp_yields_expBkg->SetMarkerColor(kBlack);
    h_bp_yields_expBkg->SetLineColor(kRed-7);
    h_bp_yields_expBkg->SetLineWidth(4);
    h_bp_yields_expBkg->Draw("same");
    
    cout<<"---- Corrected yield ----- "<<endl;
    for (int i=1; i<=h_bp_yields_expBkg->GetNbinsX(); i++){
        double exp_fit = h_bp_yields_expBkg->GetBinContent(i);
        double poly_fit = h_bp_yields_polyBkg->GetBinContent(i);
        double fractional_change = ((poly_fit - exp_fit)/exp_fit)*100.0;
        
        cout<<ptBins[i-1]<<" < pT < "<<ptBins[i]<<'\t'<<h_bp_yields_expBkg->GetBinContent(i)<<'\t'<<h_ptEffAcc->GetBinContent(i)<<'\t'<<h_bp_yields_polyBkg->GetBinContent(i)<<'\t'<<fractional_change<<endl;
    }
    cout<<"-------------------------- "<<endl;
    
    /*c2->cd();
     
    pp1_1 = new TPad("p1_1","p1_1",0,0.0,1,0.34,0,0,0);
    pp1_1->SetBottomMargin(0.2);
    pp1_1->SetTopMargin(0.04);
    pp1_1->SetLeftMargin(0.15);
    pp1_1->SetRightMargin(0.035);
     
    pp1_1->Draw();
    pp1_1->cd();
    pp1_1->SetNumber(2);
     
    const int n = 9;
    double x[n] = {15., 25., 35., 45., 55., 65., 75., 85., 95.};
    double ex[n] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
    TGraphErrors *grDataFONLL = new TGraphErrors(n,x,y,ex,ey);
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
    grDataFONLL->Draw("ALP");*/
#endif
    
    /*for (int i=1; i<=h_bp_yields->GetNbinsX(); i++) {
         ey[i-1] = sqrt( h_bp_yields->GetBinError(i)/h_bp_yields->GetBinContent(i) + FONLL_13TeVE[i-1]/FONLL_13TeV[i-1]);
        y[i-1] = h_bp_yields->GetBinContent(i)/FONLL_13TeV[i-1];
        cout<<i<<'\t'<<h_bp_yields->GetBinContent(i)<<'\t'<<FONLL_13TeV[i-1]<<'\t'<<y[i-1]<<'t'<<ey[i-1]<<endl;
    }*/
}