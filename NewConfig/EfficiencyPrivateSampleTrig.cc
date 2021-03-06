#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TEfficiency.h>
#include "../Bfinder/interface/format.h"
#include "plotDressing.h"
#include <iostream>

#define MUON_MASS    0.10565837
#define PION_MASS    0.13957018
#define KAON_MASS    0.493677
#define KSHORT_MASS  0.497614
#define KSTAR_MASS   0.89594
#define PHI_MASS     1.019455
#define JPSI_MASS    3.096916
#define PSI2S_MASS   3.686109
#define PROTON_MASS  0.938272046
#define LAMBDA_MASS  1.115683
#define BP_MASS      5.27926
#define B0_MASS      5.27958
#define BS_MASS      5.36677
#define BC_MASS      6.2756
#define LAMBDAB_MASS 5.6195

// HLT paths to be booked
enum {
    HLT_Dimuon0_Jpsi_Muon_v1,            // 0
    HLT_Dimuon10_Jpsi_Barrel_v1,         // 1
    HLT_Dimuon16_Jpsi_v1,                // 2
    HLT_Dimuon20_Jpsi_v1,                // 3
    HLT_DoubleMu4_3_Jpsi_Displaced_v1,   // 4
    HLT_DoubleMu4_JpsiTrk_Displaced_v1,  // 5
    N_HLT_BOOKINGS
};

const char HLT_paths[N_HLT_BOOKINGS][64] = {
    "HLT_Dimuon0_Jpsi_Muon_v1",          // 0
    "HLT_Dimuon10_Jpsi_Barrel_v1",       // 1
    "HLT_Dimuon16_Jpsi_v1",              // 2
    "HLT_Dimuon20_Jpsi_v1",              // 3
    "HLT_DoubleMu4_3_Jpsi_Displaced_v1", // 4
    "HLT_DoubleMu4_JpsiTrk_Displaced_v1" // 5
};

TH1D *binomialError(TH1D* num, TH1D* den)
{
    for (int i=1; i<=num->GetNbinsX(); i++) {
        double p = num->GetBinContent(i);
        double N = den->GetBinContent(i);
        num->SetBinError(i,sqrt(p*(1.-p)/N)); // binominal error
    }
    
    return num;
}

TH1D *dressing(TH1D* h)
{
    h->SetTitle("");
    h->GetXaxis()->SetTitle("p_{T} [GeV]");
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(1.09);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetLabelSize(0.055);
    h->GetXaxis()->SetTitleFont(42);
    h->GetYaxis()->SetTitle("Efficiency");
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetTitleOffset(1.14);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleFont(42);
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelSize(0.055);
    
    return h;
}

void EfficiencyPrivateSampleTrig()
{
    TChain *root = new TChain("demo/root");
    TChain *HltTree = new TChain("hltanalysis/HltTree");
    
    root->Add("/tmp/msharma/*.root");
    HltTree->Add("/tmp/msharma/*.root");
    
    int n_entries = root->GetEntries();
    if (n_entries!=HltTree->GetEntries()) {
        printf("Error: # of entries are different in two main trees.\n");
        return;
    }
    printf("Going to process %d entries.\n",n_entries);
    
    //-----------------------------------------------------------------
    // setting memory addresses to the branches
    // create new output trees
    
    EvtInfoBranches *EvtInfo = new EvtInfoBranches;
    GenInfoBranches *GenInfo = new GenInfoBranches;
    BInfoBranches *BInfo = new BInfoBranches;
    VtxInfoBranches *VtxInfo = new VtxInfoBranches;
    MuonInfoBranches *MuonInfo = new MuonInfoBranches;
    TrackInfoBranches *TrackInfo = new TrackInfoBranches;
    
    GenInfo->setbranchadd(root);
    EvtInfo->setbranchadd(root);
    VtxInfo->setbranchadd(root);
    MuonInfo->setbranchadd(root);
    TrackInfo->setbranchadd(root);
    BInfo->setbranchadd(root);
    
    ULong64_t HltTree_Event;
    int HltTree_Run;
    int HLT_book[N_HLT_BOOKINGS];
    HltTree->SetBranchAddress("Event",&HltTree_Event);
    HltTree->SetBranchAddress("Run",&HltTree_Run);
    for (int i=0;i<N_HLT_BOOKINGS;i++)
        HltTree->SetBranchAddress(HLT_paths[i],&HLT_book[i]);
    
    if (HltTree->GetListOfBranches()->FindObject("HLT_DoubleMu4_JpsiTrk_Displaced_v1"))
        cout<<"SUCCESS!!"<<endl;
    
    TFile *fout = new TFile("efficiency.root","recreate");
    
    double ptBins[] = {10., 15., 20., 25., 30., 40., 50., 60., 70., 80., 100.};
    int ptBinsSize = (sizeof(ptBins)/sizeof(*ptBins))-1;
    cout<<"pT bin size = "<<ptBinsSize<<'\t'<<ptBins[0]<<endl;
    
    TH1D* h_bppt_genLevel = new TH1D("h_bppt_genLevel","generator B+ p_{T}",40,10.,100.);
    TH1D* h_bpy_genLevel = new TH1D("h_bpy_genLevel","generator B+ rapidity",40,-2.4,2.4);
    TH1D* h_bppt_GL_coarseBins = new TH1D("h_bppt_GL_coarseBins", "gen-level B+ pT coarse bins", ptBinsSize, &ptBins[0]);
    
    TH1D* h_bppt_recoLevel = new TH1D("h_bppt_recoLevel","reco-level B+ p_{T}",40,10.,100.);
    TH1D* h_bpy_recoLevel = new TH1D("h_bpy_recoLevel","reco-level B+ rapidity",40,-2.4,2.4);
    TH1D* h_bppt_RL_coarseBins = new TH1D("h_bppt_RL_coarseBins", "reco-level B+ pT coarse bins", ptBinsSize, &ptBins[0]);
    
    TH2D* h_bpt_y_genlevel = new TH2D("h_bpt_y_genlevel", "pT vs y", 80,-2.4,2.4, 80,10.,100.);
    TH2D* h_bpt_y_recolevel = new TH2D("h_bpt_y_recolevel", "pT vs y", 80,-2.4,2.4, 80,10.,100.);
    
    TH1D* h_bpfilter_pt = new TH1D("h_bpfilter_pt","Filter denominator",40,10.,100.);
    TH1D* h_bpfilter_y = new TH1D("h_bpfilter_y","Filter denominator",40,-2.4,2.4);
    
    TH1D *h_bp_pt_effPrime_BE = new TH1D("h_bp_pt_effPrime_BE","Filter denominator",40,10.,100.);
    
    TH1D* h_cosang_sig = new TH1D("h_cosang_sig", "Signal cosine alpha", 100, 0., 1.);
    TH1D* h_vtxProb = new TH1D("h_vtxProb", "Vertex probability distribution", 100, 0., 1.);
    TH1D* h_decaylengthSig = new TH1D("h_decaylengthSig", "Decay length significance distribution", 100, 0., 10.);
    
    TH1D* h_Bmass_gen = new TH1D("h_Bmass_gen", "Gen-level J/#psi K+ inavriant mass distribution", 100, 5.0, 6.0);
    TH1D* h_Bmass = new TH1D("h_Bmass", "J/#psi K+ inavriant mass distribution", 100, 5.0, 6.0);
    TH1D* h_dR = new TH1D("h_dR", "Opening angle", 100, 0., 1.0);
    TH1D* h_ct = new TH1D("h_ct", "Life time distribution", 100, -0.1, 0.4);
    
    TH1D *bPt_WOtrig = new TH1D("bPt_WOtrig", "B-meson p_{T} distribution", ptBinsSize, &ptBins[0]);
    TH1D *bPt_Dimuon4trig = new TH1D("bPt_Dimuon4trig", "B-meson p_{T} distribution", ptBinsSize, &ptBins[0]);
    
    for (int evt=0; evt<n_entries; evt++) {
        if (evt%5000==0 || evt==n_entries-1) printf("processing %d/%d (%.2f%%).\n",evt,n_entries-1,(double)evt/(double)(n_entries-1)*100.);
        
        root->GetEntry(evt);
        HltTree->GetEntry(evt);
        
        // Vector of TVector3 for gen particle
        vector<TVector3> vect_gen;
        // Look for indices of the whole decay tree
        //-----------------------------------------------------------------
        // Calculating dR
        for (int idx = 0; idx < GenInfo->size; idx++) {
            if (abs(GenInfo->pdgId[idx])==521) {
                
                int idx_bp   = idx;
                TVector3 v3_bGen;
                v3_bGen.SetPtEtaPhi(GenInfo->pt[idx_bp],GenInfo->eta[idx_bp],GenInfo->phi[idx_bp]);
                vect_gen.push_back(v3_bGen);
                
                for (int bidx = 0; bidx < BInfo->size; bidx++) {
                    if (!(BInfo->type[bidx] == 1)) continue; //Only B+
                    TVector3 v3_bReco;
                    v3_bReco.SetPtEtaPhi(BInfo->pt[bidx],BInfo->eta[bidx],BInfo->phi[bidx]);
                    
                    if (BInfo->mass[bidx] >=5.16 && BInfo->mass[bidx] <= 5.365){
                        double dR = v3_bGen.Angle(v3_bReco);
                        h_dR->Fill(dR);
                    }
                }
            }
        }
        
        //-----------------------------------------------------------------
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
                
                bool muon1Acc = ((fabs(GenInfo->eta[idx_mu1])<1.3 && GenInfo->pt[idx_mu1]>3.3) ||
                                 (fabs(GenInfo->eta[idx_mu1])>1.3 && fabs(GenInfo->eta[idx_mu1])<2.2 && v4_mu1.P()>2.9) ||
                                 (fabs(GenInfo->eta[idx_mu1])>2.2 && fabs(GenInfo->eta[idx_mu1])<2.4 && GenInfo->pt[idx_mu1]>0.8));
                bool muon2Acc = ((fabs(GenInfo->eta[idx_mu2])<1.3 && GenInfo->pt[idx_mu2]>3.3) ||
                                 (fabs(GenInfo->eta[idx_mu2])>1.3 && fabs(GenInfo->eta[idx_mu2])<2.2 && v4_mu2.P()>2.9) ||
                                 (fabs(GenInfo->eta[idx_mu2])>2.2 && fabs(GenInfo->eta[idx_mu2])<2.4 && GenInfo->pt[idx_mu2]>0.8));

                bool kplusAcc = (GenInfo->pt[idx_kp]>0.9 && fabs(GenInfo->eta[idx_kp])<2.4);
                
                bool muon1Filter = fabs(GenInfo->eta[idx_mu1])<2.4;
                bool muon2Filter = fabs(GenInfo->eta[idx_mu2])<2.4;
                
                bool bRapFilter = fabs(v4_bp.Rapidity())<2.4;
                
                h_Bmass_gen->Fill(GenInfo->mass[idx_bp]);
                
                if (muon1Filter && muon2Filter && bRapFilter) {
                    h_bpfilter_pt->Fill(GenInfo->pt[idx_bp]);
                    h_bpfilter_y->Fill(v4_bp.Rapidity());
                }
                
                if (muon1Acc && muon2Acc && kplusAcc && GenInfo->mass[idx_bp] >=5.16 && GenInfo->mass[idx_bp] <= 5.365){
                    h_bppt_genLevel->Fill(GenInfo->pt[idx_bp]);
                    h_bpy_genLevel->Fill(v4_bp.Rapidity());
                    h_bpt_y_genlevel->Fill(v4_bp.Rapidity(), GenInfo->pt[idx_bp]);
                    h_bppt_GL_coarseBins->Fill(GenInfo->pt[idx_bp]);
                }

            }
        }
    
        int nMult_gen = (int)vect_gen.size();
    // Start of BInfo loop
    for (int bidx = 0; bidx < BInfo->size; bidx++) {
       
        if (!(BInfo->type[bidx] == 1)) continue; //Only B+
        //-----------------------------------------------------------------
        int ujidx = BInfo->rfuj_index[bidx];
        int tk1idx = BInfo->rftk1_index[bidx];
        int tk2idx = BInfo->rftk2_index[bidx];
        int mu1idx = BInfo->uj_rfmu1_index[ujidx];
        int mu2idx = BInfo->uj_rfmu2_index[ujidx];
        
        //-----------------------------------------------------------------
        TVector3 v3_muon1, v3_muon2;
        v3_muon1.SetPtEtaPhi(MuonInfo->pt[mu1idx],MuonInfo->eta[mu1idx],MuonInfo->phi[mu1idx]);
        v3_muon2.SetPtEtaPhi(MuonInfo->pt[mu2idx],MuonInfo->eta[mu2idx],MuonInfo->phi[mu2idx]);
        
        // Basic muon selections
        if (MuonInfo->pt[mu1idx]<=4.0) continue;
        if (MuonInfo->pt[mu2idx]<=4.0) continue;
        if (fabs(MuonInfo->eta[mu1idx])>=2.4) continue;
        if (fabs(MuonInfo->eta[mu2idx])>=2.4) continue;
        if (!MuonInfo->SoftMuID[mu1idx]) continue;
        if (!MuonInfo->SoftMuID[mu2idx]) continue;
    
        //-----------------------------------------------------------------
        // Basic track selections
        if (BInfo->type[bidx]==1) { // k, pi
            if (TrackInfo->pt[tk1idx]<=0.8) continue;
            if (fabs(TrackInfo->eta[tk1idx])>=2.5) continue;
            if (TrackInfo->chi2[tk1idx]/TrackInfo->ndf[tk1idx]>=5.) continue;
        }else { // others (2 tracks)
            if (TrackInfo->pt[tk1idx]<=0.7) continue;
            if (TrackInfo->pt[tk2idx]<=0.7) continue;
            if (fabs(TrackInfo->eta[tk1idx])>=2.5) continue;
            if (fabs(TrackInfo->eta[tk2idx])>=2.5) continue;
            if (TrackInfo->chi2[tk1idx]/TrackInfo->ndf[tk1idx]>=5.) continue;
            if (TrackInfo->chi2[tk2idx]/TrackInfo->ndf[tk2idx]>=5.) continue;
        }
        
        //-----------------------------------------------------------------
        // J/psi cut
        if (fabs(BInfo->uj_mass[ujidx]-JPSI_MASS)>=0.150) continue;
        if (BInfo->uj_pt[ujidx]<=8.0) continue;

        //-----------------------------------------------------------------
        // Find the best pointing PV
        //
        // keep the selecton code but PV is replaced with BS in the end.
        
        TVector3 bvtx(BInfo->vtxX[bidx],BInfo->vtxY[bidx],BInfo->vtxZ[bidx]);
        TVector3 bvtx_err(BInfo->vtxXE[bidx],BInfo->vtxYE[bidx],BInfo->vtxZE[bidx]);
        TVector3 bmom(BInfo->px[bidx],BInfo->py[bidx],BInfo->pz[bidx]);
        
        int vidx = -1;
        double max_cosang = -1.;
        for (int idx = 0; idx < VtxInfo->Size; idx++) {
            TVector3 vtx(VtxInfo->x[idx],VtxInfo->y[idx],VtxInfo->z[idx]);
            
            double cosang = bmom.Dot(bvtx-vtx)/(bmom.Mag()*(bvtx-vtx).Mag());
            if (cosang>max_cosang) {
                vidx = idx;
                max_cosang = cosang;
            }
        }
        
        if (vidx==-1) {
            printf("Error: no PV found. Run: %d, Event: %d.\n",EvtInfo->RunNo,EvtInfo->EvtNo);
            continue;
        }
        
        TVector3 BS(EvtInfo->PVx,EvtInfo->PVy,EvtInfo->PVz);
        TVector3 BS_err(EvtInfo->PVxE,EvtInfo->PVyE,EvtInfo->PVzE);
        TVector3 PV(VtxInfo->x[vidx],VtxInfo->y[vidx],VtxInfo->z[vidx]);
        TVector3 PV_err(VtxInfo->xE[vidx],VtxInfo->yE[vidx],VtxInfo->zE[vidx]);
        
        // Always use the beamspot
        PV = BS; PV_err = BS_err;
    
        TLorentzVector v4_uj;
        v4_uj.SetPtEtaPhiM(BInfo->uj_pt[ujidx],BInfo->uj_eta[ujidx],BInfo->uj_phi[ujidx],BInfo->uj_mass[ujidx]);
        
        TLorentzVector v4_b;
        v4_b.SetPtEtaPhiM(BInfo->pt[bidx],BInfo->eta[bidx],BInfo->phi[bidx],BInfo->mass[bidx]);
        
        TLorentzVector v4_tktk;
        TVector3 tktkvtx, tktkvtx_err;
        
        if (BInfo->type[bidx]!=1 && BInfo->type[bidx]!=2) { // other then K+, pi
            v4_tktk.SetPtEtaPhiM(BInfo->tktk_pt[bidx],BInfo->tktk_eta[bidx],BInfo->tktk_phi[bidx],BInfo->tktk_mass[bidx]);
            tktkvtx.SetXYZ(BInfo->tktk_vtxX[bidx],BInfo->tktk_vtxY[bidx],BInfo->tktk_vtxZ[bidx]);
            tktkvtx_err.SetXYZ(BInfo->tktk_vtxXE[bidx],BInfo->tktk_vtxYE[bidx],BInfo->tktk_vtxZE[bidx]);
        }
        
        TVector3 v3_bReco;
        v3_bReco.SetPtEtaPhi(BInfo->pt[bidx],BInfo->eta[bidx],BInfo->phi[bidx]);
        // Start to fill the B hadron information
        for (int ngen=0; ngen<nMult_gen; ngen++)
        {
            //-----------------------------------------------------------------
            // Topologocal cuts for B+
            double cosalpha2d = bmom.XYvector()*(bvtx-PV).XYvector()/(bmom.Perp()*(bvtx-PV).Perp());
            double lxy = (bvtx-PV).Perp();
            double errxy = sqrt(bvtx_err.Perp2()+PV_err.Perp2());
            double lxy_significance = lxy/errxy;
            double vtxprob = TMath::Prob(BInfo->vtxchi2[bidx],BInfo->vtxdof[bidx]);
            
            // ct calculations
            TVector3 v_l = bvtx-PV, v_lerr2;
            TVector3 v_p = v4_b.Vect();
            double ct = v_l.XYvector()*v_p.XYvector()*BP_MASS/v_p.Perp2();
            
            //-----------------------------------------------------------------
            TVector3 pvector_gen = (vect_gen)[ngen];
            
            double dR = pvector_gen.Angle(v3_bReco);
            if (dR > 0.011) continue;
            
            h_cosang_sig->Fill(cosalpha2d);
            h_decaylengthSig->Fill(lxy_significance);
            h_vtxProb->Fill(vtxprob);
            h_ct->Fill(ct);
            
            if (cosalpha2d<=0.99) continue;
            if (lxy_significance <=3.0) continue;
            if (vtxprob<=0.1) continue;
            
            h_Bmass->Fill(BInfo->mass[bidx]);
            
            if (BInfo->mass[bidx] >= 5.16 && BInfo->mass[bidx] <= 5.365)
                bPt_WOtrig->Fill(BInfo->pt[bidx]);
            
          //  for (int i=0;i<N_HLT_BOOKINGS;i++)
          //      cout<<HLT_book[i]<<endl;
            
            if (HLT_book[5]!=1) continue;
                bPt_Dimuon4trig->Fill(BInfo->pt[bidx]);
            
            if (BInfo->mass[bidx] >=5.16 && BInfo->mass[bidx] <= 5.365){
                h_bppt_recoLevel->Fill(BInfo->pt[bidx]);
                h_bpy_recoLevel->Fill(v4_b.Rapidity());
                h_bpt_y_recolevel->Fill(v4_b.Rapidity(), BInfo->pt[bidx]);
                h_bppt_RL_coarseBins->Fill(BInfo->pt[bidx]);
                }
            }
        }// End of BInfo loop
   } // end of evt loop
    
    TH1D *h_bp_pt_eff = (TH1D*)h_bppt_recoLevel->Clone("h_bp_pt_eff");
    h_bp_pt_eff->Divide(h_bppt_genLevel);
    binomialError(h_bp_pt_eff, h_bppt_genLevel);
    
    TH1D *h_bp_pt_eff_coarseBins = (TH1D*)h_bppt_RL_coarseBins->Clone("h_bp_pt_eff_coarseBins");
    h_bp_pt_eff_coarseBins->Divide(h_bppt_GL_coarseBins);
    binomialError(h_bp_pt_eff_coarseBins, h_bppt_GL_coarseBins);
    
    TH1D *h_bp_y_eff = (TH1D*)h_bpy_recoLevel->Clone("h_bp_y_eff");
    h_bp_y_eff->Divide(h_bpy_genLevel);
    binomialError(h_bp_y_eff, h_bpy_genLevel);
    
    TH1D *h_bppt_y_eff = (TH1D*)h_bpt_y_recolevel->Clone("h_bppt_y_eff");
    h_bppt_y_eff->Divide(h_bpt_y_genlevel);
    
    TH1D *hbPt_TrigEff = (TH1D*)bPt_Dimuon4trig->Clone("hbPt_TrigEff");
    hbPt_TrigEff->Divide(bPt_WOtrig);
    binomialError(hbPt_TrigEff, bPt_WOtrig);
    
    TEfficiency *pEff = 0;
    if(TEfficiency::CheckConsistency(*bPt_Dimuon4trig,*bPt_WOtrig))
    {
        pEff = new TEfficiency(*bPt_Dimuon4trig,*bPt_WOtrig);
        //this will attach the TEfficiency object to the current directory
        pEff->SetDirectory(gDirectory);
        //now all objects in gDirectory will be written to "myfile.root"
       // pFile->Write();
    }

    //-------------------------------------------------------------
    // Taking the pre-filter into account in the official MC sample
    
    TH1D *h_bp_pt_effPrime = (TH1D*)h_bppt_recoLevel->Clone("h_bp_pt_effPrime");
    h_bp_pt_effPrime->Divide(h_bpfilter_pt);
    binomialError(h_bp_pt_effPrime, h_bpfilter_pt);
    dressing(h_bp_pt_effPrime);
    
    fout->Write();
    fout->Close();
    
    delete GenInfo;
    delete EvtInfo;
    delete VtxInfo;
    delete MuonInfo;
    delete TrackInfo;
    delete BInfo;
}


