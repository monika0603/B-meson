
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

class ReducedBranches{
    public:
    int     run;
    int     event;

    int     type;  // B hadron information
    double  mass;
    double  pt;
    double  eta;
    double  phi;
    double  y;
    double  genmass;
    double  genpt;
    double  geneta;
    double  genphi;
    double  geny;
    double  vx;
    double  vy;
    double  vz;
    double  lxy;
    double  lxyz;
    double  errxy;
    double  errxyz;
    double  vtxprob;
    double  cosalpha2d;
    double  cosalpha3d;
    double  ctau2d;
    double  ctau3d;
    double  ctau2derr;
    double  ctau3derr;
    
    double  ujmass; // dimuon information
    double  ujpt;
    double  ujeta;
    double  ujphi;
    double  ujy;
    double  ujvtxprob;
    
    double  tktkmass; // ditrack information
    double  tktkpt;
    double  tktketa;
    double  tktkphi;
    double  tktky;
    double  tktkvtxprob;
    double  tktklxy;
    double  tktklxyz;
    double  tktkerrxy;
    double  tktkerrxyz;
    double  tktkblxy;
    double  tktkblxyz;
    double  tktkberrxy;
    double  tktkberrxyz;
    
    int     mu1idx;
    double  mu1pt;
    double  mu1eta;
    double  mu1phi;
    int     mu2idx;
    double  mu2pt;
    double  mu2eta;
    double  mu2phi;

    int     tk1idx;
    double  tk1pt;
    double  tk1eta;
    double  tk1phi;
    int     tk2idx;
    double  tk2pt;
    double  tk2eta;
    double  tk2phi;
    
    int nhltbook; // triggers
    int hltbook[N_HLT_BOOKINGS];
    
    std::vector<std::vector<double > > *MuTrgObjPt; //Trigger matching
    std::vector<std::vector<double > > *MuTrgObjEta;
    std::vector<std::vector<double > > *MuTrgObjPhi;
    
    
    void regTree(TTree *root){
        root->Branch("run",&run,"run/I");
        root->Branch("event",&event,"event/I");
        root->Branch("type",&type,"type/I");
        root->Branch("mass",&mass,"mass/D");
        root->Branch("pt",&pt,"pt/D");
        root->Branch("eta",&eta,"eta/D");
        root->Branch("phi",&phi,"phi/D");
        root->Branch("y",&y,"y/D");
        root->Branch("genmass",&genmass,"genmass/D");
        root->Branch("genpt",&genpt,"genpt/D");
        root->Branch("geneta",&geneta,"geneta/D");
        root->Branch("genphi",&genphi,"genphi/D");
        root->Branch("geny",&geny,"geny/D");
        root->Branch("vx",&vx,"vx/D");
        root->Branch("vy",&vy,"vy/D");
        root->Branch("vz",&vz,"vz/D");
        root->Branch("lxy",&lxy,"lxy/D");
        root->Branch("lxyz",&lxyz,"lxyz/D");
        root->Branch("errxy",&errxy,"errxy/D");
        root->Branch("errxyz",&errxyz,"errxyz/D");
        root->Branch("vtxprob",&vtxprob,"vtxprob/D");
        root->Branch("cosalpha2d",&cosalpha2d,"cosalpha2d/D");
        root->Branch("cosalpha3d",&cosalpha3d,"cosalpha3d/D");
        root->Branch("ctau2d",&ctau2d,"ctau2d/D");
        root->Branch("ctau3d",&ctau3d,"ctau3d/D");
        root->Branch("ctau2derr",&ctau2derr,"ctau2derr/D");
        root->Branch("ctau3derr",&ctau3derr,"ctau3derr/D");
        root->Branch("ujmass",&ujmass,"ujmass/D");
        root->Branch("ujpt",&ujpt,"ujpt/D");
        root->Branch("ujeta",&ujeta,"ujeta/D");
        root->Branch("ujphi",&ujphi,"ujphi/D");
        root->Branch("ujy",&ujy,"ujy/D");
        root->Branch("ujvtxprob",&ujvtxprob,"ujvtxprob/D");
        root->Branch("tktkmass",&tktkmass,"tktkmass/D");
        root->Branch("tktkpt",&tktkpt,"tktkpt/D");
        root->Branch("tktketa",&tktketa,"tktketa/D");
        root->Branch("tktkphi",&tktkphi,"tktkphi/D");
        root->Branch("tktky",&tktky,"tktky/D");
        root->Branch("tktkvtxprob",&tktkvtxprob,"tktkvtxprob/D");
        root->Branch("tktklxy",&tktklxy,"tktklxy/D");
        root->Branch("tktklxyz",&tktklxyz,"tktklxyz/D");
        root->Branch("tktkerrxy",&tktkerrxy,"tktkerrxy/D");
        root->Branch("tktkerrxyz",&tktkerrxyz,"tktkerrxyz/D");
        root->Branch("tktkblxy",&tktkblxy,"tktkblxy/D");
        root->Branch("tktkblxyz",&tktkblxyz,"tktkblxyz/D");
        root->Branch("tktkberrxy",&tktkberrxy,"tktkberrxy/D");
        root->Branch("tktkberrxyz",&tktkberrxyz,"tktkberrxyz/D");
        root->Branch("mu1idx",&mu1idx,"mu1idx/I");
        root->Branch("mu1pt",&mu1pt,"mu1pt/D");
        root->Branch("mu1eta",&mu1eta,"mu1eta/D");
        root->Branch("mu1phi",&mu1phi,"mu1phi/D");
        root->Branch("mu2idx",&mu2idx,"mu2idx/I");
        root->Branch("mu2pt",&mu2pt,"mu2pt/D");
        root->Branch("mu2eta",&mu2eta,"mu2eta/D");
        root->Branch("mu2phi",&mu2phi,"mu2phi/D");
        root->Branch("tk1idx",&tk1idx,"tk1idx/I");
        root->Branch("tk1pt",&tk1pt,"tk1pt/D");
        root->Branch("tk1eta",&tk1eta,"tk1eta/D");
        root->Branch("tk1phi",&tk1phi,"tk1phi/D");
        root->Branch("tk2idx",&tk2idx,"tk2idx/I");
        root->Branch("tk2pt",&tk2pt,"tk2pt/D");
        root->Branch("tk2eta",&tk2eta,"tk2eta/D");
        root->Branch("tk2phi",&tk2phi,"tk2phi/D");
        root->Branch("nhltbook",&nhltbook,"nhltbook/I");
        root->Branch("hltbook",hltbook,"hltbook[nhltbook]/I");
        root->Branch("MuTrgObjEta","std::vector<std::vector<double>>",&MuTrgObjEta);
        root->Branch("MuTrgObjPhi","std::vector<std::vector<double>>",&MuTrgObjPhi);
        root->Branch("MuTrgObjPt","std::vector<std::vector<double>>",&MuTrgObjPt);
    }
    void setbranchadd(TTree *root){
        root->SetBranchAddress("run",&run);
        root->SetBranchAddress("event",&event);
        root->SetBranchAddress("type",&type);
        root->SetBranchAddress("mass",&mass);
        root->SetBranchAddress("pt",&pt);
        root->SetBranchAddress("eta",&eta);
        root->SetBranchAddress("phi",&phi);
        root->SetBranchAddress("y",&y);
        root->SetBranchAddress("genmass",&genmass);
        root->SetBranchAddress("genpt",&genpt);
        root->SetBranchAddress("geneta",&geneta);
        root->SetBranchAddress("genphi",&genphi);
        root->SetBranchAddress("geny",&geny);
        root->SetBranchAddress("vx",&vx);
        root->SetBranchAddress("vy",&vy);
        root->SetBranchAddress("vz",&vz);
        root->SetBranchAddress("lxy",&lxy);
        root->SetBranchAddress("lxyz",&lxyz);
        root->SetBranchAddress("errxy",&errxy);
        root->SetBranchAddress("errxyz",&errxyz);
        root->SetBranchAddress("vtxprob",&vtxprob);
        root->SetBranchAddress("cosalpha2d",&cosalpha2d);
        root->SetBranchAddress("cosalpha3d",&cosalpha3d);
        root->SetBranchAddress("ctau2d",&ctau2d);
        root->SetBranchAddress("ctau3d",&ctau3d);
        root->SetBranchAddress("ctau2derr",&ctau2derr);
        root->SetBranchAddress("ctau3derr",&ctau3derr);
        root->SetBranchAddress("ujmass",&ujmass);
        root->SetBranchAddress("ujpt",&ujpt);
        root->SetBranchAddress("ujeta",&ujeta);
        root->SetBranchAddress("ujphi",&ujphi);
        root->SetBranchAddress("ujy",&ujy);
        root->SetBranchAddress("ujvtxprob",&ujvtxprob);
        root->SetBranchAddress("tktkmass",&tktkmass);
        root->SetBranchAddress("tktkpt",&tktkpt);
        root->SetBranchAddress("tktketa",&tktketa);
        root->SetBranchAddress("tktkphi",&tktkphi);
        root->SetBranchAddress("tktky",&tktky);
        root->SetBranchAddress("tktkvtxprob",&tktkvtxprob);
        root->SetBranchAddress("tktklxy",&tktklxy);
        root->SetBranchAddress("tktklxyz",&tktklxyz);
        root->SetBranchAddress("tktkerrxy",&tktkerrxy);
        root->SetBranchAddress("tktkerrxyz",&tktkerrxyz);
        root->SetBranchAddress("tktkblxy",&tktkblxy);
        root->SetBranchAddress("tktkblxyz",&tktkblxyz);
        root->SetBranchAddress("tktkberrxy",&tktkberrxy);
        root->SetBranchAddress("tktkberrxyz",&tktkberrxyz);
        root->SetBranchAddress("mu1idx",&mu1idx);
        root->SetBranchAddress("mu1pt",&mu1pt);
        root->SetBranchAddress("mu1eta",&mu1eta);
        root->SetBranchAddress("mu1phi",&mu1phi);
        root->SetBranchAddress("mu2idx",&mu2idx);
        root->SetBranchAddress("mu2pt",&mu2pt);
        root->SetBranchAddress("mu2eta",&mu2eta);
        root->SetBranchAddress("mu2phi",&mu2phi);
        root->SetBranchAddress("tk1idx",&tk1idx);
        root->SetBranchAddress("tk1pt",&tk1pt);
        root->SetBranchAddress("tk1eta",&tk1eta);
        root->SetBranchAddress("tk1phi",&tk1phi);
        root->SetBranchAddress("tk2idx",&tk2idx);
        root->SetBranchAddress("tk2pt",&tk2pt);
        root->SetBranchAddress("tk2eta",&tk2eta);
        root->SetBranchAddress("tk2phi",&tk2phi);
        root->SetBranchAddress("nhltbook",&nhltbook);
        root->SetBranchAddress("hltbook",hltbook);
        
        MuTrgObjPt= new std::vector<std::vector<double > >();                                    \
        MuTrgObjEta= new std::vector<std::vector<double > >();
        MuTrgObjPhi= new std::vector<std::vector<double > >();
        
        root->SetBranchAddress("MuTrgObjEta",&MuTrgObjEta);
        root->SetBranchAddress("MuTrgObjPhi",&MuTrgObjPhi);
        root->SetBranchAddress("MuTrgObjPt",&MuTrgObjPt);
    }
};
