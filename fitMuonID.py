import FWCore.ParameterSet.Config as cms

### USAGE:
###    cmsRun fitMuonID.py <scenario>
### scenarios:
###   - data_all (default)
###   - signal_mc

import sys
args = sys.argv[1:]
if (sys.argv[0] == "cmsRun"): args =sys.argv[2:]
scenario = "data_all"
if len(args) > 0: scenario = args[0]
print "Will run scenario ", scenario

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

Template = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-muon Mass", "2.9", "3.3", "GeV/c^{2}"), #2.8-3.35
        p  = cms.vstring("muon p", "0", "1000", "GeV/c"),
        pt = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
        eta = cms.vstring("muon #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("muon |#eta|", "0", "2.5", ""),
        tag_pt = cms.vstring("Tag p_{T}",    "0", "1000", "GeV/c"),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
        tag_nVerticesDA = cms.vstring("Number of vertices", "0", "999", ""),
        pair_dphiVtxTimesQ = cms.vstring("q1 * (#phi1-#phi2)", "-6", "6", ""),
        pair_drM1   = cms.vstring("#Delta R(Station 1)", "-99999", "999999", "rad"),
        pair_distM1 = cms.vstring("dist(Station 1)", "-99999", "999999", "cm"),
        pair_dz = cms.vstring("dz","-5","5",""),
        pair_probeMultiplicity = cms.vstring("multiplicity","0","99",""),
        dB = cms.vstring("dB", "-1000", "1000", ""),
        dzPV = cms.vstring("dzPV", "-1000", "1000", ""),
        dxyBS = cms.vstring("dxyBS", "-1000", "1000", ""),
        tkValidHits = cms.vstring("track.numberOfValidHits", "-1", "999", ""),
        tkTrackerLay = cms.vstring("track.hitPattern.trackerLayersWithMeasurement", "-1", "999", ""),
        tkValidPixelHits = cms.vstring("track.hitPattern.numberOfValidPixelHits", "-1", "999", ""),
        tkPixelLay = cms.vstring("track.hitPattern.pixelLayersWithMeasurement", "-1", "999", ""),
        tkChi2 = cms.vstring("track.normalizedChi2", "-1", "999", ""),
        numberOfMatchedStations = cms.vstring("numberOfMatchedStations", "-1", "99", ""),
        glbChi2 = cms.vstring("global.normalizedChi2", "-9999", "9999", ""),
        glbValidMuHits = cms.vstring("globalTrack.numberOfValidMuonHits", "-1", "9999", ""),
        caloComp = cms.vstring("caloCompatibility","-1","5",""),
    ),

    Categories = cms.PSet(
        TM   = cms.vstring("Tracker muon", "dummy[pass=1,fail=0]"),
        TMA   = cms.vstring("Tracker muon", "dummy[pass=1,fail=0]"),
        Glb   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        VBTF  = cms.vstring("VBTFLike", "dummy[pass=1,fail=0]"),
        TMOST = cms.vstring("TMOneStationTight", "dummy[pass=1,fail=0]"),
        PF = cms.vstring("PF", "dummy[pass=1,fail=0]"),
        Track_HP = cms.vstring("Track_HP", "dummy[pass=1,fail=0]"),
        Tight2012 = cms.vstring("Tight Muon", "dummy[pass=1,fail=0]"),
        Mu5_Track2_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu7_Track7_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu5_Track2_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu7_Track7_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        Mu5_Track0_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu3_Track3_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        Mu5_Track5_Jpsi_TK = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu5_Track0_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu3_Track3_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
        tag_Mu5_Track5_Jpsi_MU = cms.vstring("ProbeTrigger_Track0", "dummy[pass=1,fail=0]"),
    ),


    Expressions = cms.PSet(
     LooseVar = cms.vstring("LooseVar", "PF==1 && (Glb==1 || TM==1) ", "PF", "Glb", "TM"),
     SoftVar = cms.vstring("SoftVar", "TMOST ==1 && tkTrackerLay > 5 && tkPixelLay > 1 && tkChi2 < 1.8 && abs(dzPV) < 30 && abs(dB) < 3", "TMOST","tkTrackerLay", "tkPixelLay", "tkChi2", "dzPV", "dB"),
     #without chi2 cut
     #SoftVar = cms.vstring("SoftVar", "TMOST ==1 && tkTrackerLay > 5 && tkPixelLay > 1 && abs(dzPV) < 30 && abs(dB) < 3", "TMOST","tkTrackerLay", "tkPixelLay", "dzPV", "dB"),
     #changed cuts
     #SoftVar = cms.vstring("SoftVar", "TMOST ==1", "TMOST"),
     #SoftVar = cms.vstring("SoftVar", "TMOST ==1 && abs(dzPV) < 20", "TMOST", "dzPV"),
     #SoftVar = cms.vstring("SoftVar", "TMOST ==1 && abs(dzPV) < 20 && abs(dB) < 0.3", "TMOST", "dzPV", "dB"),
     #SoftVar = cms.vstring("SoftVar", "TMOST ==1 && abs(dzPV) < 20 && abs(dB) < 0.3 && tkPixelLay > 0", "TMOST", "dzPV", "dB", "tkPixelLay"),
     #SoftVar = cms.vstring("SoftVar", "TMOST ==1 && abs(dzPV) < 20 && abs(dB) < 0.3 && tkPixelLay > 0 && tkTrackerLay > 5", "TMOST", "dzPV", "dB", "tkPixelLay", "tkTrackerLay"),
     #newID
     newSoftVar = cms.vstring("newSoftVar", "TMOST ==1 && tkTrackerLay > 5 && tkPixelLay > 0 && abs(dzPV) < 20 && abs(dB) < 0.3 && Track_HP == 1", "TMOST","tkTrackerLay", "tkPixelLay", "dzPV", "dB", "Track_HP"),
     ),

    Cuts = cms.PSet(
     Loose2012 = cms.vstring("Loose", "LooseVar", "0.5"),
     Soft2012 = cms.vstring("Soft", "SoftVar", "0.5"),
     newSoft2012 = cms.vstring("newSoft", "newSoftVar", "0.5"),
    ),

    PDFs = cms.PSet(
        gaussPlusExpo = cms.vstring(
            "CBShape::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.06], alpha[3., 0.5, 5.], n[1, 0.1, 100.])",
            #"Chebychev::backgroundPass(mass, {cPass[0,-0.5,0.5], cPass2[0,-0.5,0.5]})",
            #"Chebychev::backgroundFail(mass, {cFail[0,-0.5,0.5], cFail2[0,-0.5,0.5]})",
            #"Gaussian::signal(mass, mean[3.1,3.0,3.2], sigma[0.05,0.02,0.1])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        )
    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),

    Efficiencies = cms.PSet(), # will be filled later
)

# pick muons that bend apart from each other
SEPARATED = cms.PSet(pair_drM1 = cms.vdouble(0.5,10),
                     pair_probeMultiplicity = cms.vdouble(0.5,1.5),
                     #pair_dz = cms.vdouble(-0.1,0.1),
                     #caloComp = cms.vdouble(0.5,5),
                     )
#                     pair_distM1 = cms.vdouble(200,1000))
#SEPARATED = cms.PSet(pair_dphiVtxTimesQ = cms.vdouble(-3.14,0), #seagulls
#                     pair_drM1 = cms.vdouble(0.5,10),
#                     pair_distM1 = cms.vdouble(200,1000))
#SEPARATED = cms.PSet(pair_dphiVtxTimesQ = cms.vdouble(0,3.14),) #cowboys
#                     pair_drM1 = cms.vdouble(0.5,10),
#                     pair_distM1 = cms.vdouble(200,1000))

PT_ETA_BINS = cms.PSet(SEPARATED,
                       #pt     = cms.vdouble(2, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 9.0, 11.0, 14.0, 17.0, 20.0, 25., 30., 35., 40.),
                       pt = cms.vdouble(2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0),
                       abseta = cms.vdouble(0.0,0.9)
                       #pt = cms.vdouble(2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 20.0),
                       #abseta = cms.vdouble(0.9,1.2)
                       #abseta = cms.vdouble(1.2,2.1)
                       #abseta = cms.vdouble(0.0, 0.9, 1.2, 2.1)
)

VTX_BINS = cms.PSet(SEPARATED,
                    abseta = cms.vdouble(0.0, 2.1),
                    pt     = cms.vdouble(8.0, 20.0),
                    #pt     = cms.vdouble(7.0, 20.0),
                    tag_nVertices = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5)
                    #tag_nVerticesDA = cms.vdouble(0.5,2.5,4.5,6.5,8.5,10.5,12.5,14.5,16.5,18.5,20.5,22.5,24.5,26.5,28.5,30.5)
)

ETA_BINS = cms.PSet(SEPARATED,
                    pt     = cms.vdouble(8.0, 20.0),
                    #pt     = cms.vdouble(7.0, 20.0),
                    eta = cms.vdouble(-2.1,-1.6,-1.2,-0.9,-0.6,-0.3,-0.2,0.2,0.3,0.6,0.9,1.2,1.6,2.1,),
                    )

PLATEAU_ABSETA = cms.PSet(SEPARATED,
                          abseta = cms.vdouble(0.0, 0.9, 1.2, 2.1),
                          pt     = cms.vdouble(8.0, 20.0),
                          #pt     = cms.vdouble(4.0, 8.0),
                          #pt     = cms.vdouble(7.0, 20.0),
                          #pt     = cms.vdouble(20.0, 30.0),
)

PT_ABSETA_WIDE = cms.PSet(SEPARATED,
    abseta = cms.vdouble(0.0, 1.2, 2.4),
    pt     = cms.vdouble(5.0, 7.0, 20.0),
)


PREFIX="root://eoscms//eos/cms/store/caf/user/gpetrucc/TnP/V5/"
process.TnP_MuonID = Template.clone(
    InputFileNames = cms.vstring(
        PREFIX+'tnpJPsi_Run2012A.root',
        PREFIX+'tnpJPsi_Run2012B.root',
        PREFIX+'tnpJPsi_Run2012C.root',
        PREFIX+'tnpJPsi_Run2012D.root',
    ),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string("TnP_MuonID_%s.root" % scenario),
    Efficiencies = cms.PSet(),
)

IDS = ["newSoft2012"]
#IDS = ["Soft2012"]
#IDS = ["Loose2012"]
#IDS = ["Loose2012", "Soft2012", "newSoft2012"]
#IDS = [ "Glb", "TMOST", "VBTF", "PF" ]
TRIGS = [ (2,'Mu5_Track2'), (7,'Mu7_Track7') ]

if "mc" in scenario:
     process.TnP_MuonID.InputFileNames = [PREFIX+'tnpJPsi_MC53X.root']

#ALLBINS =  [("plateau_abseta",PLATEAU_ABSETA)]
ALLBINS =  [("pt_abseta",PT_ETA_BINS)]
#ALLBINS =  [("vtx",VTX_BINS)]
#ALLBINS =  [("plateau_abseta",PLATEAU_ABSETA), ("vtx",VTX_BINS), ("eta",ETA_BINS)]
#ALLBINS =  [("pt_abseta",PT_ETA_BINS), ("vtx",VTX_BINS), ("eta",ETA_BINS)]
#ALLBINS =  [("pt_abseta",PT_ETA_BINS), ("vtx",VTX_BINS), ("plateau",PLATEAU_ABSETA)]
#ALLBINS += [("pt_abseta_wide",PT_ABSETA_WIDE)]

for ID in IDS:
     if len(args) > 1 and args[1] in IDS and ID != args[1]: continue
     for X,B in ALLBINS:
          if len(args) > 2 and X not in args[2:]: continue
          module = process.TnP_MuonID.clone(OutputFileName = cms.string(#"TnP_MuonID_%s_%s_%s.root" % (scenario, ID, X)))
               "/scratch/ikratsch/TnP2012/MuonPOG/official6March2014/changedMass/TnP_MuonID_%s_%s_%s0_multiplicity_noTrigger.root" % (scenario, ID, X)))
          for PTMIN, TRIG in TRIGS:
               TRIGLABEL=""
               if "pt_" in X:
                    TRIGLABEL="_"+TRIG
               else:
                    if TRIG != "Mu7_Track7": continue # use only one trigger except for turn-on
                    #if TRIG != "Mu5_Track2": continue # use only one trigger except for turn-on
               DEN = B.clone()
               if hasattr(DEN, "pt"):
                    DEN.pt = cms.vdouble(*[i for i in B.pt if i >= PTMIN])
                    if len(DEN.pt) == 1: DEN.pt = cms.vdouble(PTMIN, DEN.pt[0])
               #setattr(DEN, "tag_%s_Jpsi_MU" % TRIG, cms.vstring("pass"))
               #setattr(DEN,     "%s_Jpsi_TK" % TRIG, cms.vstring("pass"))
               #setattr(DEN, "TM", cms.vstring("pass"))
               #if "calomu" in scenario: DEN.Calo = cms.vstring("pass")
               setattr(module.Efficiencies, ID+"_"+X+TRIGLABEL, cms.PSet(
                   EfficiencyCategoryAndState = cms.vstring(ID,"above"),
                   UnbinnedVariables = cms.vstring("mass"),
                   BinnedVariables = DEN,
                   BinToPDFmap = cms.vstring("gaussPlusExpo")
               ))
               i#if "plateau" in X: module.SaveWorkspace = True
               ## mc efficiency, if scenario is mc
               if "mc" in scenario:
                    setattr(module.Efficiencies, ID+"_"+X+TRIGLABEL+"_mcTrue", cms.PSet(
                        EfficiencyCategoryAndState = cms.vstring(ID,"above"),
                        UnbinnedVariables = cms.vstring("mass"),
                        BinnedVariables = DEN.clone(mcTrue = cms.vstring("true"))
                    ))
          setattr(process, "TnP_MuonID_"+ID+"_"+X, module)
          setattr(process, "run_"+ID+"_"+X, cms.Path(module))

