from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.requestName = 'TagAndProbeFirstAttempt'
config.General.workArea = ''
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'stepITagAndProbeTreeProducer_JPsiMuMu.py'
config.JobType.outputFiles = ['tnpJpsi_MC.root']

config.section_('Data')
config.Data.inputDataset = '/JpsiToMuMu_JPsiPt7WithFSR_13TeV-pythia6-evtgen/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v2/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
#config.Data.lumiMask = 'Cert_181530-183126_HI7TeV_PromptReco_Collisions11_JSON.txt'
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_US_Vanderbilt'