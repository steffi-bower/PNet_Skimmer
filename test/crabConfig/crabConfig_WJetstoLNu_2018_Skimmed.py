from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'MiniAOD_WJetstoLNu_2018'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../rerunTauRecoOnMiniAOD_WithClean_Custom.py'
config.JobType.numCores = 4
config.JobType.maxMemoryMB = 10000

config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset ='/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.outLFNDirBase = '/store/user/nbower/Events/WJetstoLNu_2018_Skimmed/'
config.Data.publication = True
config.Data.outputDatasetTag = 'WJetstoLNu_2018_Skimmed'
#config.Data.ignoreLocality = True
#config.Site.whitelist = ["T*_FR_*", "T*_DE_*", "T*_CH_*"]
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.blacklist = ['T3_KR_KNU', 'T3_FR_IPNL', 'T2_TR_METU', 'T2_TW_NCHC', 'T2_BE_IIHE', 'T3_US_Baylor']
config.Site.blacklist = ['T3_KR_KNU', 'T3_FR_IPNL', 'T2_TR_METU', 'T2_TW_NCHC', 'T2_BE_IIHE', 'T3_US_Baylor','T2_US_Purdue']
