from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '$NAME$'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../rerunTauRecoOnMiniAOD_WithClean_Custom.py'
config.JobType.numCores = 4
config.JobType.maxMemoryMB = 10000

config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset ='$INPUTDS$'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.outLFNDirBase = '/store/user/nbower/Events/'
config.Data.publication = True
config.Data.outputDatasetTag = '$NAME$'
#config.Data.ignoreLocality = True                                                                                                                                                                                                                                                                                                                                        
#config.Site.whitelist = ["T*_FR_*", "T*_DE_*", "T*_CH_*"]                                                                                                                                                                                                                                                                                                                
config.Site.storageSite = 'T2_US_Florida'
#config.Site.blacklist = ['T3_KR_KNU', 'T3_FR_IPNL', 'T2_TR_METU', 'T2_TW_NCHC', 'T2_BE_IIHE', 'T3_US_Baylor']                                                                                                                                                                                                                                                            
config.Site.blacklist = ['T3_KR_KNU', 'T3_FR_IPNL', 'T2_TR_METU', 'T2_TW_NCHC', 'T2_BE_IIHE', 'T3_US_Baylor','T2_US_Purdue']
