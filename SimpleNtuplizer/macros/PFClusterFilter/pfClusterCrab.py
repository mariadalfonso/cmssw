from CRABClient.UserUtilities import config
config = config()

# 20 June: First test of saturation variables
config.General.requestName              = 'HggClusters_v2'

config.General.workArea                 = './'
config.General.transferOutputs          = True
config.General.transferLogs             = False

config.JobType.pluginName               = 'Analysis'
config.JobType.psetName                 = 'pfClusterForCalibration.py'

config.Data.inputDataset                = '/GluGluHToGG_M-125_13TeV_powheg_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'
config.Data.secondaryInputDataset       = '/GluGluHToGG_M-125_13TeV_powheg_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/GEN-SIM-RAW'
config.Data.inputDBS                    = 'global'
config.Data.splitting                   = 'FileBased'
config.Data.unitsPerJob                 = 1
config.Data.publication                 = True

config.Site.storageSite                 = 'T3_US_FNALLPC'


