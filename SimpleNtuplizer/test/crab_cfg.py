from CRABClient.UserUtilities import config

config = config()

config.General.workArea                 = './'
config.General.transferOutputs          = True
config.General.transferLogs             = True

config.JobType.pluginName               = 'Analysis'
config.JobType.psetName                 = 'runAnalyzer.py'
config.JobType.allowUndistributedCMSSW  = True

config.Data.inputDBS                    = 'global'
config.Data.splitting                   = 'FileBased'
config.Data.unitsPerJob                 = 1

config.Site.storageSite                 = 'T3_US_FNALLPC'

if __name__ == '__main__':

    datasets = {'ElectronLowPt' : '/DoubleElectron_FlatPt-1To300/RunIISpring16DR80-PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
                'ElectronHighPt' : '/DoubleElectron_FlatPt-300To6500/RunIISpring16DR80-PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
                'PhotonLowPt' : '/DoublePhoton_FlatPt-5To300/RunIISpring16DR80-PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
                'PhotonHighPt' : '/DoublePhoton_FlatPt-300To6500/RunIISpring16DR80-PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'}

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    for name, dataset in datasets.iteritems():
        config.General.requestName = name
        config.Data.inputDataset = dataset
        submit(config)
