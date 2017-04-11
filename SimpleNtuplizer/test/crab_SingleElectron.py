from CRABClient.UserUtilities import config

config = config()

config.General.workArea                 = './'
config.General.transferOutputs          = True
config.General.transferLogs             = True

config.JobType.pluginName               = 'Analysis'
config.JobType.psetName                 = 'runAnalyzer.py'
config.JobType.allowUndistributedCMSSW  = True
config.JobType.sendExternalFolder       = True

config.Data.inputDBS                    = 'global'
config.Data.splitting                   = 'LumiBased'
config.Data.unitsPerJob                 = 50
config.Data.lumiMask                    = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.publication                 = False

config.Site.storageSite                 = 'T3_US_FNALLPC'

if __name__ == '__main__':

    datasets = {'2016rereco_RunBv2':'/SingleElectron/Run2016B-23Sep2016-v2/AOD','2016rereco_RunBv3':'/SingleElectron/Run2016B-23Sep2016-v3/AOD','2016rereco_RunC':'/SingleElectron/Run2016C-23Sep2016-v1/AOD','2016rereco_RunD':'/SingleElectron/Run2016D-23Sep2016-v1/AOD','2016rereco_RunE':'/SingleElectron/Run2016E-23Sep2016-v1/AOD','2016rereco_RunF':'/SingleElectron/Run2016F-23Sep2016-v1/AOD','2016rereco_RunG':'/SingleElectron/Run2016G-23Sep2016-v1/AOD','2016prompt_RunHv1':'/SingleElectron/Run2016H-PromptReco-v1/AOD','2016prompt_RunHv2':'/SingleElectron/Run2016H-PromptReco-v2/AOD','2016prompt_RunHv3':'/SingleElectron/Run2016H-PromptReco-v3/AOD'}

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
