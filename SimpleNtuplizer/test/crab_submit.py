from CRABClient.UserUtilities import config

config = config()

#config.General.workArea                 = './'
config.General.workArea                 =  'crab_PFcluster'

config.General.transferLogs             = False

config.JobType.pluginName               = 'Analysis'
config.JobType.psetName                 = 'runAnalyzer.py'
config.JobType.allowUndistributedCMSSW  = True

#config.Data.inputDBS                    = 'global'
config.Data.inputDBS                    = 'phys03'
config.Data.splitting                   = 'FileBased'
config.Data.unitsPerJob                 = 100

#config.Site.storageSite                 = 'T3_US_FNALLPC'
config.Site.storageSite                 = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/phys_egamma/PFClusterCalibration/MC18_V2/FlatTrees'

if __name__ == '__main__':

    datasets = {#'PFclusterPU0p01To5' : '/DoublePhotonNoMaterial_FlatPt-0p01To5/rcoelhol-crab_PFCluster_PU_0p01To5-36e39db79e677743e80137e00b85de2e/USER'
        'PFClusterPU0p01To5' : '/DoublePhotonNoMaterial_FlatPt-0p01To5/shilpi-crab_DoublePhotonNoMaterial_FlatPt-0p01To5-37b8bb60b9584eb61e5ec36b3b302521/USER'
#        ,'PFClusterPU5To300'  : '/DoublePhotonNoMaterial_FlatPt-5To300/shilpi-crab_DoublePhotonNoMaterial_FlatPt-5To300-37b8bb60b9584eb61e5ec36b3b302521/USER'
#        'PFClusterPU5To300'  : '/DoublePhotonNoMaterial_FlatPt-5To300/shilpi-crab_DoublePhotonNoMaterial_FlatPt-5To300-37b8bb60b9584eb61e5ec36b3b302521/USER'
        }

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
