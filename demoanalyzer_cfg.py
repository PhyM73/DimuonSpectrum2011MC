import FWCore.ParameterSet.Config as cms
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
process = cms.Process("Demo")

# intialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))
# **********************************************************************
# set the maximum number of events to be processed                     *
#    this number (argument of int32) is to be modified by the user     *
#    according to need and wish                                        *
#    default is preset to 10000 events                                 *
# **********************************************************************
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))

# set the number of events to be skipped (if any) at end of file below

# ****************************************************************************
# define the input data set here by inserting the appropriate .txt file list *
# ****************************************************************************
import FWCore.Utilities.FileUtils as FileUtils

#
# ****************************************************************
# load the data set                                              *
# useful datasets are SingleMu and DoubleMu (default)            *
# To run over all data subsets, replace '10000' by '10001' etc.  *
# consecutively (make sure you save the output before rerunning) *
# and add up the histograms using root tools.                    *
# ****************************************************************
#
# *** MonteCarlo data set ***
mcfilelist = []
# files2011data = FileUtils.loadListFromFile ('datasets/CMS_Run2011A_DoubleMu_AOD_12Oct2013-v1_10000_file_index.txt')
datasets = FileUtils.os.walk(r"./data")
for path, dir_list, file_list in datasets:
    for indexfile in file_list:
        print(indexfile)
        mcfilelist.extend(FileUtils.loadListFromFile(FileUtils.os.path.join(path, indexfile)))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(*mcfilelist))

# *************************************************
# number of events to be skipped (0 by default)   *
# *************************************************
process.source.skipEvents = cms.untracked.uint32(0)

process.demo = cms.EDAnalyzer('DimuonSpectrum2011MC')
# ***********************************************************
# output file name                                          *
# default is DoubleMuMC.root                                *
# change this according to your wish                        *
# ***********************************************************
process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string('DoubleMuMC.root'))

process.p = cms.Path(process.demo)
