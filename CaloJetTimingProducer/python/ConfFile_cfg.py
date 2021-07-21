import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/m/mcitron/triggerOutput.root'
    )
)

process.demo = cms.EDProducer('CaloJetTimingProducer'
)

process.caloTimingTestFilter = cms.EDFilter('CaloJetTimingFilter',
    saveTags = cms.bool( True ),
)
process.output = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "/afs/cern.ch/work/m/mcitron/tempOutput.root" ),
    )
process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')

process.p = cms.Path(process.demo)
process.Output = cms.EndPath(process.output)
