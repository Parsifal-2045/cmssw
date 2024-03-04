import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *

from DPGAnalysis.MuonTools.nano_mu_local_reco_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi import *


muDPGNanoProducerHLT = cms.Sequence(muLocalRecoProducers
                                    + prunedGenParticles
                                    + finalGenParticles
                                    + genParticleTable)

def muDPGNanoCustomize(process) :

     if hasattr(process, "NANOAODSIMoutput"):
          process.prunedGenParticles.src = "genParticles"
          process.NANOAODSIMoutput.outputCommands.append("keep nanoaodFlatTable_*Table*_*_*")
          process.NANOAODSIMoutput.outputCommands.append("drop edmTriggerResults_*_*_*")

     return process
