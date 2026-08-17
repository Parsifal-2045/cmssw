import FWCore.ParameterSet.Config as cms

from Configuration.ProcessModifiers.phase2CAStubs_cff import phase2CAStubs
from Configuration.ProcessModifiers.pixelTrackMask_cff import pixelTrackMask

# Modifier carrying the ngt scouting menu's own customisations; use for toModify/toReplaceWith.
ngtScoutingBase = cms.Modifier()

# NGT scouting menu: runs pixel CA with OT stubs over both iterations (prompt+displaced),
# uses the resulting pixel tracks as general tracks.
ngtScouting = cms.ModifierChain(phase2CAStubs, pixelTrackMask, ngtScoutingBase)

# Single-iteration variant: stub-seeded pixel CA runs prompt iteration only (no hit masking).
ngtScoutingSingleIter = cms.ModifierChain(phase2CAStubs, ngtScoutingBase)
