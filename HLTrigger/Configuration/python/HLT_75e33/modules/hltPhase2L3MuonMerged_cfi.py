import FWCore.ParameterSet.Config as cms

from HLTrigger.Configuration.HLT_75e33.modules.hltL2MuonSeedsFromL1TkMuon_cfi import PHASE2_TAG

from HLTrigger.Configuration.HLT_75e33.modules.Phase2HLTMuonSelectorForL3_cfi import L3IOFIRST

if PHASE2_TAG:
    TracksToMerge = cms.VInputTag()
    if L3IOFIRST:
        TracksToMerge = cms.VInputTag("hltPhase2L3OIMuonTrackSelectionHighPurity", cms.InputTag("phase2L3FilteredObjects","L3IOTracksFiltered"))
    else:
        TracksToMerge = cms.VInputTag(cms.InputTag("phase2L3FilteredObjects","L3OITracksFiltered"), "hltIter2Phase2L3FromL1TkMuonMerged")

    hltPhase2L3MuonMerged = cms.EDProducer("TrackListMerger",
        Epsilon = cms.double(-0.001),
        FoundHitBonus = cms.double(5.0),
        LostHitPenalty = cms.double(20.0),
        MaxNormalizedChisq = cms.double(1000.0),
        MinFound = cms.int32(3),
        MinPT = cms.double(0.05),
        ShareFrac = cms.double(0.19),
        TrackProducers = TracksToMerge,
        allowFirstHitShare = cms.bool(True),
        copyExtras = cms.untracked.bool(True),
        copyMVA = cms.bool(False),
        hasSelector = cms.vint32(0, 0),
        indivShareFrac = cms.vdouble(1.0, 1.0),
        newQuality = cms.string('confirmed'),
        selectedTrackQuals = TracksToMerge,
        setsToMerge = cms.VPSet(cms.PSet(
            pQual = cms.bool(False),
            tLists = cms.vint32(0, 1)
        )),
        trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
        writeOnlyTrkQuals = cms.bool(False)
    )

else:
    hltPhase2L3MuonMerged = cms.EDProducer("TrackListMerger",
        Epsilon = cms.double(-0.001),
        FoundHitBonus = cms.double(5.0),
        LostHitPenalty = cms.double(20.0),
        MaxNormalizedChisq = cms.double(1000.0),
        MinFound = cms.int32(3),
        MinPT = cms.double(0.05),
        ShareFrac = cms.double(0.19),
        TrackProducers = cms.VInputTag("hltPhase2L3OIMuonTrackSelectionHighPurity", "hltIter2Phase2L3FromL1TkMuonMerged"),
        allowFirstHitShare = cms.bool(True),
        copyExtras = cms.untracked.bool(True),
        copyMVA = cms.bool(False),
        hasSelector = cms.vint32(0, 0),
        indivShareFrac = cms.vdouble(1.0, 1.0),
        newQuality = cms.string('confirmed'),
        selectedTrackQuals = cms.VInputTag("hltPhase2L3OIMuonTrackSelectionHighPurity", "hltIter2Phase2L3FromL1TkMuonMerged"),
        setsToMerge = cms.VPSet(cms.PSet(
            pQual = cms.bool(False),
            tLists = cms.vint32(0, 1)
        )),
        trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
        writeOnlyTrkQuals = cms.bool(False)
    )
