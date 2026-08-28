import FWCore.ParameterSet.Config as cms

import collections

# Release the pixel digi and cluster SoA of the stub-seeded pixel CA once their only
# consumer (SiPixelRecHitAlpaka) has run.
#
# Only the DEVICE branches go into canDeleteEarly: the asynchronous backend's host copies
# share the host spelling, so listing them causes an UnusedProductsForCanDeleteEarly warning.
# The host copies are still registered as holding references to the device branches so
# the framework's indirect reader count covers the device branch.
#
# SAFETY: an early delete is safe only when consumers run on the producer's queue. That holds
# for SiPixelRecHitAlpaka (it calls get() before queue()). The stub merger, hit masker, SoA
# merger and CAHitNtuplet all take the queue before their first get -- their branches must
# NOT be added here (early deletion would make the CA read junk).

_CLUSTER_PRODUCER = "SiPixelPhase2DigiToCluster@alpaka"

# (device friendly class name, host-copy friendly class name) for each product of that module.
# The device spelling is the CUDA one; another asynchronous backend needs its own row.
_PRODUCTS = (
    ("alpakaDevCudaRtSiPixelDigisDeviceedmDeviceProduct", "SiPixelDigisHost"),
    ("alpakaDevCudaRtSiPixelClustersDeviceedmDeviceProduct", "SiPixelClustersHost"),
)


def customiseEarlyDeleteForPixelStubsCA(process, products):
    references = collections.defaultdict(list)

    def _branchName(productType, moduleLabel, instanceLabel=""):
        return "%s_%s_%s_%s" % (productType, moduleLabel, instanceLabel, process.name_())

    for name, module in process.producers_().items():
        if module.type_() != _CLUSTER_PRODUCER:
            continue
        for deviceType, hostType in _PRODUCTS:
            deviceBranch = _branchName(deviceType, name)
            hostBranch = _branchName(hostType, name)
            products[name].append(deviceBranch)
            # the host copy is produced by a transform that READS the device branch; its consumers
            # count as readers of the device branch through this reference
            references[hostBranch] = [deviceBranch]

    return (products, references)
