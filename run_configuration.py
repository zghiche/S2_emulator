import os
import yaml
from Toolkit import *

import cppyy
cppyy.add_include_path(os.path.abspath(''))

cppyy.load_library('lib_configuration.so')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h')

emulator = cppyy.gbl.l1thgcfirmware

def get_input_data():
    LinksInData = cppyy.gbl.std.vector['std::unique_ptr<l1thgcfirmware::HGCalLinkTriggerCell>']()
    HGCalLinkTriggerCell = emulator.HGCalLinkTriggerCell
    
    for i in range(54432):
        LinksInData.push_back(cppyy.gbl.std.make_unique[HGCalLinkTriggerCell]())
        if i == 264:  
            LinksInData[-1].data_ = 99
            print(LinksInData[-1].data_)
    
    return LinksInData

def run_algorithm(config):
    ''' Calling the emulator algorithm in all its steps '''

    linkUnpacking_ = emulator.HGCalLinkUnpacking(config)
    seeding_       = emulator.HGCalHistoSeeding(config)
    clustering_    = emulator.HGCalHistoClustering(config)

    unpackedTCs = emulator.HGCalTriggerCellSAPtrCollection()
    linkUnpacking_.runLinkUnpacking(get_input_data(), unpackedTCs);

    histogram = emulator.HGCalHistogramCellSAPtrCollection()
    seeding_.runSeeding(unpackedTCs, histogram)

    clusters = emulator.HGCalClusterSAPtrCollection()
    clustering_.runClustering(unpackedTCs, histogram, clusters)


if __name__ == '__main__':

    define_map(params)

    config = emulator.ClusterAlgoConfig(**params)
    config.printConfiguration()

    run_algorithm(config)
