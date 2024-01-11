import os
import yaml
import Toolkit as tool
import cppyy

cppyy.add_include_path(os.path.abspath(''))

cppyy.load_library('lib_configuration.so')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h')

from cppyy.gbl import l1thgcfirmware, std

with open('config.yaml', "r") as afile:
    params = yaml.safe_load(afile)

def get_input_data():
    LinksInData = std.vector['std::unique_ptr<l1thgcfirmware::HGCalLinkTriggerCell>']()
    HGCalLinkTriggerCell = l1thgcfirmware.HGCalLinkTriggerCell
    
    for i in range(54432):
        LinksInData.push_back(std.make_unique[HGCalLinkTriggerCell]())
        if i == 264:  
            LinksInData[-1].data_ = 99
            print(LinksInData[-1].data_)
    
    return LinksInData

def run_algorithm(config):
    ''' Calling the emulator algorithm in all its steps '''

    linkUnpacking_ = l1thgcfirmware.HGCalLinkUnpacking(config)
    seeding_       = l1thgcfirmware.HGCalHistoSeeding(config)
    clustering_    = l1thgcfirmware.HGCalHistoClustering(config)

    unpackedTCs = l1thgcfirmware.HGCalTriggerCellSAPtrCollection()
    linkUnpacking_.runLinkUnpacking(get_input_data(), unpackedTCs);

    histogram = l1thgcfirmware.HGCalHistogramCellSAPtrCollection()
    seeding_.runSeeding(unpackedTCs, histogram)

    clusters = l1thgcfirmware.HGCalClusterSAPtrCollection()
    clustering_.runClustering(unpackedTCs, histogram, clusters)


if __name__ == '__main__':
    tool.process_event()
    exit()

    tool.define_map(params)
    config = l1thgcfirmware.ClusterAlgoConfig(**params)
    config.printConfiguration()

    run_algorithm(config)
