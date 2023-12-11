import os
import cppyy
cppyy.add_include_path(os.path.abspath(''))

cppyy.load_library('lib_configuration.so')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h')

import yaml
with open('config.yaml', "r") as afile:
    params = yaml.safe_load(afile)

def define_map(params):
    enum = cppyy.gbl.l1thgcfirmware.Step
    map_custom = cppyy.gbl.std.map[enum, "unsigned int"]()
    
    for d in params['stepLatency']:
      for latency in d.keys():
        map_custom[int(vars(enum)[latency])] = d[latency]

    params['stepLatency'] = map_custom

def get_input_data():
    LinksInData = cppyy.gbl.std.vector['std::unique_ptr<l1thgcfirmware::HGCalLinkTriggerCell>']()
    HGCalLinkTriggerCell = cppyy.gbl.l1thgcfirmware.HGCalLinkTriggerCell
    
    for i in range(54432):
        LinksInData.push_back(cppyy.gbl.std.make_unique[HGCalLinkTriggerCell]())
        if i == 264:  
            LinksInData[-1].data_ = 99
            print(LinksInData[-1].data_)
    
    return LinksInData

if __name__ == '__main__':

    define_map(params)
    data = get_input_data()

    print('setting parameters from yaml file')
    ClusterAlgoConfig = cppyy.gbl.l1thgcfirmware.ClusterAlgoConfig
    config = ClusterAlgoConfig(**params)
    # config.printConfiguration()

    print('calling the algorithm')
    HGCalHistoClusteringImplSA = cppyy.gbl.l1thgcfirmware.HGCalHistoClusteringImplSA
    theAlgo = HGCalHistoClusteringImplSA(config)
    theAlgo.runAlgorithm(data)
