import os
import numpy as np
import matplotlib.pyplot as plt
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
    
    data_input = tool.data_packer()
    for link_frame in range(84*216):
        LinksInData.push_back(std.make_unique[HGCalLinkTriggerCell]())
        if link_frame in data_input.keys():
            LinksInData[-1].data_     = data_input[link_frame][0]
            LinksInData[-1].r_over_z_ = data_input[link_frame][1]
            LinksInData[-1].phi_      = data_input[link_frame][2]
    
    return LinksInData

def create_plot(objects, step):
    heatmap_data = np.zeros((64, 124)) 
    
    for object in objects:
      if (step==0) and (object.energy()>0): heatmap_data[int((object.rOverZ()-440)/64)-1, object.index()-1] += (object.energy())/100
      # print("Energy : ", tc.energy(), "Phi", tc.phi(), "R/Z", tc.rOverZ(), "Column", tc.index())
      
      if (step==1) and (object.X()>0):      heatmap_data[object.sortKey()-1, object.index()-1] += (object.S())/100
      # print("Smeared Energy : ", bin.S(), "R/Z bin", bin.sortKey(), "col", bin.index())
    
    title = 'Unpacked Energy' if step==0 else 'Smeared Energy'
    plt.imshow(heatmap_data, cmap='viridis', aspect='auto')
    plt.xlabel('Column')
    plt.ylabel('R/Z Bin')
    plt.title('Heatmap ' + title)
    plt.colorbar(label = title)
    plt.savefig('plots/'+title.replace(' ', '_')+'.pdf')
    plt.clf()


def run_algorithm(config):
    ''' Calling the emulator algorithm in all its steps '''

    linkUnpacking_ = l1thgcfirmware.HGCalLinkUnpacking(config)
    seeding_       = l1thgcfirmware.HGCalHistoSeeding(config)
    clustering_    = l1thgcfirmware.HGCalHistoClustering(config)

    unpackedTCs = l1thgcfirmware.HGCalTriggerCellSAPtrCollection()
    linkUnpacking_.runLinkUnpacking(get_input_data(), unpackedTCs);
    # create_plot(unpackedTCs, 0)

    histogram = l1thgcfirmware.HGCalHistogramCellSAPtrCollection()
    seeding_.runSeeding(unpackedTCs, histogram)
    # create_plot(histogram, 1)

    clusters = l1thgcfirmware.HGCalClusterSAPtrCollection()
    clustering_.runClustering(unpackedTCs, histogram, clusters)

    #for cluster in clusters:
    #  print('ciao')

if __name__ == '__main__':
    tool.define_map(params)
    config = l1thgcfirmware.ClusterAlgoConfig(**params)
    # config.printConfiguration()

    run_algorithm(config)
