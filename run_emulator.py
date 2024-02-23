import os
import argparse
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

def LinksInData(data_packed):
    LinksInData = std.vector['std::unique_ptr<l1thgcfirmware::HGCalLinkTriggerCell>']()
    HGCalLinkTriggerCell = l1thgcfirmware.HGCalLinkTriggerCell
    
    # number of links = NChannels*Nframes
    for link_frame in range(84*108):
        LinksInData.push_back(std.make_unique[HGCalLinkTriggerCell]())
        if link_frame in data_packed.keys():
            LinksInData[-1].data_     = data_packed[link_frame][0]
            LinksInData[-1].r_over_z_ = data_packed[link_frame][1]
            LinksInData[-1].phi_      = data_packed[link_frame][2]
    
    return LinksInData

def create_plot(objects, step):
    heatmap_data = np.zeros((64, 124)) 
   
    for object in objects:
      if (step==0) and (object.energy()>0): 
         heatmap_data[int((object.rOverZ()-440)/64)-1, int(124/1944*object.phi())-1] += (object.energy())/100
         print("Energy : ", object.energy(), "Phi", object.phi(), "R/Z", int((object.rOverZ()-440)/64), "Column", object.index())
      
      if (step==1) and (object.S()>0):
         heatmap_data[object.sortKey()-1, object.index()-1] += (object.S())/100
         print("Smeared Energy : ", object.S(), "R/Z bin", object.sortKey(), "col", object.index())
    
    title = 'Unpacked Energy' if step==0 else 'Smeared Energy'
    plt.imshow(heatmap_data, cmap='viridis', aspect='auto')
    plt.xlabel('Column')
    plt.ylabel('R/Z Bin')
    plt.title('Heatmap ' + title)
    plt.colorbar(label = title)
    plt.savefig('plots/'+title.replace(' ', '_')+'.pdf')
    plt.clf()


def run_algorithm(config, data_links):
    ''' Calling the emulator algorithm in all its steps '''

    linkUnpacking_ = l1thgcfirmware.HGCalLinkUnpacking(config)
    seeding_       = l1thgcfirmware.HGCalHistoSeeding(config)
    clustering_    = l1thgcfirmware.HGCalHistoClustering(config)

    unpackedTCs = l1thgcfirmware.HGCalTriggerCellSAPtrCollection()
    linkUnpacking_.runLinkUnpacking(data_links, unpackedTCs);
    # create_plot(unpackedTCs, 0)

    histogram = l1thgcfirmware.HGCalHistogramCellSAPtrCollection()
    seeding_.runSeeding(unpackedTCs, histogram)
    # create_plot(histogram, 1)

    # clusters = l1thgcfirmware.HGCalClusterSAPtrCollection()
    # clustering_.runClustering(unpackedTCs, histogram, clusters)

    # print(clusters)
    # for cluster in clusters:
    # print('ciao')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stage-2 Emulator Parameters')
    parser.add_argument('-n', type=int, help='Provide the number of events')
    args = parser.parse_args()

    tool.define_map(params)
    config = l1thgcfirmware.ClusterAlgoConfig(**params)
    # config.printConfiguration()

    ds = tool.provide_events(args.n)
    for event in ds:
      data_links = LinksInData(tool.data_packer(event))
      run_algorithm(config, data_links)
