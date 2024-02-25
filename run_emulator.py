import os
import argparse
import yaml
import cppyy

import Toolkit as tool
import data_handle.plot_tools as plot

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

def run_algorithm(config, data_links, args, event):
    ''' Calling the emulator algorithm in all its steps '''

    linkUnpacking_ = l1thgcfirmware.HGCalLinkUnpacking(config)
    seeding_       = l1thgcfirmware.HGCalHistoSeeding(config)
    clustering_    = l1thgcfirmware.HGCalHistoClustering(config)

    unpackedTCs = l1thgcfirmware.HGCalTriggerCellSAPtrCollection()
    linkUnpacking_.runLinkUnpacking(data_links, unpackedTCs);
    if args.plot: plot.create_plot(unpackedTCs, 'post_unpacking', event, args.col)

    histogram = l1thgcfirmware.HGCalHistogramCellSAPtrCollection()
    seeding_.runSeeding(unpackedTCs, histogram)
    # if args.plot: tool.create_plot(histogram, 'post_seeding', event)

    # clusters = l1thgcfirmware.HGCalClusterSAPtrCollection()
    # clustering_.runClustering(unpackedTCs, histogram, clusters)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stage-2 Emulator Parameters')
    parser.add_argument('-n', type=int, default=1, help='Provide the number of events')
    parser.add_argument('--plot', action='store_true', help='Create plots for each event')
    parser.add_argument('--col',  action='store_true', help='Create plots using column numbers')
    args = parser.parse_args()

    tool.define_map(params)
    config = l1thgcfirmware.ClusterAlgoConfig(**params)
    # config.printConfiguration()

    ds = tool.provide_events(args.n)
    for ds_event in ds:
      ds_TCs, ds_gen = ds_event[0], ds_event[1]
      print('Processing event {}. (\u03B7, \u03C6) = {:.2f}, {:.2f}'.format(ds_gen.event,
            round(ds_gen.good_genpart_exeta[0], 2), round(ds_gen.good_genpart_exphi[0], 2)))

      data_links = LinksInData(tool.data_packer(ds_TCs, ds_gen, args))
      run_algorithm(config, data_links, args, ds_gen.event)

