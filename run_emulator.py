import os
import argparse
import cppyy
cppyy.add_include_path(os.path.abspath(''))
cppyy.load_library('lib_configuration.so')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h')
cppyy.include('L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h')

from cppyy.gbl import l1thgcfirmware
import data_handle.tools as tool
import data_handle.plot_tools as plot
from   data_handle.event import provide_events
import data_handle.geometry as geometry
import numpy as np

def run_algorithm(config, event, args, result):
    ''' Calling the emulator algorithm in all its steps '''

    linkUnpacking_ = l1thgcfirmware.HGCalLinkUnpacking(config)
    seeding_       = l1thgcfirmware.HGCalHistoSeeding(config)
    clustering_    = l1thgcfirmware.HGCalHistoClustering(config)

    unpackedTCs = l1thgcfirmware.HGCalTriggerCellSAPtrCollection()
    linkUnpacking_.runLinkUnpacking(event.data_packer, unpackedTCs);
    if args.plot: result.append(plot.create_plot(unpackedTCs, 'unpacking', event, args))

    histogram = l1thgcfirmware.HGCalHistogramCellSAPtrCollection()
    seeding_.runSeeding(unpackedTCs, histogram)
    if args.plot: result.append(plot.create_plot(histogram, 'seeding', event, args))

    clusters = l1thgcfirmware.HGCalClusterSAPtrCollection()
    clustering_.runClustering(unpackedTCs, histogram, clusters)
    if args.plot: result.append(plot.create_plot(histogram, 'clustering', event, args, clusters))
    
if __name__ == '__main__':
    ''' python run_emulator.py -n 2 --pileup PU0 --particles photons '''

    parser = argparse.ArgumentParser(description='Stage-2 Emulator Parameters')
    parser.add_argument('-n',          type=int, default=1,         help='Provide the number of events')
    parser.add_argument('--particles', type=str, default='photons', help='Choose the particle sample')
    parser.add_argument('--pileup',    type=str, default='PU0',     help='Choose the pileup - PU0 or PU200')
    parser.add_argument('--plot',        action='store_true', help='Create plots')
    parser.add_argument('--col',         action='store_true', help='Create plots using column numbers')
    parser.add_argument('--phi',         action='store_true', help='Create plots using phi coordinates')
    parser.add_argument('--performance', action='store_true', help='Create performance plots: distance gen_particle/max_TC')
    parser.add_argument('--thr_seed',    action='store_true', help='Create efficiency plots post seeding')
    parser.add_argument('--cl_energy',   action='store_true', help='Create plot of gen_pt vs recontructed energy')
    args = parser.parse_args()

    params = tool.define_map()
    config = l1thgcfirmware.ClusterAlgoConfig(**params)

    results = []
    events = provide_events(args.n, args.particles, args.pileup)
    xml_data, xml_MB = geometry.read_xml(), geometry.MB_geometry()
    for idx, event in enumerate(events):
      if idx % 50 == 0: print('Processing event', idx)
      if args.n <= 20: print('Processing event {}. (\u03B7, \u03C6) = {:.2f}, {:.2f}. pT = {:.2f} GeV'.format(
                              event.event, event.eta_gen, event.phi_gen, event.pT_gen))

      # if (event.pT_gen < 10): continue
      event._data_packer(args, xml_data, xml_MB)
      for thr_a in params['thresholdMaximaParam_a']:
        config.setThresholdMaximaConstants(params['cRows'], int(thr_a/event.LSB), 0, 0)
        run_algorithm(config, event, args, results)

    if args.performance: plot.produce_plots(results, args)
    if args.thr_seed:    plot.plot_seeds(results, args)
    if args.cl_energy:   plot.plot_cluster_energy(results, args)
