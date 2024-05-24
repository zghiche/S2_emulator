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
from data_handle.event import provide_events

def run_algorithm(config, event, args, shift):
    ''' Calling the emulator algorithm in all its steps '''

    linkUnpacking_ = l1thgcfirmware.HGCalLinkUnpacking(config)
    seeding_       = l1thgcfirmware.HGCalHistoSeeding(config)
    clustering_    = l1thgcfirmware.HGCalHistoClustering(config)

    unpackedTCs = l1thgcfirmware.HGCalTriggerCellSAPtrCollection()
    linkUnpacking_.runLinkUnpacking(event.data_packer, unpackedTCs);
    if args.plot: shift.append(plot.create_plot(unpackedTCs, 'post_unpacking', event, args))

    # histogram = l1thgcfirmware.HGCalHistogramCellSAPtrCollection()
    # seeding_.runSeeding(unpackedTCs, histogram)

    # clusters = l1thgcfirmware.HGCalClusterSAPtrCollection()
    # clustering_.runClustering(unpackedTCs, histogram, clusters)

if __name__ == '__main__':
    ''' python run_emulator.py --plot -n 2 '''

    parser = argparse.ArgumentParser(description='Stage-2 Emulator Parameters')
    parser.add_argument('-n', type=int, default=1, help='Provide the number of events')
    parser.add_argument('--plot',        action='store_true', help='Create plots for each event')
    parser.add_argument('--col',         action='store_true', help='Create plots using column numbers')
    parser.add_argument('--performance', action='store_true', help='Create plots calculating the distance gen_particle/max_TC')
    args = parser.parse_args()

    params = tool.define_map()
    config = l1thgcfirmware.ClusterAlgoConfig(**params)
    # config.printConfiguration()

    shift_pre, shift_post = [], []
    events = provide_events(args.n)
    with open('output.txt', 'a') as file:   
      #write the first comments for the EMP file
      num_columns=56
      file.write(f"ID: x1 \n")
      file.write(f"Metadata: (strobe,) start of orbit, start of packet, end of packet, valid \n \n")
      column_str = '          '.join(str(j).zfill(3).rjust(20) for j in range(num_columns))
      file.write(f"      Link {column_str}\n")
      file.close()
    nevents = None 
    for idx, event in enumerate(events):
      if idx % 50 == 0: print('Processing event', idx)
      print('Processing event {}. (\u03B7, \u03C6) = {:.2f}, {:.2f}. pT = {:.2f} GeV'.format(
            event.event, event.eta_gen, event.phi_gen, event.pT_gen))
      nevents=idx+1

      event._data_packer(args, shift_pre,idx,nevents)
      run_algorithm(config, event, args, shift_post)

    if args.performance: plot.produce_plots(shift_pre, shift_post)


    # with open('output.txt', 'a') as file:  
    #   #metadata = [1101 if i == 0 else 11 if i == 107 else 1 for i in range(num_rows)],below only for the last line of the EMP file
    #   metadata = 11 
    # # Iterate over each row and write the data for each link
    #   LinksInData_2d = [0] * 56       
    #   row_str = ' '.join(str(metadata).zfill(4) + " " + str(f'{LinksInData_2d[j]:016x}' ) + " " for j in range(56))
    #   # Write the last row to the file
    #   frame = 108*nevents
    #   file.write(f"Frame " + str(frame).zfill(4).rjust(4) + "    " + row_str + '\n')            
    #   file.close()