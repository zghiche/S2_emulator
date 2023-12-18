import cppyy
from cppyy.gbl import std
import numpy as np
import awkward as ak
import uproot

def define_map(params):
    enum = cppyy.gbl.l1thgcfirmware.Step
    map_custom = std.map[enum, "unsigned int"]()

    for d in params['stepLatency']:
      for latency in d.keys():
        map_custom[int(vars(enum)[latency])] = d[latency]

    params['stepLatency'] = map_custom

def unpack_values(lut_out):
    print(lut_out)
    R_over_Z = (lut_out >> 0) & 0xFFF
    Phi = (lut_out >> 12) & 0xFFF
    print(Phi)
    Layer = (lut_out >> 24) & 0x3F
    index = (lut_out >> 30) & 0x1FF
    exit()
    return R_over_Z, Phi, Layer, index

def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    values_list = []
    for line in lines:
        lut_out = hex(line.strip())
        valid = ( lut_out >> 39 ) & 0x1;

        if valid: 
            unpacked_values = unpack_values(lut_out)
            values_list.append(unpacked_values)

    values = np.array(values_list)
    return values

def r_over_z(entry):
    x = entry['good_tc_x'][0]
    y = entry['good_tc_y'][0]
    z = entry['good_tc_z'][0]
    r = np.sqrt(np.sum(x**2 + y**2))
   
    return np.array(r/z)

def merge_arrays(event, mif_data):
    for tc in event:
        for row_index, coords in enumerate(mif_data):
            print(coords[:3], tc[:3])
            if np.all(coords[:3] == tc[:3]): 
                print(mif_data[row_index])
                continue


def process_event():
    mif_data = process_file('S2.CombinedTD.Balanced60.MixedTypes.NoSplit.mif')
    
    ds = reading_input_file()
    for event in ds.event:
        print('Processing event', event)
        R_over_Z = r_over_z(ds[ds.event == event]).astype(int)
        Phi      = np.array(ds[ds.event == event]['good_tc_phi'])[0].astype(int)
        Plane    = np.array(ds[ds.event == event]['good_tc_layer'])[0]
        energy   = np.array(ds[ds.event == event]['good_tc_energy'])[0]

        event = np.stack((R_over_Z, Phi, Plane, energy), axis=-1)    
        # LiksInData = merge_arrays(event, mif_data)

def reading_input_file():
    filepath = '../../building_ROI/skim_small_photons_0PU_bc_stc_hadd.root'
    name_tree = "FloatingpointMixedbcstcrealsig4DummyHistomaxxydr015GenmatchGenclustersntuple/HGCalTriggerNtuple"

    branches = ['event',
    'good_tc_x',
    'good_tc_y',
    'good_tc_z',
    'good_tc_phi',
    'good_tc_layer',
    'good_tc_cellu',
    'good_tc_cellv',
    'good_tc_energy',
    'good_tc_mipPt',
    'good_tc_cluster_id']

    tree  = uproot.open(filepath)[name_tree]
    return tree.arrays(branches, entry_stop=1, library='ak')

