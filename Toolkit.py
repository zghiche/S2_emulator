import cppyy
from cppyy.gbl import std
import numpy as np
import awkward as ak
import uproot
import xml.etree.ElementTree as ET
from xml.dom import minidom
from collections import defaultdict

def get_channel_frame():
    tree = ET.parse('../../mapping/FE2BE_mapping/checks_FE_BE/xml/S1toChannels.SeparateTD.120.SingleTypes.NoSplit.xml')
    root = tree.getroot()

    reversed_data = defaultdict(lambda: defaultdict(list))

    for s1_element in root.findall('.//S1'):
        s1_id = s1_element.get('id')
        for channel_element in s1_element.findall('.//Channel'):
            channel_id = channel_element.get('id')
            for frame_element in channel_element.findall('.//Frame'):
                if all(attr in frame_element.attrib for attr in ['id', 'column', 'Module']):
                    frame_id = frame_element.get('id')
                    column = frame_element.get('column')
                    module = frame_element.get('Module')
                    reversed_data[module][column].append({'frame_id': frame_id, 'channel_id': channel_id, 's1_id': s1_id})

    return reversed_data

def define_map(params):
    enum = cppyy.gbl.l1thgcfirmware.Step
    map_custom = std.map[enum, "unsigned int"]()

    for d in params['stepLatency']:
      for latency in d.keys():
        map_custom[int(vars(enum)[latency])] = d[latency]

    params['stepLatency'] = map_custom

def unpack_values(lut_out):
    R_over_Z = (lut_out >> 0) & 0xFFF
    Phi = (lut_out >> 12) & 0xFFF
    Layer = (lut_out >> 24) & 0x3F
    index = (lut_out >> 30) & 0x1FF
    return R_over_Z, Phi, Layer, index

def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    values_list = []
    for line in lines:
        lut_out = int(line.strip(), 16)#, "040b")
        print(lut_out)
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

def prepare_geometry_txt():
    ''' old but working version, it reads the geometry 
        file from Pedro's txt '''
    geometry_file = '../../mapping/FE_mapping/geometries/v15.3/geometry.hgcal.txt'
    
    # 'plane','u','v','hash'
    return np.loadtxt(geometry_file, delimiter=' ', usecols=(0,1,2,37), skiprows=1)

def reading_input_file():
    filepath = '../../building_ROI/skim_small_photons_0PU_bc_stc_hadd.root'
    name_tree = "FloatingpointMixedbcstcrealsig4DummyHistomaxxydr015GenmatchGenclustersntuple/HGCalTriggerNtuple"

    branches = ['event',
    'good_tc_x',
    'good_tc_y',
    'good_tc_z',
    'good_tc_phi',
    'good_tc_layer',
    'good_tc_waferu',
    'good_tc_waferv',
    'good_tc_energy',
    'good_tc_mipPt',
    'good_tc_cluster_id']

    tree  = uproot.open(filepath)[name_tree]
    return tree.arrays(branches, entry_start=19, entry_stop=20, library='ak')

def get_module_hash(conversion, plane, u, v):
    filtered_rows = conversion[(conversion[:, 0] == plane) & 
                               (conversion[:, 1] == u) &
                               (conversion[:, 2] == v)]
    return filtered_rows[:, 3]

def process_event():
    xml_data = get_channel_frame()
    module_conversion = prepare_geometry_txt()

    ds = reading_input_file()
    for event in ds.event:
        print('Processing event', event)
        ds_event = ds[ds.event == event]
        for tc_idx in range(len(ds_event.good_tc_x[0])):
            module = get_module_hash(module_conversion,
                                     ds_event.good_tc_layer[0][tc_idx],
                                     ds_event.good_tc_waferu[0][tc_idx],
                                     ds_event.good_tc_waferv[0][tc_idx])
            column = int((3*ds_event.good_tc_phi[0][tc_idx]*84)/(2*3.14))
            print(module, column)
            exit()
            print(xml_data[module][column])

    # mif_data = process_file('S2.CombinedTD.Balanced60.MixedTypes.NoSplit.mif')
    # 
    # print(mif_data)
    # ds = reading_input_file()
    # for event in ds.event:
    #     print('Processing event', event)
    #     R_over_Z = r_over_z(ds[ds.event == event]).astype(int)
    #     Phi      = np.array(ds[ds.event == event]['good_tc_phi'])[0].astype(int)
    #     Plane    = np.array(ds[ds.event == event]['good_tc_layer'])[0]
    #     energy   = np.array(ds[ds.event == event]['good_tc_energy'])[0]

    #     event = np.stack((R_over_Z, Phi, Plane, energy), axis=-1)    
    #     # LiksInData = merge_arrays(event, mif_data)

