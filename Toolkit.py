import cppyy
from cppyy.gbl import std
import numpy as np
import awkward as ak
import uproot
import math
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
                    glob_channel = 12*int(s1_id[1:])+int(channel_id)
                    reversed_data[module][column].append({'frame_id': frame_id, 'channel_id': channel_id, 'glob_channel': glob_channel, 's1_id': s1_id})

    return reversed_data

def define_map(params):
    enum = cppyy.gbl.l1thgcfirmware.Step
    map_custom = std.map[enum, "unsigned int"]()

    for d in params['stepLatency']:
      for latency in d.keys():
        map_custom[int(vars(enum)[latency])] = d[latency]

    params['stepLatency'] = map_custom

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
    event, phi_gen = 0, -1
    while not 0 < phi_gen < 2.02:
        event += 1
        phi_gen = tree.arrays(['event','good_genpart_exphi'], entry_start=event, 
                             entry_stop=event+1, library='ak').good_genpart_exphi

    data = tree.arrays(branches, entry_start=event, entry_stop=event+1, library='ak')
    sorted_indices = ak.argsort(data['good_tc_energy'], ascending=False)
    for key in data.fields:
        if key != 'event':
            data[key] = data[key][sorted_indices]
    return data

def get_module_hash(conversion, plane, u, v):
    filtered_rows = conversion[(conversion[:, 0] == plane) & 
                               (conversion[:, 1] == u) &
                               (conversion[:, 2] == v)]
    return int(filtered_rows[0][3])

def get_frame_channel(xml_data, module, column):
    return xml_data[str(module)][str(column)]

def compress_value(value, exponent_bits=4, mantissa_bits=3, truncation_bits=0):
    saturation_code = (1 << (exponent_bits + mantissa_bits)) - 1
    saturation_value = ((1 << (mantissa_bits + truncation_bits + 1)) - 1) << ((1 << exponent_bits) - 2)

    if value > saturation_value:
        return saturation_code

    bitlen = 0
    shifted_value = value >> truncation_bits
    valcopy = shifted_value
    while valcopy != 0:
        valcopy >>= 1
        bitlen += 1

    if bitlen <= mantissa_bits:
        compressed_code = shifted_value
        compressed_value = shifted_value << truncation_bits
        return compressed_code

    # Build exponent and mantissa
    exponent = bitlen - mantissa_bits
    mantissa = (shifted_value >> (exponent - 1)) & ~(1 << mantissa_bits)
   
    # Assemble floating-point
    compressed_value = ((1 << mantissa_bits) | mantissa) << (exponent - 1)
    compressed_code = (mantissa << exponent_bits) | exponent #(exponent << mantissa_bits) | mantissa
    return compressed_code

def create_link(data_in):
    return (data_in[2] << 30) | (data_in[1] << 15) | (data_in[0])

def process_event():
    Nchannels = 84 * 3
    LSB = 1/100 # 10 MeV
    xml_data = get_channel_frame()
    module_conversion = prepare_geometry_txt()

    data_input = {}
    ds = reading_input_file()
    for event in ds.event:
        print('Processing event', event)
        ds_event = ds[ds.event == event]
        data_in_frames = {}
        data_index = {}
        for tc_idx in range(len(ds_event.good_tc_x[0])):
            module = get_module_hash(module_conversion,
                                     ds_event.good_tc_layer[0][tc_idx],
                                     ds_event.good_tc_waferu[0][tc_idx],
                                     ds_event.good_tc_waferv[0][tc_idx])
            column = int((3*(ds_event.good_tc_phi[0][tc_idx])*84)/(2*3.14))
            data = get_frame_channel(xml_data, module, column)
            if not data: continue
            
            if not (module, column) in data_index.keys():
                data_index[(module, column)] = len(data)-1
            elif data_index[(module, column)] != 0: 
                 data_index[(module, column)] -= 1
            else: continue

            # there is room for this TC
            frame, channel, glob_channel = (int(data[data_index[(module, column)]]['frame_id']), 
                                            int(data[data_index[(module, column)]]['channel_id']), 
                                            int(data[data_index[(module, column)]]['glob_channel']))

            # print("TC having energy ", ds_event.good_tc_energy[0][tc_idx], "has frame and channel ", 
            #       frame, glob_channel, "module, column", module, column)
            
            # packing of the data
            n_link = math.floor(glob_channel/3)
            if frame not in data_in_frames.keys():
                data_in_frames[frame] = {}
            if n_link not in data_in_frames[frame].keys():
                data_in_frames[frame][n_link] = [0]*3
            
            data_in_frames[frame][n_link][channel%3] = compress_value(int(ds_event.good_tc_energy[0][tc_idx]/LSB))
        break # only one event
  
    for frame in data_in_frames.keys():
        for n_link in data_in_frames[frame].keys():
            link_data = create_link(data_in_frames[frame][n_link])
            data_input[84*frame+n_link] = link_data
    
    return data_input
