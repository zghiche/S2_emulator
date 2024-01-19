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
    return tree.arrays(branches, entry_start=1, entry_stop=2, library='ak')

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

def process_event():
    Nchannels = 84 * 3
    xml_data = get_channel_frame()
    module_conversion = prepare_geometry_txt()

    ds = reading_input_file()
    for event in ds.event:
        print('Processing event', event)
        ds_event = ds[ds.event == event]
        data_input = {}
        for tc_idx in range(len(ds_event.good_tc_x[0])):
            module = get_module_hash(module_conversion,
                                     ds_event.good_tc_layer[0][tc_idx],
                                     ds_event.good_tc_waferu[0][tc_idx],
                                     ds_event.good_tc_waferv[0][tc_idx])
            column = int((3*(ds_event.good_tc_phi[0][tc_idx]-1.04)*84)/(2*3.14))
            data = get_frame_channel(xml_data, module, column)
            if not data: continue
            
            frame, channel = int(data[0]['frame_id']), int(data[0]['channel_id'])
            data_input[Nchannels*frame+channel] = compress_value(int(10*ds_event.good_tc_energy[0][tc_idx]))
        break
    return data_input
