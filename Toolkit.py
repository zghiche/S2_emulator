import cppyy
from cppyy.gbl import std
import math
import numpy as np
import awkward as ak

from data_handle.geometry import provide_events, read_geometry_txt, read_xml

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

def define_map(params):
    enum = cppyy.gbl.l1thgcfirmware.Step
    map_custom = std.map[enum, "unsigned int"]()

    for d in params['stepLatency']:
      for latency in d.keys():
        map_custom[int(vars(enum)[latency])] = d[latency]

    params['stepLatency'] = map_custom

def get_module_hash(conversion, plane, u, v):
    filtered_rows = conversion[(conversion[:, 0] == plane) & 
                               (conversion[:, 1] == u) &
                               (conversion[:, 2] == v)]
    return int(filtered_rows[0][3])

def get_TC_allocation(xml_data, module):
    return xml_data[module]

def create_link(data_in):
    energy   = (data_in[2][0] << 30) | (data_in[1][0] << 15) | (data_in[0][0])
    r_over_z = (data_in[2][1] << 30) | (data_in[1][1] << 15) | (data_in[0][1])
    phi      = (data_in[2][2] << 30) | (data_in[1][2] << 15) | (data_in[0][2])
    return [energy, r_over_z, phi]

def process_event(ds_event, geometry, xml_data):
    LSB = 1/100 # 10 MeV
    LSB_phi = np.pi/1944
    LSB_r_z = 0.7/4096
    data_TCs = {}
    
    for module_idx in range(len(ds_event.good_tc_x)):
        module = get_module_hash(geometry,
                                 ds_event.good_tc_layer[module_idx][0],
                                 ds_event.good_tc_waferu[module_idx][0],
                                 ds_event.good_tc_waferv[module_idx][0])
        xml_alloc = get_TC_allocation(xml_data, module)
        
        # calculating the number of TC that ca be allocated / module
        # apply the cut to simulate the BC algorithm
        n_TCs = xml_alloc[-1]['index']  # dangerous

        # simulating the phi sorting by the S1 FPGA
        mod_phi = ds_event.good_tc_phi[module_idx][:n_TCs]
        mod_energy = ds_event.good_tc_mipPt[module_idx][:n_TCs][ak.argsort(mod_phi)]
        mod_r_over_z = ds_event.r_over_z[module_idx][:n_TCs][ak.argsort(mod_phi)]
        mod_phi = ak.sort(mod_phi)       

        print('Analysing module', module)
        for tc_idx, TC_xml in enumerate(xml_alloc):
            if tc_idx > len(mod_energy)-1: break
            n_link = math.floor(TC_xml['glb_channel']/3)

            if TC_xml['frame'] not in data_TCs.keys():
                data_TCs[TC_xml['frame']] = {}
            if n_link not in data_TCs[TC_xml['frame']].keys():
                data_TCs[TC_xml['frame']][n_link] = [[0]*3]*3

            data_TCs[TC_xml['frame']][n_link][TC_xml['channel']%3] = [
                compress_value(int(mod_energy[tc_idx]/LSB)),
                compress_value(int(mod_r_over_z[tc_idx]/LSB_r_z)),
                compress_value(int(mod_phi[tc_idx]/LSB_phi))
                ]
    return data_TCs
    
def data_packer():
    xml_data = read_xml()
    geometry = read_geometry_txt()

    ds = provide_events()
    print('Processing event..')
    data_TCs = process_event(ds[0], geometry, xml_data)

    # packing data in links 
    data_links = {}
    for frame in data_TCs.keys():
        for n_link in data_TCs[frame].keys():
            link_data = create_link(data_TCs[frame][n_link])
            data_links[84*frame+n_link] = link_data
    
    return data_links
