import cppyy
from cppyy.gbl import std
import math

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

def get_frame_channel(xml_data, module, column):
    return xml_data[module][column]

def create_link(data_in):
    return (data_in[2] << 30) | (data_in[1] << 15) | (data_in[0])

def process_event(ds_event, geometry, xml_data):
    LSB = 1/100 # 10 MeV
    data_TCs = {}
    TC_index = {}
    for tc_idx in range(len(ds_event.good_tc_x[0])):
        module = get_module_hash(geometry,
                                 ds_event.good_tc_layer[0][tc_idx],
                                 ds_event.good_tc_waferu[0][tc_idx],
                                 ds_event.good_tc_waferv[0][tc_idx])
        column = int((3*(ds_event.good_tc_phi[0][tc_idx])*84)/(2*3.14))
    
        data = get_frame_channel(xml_data, module, column)
        if not data: continue
        
        if not (module, column) in TC_index.keys():
             TC_index[(module,column)] = len(data)-1
        elif TC_index[(module,column)] != 0: 
             TC_index[(module,column)] -= 1
        else: continue

        # allocating space for the current TC
        allocation = data[TC_index[(module,column)]]
        frame, channel, glb_channel = allocation['frame'], allocation['channel'], allocation['glb_channel']
 
        # print("TC having energy ", ds_event.good_tc_energy[0][tc_idx], "has frame and channel ", 
        #      frame, glb_channel, "module, column", module, column)
       
        # TC data packing
        n_link = math.floor(glb_channel/3)
        if frame not in data_TCs.keys():
            data_TCs[frame] = {}
        if n_link not in data_TCs[frame].keys():
            data_TCs[frame][n_link] = [0]*3
        
        data_TCs[frame][n_link][channel%3] = compress_value(int(ds_event.good_tc_energy[0][tc_idx]/LSB))
    return data_TCs

def data_packer():
    xml_data = read_xml()
    geometry = read_geometry_txt()

    ds = provide_events()
    for event in ds.event:
        print('Processing event', event)
        ds_event = ds[ds.event == event]
        data_TCs = process_event(ds_event, geometry, xml_data)
        break # only one event
  
    # packing data in links 
    data_links = {}
    for frame in data_TCs.keys():
        for n_link in data_TCs[frame].keys():
            link_data = create_link(data_TCs[frame][n_link])
            data_links[84*frame+n_link] = link_data
    
    return data_links
