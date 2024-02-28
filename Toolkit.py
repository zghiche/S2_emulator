import cppyy
from cppyy.gbl import std
import math
import numpy as np
import awkward as ak

from data_handle.geometry import provide_events, read_xml
import data_handle.plot_tools as plot 

def compress_value(value, exponent_bits=4, mantissa_bits=3, truncation_bits=0):
    saturation_code = (1 << (exponent_bits + mantissa_bits)) - 1
    saturation_value = ((1 << (mantissa_bits + truncation_bits + 1)) - 1) << ((1 << exponent_bits) - 2)

    if value > saturation_value:
        return saturation_value, saturation_code

    bitlen = 0
    shifted_value = int(value) >> truncation_bits
    valcopy = shifted_value
    while valcopy != 0:
        valcopy >>= 1
        bitlen += 1

    if bitlen <= mantissa_bits:
        compressed_code = shifted_value
        compressed_value = shifted_value << truncation_bits
        return compressed_value, compressed_code

    # Build exponent and mantissa
    exponent = bitlen - mantissa_bits
    mantissa = (shifted_value >> (exponent - 1)) & ~(1 << mantissa_bits)
   
    # Assemble floating-point
    compressed_value = ((1 << mantissa_bits) | mantissa) << (exponent - 1)
    compressed_code = (mantissa << exponent_bits) | exponent #(exponent << mantissa_bits) | mantissa
    return compressed_value, compressed_code

def define_map(params):
    enum = cppyy.gbl.l1thgcfirmware.Step
    map_custom = std.map[enum, "unsigned int"]()

    for d in params['stepLatency']:
      for latency in d.keys():
        map_custom[int(vars(enum)[latency])] = d[latency]

    params['stepLatency'] = map_custom

def ObjectType( object_type ):
    return ( ( object_type & 0xF ) << 22 )

def get_module_id(plane, u, v):
    # CMSSW to our u v convention u'=v-u, v'=v
    print('Analysing module ', plane, v-u, v)
    if plane & ~0x3F : raise Exception( "Invalid plane" )
    if v-u & ~0xF : raise Exception( "Invalid u" )
    if v   & ~0xF : raise Exception( "Invalid v" )
    return hex(0x60000000 | ObjectType(0) | ((plane & 0x3F) << 16) | ((v-u & 0xF) << 12) | ((v & 0xF) << 8))
 
def get_TC_allocation(xml_data, module):
    return xml_data[module]

def create_link(data_in):
    energy   = (data_in[2][0] << 30) | (data_in[1][0] << 15) | (data_in[0][0])
    r_over_z = (data_in[2][1] << 30) | (data_in[1][1] << 15) | (data_in[0][1])
    phi      = (data_in[2][2] << 30) | (data_in[1][2] << 15) | (data_in[0][2])
    return [energy, r_over_z, phi]

def process_event(ds_TCs, ds_gen, args, shift):
    LSB = 1/100000 # 10 keV
    LSB_r_z = 0.7/4096
    LSB_phi = np.pi/1944
    offset_phi = -0.3
    data_TCs = {}
    xml_data = read_xml()
   
    if args.plot: data_heatmap = []
    for module_idx in range(len(ds_TCs.good_tc_x)):
        module = get_module_id(ds_TCs.good_tc_layer[module_idx][0],
                               ds_TCs.good_tc_waferu[module_idx][0],
                               ds_TCs.good_tc_waferv[module_idx][0])
        if not module: continue
        xml_alloc = get_TC_allocation(xml_data, module)

        # calculating the number of TC that ca be allocated / module
        n_TCs = xml_alloc[-1]['index']  # dangerous
        columns = [frame['column'] for frame in xml_alloc]

        # simulating the BC algorithm (ECON-T) and the phi sorting in the S1 FPGA
        mod_phi = ds_TCs.good_tc_phi[module_idx][:n_TCs+1]
        mod_energy = ds_TCs.good_tc_mipPt[module_idx][:n_TCs+1][ak.argsort(mod_phi)]
        mod_r_over_z = ds_TCs.r_over_z[module_idx][:n_TCs+1][ak.argsort(mod_phi)]
        mod_phi = ak.sort(mod_phi)       

        for tc_idx, TC_xml in enumerate(xml_alloc):
            if tc_idx > len(mod_energy)-1: break
            n_link = math.floor(TC_xml['glb_channel']/3)

            if TC_xml['frame'] not in data_TCs.keys():
                data_TCs[TC_xml['frame']] = {}
            if n_link not in data_TCs[TC_xml['frame']].keys():
                data_TCs[TC_xml['frame']][n_link] = [[0]*3]*3

            value_energy, code_energy = compress_value(mod_energy[tc_idx]/LSB)
            value_r_z = int(mod_r_over_z[tc_idx]/LSB_r_z) & 0xFFF # 12 bits
            value_phi = int((mod_phi[tc_idx]-offset_phi)/LSB_phi) & 0xFFF # 12 bits

            if args.plot:
                data_heatmap.append({
                    'rOverZ': value_r_z,
                    'phi'   : value_phi,
                    'column': columns[tc_idx],
                    'energy': value_energy
                })

            data_TCs[TC_xml['frame']][n_link][TC_xml['channel']%3] = [
                code_energy, value_r_z, value_phi
                ]
   
    if args.plot: shift.append(plot.create_plot_py(data_heatmap, ds_gen, args))
    return data_TCs
    
def data_packer(ds_TCs, ds_gen, args, shift):
    data_TCs = process_event(ds_TCs, ds_gen, args, shift)

    # packing data in links 
    data_links = {}
    for frame in data_TCs.keys():
        for n_link in data_TCs[frame].keys():
            link_data = create_link(data_TCs[frame][n_link])
            data_links[84*frame+n_link] = link_data
    
    return data_links
