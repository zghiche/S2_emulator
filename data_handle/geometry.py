import numpy as np
import awkward as ak
import uproot
import xml.etree.ElementTree as ET
from collections import defaultdict

# def read_geometry_txt():
#     ''' reads the geometry file from Pedro's txt '''
#     
#     geometry_file = '../../mapping/FE_mapping/geometries/v15.3/geometry.hgcal.txt'
#     return np.loadtxt(geometry_file, delimiter=' ', usecols=(0,1,2,37), skiprows=1) #'plane','u','v','hash'

def read_xml():
    tree = ET.parse('config_files/S1.ChannelAllocation.xml')
    root = tree.getroot()
    reversed_data = defaultdict(list)

    S1_index = 0
    for s1_element in root.findall('.//S1'):
        for channel_element in s1_element.findall('.//Channel'):
            channel = int(channel_element.get('aux-id'))
            for frame_element in channel_element.findall('.//Frame'):
                if all(attr in frame_element.attrib for attr in ['id', 'column', 'Module']):
                    frame  = int(frame_element.get('id'))
                    column = int(frame_element.get('column'))
                    module = hex(int(frame_element.get('Module'),16))
                    glb_channel = 12*S1_index+channel
                    index  = int(frame_element.get('index'))
                    reversed_data[module].append({'column'     : column,
                                                  'frame'      : frame, 
                                                  'channel'    : channel, 
                                                  'glb_channel': glb_channel,
                                                  'index'      : index})
        S1_index += 1
    return reversed_data

def apply_sort(df, counts, axis):
    for field in df.fields:
        df[field] = ak.unflatten(df[field], counts, axis)
    return df

def provide_event(tree, event):
    branches_tc = [
        'good_tc_x', 'good_tc_y', 'good_tc_z',
        'good_tc_phi', 'good_tc_layer',
        'good_tc_waferu', 'good_tc_waferv',
        'good_tc_energy', 'good_tc_mipPt'
    ]
    
    branches_gen = [
        'event', 'good_genpart_exeta', 'good_genpart_exphi', 'good_genpart_energy'
    ]

    data = tree.arrays(branches_tc, entry_start=event, entry_stop=event+1, library='ak')
    data_gen = tree.arrays(branches_gen, entry_start=event, entry_stop=event+1, library='ak')[0]
    data['r_over_z'] = np.sqrt(data.good_tc_x**2 + data.good_tc_y**2)/data.good_tc_z

    # sorting by modules  
    sorted_waferu = data[ak.argsort(data['good_tc_waferu'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferu.good_tc_waferu), axis=None)
    sorted_df = apply_sort(sorted_waferu, counts, 1)

    sorted_waferv = sorted_df[ak.argsort(sorted_df['good_tc_waferv'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferv.good_tc_waferv), axis=None)
    sorted_df = apply_sort(sorted_waferv, counts, 2)
    sorted_df = ak.flatten(sorted_df, axis=2)
    
    # sorting by transverse energy, simulating the ECONT_T
    sorted_df = sorted_df[ak.argsort(sorted_df['good_tc_mipPt'], ascending=False)][0]
    return [sorted_df, data_gen]
    
def provide_events(n=1):
    filepath = '/data_CMS/cms/ehle/L1HGCAL/PU0/photons/skims/skim_tightTC_dRxy_hadd.root' 
    name_tree = "FloatingpointMixedbcstcrealsig4DummyHistomaxxydr015GenmatchGenclustersntuple/HGCalTriggerNtuple"

    tree  = uproot.open(filepath)[name_tree]
    event, events = 0, []
    for n_ev in range(n):
        phi_gen, eta_gen = -1, -1
        event_found = False
    
        while not 0.2 < phi_gen < 1.8 or not 1.6 < eta_gen < 2.8:
            event += 1
            tree_ev = tree.arrays(['event', 'good_genpart_exphi', 'good_genpart_exeta'], entry_start=event,
                                  entry_stop=event+1, library='ak')
            
            phi_gen, eta_gen = tree_ev.good_genpart_exphi[0], tree_ev.good_genpart_exeta[0]
    
            if 0.2 < phi_gen < 1.8 and 1.6 < eta_gen < 2.8:
                event_found = True
                break
    
        if event_found:
            events.append(event)

    events_ds = [provide_event(tree, ev) for ev in events]
    return events_ds
