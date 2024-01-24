import numpy as np
import awkward as ak
import uproot
import xml.etree.ElementTree as ET
from collections import defaultdict

def read_geometry_txt():
    ''' reads the geometry file from Pedro's txt '''
    
    geometry_file = '../../mapping/FE_mapping/geometries/v15.3/geometry.hgcal.txt'
    return np.loadtxt(geometry_file, delimiter=' ', usecols=(0,1,2,37), skiprows=1) #'plane','u','v','hash'

def read_xml():
    tree = ET.parse('../../mapping/FE2BE_mapping/checks_FE_BE/xml/S1toChannels.SeparateTD.120.SingleTypes.NoSplit.xml')
    root = tree.getroot()

    reversed_data = defaultdict(lambda: defaultdict(list))

    for s1_element in root.findall('.//S1'):
        s1_id = s1_element.get('id')
        for channel_element in s1_element.findall('.//Channel'):
            channel = int(channel_element.get('id'))
            for frame_element in channel_element.findall('.//Frame'):
                if all(attr in frame_element.attrib for attr in ['id', 'column', 'Module']):
                    frame  = int(frame_element.get('id'))
                    column = int(frame_element.get('column'))
                    module = int(frame_element.get('Module'))
                    glb_channel = 12*int(s1_id[1:])+channel

                    reversed_data[module][column].append({'frame'      : frame, 
                                                          'channel'    : channel, 
                                                          'glb_channel': glb_channel})

    return reversed_data

def provide_events():
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
