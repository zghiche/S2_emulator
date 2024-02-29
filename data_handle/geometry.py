import xml.etree.ElementTree as ET
from collections import defaultdict

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

