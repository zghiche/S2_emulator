import xml.etree.ElementTree as ET
from collections import defaultdict
import matplotlib.pyplot as plt

def ObjectType( object_type ):
  return ( ( object_type & 0xF ) << 22 )

def extract_S1_values(s1_output):
    s1_output_index = s1_output & 0x7F # range [0, 18*6-1]
    S1id = s1_output & ~(ObjectType(2) | ((s1_output_index & 0x7F) << 0))
    return hex((S1id>>16)&0x3F), s1_output_index

def extract_S2_values(s2_input):
    s2_input_index = s2_input & 0x7F # range [0, 18*6-1]
    S2id = s2_input & ~(ObjectType(1) | (( s2_input_index & 0x7F ) << 0 ))
    return hex(S2id), s2_input_index

def read_backend_mapping():
    tree = ET.parse('config_files/BackendMapping.xml')
    root = tree.getroot()
    reversed_data = []

    for endcap in root.findall('.//Endcap'):
        for fibre in endcap.findall('.//fibre'):
            fibre_id  = fibre.get('id')
            fibre_ref = fibre.get('ref-fibre')
            S1id, s1_output_index = extract_S1_values(int(fibre.get('s1-output'),16))
            S2id, s2_input_index  = extract_S2_values(int(fibre.get('s2-input'),16))
            reversed_data.append({'S1id'            : S1id,
                                  's1_output_index' : s1_output_index, 
                                  'S2id'            : S2id, 
                                  's2_input_index'  : s2_input_index,
                                  'fibre_ref'       : fibre_ref})
        break
    return reversed_data

if __name__ == '__main__':
    xml = read_backend_mapping()

    # print(set([x['S1id'] for x in xml]))
    S2id = [x['S2id'] for x in xml]
    S1id = [x['S1id'] for x in xml]
    S2_devices = list(set(S2id))
    S1_devices = list(set(S1id))
    s1_output = [x['s1_output_index'] for x in xml if x['S2id'] in S2_devices]
    s2_input  = [x['s2_input_index']  for x in xml if x['S2id'] in S2_devices]
    S2id = [x['S2id'] for x in xml if x['S2id'] in S2_devices]

    print(len(S2_devices))
    cmap = plt.cm.get_cmap('rainbow', len(S1_devices))
    color_map = {category: cmap(i) for i, category in enumerate(S1_devices)}
    colors = [color_map[string] for string in [x['S1id'] for x in xml if x['S2id'] in S2_devices]]
    
    plt.scatter(s1_output, s2_input, s=3, c=colors)
    
    plt.xlabel('s1_output')
    plt.ylabel('s2_input')
    plt.title('Scatter Plot')
    plt.savefig('s1_vs_s2_fibres.pdf')
