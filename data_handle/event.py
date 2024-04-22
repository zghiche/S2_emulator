import numpy as np
import awkward as ak
import uproot
import math
import cppyy
from cppyy.gbl import l1thgcfirmware, std

import data_handle.plot_tools as plot
from data_handle.tools import compress_value
        
cppyy.cppdef("""
#import "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkTriggerCell_SA.h"
std::vector<std::unique_ptr<l1thgcfirmware::HGCalLinkTriggerCell>> fill_Links(std::map<int, std::vector<long int>>& data_links) {
std::vector<std::unique_ptr<l1thgcfirmware::HGCalLinkTriggerCell>> LinksInData;
for (int link_frame = 0; link_frame < 84 * 108; ++link_frame) {
    LinksInData.push_back(std::make_unique<l1thgcfirmware::HGCalLinkTriggerCell>());
    if (data_links.find(link_frame) != data_links.end()) {
        LinksInData.back()->data_     = data_links[link_frame][0];
        LinksInData.back()->r_over_z_ = data_links[link_frame][1];
        LinksInData.back()->phi_      = data_links[link_frame][2];
    }
};
return LinksInData;
}""")

cppyy.cppdef("""
std::vector<long int> create_link(const std::map<int, std::vector<long int>>& data_in) {
    int64_t energy = 0, r_over_z = 0, phi = 0;
    for (const auto& pair : data_in) {
        int channel = pair.first;
        energy   = energy | (pair.second[0] << (channel*15));
        r_over_z = r_over_z | (pair.second[1] << (channel*15));
        phi      = phi | (pair.second[2] << (channel*15));
    }
    return {energy, r_over_z, phi};
}

std::map<int, std::vector<long int>> pack_Links(std::map<int, std::map<int, std::map<int, std::vector<long int>>>>& data_TCs) {
std::map<int, std::vector<long int>> data_links;
// Iterate over data_TCs
for (const auto& frame_pair : data_TCs) {
    int frame = frame_pair.first;
    for (const auto& n_link_pair : frame_pair.second) {
        int n_link = n_link_pair.first;
        // Create link data
        // std::cout << n_link_pair.first << std::endl;
        std::vector<long int> link_data = create_link(n_link_pair.second);
        // Insert into data_links
        data_links[84 * frame + n_link] = link_data;
    }
}
return data_links;
}""")

class EventData():
    def __init__(self, ds_TCs, gen):
        self.ds_TCs  = ds_TCs
        self.gen     = gen
        self.event   = gen.event
        self.eta_gen = gen.good_genpart_exeta[0]
        self.phi_gen = gen.good_genpart_exphi[0]
        self.pT_gen  = self._compute_pt(self.eta_gen,
                             gen.good_genpart_energy[0])

        self.data_packer = None

    def _compute_pt(self, eta, energy):
        return energy/np.cosh(eta)

    def ObjectType(self, object_type):
        return ((object_type & 0xF) << 22)
    
    def get_module_id(self, plane, u, v):
        # CMSSW to our u v convention u'=v-u, v'=v
        # print('Analysing module ', plane, v-u, v)
        if plane & ~0x3F : return 0 # raise Exception( "Invalid plane" )
        if v-u & ~0xF : return 0 # raise Exception( "Invalid u" )
        if v   & ~0xF : return 0 # raise Exception( "Invalid v" )
        return hex(0x60000000 | self.ObjectType(0) | ((plane & 0x3F) << 16) | ((v-u & 0xF) << 12) | ((v & 0xF) << 8))
    
    def get_TC_allocation(self, xml_data, module):
        return xml_data[module]
     
    def _process_event(self, args, xml, shift):
        LSB = 1/10000 # 100 keV
        LSB_r_z = 0.7/4096
        LSB_phi = np.pi/1944
        offset_phi = -0.3
        data_TCs = std.map[int,std.map[int,std.map[int,'std::vector<long int>']]]()

        # if args.plot: data_heatmap = []
        for module_idx in range(len(self.ds_TCs.good_tc_x)):
            module = self.get_module_id(self.ds_TCs.good_tc_layer[module_idx][0],
                                        self.ds_TCs.good_tc_waferu[module_idx][0],
                                        self.ds_TCs.good_tc_waferv[module_idx][0])
            xml_alloc = self.get_TC_allocation(xml, module)
            if not xml_alloc: continue
    
            # calculating the number of TC that ca be allocated / module
            n_TCs = xml_alloc[-1]['index']  # dangerous
            columns = [frame['column'] for frame in xml_alloc]
   
            # simulating the BC algorithm (ECON-T) and the phi sorting in the S1 FPGA
            mod_phi = self.ds_TCs.good_tc_phi[module_idx][:n_TCs+1]
            mod_energy = self.ds_TCs.good_tc_pt[module_idx][:n_TCs+1][ak.argsort(mod_phi)]
            mod_r_over_z = self.ds_TCs.r_over_z[module_idx][:n_TCs+1][ak.argsort(mod_phi)]
            mod_phi = ak.sort(mod_phi)

            # assigning each TCs to a columns
            xml_alloc = sorted(xml_alloc, key=lambda x: x['column'])

            for tc_idx, TC_xml in enumerate(xml_alloc):
                if tc_idx > len(mod_energy)-1: break
                n_link = TC_xml['n_link']
    
                value_energy, code_energy = compress_value(mod_energy[tc_idx]/LSB)
                value_r_z = int(mod_r_over_z[tc_idx]/LSB_r_z) & 0xFFF # 12 bits
                value_phi = int((mod_phi[tc_idx]-offset_phi)/LSB_phi) & 0xFFF # 12 bits
                # if args.plot:
                #     data_heatmap.append({
                #         'rOverZ': value_r_z,
                #         'phi'   : value_phi,
                #         'column': columns[tc_idx],
                #         'energy': value_energy
                #     })

                data_TCs[TC_xml['frame']][n_link][TC_xml['channel']%3] = [ 
                    code_energy, value_r_z, value_phi
                    ]

        # if args.plot: shift.append(plot.create_plot_py(data_heatmap, self, args))
        return data_TCs
    
    def _data_packer(self, args, xml, shift):
        data_TCs = self._process_event(args, xml, shift)
        data_links = cppyy.gbl.pack_Links(data_TCs)

        # filling data into emulator c++ variables
        self.data_packer = cppyy.gbl.fill_Links(data_links) 


#######################################################################################
############################### PROVIDE EVENTS ########################################
#######################################################################################

def apply_sort(df, counts, axis):
    for field in df.fields:
        df[field] = ak.unflatten(df[field], counts, axis)
    return df

def provide_event(ev, gen):
    ev['r_over_z'] = np.sqrt(ev.good_tc_x**2 + ev.good_tc_y**2)/ev.good_tc_z

    # sorting by modules  
    sorted_waferu = ev[ak.argsort(ev['good_tc_waferu'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferu.good_tc_waferu), axis=None)
    sorted_df = apply_sort(sorted_waferu, counts, 1)

    sorted_waferv = sorted_df[ak.argsort(sorted_df['good_tc_waferv'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferv.good_tc_waferv), axis=None)
    sorted_df = apply_sort(sorted_waferv, counts, 2)

    sorted_layer = sorted_df[ak.argsort(sorted_df['good_tc_layer'])]
    counts = ak.flatten(ak.run_lengths(sorted_layer.good_tc_layer), axis=None)
    sorted_df = apply_sort(sorted_layer, counts, 3)
    sorted_df = ak.flatten(sorted_df, axis=3)
    sorted_df = ak.flatten(sorted_df, axis=2)

    # sorting by transverse energy, simulating the ECONT_T
    sorted_df = sorted_df[ak.argsort(sorted_df['good_tc_pt'], ascending=False)][0]
    return EventData(sorted_df, gen)

def provide_events(n=1):
    filepath = '/data_CMS/cms/mchiusi/L1HGCAL/skim_photons_emulator_02to18phi.root'
    name_tree = "FloatingpointMixedbcstcrealsig4DummyHistomaxxydr015GenmatchGenclustersntuple/HGCalTriggerNtuple"

    branches_tc = [
        'good_tc_x', 'good_tc_y', 'good_tc_z',
        'good_tc_phi', 'good_tc_layer',
        'good_tc_waferu', 'good_tc_waferv',
        'good_tc_energy', 'good_tc_pt'
    ]

    branches_gen = [
        'event', 'good_genpart_exeta', 'good_genpart_exphi', 'good_genpart_energy'
    ]

    tree  = uproot.open(filepath)[name_tree]
    events_ds = []
    for ev in range(n):
        data = tree.arrays(branches_tc, entry_start=ev, entry_stop=ev+1, library='ak')
        data_gen = tree.arrays(branches_gen, entry_start=ev, entry_stop=ev+1, library='ak')[0]
        events_ds.append(provide_event(data, data_gen))
    return events_ds
