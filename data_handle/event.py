import numpy as np
import awkward as ak
import uproot
import math
import cppyy
from cppyy.gbl import l1thgcfirmware, std

import data_handle.plot_tools as plot
from data_handle.tools import compress_value, printProgressBar
        
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
        std::vector<long int> link_data = create_link(n_link_pair.second);
        // Insert into data_links
        data_links[84 * frame + n_link] = link_data;
    }
}
return data_links;
}""")

class EventData():
    def __init__(self, ds_si, ds_sci, gen):
        self.ds_si  = ds_si
        self.ds_sci  = ds_sci
        self.gen     = gen
        self.event   = gen.event
        self.eta_gen = gen.good_genpart_exeta[0]
        self.phi_gen = gen.good_genpart_exphi[0]
        self.pT_gen  = self._compute_pt(self.eta_gen,
                             gen.good_genpart_energy[0])

        self.data_packer = None
        self.LSB = 1/10000 # 100 keV
        self.LSB_r_z = 0.7/4096
        self.LSB_phi = np.pi/1944
        self.offset_phi = -0.3

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
    
    def get_MB_id(self, plane, v, MB):
        return MB[int(plane)][int(v)]
    
    def get_TC_allocation(self, xml_data, module):
        return xml_data[module]
    
    def _process_module(self, ds_TCs, idx, xml_alloc, data_TCs):
        n_TCs = xml_alloc[-1]['index']  # dangerous
        columns = [frame['column'] for frame in xml_alloc]
   
        # simulating the BC algorithm (ECON-T) and the phi sorting in the S1 FPGA
        mod_phi = ds_TCs.good_tc_phi[idx][:n_TCs+1]
        mod_energy = ds_TCs.good_tc_pt[idx][:n_TCs+1][ak.argsort(mod_phi)]
        mod_r_over_z = ds_TCs.r_over_z[idx][:n_TCs+1][ak.argsort(mod_phi)]
        mod_phi = ak.sort(mod_phi)

        # assigning each TCs to a columns
        xml_alloc = sorted(xml_alloc, key=lambda x: x['column'])
        
        for tc_idx, TC_xml in enumerate(xml_alloc):
            if tc_idx > len(mod_energy)-1: break
            n_link = TC_xml['n_link']
    
            value_energy, code_energy = compress_value(mod_energy[tc_idx]/self.LSB)
            value_r_z = int(mod_r_over_z[tc_idx]/self.LSB_r_z) & 0xFFF # 12 bits
            value_phi = int((mod_phi[tc_idx]-self.offset_phi)/self.LSB_phi) & 0xFFF # 12 bits

            data_TCs[TC_xml['frame']][n_link][TC_xml['channel']%3] = [ 
                code_energy, value_r_z, value_phi
                ]
 
    def _process_event(self, args, xml, MB_conv):
        data_TCs = std.map[int,std.map[int,std.map[int,'std::vector<long int>']]]()

        for module_idx in range(len(self.ds_si.good_tc_layer)):
            module = self.get_module_id(self.ds_si.good_tc_layer[module_idx][0],
                                        self.ds_si.good_tc_waferu[module_idx][0],
                                        self.ds_si.good_tc_waferv[module_idx][0])
            xml_alloc = self.get_TC_allocation(xml[0], module)
            if xml_alloc: self._process_module(self.ds_si, module_idx, xml_alloc, data_TCs)

        for MB_idx in range(len(self.ds_sci.good_tc_layer)):
            if self.ds_sci.MB_v[MB_idx][0] > 11: continue
            if self.ds_sci.good_tc_layer[MB_idx][0] > 47: continue
            MB = self.get_MB_id(self.ds_sci.good_tc_layer[MB_idx][0],
                                self.ds_sci.MB_v[MB_idx][0], MB_conv)
            xml_alloc = self.get_TC_allocation(xml[1], MB)
            if xml_alloc: 
                self._process_module(self.ds_sci, MB_idx, xml_alloc, data_TCs)
        return data_TCs
    
    def _data_packer(self, args, xml, xml_MB):
        data_TCs = self._process_event(args, xml, xml_MB)
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
    ev['MB_v'] = np.floor(ev.good_tc_cellv/4)
    ev = ev[[x for x in ak.fields(ev) if not x in ["good_tc_x","good_tc_y","good_tc_z","good_tc_cellv"]]]
    
    # dividing silicon and scintillators
    sci = ev[ev['good_tc_subdet'] == 10]
    si  = ev[ev['good_tc_subdet'] != 10]

    # sorting by modules  
    sorted_waferu = si[ak.argsort(si['good_tc_waferu'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferu.good_tc_waferu), axis=None)
    sorted_si = apply_sort(sorted_waferu, counts, 1)

    sorted_waferv = sorted_si[ak.argsort(sorted_si['good_tc_waferv'])]
    counts = ak.flatten(ak.run_lengths(sorted_waferv.good_tc_waferv), axis=None)
    sorted_si = apply_sort(sorted_waferv, counts, 2)

    sorted_layer = sorted_si[ak.argsort(sorted_si['good_tc_layer'])]
    counts = ak.flatten(ak.run_lengths(sorted_layer.good_tc_layer), axis=None)
    sorted_si = apply_sort(sorted_layer, counts, 3)
    sorted_si = ak.flatten(sorted_si, axis=3)
    sorted_si = ak.flatten(sorted_si, axis=2)

    # sorting by transverse energy, simulating the ECONT_T
    sorted_si = sorted_si[ak.argsort(sorted_si['good_tc_pt'], ascending=False)][0]
    
    # sorting sci by MB (cellv) and plane
    sorted_MB = sci[ak.argsort(sci['MB_v'])]
    counts = ak.flatten(ak.run_lengths(sorted_MB.MB_v), axis=None)
    sorted_sci = apply_sort(sorted_MB, counts, 1)

    sorted_layer = sorted_sci[ak.argsort(sorted_sci['good_tc_layer'])]
    counts = ak.flatten(ak.run_lengths(sorted_layer.good_tc_layer), axis=None)
    sorted_sci = apply_sort(sorted_layer, counts, 2)
    sorted_sci = ak.flatten(sorted_sci, axis=2)

    # sorting by transverse energy, simulating the ECONT_T
    sorted_sci = sorted_sci[ak.argsort(sorted_sci['good_tc_pt'], ascending=False)][0]

    return EventData(sorted_si, sorted_sci, gen)

def provide_events(n=1):
    # filepath = '/data_CMS/cms/mchiusi/ntupleProduction/DoublePhoton_FlatPt-1To100_PU200/DoublePhoton_FlatPt-1To100_PU200_TC_STC_STC.root'
    # name_tree = 'l1tHGCalTriggerNtuplizer/HGCalTriggerNtuple'
    # filepath = '/data_CMS/cms/mchiusi/L1HGCAL/skim_photons_emulator_02to18phi_new.root'
    filepath = '/data_CMS/cms/mchiusi/L1HGCAL/skim_pions_emulator_02to18phi_200PU.root'
    name_tree = "FloatingpointMixedbcstcrealsig4DummyHistomaxxydr015GenmatchGenclustersntuple/HGCalTriggerNtuple"

    branches_tc = [
        'good_tc_x', 'good_tc_y', 'good_tc_z',
        'good_tc_phi', 'good_tc_layer', 'good_tc_cellv',
        'good_tc_waferu', 'good_tc_waferv',
        'good_tc_pt', 'good_tc_subdet'
    ]

    branches_gen = [
        'event', 'good_genpart_exeta', 'good_genpart_exphi', 'good_genpart_energy'
    ]

    tree  = uproot.open(filepath)[name_tree]
    events_ds = []
    printProgressBar(0, n, prefix='Reading '+str(n)+' events from ROOT file:', suffix='Complete', length=50)
    for ev in range(n):
      data = tree.arrays(branches_tc, entry_start=ev, entry_stop=ev+1, library='ak')
      data_gen = tree.arrays(branches_gen, entry_start=ev, entry_stop=ev+1, library='ak')[0]
      events_ds.append(provide_event(data, data_gen))
      printProgressBar(ev+1, n, prefix='Reading '+str(n)+' events from ROOT file:', suffix='Complete', length=50)
    return events_ds
