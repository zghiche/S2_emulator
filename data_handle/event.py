import numpy as np
import awkward as ak
import uproot
import math
import cppyy
from cppyy.gbl import l1thgcfirmware, std

import data_handle.plot_tools as plot
import data_handle.geometry as geometry
from data_handle.tools import compress_value

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
        print('Analysing module ', plane, v-u, v)
        if plane & ~0x3F : raise Exception( "Invalid plane" )
        if v-u & ~0xF : raise Exception( "Invalid u" )
        if v   & ~0xF : raise Exception( "Invalid v" )
        return hex(0x60000000 | self.ObjectType(0) | ((plane & 0x3F) << 16) | ((v-u & 0xF) << 12) | ((v & 0xF) << 8))
    
    def get_TC_allocation(self, xml_data, module):
        return xml_data[module]
    
    def create_link(self, data_in):
        energy   = (data_in[2][0] << 30) | (data_in[1][0] << 15) | (data_in[0][0])
        r_over_z = (data_in[2][1] << 30) | (data_in[1][1] << 15) | (data_in[0][1])
        phi      = (data_in[2][2] << 30) | (data_in[1][2] << 15) | (data_in[0][2])
        return [energy, r_over_z, phi]
   
     
    def _process_event(self, args, shift):
        LSB = 1/100000 # 10 keV
        LSB_r_z = 0.7/4096
        LSB_phi = np.pi/1944
        offset_phi = -0.3
        data_TCs = {}
        xml_data = geometry.read_xml()
    
        if args.plot: data_heatmap = []
        for module_idx in range(len(self.ds_TCs.good_tc_x)):
            module = self.get_module_id(self.ds_TCs.good_tc_layer[module_idx][0],
                                        self.ds_TCs.good_tc_waferu[module_idx][0],
                                        self.ds_TCs.good_tc_waferv[module_idx][0])
            xml_alloc = self.get_TC_allocation(xml_data, module)
            if not xml_alloc: continue
    
            # calculating the number of TC that ca be allocated / module
            n_TCs = xml_alloc[-1]['index']  # dangerous
            columns = [frame['column'] for frame in xml_alloc]
    
            # simulating the BC algorithm (ECON-T) and the phi sorting in the S1 FPGA
            mod_phi = self.ds_TCs.good_tc_phi[module_idx][:n_TCs+1]
            mod_energy = self.ds_TCs.good_tc_mipPt[module_idx][:n_TCs+1][ak.argsort(mod_phi)]
            mod_r_over_z = self.ds_TCs.r_over_z[module_idx][:n_TCs+1][ak.argsort(mod_phi)]
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
    
        if args.plot: shift.append(plot.create_plot_py(data_heatmap, self.gen, args))
        return data_TCs
    
    def _data_packer(self, args, shift):
        data_TCs = self._process_event(args, shift)
    
        # packing data in links 
        data_links = {}
        for frame in data_TCs.keys():
            for n_link in data_TCs[frame].keys():
                link_data = self.create_link(data_TCs[frame][n_link])
                data_links[84*frame+n_link] = link_data

        # filling data into emulator c++ variables
        LinksInData = std.vector['std::unique_ptr<l1thgcfirmware::HGCalLinkTriggerCell>']()
        HGCalLinkTriggerCell = l1thgcfirmware.HGCalLinkTriggerCell

        # number of links = NChannels*Nframes
        for link_frame in range(84*108):
            LinksInData.push_back(std.make_unique[HGCalLinkTriggerCell]())
            if link_frame in data_links.keys():
                LinksInData[-1].data_     = data_links[link_frame][0]
                LinksInData[-1].r_over_z_ = data_links[link_frame][1]
                LinksInData[-1].phi_      = data_links[link_frame][2]

        self.data_packer = LinksInData


#######################################################################################
############################### PROVIDE EVENTS ########################################
#######################################################################################

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
    return EventData(sorted_df, data_gen)

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
