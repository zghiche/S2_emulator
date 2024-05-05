import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import binomtest
import numpy as np
import yaml

white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0,    '#ffffff'),
    (1e-10,'#440053'),
    (0.2,  '#404388'),
    (0.4,  '#2a788e'),
    (0.6,  '#21a784'),
    (0.8,  '#78d151'),
    (1,    '#fde624'),
], N=1000)

with open('config.yaml', "r") as afile:
    cfg = yaml.safe_load(afile)

def distance(tc, gen):
    r_z_bin, phi_bin = bin2coord(tc.sortKey()+0.5, tc.index()+0.5)
    eta_bin = -np.log(np.tan((r_z_bin*gen.LSB_r_z)/2))
    return np.sqrt((eta_bin-gen.eta_gen)**2+(phi_bin-gen.phi_gen)**2)

def define_bin(r_z, phi=0):
    return int((r_z-440)/64), 23+int(124*phi/np.pi)

def bin2coord(r_z_bin, phi_bin):
    return 64*(r_z_bin)+440, np.pi*(phi_bin-23)/124

def calculate_shift(heatmap, ev):
    max_bin = np.unravel_index(np.argmax(heatmap), heatmap.shape)
    max_r_z, max_phi = bin2coord(max_bin[0]+0.5, max_bin[1]+0.5)
    r_over_z = np.tan(2*np.arctan(np.exp(-ev.eta_gen)))
    return [r_over_z-max_r_z*ev.LSB_r_z, 
            ev.phi_gen - max_phi,
            np.sum(heatmap)/ev.pT_gen]

def create_plot(objects, step, ev, args):
    heatmap, seed = np.zeros((64, 124)), []

    for bin in objects:
      if (step=='post_unpacking') and (bin.energy()>0):
        if args.phi: heatmap[define_bin(bin.rOverZ(), ev.LSB_phi*bin.phi()+ev.offset_phi)] += bin.energy()*ev.LSB
        elif args.col: heatmap[define_bin(bin.rOverZ())[0], bin.index()] += bin.energy()*ev.LSB
        # print("Energy", bin.energy(), "R/Z", int((bin.rOverZ()-440)/64), "Column", bin.index())

      if (step=='post_seeding') and (bin.S()>0):
        heatmap[bin.sortKey(), bin.index()] += (bin.S())*ev.LSB
        # print("Smeared Energy : ", bin.S(), "R/Z bin", bin.sortKey(), "col", bin.index())
        if bin.maximaOffset() == cfg['fanoutWidths'][bin.sortKey()]: seed.append([bin.sortKey(), bin.index()])    

    if len(seed) == 3 and distance(bin, ev) < 10:
        print(f'3 seeds found for event {ev.event}, (pT, \u03B7, \u03C6)=({ev.pT_gen:.0f}, {ev.eta_gen:.2f},{ev.phi_gen:.2f})') 
        create_heatmap(heatmap, step, ev, seed)
    if args.performance: return calculate_shift(heatmap, ev)
    elif args.col or args.phi: create_heatmap(heatmap, 'columns_'+step if (args.col and step!='post_seeding') else step, ev, seed)
    if args.thr_seed: return [len(seed), ev.eta_gen, ev.pT_gen]

def hgcal_limits(ev):
    plt.axhline(y=(0.476/ev.LSB_r_z-440)/64, color='red', linestyle='--')
    plt.text(3, ((0.476+0.015)/ev.LSB_r_z-440)/64, 'Layer1', color='red', fontsize=6)
    plt.axhline(y=(0.462/ev.LSB_r_z-440)/64, color='red', linestyle='--')
    plt.text(3, ((0.462-0.02)/ev.LSB_r_z-440)/64, 'Layer27', color='red', fontsize=6)
    plt.axhline(y=(0.085/ev.LSB_r_z-440)/64, color='red', linestyle='--')
    plt.text(3, ((0.085+0.015)/ev.LSB_r_z-440)/64, 'Layer1', color='red', fontsize=6)
    plt.axhline(y=(0.076/ev.LSB_r_z-440)/64, color='red', linestyle='--')
    plt.text(3, ((0.076-0.02)/ev.LSB_r_z-440)/64, 'Layer27', color='red', fontsize=6)

def create_heatmap(heatmap, title, gen, seeds=[]):
    plt.imshow(heatmap, cmap=white_viridis, origin='lower', aspect='auto')
    x_tick_labels = [int(val) for val in np.linspace(-30, 150, num=7)]
    y_tick_labels = ['{:.2f}'.format(val) for val in np.linspace(440*gen.LSB_r_z, (64**2+440)*gen.LSB_r_z, num=8)]
    plt.xticks(np.linspace(0, 123, num=7), labels=x_tick_labels)
    plt.yticks(np.linspace(1, 64,  num=8), labels=y_tick_labels)
    plt.colorbar(label='Transverse Energy [GeV]')
    plt.xlabel('\u03C6 (degrees)')
    plt.ylabel('r/z')
    plt.scatter(23+int(124*gen.phi_gen/np.pi), (np.tan(2*np.arctan(np.exp(-gen.eta_gen)))/gen.LSB_r_z-440)/64, 
                color='red', marker='x', s=50)
    for seed in seeds:
      plt.scatter(seed[1], seed[0], color='white', marker='o', s=20)
    plt.title(f'{title} - Event {gen.event} \n pT:{gen.pT_gen:.0f} GeV, \u03B7:{gen.eta_gen:.2f}, \u03C6:{gen.phi_gen:.2f}'.replace('_', ' '))
    plt.grid(linestyle='--')
    hgcal_limits(gen)
    plt.savefig(f'plots/single_events/{gen.event}_{title}.pdf')
    plt.savefig(f'plots/single_events/{gen.event}_{title}.png')
    plt.clf()

def produce_plots(shift_pre, shift_post):
    create_histo([r_z[0] for r_z in shift_pre], 'r_z', 'pre_post_unpacking',
                 [r_z[0] for r_z in shift_post])
    create_histo([phi[1] for phi in shift_pre], 'phi', 'pre_post_unpacking',
                 [phi[1] for phi in shift_post])
    create_histo([p_t[2] for p_t in shift_pre], 'p_t', 'pre_post_unpacking',
                 [p_t[2] for p_t in shift_post])

def create_histo(data, variable, title, data2=None, data3=None):
    if variable == 'p_t':
        plt.hist(data,  bins=25, alpha=0.5, label='pre unpacking')
        plt.hist(data2, bins=25, alpha=0.5, label='post unpacking')
    if variable == 'r_z': 
        plt.hist(data, bins=25, range=(-0.01,0.01), alpha=0.5, label='pre unpacking')
        plt.hist(data2, bins=25, range=(-0.01,0.01), alpha=0.5, label='post unpacking')
    if variable == 'phi': 
        plt.hist(data, bins=25, alpha=0.5, label='columns pre')
        plt.hist(data2, bins=25, alpha=0.5, label='columns post')
    plt.legend()
    xlabel = (r'$\phi_{bin}^{max energy} - \phi_{gen particle}$' if variable == 'phi' else
          r'$\frac{r}{z}_{bin}^{max energy} - \frac{r}{z}_{gen particle}$' if variable == 'r_z' else
          r'$p^{T}_{bin} / p^{T}_{gen particle}$' )
    plt.xlabel(xlabel)
    plt.ylabel('Counts')
    plt.title('Histogram '+ title + ' ' + xlabel)
    plt.savefig('plots/histogram_'+ variable + '_' + title +'.pdf')
    plt.savefig('plots/histogram_'+ variable + '_' + title +'.png')
    plt.clf()

#########################################################
################## Efficiency plots #####################
#########################################################

def compute_efficiency_plots(seeds, variable, thr, bin_n=10):
    bin_edges = np.linspace(min(variable), max(variable), num=bin_n+1)
    indices = np.digitize(variable, bin_edges) - 1
 
    eff, lo_err, up_err = {}, {}, {}
    for index in range(bin_n):
      """k is number of successes, n is number of trials"""
      bin_indices = np.where(indices == index)[0]
      seeds_bin = [seeds[i] for i in bin_indices]
      k, n = sum(1 for x in seeds_bin if x >= 1), len(seeds_bin)
      if n == 0: eff[index], lo_err[index], up_err[index] = 0, 0, 0; continue

      result = binomtest(k, n, p=k/n)
      eff[index] = k/n
      lo_err[index] = k/n - result.proportion_ci(confidence_level=0.95).low
      up_err[index] = result.proportion_ci(confidence_level=0.95).high - k/n

    plt.errorbar((bin_edges[1:] + bin_edges[:-1])/2, eff.values(), 
                 yerr=np.array(list(zip(lo_err.values(), up_err.values()))).T,
                 xerr=(bin_edges[1] - bin_edges[0])/2, fmt='o', capsize=3, label=thr, alpha=0.7) 

def produce_efficiency_plots(variable, thr_b, args):
    plt.legend()
    plt.grid()
    plt.xlabel('identified seeds' if variable=='thr' else r'$p_{T}$ [GeV]' if variable=='pT' else r'$\eta$')
    plt.ylabel('Counts' if variable=='thr' else '% of identified seeds')
    plt.title('Histogram 2D seeds '+variable)
    thresholds = '_a'+'_'.join(map(str, [int(i*10) for i in cfg['thresholdMaximaParam_a']]))+ \
                 '_b'+str(thr_b)+'_c'+str(cfg['thresholdMaximaParam_c'])
    plt.savefig('plots/histogram_seed_vs_'+variable+thresholds+'.pdf')
    plt.savefig('plots/histogram_seed_vs_'+variable+thresholds+'.png')
    plt.clf()

def plot_seeds(seeds, args):
    for thr_b in seeds.keys():
      n_params = len(cfg['thresholdMaximaParam_a'])
      seeds_list, p_t_list, eta_list, thr_list = [], [], [], []
      for index, thr in enumerate(cfg['thresholdMaximaParam_a']):
          seeds_list.append([seed[0] for idx, seed in enumerate(seeds[thr_b]) if idx%n_params == index])
          eta_list.append([eta[1] for idx, eta in enumerate(seeds[thr_b]) if idx%n_params == index])
          p_t_list.append([p_t[2] for idx, p_t in enumerate(seeds[thr_b]) if idx%n_params == index])
          thr_list.append('a:'+str(thr)+' GeV b:'+str(thr_b))
      
      plt.hist(seeds_list, bins=4, label=thr_list) 
      produce_efficiency_plots('thr', thr_b, args)

      for thr in range(len(thr_list)):
          compute_efficiency_plots(seeds_list[thr], p_t_list[thr], thr_list[thr], 7)
      produce_efficiency_plots('pT', thr_b, args)
 
      for thr in range(len(thr_list)):
          compute_efficiency_plots(seeds_list[thr], eta_list[thr], thr_list[thr], 5)
      produce_efficiency_plots('eta', thr_b, args)

## not used ##
def create_plot_py(objects, ev, args):
    heatmap = np.zeros((64, 124))

    for bin in objects: # min col = -23 from S2.ChannelAllocation.xml
        if args.phi: heatmap[define_bin(bin['rOverZ'], ev.LSB_phi*bin['phi']+ev.offset_phi)] += bin['energy']*ev.LSB
        elif args.col: heatmap[define_bin(bin['rOverZ'])[0], 23+bin['column']] += bin['energy']*ev.LSB

    if args.performance: return calculate_shift(heatmap, ev) 
    elif args.col or args.phi: create_heatmap(heatmap, 'columns_pre_unpacking' if args.col else 'pre_unpacking', ev)

