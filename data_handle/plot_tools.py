import matplotlib.pyplot as plt
from scipy.stats import binomtest
import numpy as np
import yaml

with open('config.yaml', "r") as afile:
    cfg = yaml.safe_load(afile)

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
    thresholds = '_a'+'_'.join(map(str, [int(i/1000) for i in cfg['thresholdMaximaParam_a']]))+ \
                 '_b'+str(thr_b)+'_c'+str(cfg['thresholdMaximaParam_c'])
    plt.savefig('plots/histogram_seed_vs_'+variable+thresholds+'.pdf')
    plt.clf()

def plot_seeds(seeds, args):
    for thr_b in seeds.keys():
      n_params = len(cfg['thresholdMaximaParam_a'])
      seeds_list, p_t_list, eta_list, thr_list = [], [], [], []
      for index, thr in enumerate(cfg['thresholdMaximaParam_a']):
          seeds_list.append([seed[0] for idx, seed in enumerate(seeds[thr_b]) if idx%n_params == index])
          eta_list.append([eta[1] for idx, eta in enumerate(seeds[thr_b]) if idx%n_params == index])
          p_t_list.append([p_t[2] for idx, p_t in enumerate(seeds[thr_b]) if idx%n_params == index])
          thr_list.append('a:'+str(int(thr/1000))+' b:'+str(thr_b))
      
      plt.hist(seeds_list, bins=4, label=thr_list) 
      produce_efficiency_plots('thr', thr_b, args)

      for thr in range(len(thr_list)):
          compute_efficiency_plots(seeds_list[thr], p_t_list[thr], thr_list[thr], 7)
      produce_efficiency_plots('pT', thr_b, args)
 
      for thr in range(len(thr_list)):
          compute_efficiency_plots(seeds_list[thr], eta_list[thr], thr_list[thr], 5)
      produce_efficiency_plots('eta', thr_b, args)

def produce_plots(shift_pre, shift_post):
    create_histo([r_z[0] for r_z in shift_pre], 'r_z', 'pre_post_unpacking',
                 [r_z[0] for r_z in shift_post])
    create_histo([phi[1] for phi in shift_pre], 'phi', 'pre_post_unpacking',
                 [phi[1] for phi in shift_post])
    create_histo([p_t[2] for p_t in shift_pre], 'p_t', 'pre_post_unpacking',
                 [p_t[2] for p_t in shift_post])

def define_bin(r_z, phi=0):
    return int((r_z-440)/64)-1, int(124*phi/3.14)-1

def distance(tc, gen):
    r_z_bin, phi_bin = bin2coord(tc.sortKey()-1+0.5, -23+tc.index()-1+0.5)
    eta_bin = -np.log(np.tan((r_z_bin*0.7)/(4096*2)))
    print(np.sqrt((eta_bin-gen.eta_gen)**2+(phi_bin-gen.phi_gen)**2))
    return np.sqrt((eta_bin-gen.eta_gen)**2+(phi_bin-gen.phi_gen)**2)

def bin2coord(r_z_bin, phi_bin):
    return 64*(r_z_bin+1)+440, np.pi*(phi_bin+1)/124

def calculate_shift(heatmap, ev):
    max_bin = np.unravel_index(np.argmax(heatmap), heatmap.shape)
    max_r_z, max_phi = bin2coord(max_bin[0]+0.5, max_bin[1]+0.5)
    r_over_z = np.tan(2*np.arctan(np.exp(-ev.eta_gen)))
    return [r_over_z-max_r_z*0.7/4096, 
            ev.phi_gen - max_phi+0.3 ,
            np.sum(heatmap)/ev.pT_gen]

def create_plot_py(objects, ev, args):
    heatmap = np.zeros((64, 124))

    for bin in objects: # min col = -23 from S2.ChannelAllocation.xml
        if args.col: heatmap[define_bin(bin['rOverZ'])[0], 23+bin['column']] += bin['energy']/10000
        else: heatmap[define_bin(bin['rOverZ'], (np.pi/1944)*bin['phi'])] += bin['energy']/10000

    if args.performance: return calculate_shift(heatmap, ev) 
    step = 'columns_pre_unpacking' if args.col else 'pre_unpacking'
    create_heatmap(heatmap, step, ev.event)

def create_plot(objects, step, ev, args):
    heatmap = np.zeros((64, 124))

    seed = 0
    for bin in objects:
      if (step=='post_unpacking') and (bin.energy()>0):
        if args.col:  heatmap[define_bin(bin.rOverZ())[0], bin.index()] += bin.energy()/10000
        else: heatmap[define_bin(bin.rOverZ(), (np.pi/1944)*bin.phi())] += bin.energy()/10000
        # print("Energy", bin.energy(), "R/Z", int((bin.rOverZ()-440)/64), "Column", bin.index())

      if (step=='post_seeding') and (bin.S()>0):
        heatmap[bin2coord(bin.sortKey()-1, bin.index()-1)] += (bin.S())/10000
        # print("Smeared Energy : ", bin.S(), "R/Z bin", bin.sortKey(), "col", bin.index())
        if (bin.maximaOffset() == cfg['fanoutWidths'][bin.sortKey()]) and \
           (distance(bin, ev)<10): seed += 1    
  
    if (seed == 0 and ev.pT_gen>20 and ev.eta_gen>2.4): 
        print('No seed found for event', ev.event) 
        create_heatmap(heatmap, step, ev.event)
    if args.performance: return calculate_shift(heatmap, ev)
    if args.thr_seed: return [seed, ev.eta_gen, ev.pT_gen]
    if args.col: step = 'columns_' + step
    create_heatmap(heatmap, step, ev.event)

def create_heatmap(heatmap, title, event_number):
    plt.imshow(heatmap, cmap='viridis', aspect='auto')
    plt.colorbar(label='Transverse Energy')
    plt.xlabel('Column, \u03C6 bin')
    plt.ylabel('R/Z bin')
    plt.title(f'{title} Histogram - Event {event_number}'.replace('_', ' '))

    plt.savefig(f'plots/{event_number}_{title}.pdf')
    plt.clf()

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
    plt.clf()
