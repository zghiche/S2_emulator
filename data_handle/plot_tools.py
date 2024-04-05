import matplotlib.pyplot as plt
import numpy as np
import yaml

with open('config.yaml', "r") as afile:
    cfg = yaml.safe_load(afile)

def plot_seeds(seeds):
    n_params = len(cfg['thresholdMaximaParam_a'])
    seeds_list, pt_list, eta_list, thr_list = [], [], [], []
    for index, thr in enumerate(cfg['thresholdMaximaParam_a']):
        seeds_list.append([seed[0] for idx, seed in enumerate(seeds) if idx%n_params == index])
        eta_list.append([seed[1] for idx, seed in enumerate(seeds) if idx%n_params == index])
        pt_list.append([seed[2] for idx, seed in enumerate(seeds) if idx%n_params == index])
        thr_list.append('thr:'+str(thr))
    
    plt.hist(seeds_list, bins=5, label=thr_list) 
    plt.legend()
    plt.xlabel('identified seeds')
    plt.ylabel('Counts')
    plt.title('Histogram 2D seeds')
    plt.savefig('plots/histogram_seed_vs_thr_'+str(n_params)+'.pdf')
    plt.clf()
 
    for thr in range(n_params):
        # plt.scatter(eta_list[thr], list(np.asarray(seeds_list[thr])+0.1*thr), label=thr_list[thr]) 
        plt.scatter(pt_list[thr], list(np.asarray(seeds_list[thr])+0.1*thr), label=thr_list[thr]) 
    plt.legend()
    plt.ylabel('number of identified seeds')
    plt.xlabel('$pT [GeV]$')
    plt.title('Histogram 2D seeds')
    plt.savefig('plots/histogram_seed_vs_pT_'+str(n_params)+'.pdf')

def produce_plots(shift_pre, shift_post):
    create_histo([r_z[0] for r_z in shift_pre], 'r_z', 'pre_post_unpacking',
                 [r_z[0] for r_z in shift_post])
    create_histo([phi[1] for phi in shift_pre], 'phi', 'pre_post_unpacking',
                 [phi[1] for phi in shift_post])
    create_histo([p_t[2] for p_t in shift_pre], 'p_t', 'pre_post_unpacking',
                 [p_t[2] for p_t in shift_post])

def define_bin(r_z, phi=0):
    return int((r_z-440)/64)-1, int(124*phi/3.14)-1

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
        heatmap[bin.sortKey()-1, bin.index()-1] += (bin.S())/10000
        # print("Smeared Energy : ", bin.S(), "R/Z bin", bin.sortKey(), "col", bin.index())
        if (bin.maximaOffset() == cfg['fanoutWidths'][bin.sortKey()-1]): seed += 1    
   
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
