import matplotlib.pyplot as plt
import numpy as np

def define_bin(r_z, phi=0):
    return int((r_z-440)/64)-1, int(124*phi/3.14)-1

def bin2coord(r_z_bin, phi_bin):
    return 64*(r_z_bin+1)+440, np.pi*(phi_bin+1)/124

def calculate_shift(heatmap, gen):
    max_bin = np.unravel_index(np.argmax(heatmap), heatmap.shape)
    max_r_z, max_phi = bin2coord(max_bin[0]+0.5, max_bin[1]+0.5)
    r_over_z = np.tan(2*np.arctan(np.exp(-gen.good_genpart_exeta[0])))
    print(max_r_z, max_phi)
    return [x - y for x, y in zip([r_over_z, gen.good_genpart_exphi[0]], 
                                       [max_r_z*0.7/4096, max_phi-0.3])]

def create_plot_py(objects, gen, args):
    heatmap = np.zeros((64, 124))

    for bin in objects:
        if args.col: heatmap[define_bin(bin['rOverZ'])[0], bin['column']] += bin['energy']/100000
        else: heatmap[define_bin(bin['rOverZ'], (np.pi/1944)*bin['phi'])] += bin['energy']/100000

    if args.performance: return calculate_shift(heatmap, gen) 
    step = 'columns_pre_unpacking' if args.col else 'pre_unpacking'
    create_heatmap(heatmap, step, gen.event)

def create_plot(objects, step, gen, args):
    if len(objects) <= 1: print(f"No TCs found for event {gen.event}"); return [0,0]
    heatmap = np.zeros((64, 124))

    for bin in objects:
      if (step=='post_unpacking') and (bin.energy()>0):
        if args.col: heatmap[define_bin(bin.rOverZ())[0], bin.index()] += bin.energy()/100000
        else: heatmap[define_bin(bin.rOverZ(), (np.pi/1944)*bin.phi())] += bin.energy()/100000
        # print("Energy", bin.energy(), "R/Z", int((bin.rOverZ()-440)/64), "Column", bin.index())

      if (step=='post_seeding') and (bin.S()>0):
        heatmap[bin.sortKey()-1, bin.index()-1] += (bin.S())/100000
        # print("Smeared Energy : ", bin.S(), "R/Z bin", bin.sortKey(), "col", bin.index())

    if args.performance: return calculate_shift(heatmap, gen) 
    if args.col: step = 'columns_' + step
    create_heatmap(heatmap, step, gen.event)

def create_heatmap(heatmap, title, event_number):
    plt.imshow(heatmap, cmap='viridis', aspect='auto')
    plt.colorbar(label='Transverse Energy')
    plt.xlabel('Column, \u03C6 bin')
    plt.ylabel('R/Z bin')
    plt.title(f'{title} Histogram - Event {event_number}'.replace('_', ' '))

    plt.savefig(f'plots/{event_number}_{title}.pdf')
    plt.clf()

def create_histo(data, variable, title):
    plt.hist(data, bins=40, range=(-0.1, 0.1), alpha=0.5)
    # plt.hist(data, bins=np.linspace(min(data), max(data)), alpha=0.5)
    xlabel = r'$\phi_{bin}^{max energy} - \phi_{gen particle}$' if variable=='phi' \
    else r'$\frac{r}{z}_{bin}^{max energy} - \frac{r}{z}_{gen particle}$'
    plt.xlabel(xlabel)
    plt.ylabel('Counts')
    plt.title('Histogram '+ title + ' ' + xlabel)
    plt.savefig('plots/histogram_'+ title + '_' + variable+'.pdf')
    plt.clf()
