import matplotlib.pyplot as plt
import numpy as np

def define_bin(r_z, phi=0):
    return int((r_z-440)/64)-1, int(124*phi/3.14)-1

def create_plot(objects, step, event, col):
    heatmap = np.zeros((64, 124))

    for bin in objects:
      if (step=='post_unpacking') and (bin.energy()>0):
         if col: heatmap[define_bin(bin.rOverZ())[0], bin.index()] += bin.energy()/100
         else:   heatmap[define_bin(bin.rOverZ(), (np.pi/1944)*bin.phi())] += bin.energy()/100
         # print("Energy", bin.energy(), "R/Z", int((bin.rOverZ()-440)/64), "Column", bin.index())

      if (step=='post_seeding') and (bin.S()>0):
         heatmap[bin.sortKey()-1, bin.index()-1] += (bin.S())/100
         # print("Smeared Energy : ", bin.S(), "R/Z bin", bin.sortKey(), "col", bin.index())

    if col: step = 'columns_' + step
    create_heatmap(heatmap, step, event)

def create_heatmap(heatmap, title, event_number):
    plt.imshow(heatmap, cmap='viridis', aspect='auto')
    plt.colorbar(label='Transverse Energy')
    plt.xlabel('Column, \u03C6 bin')
    plt.ylabel('R/Z bin')
    plt.title(f'{title} Histogram - Event {event_number}'.replace('_', ' '))

    plt.savefig(f'plots/{event_number}_{title}.pdf')
    plt.clf()

