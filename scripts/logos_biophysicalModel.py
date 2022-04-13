import numpy as np
import logomaker
import pandas as pd
import matplotlib.pyplot as plt

motif = pd.read_csv("{INPUTFILE}", sep="\t", header=None)
motif = motif.set_axis(['A', 'C', 'G', 'T'], axis=1, inplace=False)
motif.index = np.arange(1, len(motif) + 1)
crp_logo = logomaker.Logo(motif,
                          shade_below=.5,
                          fade_below=.5)
# style using Logo methods
crp_logo.style_spines(visible=False)
crp_logo.style_spines(spines=['left', 'bottom'], visible=True)
crp_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
crp_logo.ax.set_ylabel("Binding Energy", labelpad=-1)
crp_logo.ax.xaxis.set_ticks_position('none')
crp_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig("{OUTPUTFILE}")
