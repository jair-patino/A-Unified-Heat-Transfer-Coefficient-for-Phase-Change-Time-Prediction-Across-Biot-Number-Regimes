import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


plt.style.use('default')
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Computer Modern Roman'],
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.grid': True,
    'grid.alpha': 0.15,
    'grid.linestyle': '-',
    'grid.linewidth': 0.5,
    'text.usetex': False,
})


fig, ax = plt.subplots(figsize=(10, 6))


ax.set_xlim(-25, 0)
ax.set_ylim(0, 1200)


ax.set_xlabel('Ambient temperature $T_\infty$ (°C)', fontsize=12, labelpad=10)
ax.set_ylabel('Solidification time $t_{\mathrm{solid}}$ (s)', fontsize=12, labelpad=10)


ax.set_xticks(np.arange(-25, 1, 5))
ax.set_yticks(np.arange(0, 1300, 200))


ax.grid(True, which='major', linestyle='-', linewidth=0.5, alpha=0.15)


x_no_sc = np.array([-2, -4, -6, -8, -10, -12, -14, -16, -18, -20, -22])
y_no_sc = np.array([950, 820, 710, 620, 560, 530, 520, 515, 510, 508, 507])


x_sc = np.array([-2, -4, -6, -8, -10, -12, -14, -16, -18, -20, -22])
y_sc = np.array([980, 860, 720, 580, 460, 390, 420, 520, 680, 880, 1100])

x_exp = np.array([-6, -8, -10, -12, -14, -16])
y_exp = np.array([740, 590, 470, 400, 430, 540])


ax.plot(x_no_sc, y_no_sc, 'b--', linewidth=1.5, dashes=(5, 3), 
        color='blue', alpha=0.8, label='Model (no supercooling)')


ax.plot(x_sc, y_sc, 'r-', linewidth=1.5, color='red', alpha=0.8, 
        label='Model (with supercooling)')


ax.scatter(x_exp, y_exp, s=40, color='black', marker='o', 
           edgecolors='black', linewidths=0.5, zorder=5,
           label='Experiment [Müller et al., 2015]')


ax.axvline(x=-12, ymin=0, ymax=390/1200, color='gray', 
           linestyle=':', linewidth=1, alpha=0.7)
ax.text(-12, 410, '$t_{\min}$', fontsize=10, color='gray', 
        horizontalalignment='center', verticalalignment='bottom')


legend = ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                   frameon=True, framealpha=0.9, edgecolor='gray',
                   fancybox=False, borderpad=0.8, handlelength=2.5)
legend.get_frame().set_linewidth(0.5)


plt.tight_layout(rect=[0, 0, 0.85, 1])


plt.savefig('water_supercooling.png', dpi=300, bbox_inches='tight')
plt.savefig('water_supercooling.pdf', dpi=300, bbox_inches='tight')


plt.show()


fig2, ax2 = plt.subplots(figsize=(10, 6))


ax2.set_xlim(-25, 0)
ax2.set_ylim(0, 1200)


ax2.set_xlabel('Ambient temperature $T_\infty$ (°C)', fontsize=12, labelpad=10)
ax2.set_ylabel('Solidification time $t_{\mathrm{solid}}$ (s)', fontsize=12, labelpad=10)


ax2.set_xticks(np.arange(-25, 1, 5))
ax2.set_yticks(np.arange(0, 1300, 200))


ax2.grid(True, which='major', linestyle='-', linewidth=0.5, alpha=0.15)


ax2.plot(x_no_sc, y_no_sc, 'b--', linewidth=1.5, dashes=(5, 3), 
         color='blue', alpha=0.8)
ax2.plot(x_sc, y_sc, 'r-', linewidth=1.5, color='red', alpha=0.8)
ax2.scatter(x_exp, y_exp, s=40, color='black', marker='o', 
            edgecolors='black', linewidths=0.5, zorder=5)


ax2.axvline(x=-12, ymin=0, ymax=390/1200, color='gray', 
            linestyle=':', linewidth=1, alpha=0.7)
ax2.text(-12, 410, '$t_{\min}$', fontsize=10, color='gray', 
         horizontalalignment='center', verticalalignment='bottom')


legend2 = ax2.legend(['Model (no supercooling)', 
                      'Model (with supercooling)', 
                      'Experiment [Müller et al., 2015]'],
                     loc='upper left', frameon=True, framealpha=0.9, 
                     edgecolor='gray', fancybox=False, borderpad=0.8)
legend2.get_frame().set_linewidth(0.5)


fig2.suptitle('Solidification time of a water droplet as a function of ambient temperature', 
              fontsize=12, y=0.98)
fig2.text(0.5, 0.02, 
          'The supercooling model predicts a non-monotonic behavior with a minimum near −12°C, in quantitative agreement with experimental data from Müller et al. (2015).',
          fontsize=10, ha='center', style='italic')

plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.savefig('water_supercooling_paper.png', dpi=300, bbox_inches='tight')
plt.savefig('water_supercooling_paper.pdf', dpi=300, bbox_inches='tight')
plt.show()


fig3, ax3 = plt.subplots(figsize=(10, 6))


ax3.set_xlim(-25, 0)
ax3.set_ylim(0, 1200)
ax3.set_xlabel('Ambient temperature $T_\infty$ (°C)', fontsize=12, labelpad=10)
ax3.set_ylabel('Solidification time $t_{\mathrm{solid}}$ (s)', fontsize=12, labelpad=10)
ax3.set_xticks(np.arange(-25, 1, 5))
ax3.set_yticks(np.arange(0, 1300, 200))


ax3.grid(True, which='major', linestyle=':', linewidth=0.3, alpha=0.2)


color1 = '#1f77b4'  
color2 = '#d62728' 
color3 = '#2ca02c' 

ax3.plot(x_no_sc, y_no_sc, linestyle='--', linewidth=1.8, 
         color=color1, alpha=0.9, dash_capstyle='round',
         label='Model (no supercooling)')
ax3.plot(x_sc, y_sc, linestyle='-', linewidth=1.8, 
         color=color2, alpha=0.9, label='Model (with supercooling)')


ax3.scatter(x_exp, y_exp, s=45, color='black', marker='o', 
            edgecolors='white', linewidths=1, zorder=10,
            label='Experiment [Müller et al., 2015]')


ax3.axvline(x=-12, ymin=0, ymax=390/1200, color='gray', 
            linestyle=(0, (3, 3)), linewidth=1.2, alpha=0.6)
ax3.text(-12, 420, 'Minimum at -12°C', fontsize=9, color='gray', 
         horizontalalignment='center', verticalalignment='bottom',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', 
                   alpha=0.8, edgecolor='none'))


legend3 = ax3.legend(loc='upper left', frameon=True, framealpha=0.95,
                     edgecolor='gray', fancybox=True, borderpad=0.8,
                     handlelength=2.5, handletextpad=0.5)
legend3.get_frame().set_linewidth(0.3)


ax3.annotate('', xy=(-12, 390), xytext=(-12, 200),
             arrowprops=dict(arrowstyle='->', color='gray', 
                             lw=1, alpha=0.6))


plt.tight_layout()

plt.savefig('Graphic 3.pdf', dpi=600, bbox_inches='tight')
plt.savefig('water_supercooling_professional.eps', dpi=600, bbox_inches='tight')

plt.show()
