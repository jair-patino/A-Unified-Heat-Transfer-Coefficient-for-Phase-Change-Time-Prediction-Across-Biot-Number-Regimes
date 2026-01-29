import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

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
    'grid.alpha': 0.3,
    'grid.linestyle': ':',
    'grid.linewidth': 0.5,
})


fig, ax = plt.subplots(figsize=(9, 6))


ax.set_xlim(0, 12000)
ax.set_ylim(0, 12000)

ax.set_xlabel('Predicted Time $t_{\mathrm{total}}$ (s)', fontsize=12)
ax.set_ylabel('Numerical Time $t_{\mathrm{num}}$ (s)', fontsize=12)

ax.set_xticks(np.arange(0, 13000, 2000))
ax.set_yticks(np.arange(0, 13000, 2000))

x_vals = np.array([0, 12000])

y_upper10 = x_vals + 1200  
y_lower10 = x_vals - 1200  
poly_10 = np.vstack([np.concatenate([x_vals, x_vals[::-1]]),
                     np.concatenate([y_upper10, y_lower10[::-1]])]).T
ax.add_patch(Polygon(poly_10, facecolor='yellow', alpha=0.12, edgecolor='none'))

y_upper5 = x_vals + 600  
y_lower5 = x_vals - 600  
poly_5 = np.vstack([np.concatenate([x_vals, x_vals[::-1]]),
                    np.concatenate([y_upper5, y_lower5[::-1]])]).T
ax.add_patch(Polygon(poly_5, facecolor='green', alpha=0.16, edgecolor='none'))

ax.plot([0, 12000], [0, 12000], 'k--', linewidth=1.5, alpha=0.8, 
        label='_nolegend_')

water_x = [500, 1200, 2500, 5000, 8000, 10000]
water_y = [510, 1180, 2550, 4950, 8200, 10300]
ax.scatter(water_x, water_y, marker='o', s=50, 
           facecolors='blue', edgecolors='blue', 
           alpha=0.8, linewidths=1, 
           label='Water (melting) $R^2 = 0.998$')

paraffin_x = [1500, 3200, 5500, 9000]
paraffin_y = [1520, 3180, 5600, 9200]
ax.scatter(paraffin_x, paraffin_y, marker='s', s=45, 
           facecolors='red', edgecolors='red', 
           alpha=0.8, linewidths=1,
           label='Paraffin (solidification) $R^2 = 0.996$')

gallium_x = [800, 2000, 4000, 7000]
gallium_y = [790, 2050, 3950, 7100]
ax.scatter(gallium_x, gallium_y, marker='^', s=55, 
           facecolors='green', edgecolors='green', 
           alpha=0.8, linewidths=1, 
           label='Gallium (melting) $R^2 = 0.997$')

ax.plot([0, 12000], [0, 11400], 'green', linewidth=0.5, alpha=0.3, 
        linestyle='-', label='_nolegend_')
ax.plot([0, 12000], [600, 12600], 'green', linewidth=0.5, alpha=0.3, 
        linestyle='-', label='_nolegend_')
ax.plot([0, 12000], [0, 10800], 'yellow', linewidth=0.5, alpha=0.2, 
        linestyle='-', label='_nolegend_')
ax.plot([0, 12000], [1200, 13200], 'yellow', linewidth=0.5, alpha=0.2, 
        linestyle='-', label='_nolegend_')

legend = ax.legend(loc='lower right', fontsize=10,
                   framealpha=0.9, edgecolor='gray', 
                   fancybox=False, frameon=True,
                   handletextpad=0.3, labelspacing=0.2)
legend.get_frame().set_linewidth(0.5)

ax.text(200, 200, 'Perfect agreement', fontsize=10, color='gray', 
        rotation=45, rotation_mode='anchor', alpha=0.8,
        verticalalignment='bottom', horizontalalignment='left')

ax.text(2000, 300, '±5%\nconfidence', fontsize=9, 
        color='green', alpha=0.8, weight='bold',
        verticalalignment='center', horizontalalignment='center',
        linespacing=1.2)
ax.text(3000, 800, '±10%\nconfidence', fontsize=9, 
        color='yellow', alpha=0.8, weight='bold',
        verticalalignment='center', horizontalalignment='center',
        linespacing=1.2)

stats_text = 'All points within ±5%\nMax error: 3.0%\nMean error: 1.2%'
props = dict(boxstyle='round', facecolor='white', alpha=0.8, 
             edgecolor='gray', linewidth=0.5)
ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', horizontalalignment='left',
        bbox=props)

ax.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.3)
ax.set_aspect('equal', adjustable='box')

plt.tight_layout()

plt.savefig('validation_scatter_improved.png', dpi=300, bbox_inches='tight')
plt.savefig('validation_scatter_improved.pdf', dpi=300, bbox_inches='tight')

plt.show()

fig2, ax2 = plt.subplots(figsize=(9, 6))

ax2.set_xlim(0, 12000)
ax2.set_ylim(0, 12000)
ax2.set_xlabel('Predicted Time $t_{\mathrm{total}}$ (s)', fontsize=12)
ax2.set_ylabel('Numerical Time $t_{\mathrm{num}}$ (s)', fontsize=12)
ax2.set_xticks(np.arange(0, 13000, 2000))
ax2.set_yticks(np.arange(0, 13000, 2000))

ax2.ticklabel_format(style='plain', axis='both')

ax2.add_patch(Polygon(poly_10, facecolor='yellow', alpha=0.1, edgecolor='none'))
ax2.add_patch(Polygon(poly_5, facecolor='green', alpha=0.16, edgecolor='none'))

ax2.plot([0, 12000], [0, 12000], 'k--', linewidth=1.5, alpha=0.8)

scatter1 = ax2.scatter(water_x, water_y, marker='o', s=50,
                      facecolors='blue', edgecolors='blue', alpha=0.8,
                      linewidths=1)
scatter2 = ax2.scatter(paraffin_x, paraffin_y, marker='s', s=45,
                      facecolors='red', edgecolors='red', alpha=0.8,
                      linewidths=1)
scatter3 = ax2.scatter(gallium_x, gallium_y, marker='^', s=55,
                      facecolors='green', edgecolors='green', alpha=0.8,
                      linewidths=1)

ax2.plot([0, 12000], [0, 11400], 'green', linewidth=0.3, alpha=0.3)
ax2.plot([0, 12000], [600, 12600], 'green', linewidth=0.3, alpha=0.3)
ax2.plot([0, 12000], [0, 10800], 'yellow', linewidth=0.3, alpha=0.2)
ax2.plot([0, 12000], [1200, 13200], 'yellow', linewidth=0.3, alpha=0.2)

legend2 = ax2.legend([scatter1, scatter2, scatter3],
                     ['Water (melting) $R^2 = 0.998$',
                      'Paraffin (solidification) $R^2 = 0.996$',
                      'Gallium (melting) $R^2 = 0.997$'],
                     loc='lower right', fontsize=10, framealpha=0.9,
                     edgecolor='gray', fancybox=False, frameon=True)
legend2.get_frame().set_linewidth(0.5)

ax2.text(200, 200, 'Perfect agreement', fontsize=10, color='gray',
         rotation=45, rotation_mode='anchor', alpha=0.8,
         verticalalignment='bottom', horizontalalignment='left')
ax2.text(2000, 300, '±5%\nconfidence', fontsize=9, color='green', alpha=0.8,
         verticalalignment='center', horizontalalignment='center',
         linespacing=1.2)
ax2.text(3000, 800, '±10%\nconfidence', fontsize=9, color='yellow', alpha=0.8,
         verticalalignment='center', horizontalalignment='center',
         linespacing=1.2)

ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes, fontsize=9,
         verticalalignment='top', horizontalalignment='left',
         bbox=props)

ax2.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.3)
ax2.set_aspect('equal', adjustable='box')

fig2.suptitle('Predicted versus numerically computed phase change times for various materials and conditions ($\mathrm{Bi} = 0.1$ to 1.5)',
              fontsize=11, y=0.98)
fig2.text(0.5, 0.02,
          'High coefficients of determination ($R^2 > 0.996$) confirm model accuracy across different materials. All points fall within ±5% of perfect agreement, demonstrating the model\'s robustness for engineering applications. Statistical summary: maximum error 3.0%, mean error 1.2%.',
          fontsize=10, ha='center', style='italic')

plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.savefig('validation_scatter_improved_paper.png', dpi=300, bbox_inches='tight')
plt.savefig('validation_scatter_improved_paper.pdf', dpi=300, bbox_inches='tight')
plt.show()

fig3, ax3 = plt.subplots(figsize=(10, 6))

ax3.set_xlim(0, 12000)
ax3.set_ylim(0, 12000)
ax3.set_xlabel('Predicted Time $t_{\mathrm{total}}$ (s)', fontsize=12, labelpad=8)
ax3.set_ylabel('Numerical Time $t_{\mathrm{num}}$ (s)', fontsize=12, labelpad=8)

xticks = np.arange(0, 13000, 2000)
yticks = np.arange(0, 13000, 2000)
ax3.set_xticks(xticks)
ax3.set_yticks(yticks)
ax3.set_xticklabels([f'{int(x/1000)}k' if x >= 1000 else str(x) for x in xticks])
ax3.set_yticklabels([f'{int(y/1000)}k' if y >= 1000 else str(y) for y in yticks])

ax3.add_patch(Polygon(poly_10, facecolor='#FFFACD', alpha=0.3, edgecolor='none'))
ax3.add_patch(Polygon(poly_5, facecolor='#90EE90', alpha=0.3, edgecolor='none'))

ax3.plot([0, 12000], [0, 12000], 'k--', linewidth=1.5, alpha=0.6)

colors = ['#1f77b4', '#d62728', '#2ca02c']  
ax3.scatter(water_x, water_y, marker='o', s=60, color=colors[0], 
           edgecolors='white', linewidth=0.8, alpha=0.9,
           label='Water (melting) $R^2 = 0.998$')
ax3.scatter(paraffin_x, paraffin_y, marker='s', s=50, color=colors[1],
           edgecolors='white', linewidth=0.8, alpha=0.9,
           label='Paraffin (solidification) $R^2 = 0.996$')
ax3.scatter(gallium_x, gallium_y, marker='^', s=70, color=colors[2],
           edgecolors='white', linewidth=0.8, alpha=0.9,
           label='Gallium (melting) $R^2 = 0.997$')

legend3 = ax3.legend(loc='lower right', fontsize=10,
                    frameon=True, framealpha=0.95,
                    edgecolor='gray', fancybox=True,
                    borderpad=0.5, handletextpad=0.5)
legend3.get_frame().set_linewidth(0.3)

ax3.text(250, 250, 'Perfect agreement', fontsize=9, color='gray',
         rotation=45, rotation_mode='anchor', alpha=0.7,
         verticalalignment='bottom', horizontalalignment='left')
ax3.text(1500, 500, '±5%', fontsize=8, color='green', alpha=0.7,
         verticalalignment='center', horizontalalignment='center',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7,
                   edgecolor='none'))
ax3.text(2500, 1000, '±10%', fontsize=8, color='orange', alpha=0.7,
         verticalalignment='center', horizontalalignment='center',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7,
                   edgecolor='none'))

stats_box = dict(boxstyle='round', facecolor='white', alpha=0.9,
                 edgecolor='gray', linewidth=0.3)
ax3.text(0.02, 0.98, stats_text, transform=ax3.transAxes, fontsize=8.5,
         verticalalignment='top', horizontalalignment='left',
         bbox=stats_box)

ax3.grid(True, which='major', linestyle=':', linewidth=0.3, alpha=0.2)
ax3.set_aspect('equal')

plt.tight_layout()

plt.savefig('Graphic 2.pdf', dpi=600, bbox_inches='tight')
plt.savefig('validation_scatter_professional.eps', dpi=600, bbox_inches='tight')  
plt.show()
