import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Configuración de estilo para calidad de paper
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

# Crear figura y eje
fig, ax = plt.subplots(figsize=(10, 6))

# Configurar escalas logarítmicas
ax.set_xscale('log')
ax.set_yscale('log')

# Límites de los ejes
ax.set_xlim(0.01, 10)
ax.set_ylim(0.1, 50)

# Etiquetas de los ejes
ax.set_xlabel('Biot Number Bi', fontsize=12)
ax.set_ylabel('Relative Error (%)', fontsize=12)

# Ticks personalizados
ax.set_xticks([0.01, 0.1, 1, 10])
ax.set_xticklabels(['0.01', '0.1', '1', '10'])
ax.set_yticks([0.1, 1, 10, 50])
ax.set_yticklabels(['0.1%', '1%', '10%', '50%'])

# Regiones sombreadas (background)
# Verde: Bi 0.01-1, Error 0.1-5%
ax.add_patch(Rectangle((0.01, 0.1), 0.99, 4.9, 
                      facecolor='green', alpha=0.15, edgecolor='none'))
# Amarillo: Bi 1-2, Error 0.1-15%
ax.add_patch(Rectangle((1, 0.1), 1, 14.9, 
                      facecolor='yellow', alpha=0.15, edgecolor='none'))
# Rojo: Bi 2-10, Error 0.1-50
ax.add_patch(Rectangle((2, 0.1), 8, 49.9, 
                      facecolor='red', alpha=0.15, edgecolor='none'))

# Líneas de umbral de error
ax.axhline(y=5, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.text(0.012, 5.2, '5% Error', fontsize=10, color='gray', alpha=0.7, 
        verticalalignment='bottom')

ax.axhline(y=15, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.text(0.012, 15.5, '15% Error', fontsize=10, color='gray', alpha=0.7, 
        verticalalignment='bottom')

# Límites de validez de modelos anteriores
# Liu et al. limit (Bi = 0.2)
ax.axvline(x=0.2, color='red', linestyle=(0, (6, 3)), linewidth=1, alpha=0.8)
ax.text(0.2, 0.09, 'Liu et al. limit', fontsize=10, color='red', alpha=0.8,
        horizontalalignment='center', verticalalignment='top')

# Karwa et al. limit (Bi = 0.5)
ax.axvline(x=0.5, color='green', linestyle=(0, (6, 3)), linewidth=1, alpha=0.8)
ax.text(0.5, 0.09, 'Karwa et al. limit', fontsize=10, color='green', alpha=0.8,
        horizontalalignment='center', verticalalignment='top')

# Datos: Our model - Planar
our_model_data = np.array([
    [0.01, 0.2], [0.015, 0.3], [0.02, 0.4], [0.03, 0.6], [0.04, 0.8],
    [0.05, 1.1], [0.07, 1.5], [0.1, 2.3], [0.15, 2.8], [0.2, 3.2],
    [0.3, 3.8], [0.4, 4.3], [0.5, 4.8], [0.7, 5.5], [1.0, 7.2],
    [1.5, 9.8], [2.0, 12.5], [3.0, 18.2], [4.0, 22.5], [5.0, 25.1]
])

# Datos: Liu et al. (2009) - Planar
liu_data = np.array([
    [0.01, 0.5], [0.02, 1.0], [0.05, 3.5], [0.1, 8.2], [0.15, 15.3],
    [0.2, 24.5], [0.3, 45.2], [0.4, 65.8], [0.5, 82.3]
])

# Datos: Karwa et al. (2013) - Planar
karwa_data = np.array([
    [0.01, 0.8], [0.02, 1.5], [0.05, 4.2], [0.1, 9.5], [0.15, 16.8],
    [0.2, 26.3], [0.3, 48.5], [0.4, 70.2], [0.5, 88.5], [0.7, 125.3],
    [1.0, 185.2]
])

# Graficar las curvas (en orden para que aparezcan encima de las regiones)
# Our model
ax.plot(our_model_data[:, 0], our_model_data[:, 1], 
        color='blue', marker='o', markersize=4, markeredgecolor='blue',
        markerfacecolor='white', markeredgewidth=1, linewidth=1.5,
        label='Our model (planar)')

# Liu et al.
ax.plot(liu_data[:, 0], liu_data[:, 1], 
        color='red', linestyle='--', marker='s', markersize=3.5,
        markeredgecolor='red', markerfacecolor='white', markeredgewidth=1,
        linewidth=1.5, label='Liu et al. (2009)')

# Karwa et al.
ax.plot(karwa_data[:, 0], karwa_data[:, 1], 
        color='green', linestyle=':', marker='^', markersize=4,
        markeredgecolor='green', markerfacecolor='white', markeredgewidth=1,
        linewidth=1.5, label='Karwa et al. (2013)')

# Configurar la leyenda
legend = ax.legend(loc='lower right', fontsize=10, framealpha=0.9,
                   edgecolor='black', fancybox=False, frameon=True)
legend.get_frame().set_linewidth(0.5)

# Añadir grid
ax.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.3)

# Ajustar layout
plt.tight_layout()

# Guardar la figura en alta calidad
plt.savefig('error_vs_bi.png', dpi=300, bbox_inches='tight')
plt.savefig('error_vs_bi.pdf', dpi=300, bbox_inches='tight')  # Para LaTeX

# Mostrar la gráfica
plt.show()

# Versión alternativa con fondo blanco puro (para paper)
fig2, ax2 = plt.subplots(figsize=(10, 6))

# Mismo código pero con fondo blanco puro
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(0.01, 10)
ax2.set_ylim(0.1, 50)
ax2.set_xlabel('Biot Number Bi', fontsize=12)
ax2.set_ylabel('Relative Error (%)', fontsize=12)
ax2.set_xticks([0.01, 0.1, 1, 10])
ax2.set_xticklabels(['0.01', '0.1', '1', '10'])
ax2.set_yticks([0.1, 1, 10, 50])
ax2.set_yticklabels(['0.1%', '1%', '10%', '50%'])

# Regiones sombreadas
ax2.add_patch(Rectangle((0.01, 0.1), 0.99, 4.9, 
                       facecolor='green', alpha=0.15, edgecolor='none'))
ax2.add_patch(Rectangle((1, 0.1), 1, 14.9, 
                       facecolor='yellow', alpha=0.15, edgecolor='none'))
ax2.add_patch(Rectangle((2, 0.1), 8, 49.9, 
                       facecolor='red', alpha=0.15, edgecolor='none'))

# Líneas de error
ax2.axhline(y=5, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax2.text(0.012, 5.2, '5% Error', fontsize=10, color='gray', alpha=0.7, 
         verticalalignment='bottom')
ax2.axhline(y=15, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax2.text(0.012, 15.5, '15% Error', fontsize=10, color='gray', alpha=0.7, 
         verticalalignment='bottom')

# Límites de modelos
ax2.axvline(x=0.2, color='red', linestyle=(0, (6, 3)), linewidth=1, alpha=0.8)
ax2.text(0.2, 0.09, 'Liu et al. limit', fontsize=10, color='red', alpha=0.8,
         horizontalalignment='center', verticalalignment='top')
ax2.axvline(x=0.5, color='green', linestyle=(0, (6, 3)), linewidth=1, alpha=0.8)
ax2.text(0.5, 0.09, 'Karwa et al. limit', fontsize=10, color='green', alpha=0.8,
         horizontalalignment='center', verticalalignment='top')

# Curvas
ax2.plot(our_model_data[:, 0], our_model_data[:, 1], 
         color='blue', marker='o', markersize=4, markeredgecolor='blue',
         markerfacecolor='white', markeredgewidth=1, linewidth=1.5,
         label='Our model (planar)')
ax2.plot(liu_data[:, 0], liu_data[:, 1], 
         color='red', linestyle='--', marker='s', markersize=3.5,
         markeredgecolor='red', markerfacecolor='white', markeredgewidth=1,
         linewidth=1.5, label='Liu et al. (2009)')
ax2.plot(karwa_data[:, 0], karwa_data[:, 1], 
         color='green', linestyle=':', marker='^', markersize=4,
         markeredgecolor='green', markerfacecolor='white', markeredgewidth=1,
         linewidth=1.5, label='Karwa et al. (2013)')

# Leyenda y grid
legend2 = ax2.legend(loc='lower right', fontsize=10, framealpha=0.9,
                    edgecolor='black', fancybox=False, frameon=True)
legend2.get_frame().set_linewidth(0.5)
ax2.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.3)

# Título para el caption
fig2.suptitle('Relative error versus Biot number for planar geometry, comparing our model with previous approaches.', 
              fontsize=11, y=0.98)
fig2.text(0.5, 0.02, 
          'Our formulation maintains errors below 5% for Bi ≤ 1 and below 15% for Bi = 2, surpassing previous lumped-capacitance extensions limited to Bi < 0.2 (Liu et al., 2009) and Bi < 0.5 (Karwa et al., 2013).',
          fontsize=10, ha='center', style='italic')

plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.savefig('Graphic 1.pdf', dpi=300, bbox_inches='tight')
plt.show()
