import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 6))

ax.set_xlim(-25, 0)
ax.set_ylim(0, 1200)

ax.set_xlabel('Ambient temperature $T_\\infty$ (°C)')
ax.set_ylabel('Solidification time $t_{\\mathrm{solid}}$ (s)')

ax.set_xticks(np.arange(-25, 1, 5))
ax.set_yticks(np.arange(0, 1300, 200))

x_no_sc = np.array([-2, -4, -6, -8, -10, -12, -14, -16, -18, -20, -22])
y_no_sc = np.array([950, 820, 710, 620, 560, 530, 520, 515, 510, 508, 507])

x_sc = np.array([-2, -4, -6, -8, -10, -12, -14, -16, -18, -20, -22])
y_sc = np.array([980, 860, 720, 580, 460, 390, 420, 520, 680, 880, 1100])

x_exp = np.array([-6, -8, -10, -12, -14, -16])
y_exp = np.array([740, 590, 470, 400, 430, 540])

ax.plot(x_no_sc, y_no_sc, '--', label='Model (no supercooling)')
ax.plot(x_sc, y_sc, '-', label='Model (with supercooling)')
ax.scatter(x_exp, y_exp, label='Experiment')

ax.axvline(x=-12)
ax.text(-12, 420, '$t_{min}$', ha='center')

ax.legend()

plt.tight_layout()
plt.savefig('water_supercooling_simple.png', dpi=300)
plt.show()
