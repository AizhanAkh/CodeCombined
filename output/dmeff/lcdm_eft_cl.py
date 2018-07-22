import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/u/aizhana/Projects/CodeCombined/output/dmeff/lcdm_eft_cl.dat', '/u/aizhana/Projects/CodeCombined/output/dmeff/dmeff_-2_eft_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['lcdm_eft_cl', 'dmeff_-2_eft_cl']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'phiphi']
tex_names = ['phiphi']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 5]))

index, curve = 1, data[1]
y_axis = [u'phiphi']
tex_names = ['phiphi']
x_axis = 'l'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 5]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
plt.show()