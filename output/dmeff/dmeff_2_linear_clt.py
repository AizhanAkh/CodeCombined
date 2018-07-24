import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/u/aizhana/Projects/CodeCombined/output/dmeff/dmeff_2_linear_clt.dat', '/u/aizhana/Projects/CodeCombined/output/dmeff/dmeff_2_eft_clt.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['dmeff_2_linear_clt', 'dmeff_2_eft_clt']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = ['phiphi']
tex_names = ['phiphi']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], curve[:, 5])

index, curve = 1, data[1]
y_axis = ['phiphi']
tex_names = ['phiphi']
x_axis = 'l'
ylim = []
xlim = []
ax.plot(curve[:, 0], curve[:, 5])

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
plt.show()