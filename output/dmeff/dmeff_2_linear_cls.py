import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/u/aizhana/Projects/CodeCombined/output/dmeff/dmeff_2_linear_cls.dat', '/u/aizhana/Projects/CodeCombined/output/dmeff/dmeff_2_eft_cls.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['dmeff_2_linear_cls', 'dmeff_2_eft_cls']

fig, ax = plt.subplots()
y_axis = ['phiphi']
tex_names = ['phiphi']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()