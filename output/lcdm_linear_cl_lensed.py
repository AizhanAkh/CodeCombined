import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/u/aizhana/Projects/CodeCombined/output/lcdm_linear_cl_lensed.dat', '/u/aizhana/Projects/CodeCombined/output/dmeff/lcdm_halofit_cl_lensed.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['lcdm_linear_cl_lensed', 'lcdm_halofit_cl_lensed']

fig, ax = plt.subplots()
y_axis = [u'phiphi']
tex_names = ['phiphi']
x_axis = 'l'
ax.set_xlabel('$\ell$', fontsize=16)
plt.show()