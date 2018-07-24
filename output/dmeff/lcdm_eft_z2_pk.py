import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/u/aizhana/Projects/CodeCombined/output/dmeff/lcdm_eft_z2_pk.dat', '/u/aizhana/Projects/CodeCombined/output/dmeff/lcdm_eft_z1_pk.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['lcdm_eft_z2_pk', 'lcdm_eft_z1_pk']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = ['P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

index, curve = 1, data[1]
y_axis = ['P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('k (h/Mpc)', fontsize=16)
plt.show()