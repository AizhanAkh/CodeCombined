import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/u/aizhana/Projects/CodeCombined/output/test_z1_pk.dat', '/u/aizhana/Projects/CodeCombined/output/test_z1_pk_nl_pt.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['test_z1_pk', 'test_z1_pk_nl_pt']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'P(Mpc/h)^3']
tex_names = ['P (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))

index, curve = 1, data[1]
y_axis = [u'P_lin+P_1loop(Mpc/h)^3', u'P_CTR(Mpc/h)^3', u'Id2d2(Mpc/h)^3', u'Id2(Mpc/h)^3', u'IG2(Mpc/h)^3', u'Id2G2(Mpc/h)^3', u'IG2G2(Mpc/h)^3', u'IFG2(Mpc/h)^3']
tex_names = ['P_lin + P_1loop (Mpc/h)^3', 'P_CTR (Mpc/h)^3', 'Id2d2 (Mpc/h)^3', 'Id2 (Mpc/h)^3', 'IG2 (Mpc/h)^3', 'Id2G2 (Mpc/h)^3', 'IG2G2 (Mpc/h)^3', 'IFG2 (Mpc/h)^3']
x_axis = 'k (h/Mpc)'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 1]))
ax.loglog(curve[:, 0], abs(curve[:, 2]))
ax.loglog(curve[:, 0], abs(curve[:, 3]))
ax.loglog(curve[:, 0], abs(curve[:, 4]))
ax.loglog(curve[:, 0], abs(curve[:, 5]))
ax.loglog(curve[:, 0], abs(curve[:, 6]))
ax.loglog(curve[:, 0], abs(curve[:, 7]))
ax.loglog(curve[:, 0], abs(curve[:, 8]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('k (h/Mpc)', fontsize=16)
plt.show()