import numpy as np
import matplotlib.pyplot as plt
import illustris_python as il
import time

parttype_stars = 4

suite = 'IllustrisTNG'
simulation = 'L75n1820TNG'  # TNG100
snapnum = 99
basedir = '/fs/san11/vicente/simulations/%s/%s/output' % (suite, simulation)

header = il.groupcat.loadHeader(basedir, snapnum)
h = header['HubbleParam']
box_size = header['BoxSize']  # ckpc/h
volume_cMpc = (box_size / 1000.0 / h)**3

# Read subhalo info
mstar_in_rad = il.groupcat.loadSubhalos(
    basedir, snapnum, fields=['SubhaloMassInRadType'])[:, parttype_stars]

fig = plt.figure()
ax = fig.add_subplot(111)

# Stellar mass bins (10^8 - 10^12 Msun)
nbins = 10
mstar_min = 1e8 / 1e10 * h  # 10^10 Msun/h
mstar_max = 1e12 / 1e10 * h  # 10^10 Msun/h
mstar_bin_edges = np.logspace(np.log10(mstar_min), np.log10(mstar_max), nbins+1)
mstar_bin_centers = np.sqrt(mstar_bin_edges[:-1] * mstar_bin_edges[1:])  # geometric mean
d_mstar = mstar_bin_edges[1:] - mstar_bin_edges[:-1]  # np.diff()
d_log_mstar = np.log10(mstar_bin_edges[1:]) - np.log10(mstar_bin_edges[:-1])  # np.diff(np.log10())

hist, mstar_bin_edges = np.histogram(mstar_in_rad, bins=mstar_bin_edges)
Phi = hist / volume_cMpc / d_log_mstar
Phi_approx = np.log(10) * mstar_bin_centers * hist / volume_cMpc / d_mstar

ax.plot(mstar_bin_centers*1e10/h, Phi, label='Exact')
ax.plot(mstar_bin_centers*1e10/h, Phi_approx, label='Approx.')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Stellar Mass, ' + r'$M_{\ast}^{< r_{\ast}} \, [{\rm M}_{\odot}]$', fontsize=14)
ax.set_ylabel('Stellar Mass Function, ' + r'$\Phi \, [{\rm Mpc^{-3} dex^{-1}}]$', fontsize=14)
ax.legend(loc=3, fontsize=14, frameon=False)
fig.subplots_adjust(left=0.125, right=0.95, bottom=0.125, top=0.95)  # defaults: 0.125, 0.9, 0.1, 0.9

#plt.show()
fig.savefig('stellar_mass_function.png')
plt.close(fig)
