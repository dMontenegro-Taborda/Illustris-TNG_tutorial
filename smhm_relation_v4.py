import numpy as np
import matplotlib.pyplot as plt
import illustris_python as il
import time

parttype_stars = 4

suite = 'IllustrisTNG'
simulation = 'L75n1820TNG'
snapnum = 99
basedir = '/fs/san11/vicente/simulations/%s/%s/output' % (suite, simulation)

header = il.groupcat.loadHeader(basedir, snapnum)
h = header['HubbleParam']
Omega_m = header['Omega0']
Omega_b = 0.0486  # Pillepich et al. (2018)

# Halo mass bins (M200 = 10^10 - 10^15 Msun)
nbins = 20
m200_min = 1e10 / 1e10 * h  # 10^10 Msun/h
m200_max = 5e14 / 1e10 * h  # 10^10 Msun/h
m200_bin_edges = np.logspace(np.log10(m200_min), np.log10(m200_max), nbins+1)
m200_bin_centers = np.sqrt(m200_bin_edges[:-1] * m200_bin_edges[1:])  # geometric mean

# Read subhalo info
print('Reading subhalo info...')
start = time.time()
mstar_in_rad = il.groupcat.loadSubhalos(basedir, snapnum, fields=['SubhaloMassInRadType'])[:, parttype_stars]
print('Time: %.2f s.' % (time.time() - start))

# Read halo info
print('Reading halo info...')
start = time.time()
cat_groups = il.groupcat.loadHalos(basedir, snapnum, fields=['GroupFirstSub', 'Group_M_Crit200'])
group_first_sub = cat_groups['GroupFirstSub']
group_m200 = cat_groups['Group_M_Crit200']
print('Time: %.2f s.' % (time.time() - start))

# Map halos to their central subhalos
group_mstar_central = mstar_in_rad[group_first_sub]

# Get halos with at least one subhalo
locs_valid = group_first_sub != -1
group_m200 = group_m200[locs_valid]
group_mstar_central = group_mstar_central[locs_valid]

# Avoid division by zero
shm_ratio = -1.0 * np.ones_like(group_m200)
locs_positive = group_m200 > 0
shm_ratio[locs_positive] = group_mstar_central[locs_positive] / group_m200[locs_positive]
baryon_conv_eff = shm_ratio / (Omega_b/Omega_m)

# Calculate median and 1-sigma scatter
y_median = np.zeros(nbins)
y_16 = np.zeros(nbins)
y_84 = np.zeros(nbins)
for i in range(nbins):
    locs = (group_m200 >= m200_bin_edges[i]) & (group_m200 < m200_bin_edges[i+1])
    y_median[i] = np.median(baryon_conv_eff[locs])
    y_16[i] = np.percentile(baryon_conv_eff[locs], 16.0)
    y_84[i] = np.percentile(baryon_conv_eff[locs], 84.0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(m200_bin_centers*1e10/h, y_median, color='b')
ax.fill_between(m200_bin_centers*1e10/h, y_16, y_84, color='b', alpha=0.2)

ax.set_xscale('log')
ax.set_xlabel('Halo Mass, ' + r'$M_{200,{\rm c}} \, [{\rm M}_{\odot}]$', fontsize=14)
ax.set_ylabel(r'$M_{\ast} / M_{200, {\rm c}} {\rm (\Omega_b / \Omega_m)^{-1}}$', fontsize=14)
fig.subplots_adjust(left=0.14, right=0.95, bottom=0.125, top=0.95)  # defaults: 0.125, 0.9, 0.1, 0.9
fig.savefig('smhm_relation.png')
plt.close(fig)
