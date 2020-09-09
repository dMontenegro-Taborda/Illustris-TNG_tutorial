import numpy as np
import matplotlib.pyplot as plt
import illustris_python as il
import time
import h5py
import sys

def plot_alignment(ax, group_M200, group_jstar_central, linecolor='k', linelabel=None):
    # Apply eq. from slide 6
    numerator = np.sum(group_j200 * group_jstar_central, axis=1)
    denominator = np.sqrt(np.sum(group_j200**2, axis=1) * np.sum(group_jstar_central**2, axis=1))
    group_cos_theta = np.zeros(ngroups, dtype=np.float32)
    locs_valid = denominator > 0
    group_cos_theta[locs_valid] = numerator[locs_valid] / denominator[locs_valid]

    # Calculate median and 1-sigma scatter
    y_median = np.zeros(nbins)
    y_16 = np.zeros(nbins)
    y_84 = np.zeros(nbins)
    for i in range(nbins):
        locs = (group_M200 >= M200_bin_edges[i]) & (group_M200 < M200_bin_edges[i+1])
        y_median[i] = np.median(group_cos_theta[locs])
        y_16[i] = np.percentile(group_cos_theta[locs], 16.0)
        y_84[i] = np.percentile(group_cos_theta[locs], 84.0)

    ax.plot(M200_bin_centers*1e10/h, y_median, color=linecolor, label=linelabel)
    ax.fill_between(M200_bin_centers*1e10/h, y_16, y_84, color=linecolor, alpha=0.2)


if __name__ == '__main__':
    parttype_stars = 4
    h = 0.6774  # should read from header

    suite = 'IllustrisTNG'
    simulation = 'L75n1820TNG'
    snapnum = 99
    basedir = '/fs/san11/vicente/simulations/%s/%s/output' % (suite, simulation)

    # Halo mass bins (M200 = 10^10 - 10^15 Msun)
    nbins = 20
    M200_min = 1e10 / 1e10 * h  # 10^10 Msun/h
    M200_max = 5e14 / 1e10 * h  # 10^10 Msun/h
    M200_bin_edges = np.logspace(np.log10(M200_min), np.log10(M200_max), nbins+1)
    M200_bin_centers = np.sqrt(M200_bin_edges[:-1] * M200_bin_edges[1:])  # geometric mean

    # Read halo info
    print('Reading halo info...')
    start = time.time()
    group_first_sub = il.groupcat.loadHalos(basedir, snapnum, fields=['GroupFirstSub'])
    with h5py.File('/fs/san11/vicente/SpinParameter/IllustrisTNG/L75n1820TNG/halo_spin_1.0_r200_099.hdf5', 'r') as f:
        group_J200 = f['GroupAngMom'][:]
        group_M200 = f['GroupMass'][:]
    print('Time: %.2f s.' % (time.time() - start))
    ngroups = len(group_M200)

    # Create figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    colors = ['r', 'b']
    names = ['jstar', 'jstar_1.0_rhalfs']
    for k in range(len(names)):
        # Read subhalo info
        print('Reading subhalo info...')
        start = time.time()
        filename = '/fs/san11/vicente/AngularMomentum/output/%s/%s/%s_%03d.hdf5' % (
            suite, simulation, names[k], snapnum)
        with h5py.File(filename, 'r') as f:
            jstar_magnitude = f['jstar'][:]
            jstar_direction = f['jstar_direction'][:]
        jstar = jstar_direction * jstar_magnitude[:, np.newaxis]  # .reshape((len(jstar_magnitude), 1))
        print('Time: %.2f s.' % (time.time() - start))

        # Calculate specific halo AM. Avoid division by zero.
        group_j200 = np.zeros_like(group_J200)
        locs_positive = group_M200 > 0
        group_j200[locs_positive] = group_J200[locs_positive] / group_M200[locs_positive, np.newaxis]
        group_jstar_central = jstar[group_first_sub, :]

        plot_alignment(ax, group_M200, group_jstar_central, 
                       linecolor=colors[k], linelabel=names[k])
        print('Finished for', names[k])

    ax.plot(np.array([M200_bin_centers[0], M200_bin_centers[-1]])*1e10/h, [0, 0], 'k--')
    ax.set_xscale('log')
    ax.set_xlabel('Halo Mass, ' + r'$M_{200,{\rm c}} \, [{\rm M}_{\odot}]$', fontsize=14)
    ax.set_ylabel(r'$\cos(\theta)$', fontsize=14)
    ax.legend(loc=3, fontsize=14)
    fig.subplots_adjust(left=0.14, right=0.95, bottom=0.125, top=0.95)  # defaults: 0.125, 0.9, 0.1, 0.9
    fig.savefig('galaxy_halo_alignment.png')
    plt.close(fig)

