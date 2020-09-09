import numpy as np
import illustris_python as il
import time
import h5py
import sys

try:
    mode = sys.argv[1]  # 'write', 'read'
except:
    print('Arguments: mode')
    sys.exit()

parttype_stars = 4
h = 0.6774  # should read from header
box_size = 75000.0  # ckpc/h
M200_min = 1e12 / 1e10 * h  # 10^10 Msun/h

suite = 'IllustrisTNG'
simulation = 'L75n1820TNG'
snapnum = 99
basedir = '/fs/san11/vicente/simulations/%s/%s/output' % (suite, simulation)

if mode == 'write':
    # Read subhalo info
    print('Reading subhalo info...')
    start = time.time()
    with h5py.File('/fs/san11/vicente/AngularMomentum/output/IllustrisTNG/L75n1820TNG/jstar_099.hdf5', 'r') as f:
        jstar = f['jstar'][:]
    print('Time: %.2f s.' % (time.time() - start))

    # Read halo info
    print('Reading halo info...')
    start = time.time()
    group_first_sub = il.groupcat.loadHalos(basedir, snapnum, fields=['GroupFirstSub'])
    group_pos = il.groupcat.loadHalos(basedir, snapnum, fields=['GroupPos'])
    with h5py.File('/fs/san11/vicente/SpinParameter/IllustrisTNG/L75n1820TNG/halo_spin_1.0_r200_099.hdf5', 'r') as f:
        group_J200 = f['GroupAngMom'][:]
        group_M200 = f['GroupMass'][:]
    print('Time: %.2f s.' % (time.time() - start))

    # Calculate specific halo AM. Avoid division by zero.
    group_j200 = np.zeros_like(group_J200)
    locs_positive = group_M200 > 0
    group_j200[locs_positive] = group_J200[locs_positive] / group_M200[locs_positive, np.newaxis]
    group_jstar_central = jstar[group_first_sub]

    # Impose minimum halo mass
    locs_massive = group_M200 >= M200_min

    # Save to HDF5 file
    with h5py.File('data.hdf5', 'w') as f:
        f.create_dataset('group_pos', data=group_pos[locs_massive])
        f.create_dataset('group_j200', data=group_j200[locs_massive])
        f.create_dataset('group_jstar_central', data=group_jstar_central[locs_massive])

elif mode == 'read':
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    with h5py.File('data.hdf5', 'r') as f:
        group_pos = f['group_pos'][:]
        group_j200 = f['group_j200'][:]
        group_jstar_central = f['group_jstar_central'][:]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax = fig.gca(projection='3d')
    
    x = group_pos[:,0]; y = group_pos[:,1]; z = group_pos[:,2]
    maximun = max( max(group_j200[:,0]),max(group_j200[:,1]),max(group_j200[:,2]))
    #u = group_j200[:,0]/maximun; v = group_j200[:,1]/maximun; w = group_j200[:,2]/maximun
    u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
    w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) * np.sin(np.pi * z))
    
    
    
    ax.quiver(x, y, z, u, v, w, length = 10 , arrow_length_ratio = 0.1)
    ax.scatter(x, y, z, s=0.01)
    ax.set_xlim(0, box_size); ax.set_ylim(0, box_size); ax.set_zlim(0, box_size)
    ax.set_xlabel('x [ckpc/h]'); ax.set_xlabel('y [ckpc/h]'); ax.set_xlabel('z [ckpc/h]')
    
    plt.show()
    plt.close(fig)
