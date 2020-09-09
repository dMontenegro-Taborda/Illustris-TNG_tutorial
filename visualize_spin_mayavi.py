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
M200_min = 1e11 / 1e10 * h  # 10^10 Msun/h

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
        f.create_dataset('group_M200', data=group_M200[locs_massive])
        f.create_dataset('group_pos', data=group_pos[locs_massive])
        f.create_dataset('group_j200', data=group_j200[locs_massive])
        f.create_dataset('group_jstar_central', data=group_jstar_central[locs_massive])

elif mode == 'read':
    from mayavi import mlab

    with h5py.File('data.hdf5', 'r') as f:
        group_M200 = f['group_M200'][:]
        group_pos = f['group_pos'][:] 
        group_j200 = f['group_j200'][:]
        group_jstar_central = f['group_jstar_central'][:]

    x = group_pos[:,0]; y = group_pos[:,1]; z = group_pos[:,2]
    
    u = group_j200[:,0]/max(group_j200[:,0]); v = group_j200[:,1]/max(group_j200[:,1]); w = group_j200[:,2]/max(group_j200[:,2])
    fig = mlab.figure(size=(500,500), fgcolor=(0.,0.,0.), bgcolor=(1,1,1))
    s = (np.log10(group_M200/np.log10(max(group_M200)))) * 1
    mlab.points3d(x, y, z, s , scale_factor = 200 , mode='sphere', scale_mode = 'scalar', figure=fig)
    mlab.quiver3d(x, y, z, u, v, w, color = (0.2, 0.1, 0.2), line_width = 4.5, scale_factor = 450, scale_mode = 'scalar')
    box_data = np.array([[0, 0, 0],
                         [0, 0, box_size],
                         [0, box_size, 0],
                         [0, box_size, box_size],
                         [box_size, 0, 0],
                         [box_size, 0, box_size],
                         [box_size, box_size, 0],
                         [box_size, box_size, box_size]])
    mlab.plot3d(box_data[:,0], box_data[:,1], box_data[:,2])
