import numpy as np
import illustris_python as il
import time
import h5py
import sys
#import pyvista

def add_cube(renderer, a):
    """Show cube between (0,0,0) and (a,a,a)
    """
    # Create points
    p_data = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, +a],
            [0.0, +a, 0.0],
            [0.0, +a, +a],
            [+a, 0.0, 0.0],
            [+a, 0.0, +a],
            [+a, +a, 0.0],
            [+a, +a, +a]])

    edges = [(0,1), (1,3), (3,2), (2,0), (4,5), (5,7), (7,6), (6,4), (0,4), (1,5), (3,7), (2,6)]

    npoints = len(p_data)

    # Create a vtkPoints object and store the points in it
    points = vtk.vtkPoints()
    for i in range(npoints):
        points.InsertNextPoint(p_data[i])
     
    # Create a cell array to store the lines in and add the lines to it
    lines = vtk.vtkCellArray()
    for i in range(len(edges)):
      line = vtk.vtkLine()
      line.GetPointIds().SetId(0, edges[i][0])
      line.GetPointIds().SetId(1, edges[i][1])
      lines.InsertNextCell(line)
     
    # Create a polydata to store everything in
    linesPolyData = vtk.vtkPolyData()
     
    # Add the points to the dataset
    linesPolyData.SetPoints(points)
     
    # Add the lines to the dataset
    linesPolyData.SetLines(lines)
     
    # Setup actor and mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(linesPolyData)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    renderer.AddActor(actor)

def add_arrows():
    
    #position, direction
    # First set of points - first arrow
    points1 = np.array([[ 11271.915, 7538.686, 6245.661],
                        [ 11271.915, 7538.686, 5897.034]])

    # Second set of points - second arrow
    points2 = np.array([[ 10532.274, 9101.572, 6313.167],
                        [ 10532.274, 9101.572, 5964.539]])
    """
    points1 = position
    points2 = direction
    """
    # Create an arrow source with some attributes
    arrow_source = vtk.vtkArrowSource()
    #arrow_source.SetTipRadius(1)
    #arrow_source.SetShaftRadius(5)
    
    arrow_source.SetShaftRadius(1000)
    arrow_source.SetTipLength(9)
    
    # Create the vtkGlyph3D
    arrow_glyph = vtk.vtkGlyph3D()
    arrow_glyph.SetScaleModeToDataScalingOff()

    # Usual mapper
    arrow_mapper = vtk.vtkPolyDataMapper()
    arrow_mapper.SetInputConnection(arrow_glyph.GetOutputPort())

    # Usual actor
    arrow_actor = vtk.vtkActor()
    arrow_actor.SetMapper(arrow_mapper)

    append_poly_data = vtk.vtkAppendPolyData()

    # Loop through points1 and points2
    for p in [points1, points2]:
        
        arrow_poly_data = vtk.vtkPolyData()
        vtk_arrow_lines = vtk.vtkCellArray()

        vtk_arrow_lines.InsertNextCell(2)
        vtk_arrow_lines.InsertCellPoint(0)
        vtk_arrow_lines.InsertCellPoint(1)

        vtk_arrow_points = vtk.vtkPoints()
        # Loop through the head and tail of a single arrow
        for j in range(2):
            
            vtk_arrow_points.InsertNextPoint(p[j])

            arrow_poly_data.SetPoints(vtk_arrow_points)
            arrow_poly_data.SetLines(vtk_arrow_lines)
            append_poly_data.AddInputData(arrow_poly_data)


    arrow_glyph.SetInputConnection(append_poly_data.GetOutputPort())
    arrow_glyph.SetSourceConnection(arrow_source.GetOutputPort())
    arrow_glyph.SetScaleFactor(100)
    arrow_actor.GetProperty().SetColor(1, 1, 1)
    """
    # ------------------------------------------------------------
    # Create the RenderWindow, Renderer and both Actors
    # ------------------------------------------------------------
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # add actors
    ren.AddActor(arrow_actor)
    #ren.ResetCamera()
    #iren.Start()
    """
    # Visualize
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(arrow_poly_data)
     
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(1)

    # Set color of actor
    #actor.GetProperty().SetColor(color)
    
    renderer.AddActor(actor)


def add_points(renderer, pos, color):
    # Create the geometry of a point (the coordinate)
    points = vtk.vtkPoints()

    # Create the topology of the point (a vertex)
    vertices = vtk.vtkCellArray()

    # Add points one by one (don't know another way)
    for i in range(len(pos)):
        pid = points.InsertNextPoint(pos[i])
        vertices.InsertNextCell(1)
        vertices.InsertCellPoint(pid)

    # Create a polydata object
    points_poly_data = vtk.vtkPolyData()

    # Set the points and vertices we created as the geometry and topology of the polydata
    points_poly_data.SetPoints(points)
    points_poly_data.SetVerts(vertices)

    # Visualize
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(points_poly_data)
     
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(1)

    # Set color of actor
    actor.GetProperty().SetColor(color)
    
    renderer.AddActor(actor)

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
    import vtk

    with h5py.File('data.hdf5', 'r') as f:
        group_M200 = f['group_M200'][:]
        group_pos = f['group_pos'][:]
        group_j200 = f['group_j200'][:]
        group_jstar_central = f['group_jstar_central'][:]

    x = group_pos[:,0]; y = group_pos[:,1]; z = group_pos[:,2]

    
    # Create renderer
    renderer = vtk.vtkRenderer()
    
    #create arrows
    arrowSource = vtk.vtkArrowSource()
    

    # Draw box
    add_cube(renderer, box_size)
    # Add points
    add_points(renderer, group_pos, (1, 1, 1))
    # Add arrows
    #add_arrows(group_pos, group_j200*-1. )
    add_arrows()
 
    # Create render window and render window interactor
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    
    # Render
    renderWindow.Render()
    renderWindowInteractor.Start()
    
   

