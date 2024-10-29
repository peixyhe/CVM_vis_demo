import vtk



def write_obj(filename, vertices, triangles):
    """Write the mesh data to an OBJ file."""
    with open(filename, 'w') as f:
        for v in vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for t in triangles:
            f.write(f"f {int(t[0]) + 1} {int(t[1]) + 1} {int(t[2]) + 1}\n")


def vtk_append_files(vtk_file, append_filter):
    """Read a VTK file and add its data to the append filter."""
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(vtk_file)
    reader.Update()
    append_filter.AddInputData(reader.GetOutput())


def visualize_vtk_file(filename):
    """
    Visualize a VTK file using VTK's rendering pipeline.
    
    Args:
        filename (str): Path to the VTK file to visualize.
    
    Returns:
        None
    """
    # Create a reader for the VTK file
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()

    # Create a mapper
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputConnection(reader.GetOutputPort())

    # Create an actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create a renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1)  # Set background to white

    # Create a render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(800, 600)

    # Create a render window interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Initialize and start the interaction
    render_window.Render()
    interactor.Start()
