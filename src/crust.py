import sys
import os
import vtk
import pandas as pd

# Add the user-defined module path
sys.path.append(r"src")
import CONST  # Import constants from user-defined module



def main():
    """
    Main function to read crustal thickness data, process it, and visualize using VTK.
    """
    input_rawData_path = r"rawData\Crustal_Thickness_of_Chinese_Mainland_2014.txt"
    
    output_path = r'resultData\crust'
    os.makedirs(output_path, exist_ok=True)
    
    # Initialize VTK components
    points, depth_array = initialize_vtk_components()

    # Read and process data
    x_set, y_set = read_and_process_data(points, depth_array, input_rawData_path)

    # Build the grid cells
    ugrid = build_grid(points, depth_array, x_set, y_set)

    # Write the unstructured grid to a VTK file
    write_vtk_file(  ugrid, output_vtk_path=os.path.join(output_path, f'crust_X{CONST.Z_SCALE}.vtu')  )

    # Visualize the unstructured grid
    visualize(ugrid)



def initialize_vtk_components():
    """
    Initialize VTK points and depth array.
    """
    points = vtk.vtkPoints()
    depth_array = vtk.vtkFloatArray()
    depth_array.SetNumberOfComponents(1)
    depth_array.SetName('crustalthickness(km)')
    return points, depth_array

def read_and_process_data(points, depth_array, data_path):
    """
    Read crustal thickness data and process it.
    """
    # Initialize sets to store unique x and y values
    x_set = set()
    y_set = set()

    # Read data from the raw text file and filter based on coordinates
    with open(data_path, "r") as file:
        lines_data = file.readlines()

    for line in lines_data[2:]:
        line_data = line.strip().split()
        lon, lat, crust_thickness = map(float, line_data)
        if CONST.longitude_left <= lon <= CONST.longitude_right and CONST.latitude_down <= lat <= CONST.latitude_up:
            points.InsertNextPoint(lon, lat, -crust_thickness / CONST.angle_to_kilometers)
            depth_array.InsertNextValue(crust_thickness)
            x_set.add(lon)
            y_set.add(lat)

    return x_set, y_set

def build_grid(points, depth_array, x_set, y_set):
    """
    Build the grid cells and create an unstructured grid.
    """
    # Create an unstructured grid
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(depth_array)

    # Sort the x and y coordinates for consistent indexing
    x_list = sorted(x_set)
    y_list = sorted(y_set)
    x_num = len(x_list)
    y_num = len(y_list)

    # Create a mapping from coordinate to index
    coord_to_id = {}
    point_id = 0
    for y in y_list:
        for x in x_list:
            coord_to_id[(x, y)] = point_id
            point_id += 1

    # Build the grid cells (triangles)
    for j in range(y_num - 1):
        for i in range(x_num - 1):
            # Get the point IDs for the corners of the grid cell
            id0 = coord_to_id[(x_list[i], y_list[j])]
            id1 = coord_to_id[(x_list[i + 1], y_list[j])]
            id2 = coord_to_id[(x_list[i + 1], y_list[j + 1])]
            id3 = coord_to_id[(x_list[i], y_list[j + 1])]
            """
               id3------id2
                |        |
                |        |
               id0------id1
            """
            # Insert two triangles to form a quadrilateral
            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [id0, id1, id2])
            ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, [id0, id2, id3])

    return ugrid

def write_vtk_file(ugrid, output_vtk_path):
    """
    Write the unstructured grid to a VTK file.
    """
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(output_vtk_path)
    # writer.SetDataModeToAscii()  # Optional: write data in ASCII format
    writer.Write()

def visualize(ugrid):
    """
    Visualize the unstructured grid using VTK.
    """
    colors = vtk.vtkNamedColors()

    # Create a mapper and actor for visualization
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(ugrid)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d('Silver'))
    actor.GetProperty().SetPointSize(2)

    # Set up the renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.SetWindowName('Crustal Thickness Visualization')
    render_window.AddRenderer(renderer)
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Add the actor to the scene
    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('Salmon'))
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(30)
    renderer.GetActiveCamera().Elevation(30)

    # Render and start interaction
    render_window.Render()
    render_window_interactor.Start()

if __name__ == "__main__":
    main()
