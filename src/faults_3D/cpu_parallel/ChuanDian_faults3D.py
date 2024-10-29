import pymeshlab
import os
import vtk
import pyproj
import numpy as np
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

# Append custom library paths
sys.path.append(r"src")
import CONST

sys.path.append(r'src\faults_3D\lib')
from vis_func import (
    write_obj,
    vtk_append_files,
    visualize_vtk_file
)



# Define input and output paths
INPUT_OBJ_FILEPATH = r'resultData\faults_3D\Final_modeling_faults\remesh_scale_faults_obj'
OUTPUT_PATH = os.path.join(r'resultData\faults_3D', f'ChuanDian_faults3D_X{CONST.Z_SCALE}')

OUTPUT_VTK_PATH = os.path.join(OUTPUT_PATH, f'ChuanDian_faults3D_vtk_scale{CONST.Z_SCALE}')
OUTPUT_OBJ_PATH = os.path.join(OUTPUT_PATH, f'ChuanDian_faults3D_obj_scale{CONST.Z_SCALE}')

# Create output directories if they do not exist
os.makedirs(OUTPUT_VTK_PATH, exist_ok=True)
os.makedirs(OUTPUT_OBJ_PATH, exist_ok=True)

# Initialize UTM projection (transform UTM-xyz to longitude and latitude)
utm_proj = pyproj.Proj(proj='utm', zone=48, ellps='WGS84')

# List of faults to skip
abandon_faults = [
    'bianba-luolong fault', 'Awancang_Fault-2', 'Hanan_fault_S_2', 'Lushan_EQ_fault_3',
    'Maduo-Gande_Fault-2-1', 'Tazang_fault_6', 'Bailongjiang_fault_S', 'Tazang_fault_1',
    'baqing-leiniaoqiao fault', 'Bailongjiang_fault_N_4'
]



def write_vtk(filename, vertices, triangles, depth_attribute):
    """
    Write the mesh data to a VTK file.
    
    Parameters:
        filename (str): Path to the output VTK file.
        vertices (numpy.ndarray): Array of vertex coordinates.
        triangles (numpy.ndarray): Array of triangle indices.
        depth_attribute (numpy.ndarray): Array of depth values.
    """
    # Create an unstructured grid
    ugrid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()

    # Create a scalar array for depth
    depth = vtk.vtkFloatArray()
    depth.SetNumberOfComponents(1)
    depth.SetName('fault_depth_km')
    
    # Add vertices and depth attribute
    for v, dep in zip(vertices, depth_attribute):
        points.InsertNextPoint(v)
        depth.InsertNextValue(round(-1.0 * dep * CONST.angle_to_kilometers, 3))
    
    # Add triangle cells
    for t in triangles:
        ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, t)
    
    # Assign points and depth data to the grid
    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(depth)
    
    # Write to VTK file
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(filename)
    writer.Write()


def process_fault(fault):
    """
    Process a single fault file: load, transform, sample, reconstruct, repair, and save.
    
    Parameters:
        fault (str): Filename of the fault OBJ file.
    
    Returns:
        str or None: Path to the generated VTK file if successful, else None.
    """
    try:
        fault_name = os.path.splitext(fault)[0]
        print(f"Processing fault:  {fault} ...")

        # Skip abandoned faults
        if fault_name in abandon_faults:
            print(f"Skipping abandoned fault:  {fault_name}")
            return None

        fault_file = os.path.join(INPUT_OBJ_FILEPATH, fault)
        
        # Create MeshSet and load mesh
        ms = pymeshlab.MeshSet()
        ms.load_new_mesh(fault_file)
        
        vertices = ms.current_mesh().vertex_matrix()
        cells = ms.current_mesh().face_matrix()
        
        # Convert UTM coordinates to longitude and latitude, and adjust depth
        lon_lat_dep = []
        for xyz in vertices:
            lon, lat = utm_proj(xyz[0], xyz[1], inverse=True)
            lon_lat_dep.append([lon, lat, xyz[2] / CONST.angle_to_meters])
        lon_lat_dep = np.round(np.array(lon_lat_dep), 6)
        
        # Create a new MeshSet with converted coordinates
        ms1 = pymeshlab.MeshSet()
        ms1.add_mesh(pymeshlab.Mesh(vertex_matrix=lon_lat_dep, face_matrix=cells))
        
        # Perform initial Poisson disk sampling
        samplenum_initial = max(1, int(len(vertices) * CONST.Z_SCALE))
        ms1.generate_sampling_poisson_disk(samplenum=samplenum_initial)
        vertices1 = ms1.current_mesh().vertex_matrix()
        
        # Filter points within the specified longitude and latitude boundaries
        points_xyz = []
        for ver in vertices1:
            if (CONST.longitude_left <= ver[0] <= CONST.longitude_right) and (CONST.latitude_down <= ver[1] <= CONST.latitude_up):
                points_xyz.append(ver)
        
        if len(points_xyz) == 0:
            print(f"No points within specified boundaries for fault:  {fault_name}")
            return None
        
        # Create another MeshSet with filtered points
        ms2 = pymeshlab.MeshSet()
        ms2.add_mesh(pymeshlab.Mesh(vertex_matrix=np.array(points_xyz)))
        
        # Perform surface reconstruction using Ball Pivoting
        if fault_name in ['Anninghe_Zemuhe_Xiaojiang_Faults']:
            ms2.generate_surface_reconstruction_ball_pivoting()
        else:
            for i in range(3):
                ms2.generate_surface_reconstruction_ball_pivoting(ballradius=pymeshlab.PercentageValue(2 * i))
        
        # Perform Poisson disk sampling again
        samplenum_final = max(1, int(len(points_xyz) / 50.0))
        ms2.generate_sampling_poisson_disk(samplenum=samplenum_final)
        
        # Perform surface reconstruction again
        if fault_name in ['Anninghe_Zemuhe_Xiaojiang_Faults']:
            ms2.generate_surface_reconstruction_ball_pivoting()
        else:
            for i in range(3):
                ms2.generate_surface_reconstruction_ball_pivoting(ballradius=pymeshlab.PercentageValue(2 * i))
        
        # Mesh repair: remove null faces, duplicate faces, duplicate vertices, and unreferenced vertices
        ms2.meshing_remove_null_faces()
        ms2.meshing_remove_duplicate_faces()
        ms2.meshing_remove_duplicate_vertices()
        ms2.meshing_remove_unreferenced_vertices()
        
        # Retrieve processed vertices and faces
        vertices2 = ms2.current_mesh().vertex_matrix()
        faces2 = ms2.current_mesh().face_matrix()
        
        # Adjust depth
        vertices2[:, 2] *= (-1.0)
        
        # Write OBJ file
        out_obj_file_path = os.path.join(OUTPUT_OBJ_PATH, f"{fault_name}.obj")
        write_obj(out_obj_file_path, vertices2, faces2)
        
        # Write VTK file
        out_vtk_file_path = os.path.join(OUTPUT_VTK_PATH, f"{fault_name}.vtk")
        write_vtk(out_vtk_file_path, vertices2, faces2, vertices2[:, -1])
        
        return out_vtk_file_path  # Return the path for later appending
    except Exception as e:
        print(f"Error processing fault {fault}:  {e}")
        return None


def main():
    """
    Main function to parallelize fault processing and combine VTK files.
    """
    # List all OBJ files in the input directory
    faults = [f for f in os.listdir(INPUT_OBJ_FILEPATH) if f.lower().endswith('.obj')]
    
    vtk_file_paths = []
    
    # Determine the number of parallel workers (half the CPU cores to balance performance and memory usage)
    cpu_cores = os.cpu_count() or 1
    max_workers = max(1, cpu_cores // 2)
    
    print(f"Number of parallel processes:  {max_workers}")
    
    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all fault processing tasks
        future_to_fault = {executor.submit(process_fault, fault): fault for fault in faults}
        
        # Iterate over completed tasks as they finish
        for future in as_completed(future_to_fault):
            fault = future_to_fault[future]
            try:
                vtk_file = future.result()
                if vtk_file:
                    vtk_file_paths.append(vtk_file)
            except Exception as exc:
                print(f"Fault {fault} generated an exception:  {exc}")
    
    if not vtk_file_paths:
        print("No VTK files were generated.")
        return
    
    # Create a vtkAppendFilter object to combine multiple datasets
    append_filter = vtk.vtkAppendFilter()
    
    # Append all VTK files to the filter
    for vtk_file in vtk_file_paths:
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtk_file)
        reader.Update()
        append_filter.AddInputData(reader.GetOutput())
    
    # Update the filter to perform the append operation
    append_filter.Update()
    
    # Write the combined dataset to a new VTK file
    combined_vtk_file = os.path.join(OUTPUT_PATH, f"ChuanDian_faults_3D_X{CONST.Z_SCALE}.vtk")
    writer2 = vtk.vtkDataSetWriter()
    writer2.SetFileName(combined_vtk_file)
    writer2.SetInputData(append_filter.GetOutput())
    writer2.Write()
    
    print(f"Combined VTK file saved to: {combined_vtk_file}")
    
    # Visualize the combined VTK file
    visualize_vtk_file(combined_vtk_file)


if __name__ == "__main__":
    main()
