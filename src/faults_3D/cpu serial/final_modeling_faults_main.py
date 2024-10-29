import pymeshlab
import os
import vtk
import numpy as np
from scipy.spatial import cKDTree
import random
import sys

# Append custom library path
sys.path.append(r'src\faults_3D\lib')
from vis_func import (
    write_obj,
    vtk_append_files,
    visualize_vtk_file
)



# Define input and output directories
input_OBJ_filePath = r'resultData\faults_3D\remesh_faults\remesh_faults_obj'

output_path = r'resultData\faults_3D\Final_modeling_faults'
result_VTK_filePath = os.path.join(output_path, r'remesh_scale_faults_vtk')
result_OBJ_filePath = os.path.join(output_path, r'remesh_scale_faults_obj')

# Create output directories if they do not exist
os.makedirs(result_OBJ_filePath, exist_ok=True)
os.makedirs(result_VTK_filePath, exist_ok=True)

# Scaling factor for the Z-axis
scale = 25.0



def main():
    # Create a vtkAppendFilter object to combine multiple datasets
    append_filter = vtk.vtkAppendFilter()

    # Iterate through each OBJ file in the input directory
    for fault in os.listdir(input_OBJ_filePath):
        fault_name = os.path.splitext(fault)[0]
        print(f"Processing fault: {fault}")
        fault_file = os.path.join(input_OBJ_filePath, fault)
        
        # Load the mesh using PyMeshLab
        ms = pymeshlab.MeshSet()               
        ms.load_new_mesh(fault_file)
        
        vertices = ms.current_mesh().vertex_matrix()
        cells = ms.current_mesh().face_matrix()
        
        # Scale the Z-axis of the vertices
        scaled_vertices = np.array(vertices)
        scaled_vertices[:, 2] *= scale

        # Create a new MeshSet with scaled vertices
        ms_scaled = pymeshlab.MeshSet()
        ms_scaled.add_mesh(pymeshlab.Mesh(vertex_matrix=scaled_vertices, face_matrix=cells))
        
        # Generate Poisson disk sampling
        ms_scaled.generate_sampling_poisson_disk(samplenum=int(len(vertices) * scale))
        
        # Perform surface reconstruction using Ball Pivoting
        for i in range(3):
            ms_scaled.generate_surface_reconstruction_ball_pivoting(ballradius=pymeshlab.PercentageValue(2 * i))
        
        # Apply scale-dependent Laplacian smoothing
        ms_scaled.apply_coord_laplacian_smoothing_scale_dependent()
        try:
            # Apply two-step and Taubin smoothing
            ms_scaled.apply_coord_two_steps_smoothing()
            ms_scaled.apply_coord_taubin_smoothing()
        except:
            # Fallback to Taubin smoothing if two-step smoothing fails
            ms_scaled.apply_coord_taubin_smoothing()
        
        # Mesh Repair Step 1
        ms_scaled.meshing_remove_null_faces()               # Remove faces with zero area
        ms_scaled.meshing_remove_duplicate_faces()          # Remove duplicate faces
        ms_scaled.meshing_remove_duplicate_vertices()       # Remove duplicate vertices   
        ms_scaled.meshing_remove_unreferenced_vertices()    # Remove unreferenced vertices
        
        # Retrieve processed vertices and faces
        processed_vertices = ms_scaled.current_mesh().vertex_matrix()
        processed_faces = ms_scaled.current_mesh().face_matrix()
        
        # Scale the Z-axis back to original
        scaled_back_vertices = np.array(processed_vertices)
        scaled_back_vertices[:, 2] /= scale
        
        # Create another MeshSet with scaled-back vertices
        ms_final = pymeshlab.MeshSet()
        ms_final.add_mesh(pymeshlab.Mesh(vertex_matrix=scaled_back_vertices, face_matrix=processed_faces))
        
        # Generate Poisson disk sampling with original number of vertices
        ms_final.generate_sampling_poisson_disk(samplenum=len(vertices))
        
        # Perform surface reconstruction using Ball Pivoting
        for j in range(3):
            ms_final.generate_surface_reconstruction_ball_pivoting(ballradius=pymeshlab.PercentageValue(2 * j))
        
        try:
            # Apply two-step smoothing
            ms_final.apply_coord_two_steps_smoothing()
            # ms_final.apply_coord_taubin_smoothing()
        except:
            # Fallback to scale-dependent Laplacian smoothing if two-step smoothing fails
            ms_final.apply_coord_laplacian_smoothing_scale_dependent()
        
        # Mesh Repair Step 2
        ms_final.meshing_remove_null_faces()               # Remove faces with zero area
        ms_final.meshing_remove_duplicate_faces()          # Remove duplicate faces
        ms_final.meshing_remove_duplicate_vertices()       # Remove duplicate vertices   
        ms_final.meshing_remove_unreferenced_vertices()    # Remove unreferenced vertices

        # Retrieve final processed vertices and faces
        final_vertices = ms_final.current_mesh().vertex_matrix()
        final_faces = ms_final.current_mesh().face_matrix()
        
        # Write the final mesh to a VTK file
        out_vtk_file_path = os.path.join(result_VTK_filePath, f'{fault_name}.vtk')
        write_vtk(out_vtk_file_path, final_vertices, final_faces, final_vertices[:, -1])
        vtk_append_files(out_vtk_file_path, append_filter)
        
        # Write the final mesh to an OBJ file
        out_obj_file_path = os.path.join(result_OBJ_filePath, f'{fault_name}.obj')
        write_obj(out_obj_file_path, final_vertices, final_faces)

    # Combine all VTK datasets and write to a single VTK file
    append_filter.Update()
    combined_vtk_file = os.path.normpath(os.path.join(output_path, "final_faults_3D.vtk"))
    writer2 = vtk.vtkDataSetWriter()
    writer2.SetFileName(combined_vtk_file)
    writer2.SetInputData(append_filter.GetOutput())
    writer2.Write()
    
    # Visualize the combined VTK file
    visualize_vtk_file(combined_vtk_file)

def write_vtk(filename, vertices, triangles, depth_attribute):
    """Write the mesh data to a VTK file."""
    
    # Create an unstructured grid
    ugrid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()

    # Create a scalar array for depth
    depth = vtk.vtkFloatArray()
    depth.SetNumberOfComponents(1)
    depth.SetName('fault_depth_km')
    
    # Insert points and depth values
    for v, dep in zip(vertices, depth_attribute):
        points.InsertNextPoint(v)
        depth.InsertNextValue(round(dep * 0.001, 3))  # Convert depth to km
    
    # Insert triangles into the unstructured grid
    for t in triangles:
        ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, t)
    
    # Assign points and depth data to the grid
    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(depth)
    
    # Write the unstructured grid to a VTK file
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(filename)
    writer.Update()

if __name__ == "__main__":
    main()
