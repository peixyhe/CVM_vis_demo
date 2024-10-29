import pymeshlab
import os
import vtk
import pyproj
import numpy as np

import sys
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
os.makedirs(OUTPUT_VTK_PATH, exist_ok=True)
os.makedirs(OUTPUT_OBJ_PATH, exist_ok=True)

# Initialize UTM projection (transform UTM-xyz to lon, lat)
utm_proj = pyproj.Proj(proj='utm', zone=48, ellps='WGS84')

# List of faults to abandon
abandon_faults = [
    'bianba-luolong fault', 'Awancang_Fault-2', 'Hanan_fault_S_2', 'Lushan_EQ_fault_3', 'Maduo-Gande_Fault-2-1', 
    'Tazang_fault_6', 'Bailongjiang_fault_S', 'Tazang_fault_1', 'baqing-leiniaoqiao fault', 'Bailongjiang_fault_N_4'
]



def main():
    # Create a vtkAppendFilter object to combine multiple datasets
    append_filter = vtk.vtkAppendFilter()

    for fault in os.listdir(INPUT_OBJ_FILEPATH):
        fault_name = fault.split('.')[0]
        print(f"Processing fault: {fault}")

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
        
        # Perform Poisson disk sampling
        ms1.generate_sampling_poisson_disk(samplenum=int(len(vertices) * CONST.Z_SCALE))
        vertices1 = ms1.current_mesh().vertex_matrix()
        
        # Filter points within the specified longitude and latitude boundaries
        points_xyz = []
        for ver in vertices1:
            if (CONST.longitude_left <= ver[0] <= CONST.longitude_right) and (CONST.latitude_down <= ver[1] <= CONST.latitude_up):
                points_xyz.append(ver)
        
        if len(points_xyz) > 0:
            ms2 = pymeshlab.MeshSet()
            ms2.add_mesh(pymeshlab.Mesh(vertex_matrix=points_xyz))
            
            # Perform surface reconstruction using Ball Pivoting
            if fault_name in ['Anninghe_Zemuhe_Xiaojiang_Faults']:
                ms2.generate_surface_reconstruction_ball_pivoting()
            else:
                for i in range(3):
                    ms2.generate_surface_reconstruction_ball_pivoting(ballradius=pymeshlab.PercentageValue(2 * i))
            
            # Perform Poisson disk sampling again
            ms2.generate_sampling_poisson_disk(samplenum=int(len(points_xyz) / 50.0))
            
            # Perform surface reconstruction again
            if fault_name in ['Anninghe_Zemuhe_Xiaojiang_Faults']:
                ms2.generate_surface_reconstruction_ball_pivoting()
            else:
                for i in range(3):
                    ms2.generate_surface_reconstruction_ball_pivoting(ballradius=pymeshlab.PercentageValue(2 * i))
                    
            # Mesh repair: remove null faces, duplicate faces, duplicate vertices, and unreferenced vertices
            ms2.meshing_remove_null_faces()               # Remove faces with zero area
            ms2.meshing_remove_duplicate_faces()          # Remove duplicate faces
            ms2.meshing_remove_duplicate_vertices()       # Remove duplicate vertices   
            ms2.meshing_remove_unreferenced_vertices()    # Remove vertices not referenced by any face
            
            if fault_name not in abandon_faults:
                vertices2 = ms2.current_mesh().vertex_matrix()
                vertices2[:, 2] *= (-1.0)  # Adjust depth
                faces2 = ms2.current_mesh().face_matrix()
                
                # Write OBJ file
                out_obj_file_path = os.path.join(OUTPUT_OBJ_PATH, f"{fault_name}.obj")
                write_obj(out_obj_file_path, vertices2, faces2)
                
                # Write VTK file
                out_vtk_file_path = os.path.join(OUTPUT_VTK_PATH, f"{fault_name}.vtk")
                write_vtk(out_vtk_file_path, vertices2, faces2, vertices2[:, -1])
                
                # Append VTK file to the combined filter
                vtk_append_files(out_vtk_file_path, append_filter)

    # Update and write the combined dataset to a new VTK file
    append_filter.Update()
    
    combined_vtk_file = os.path.join(OUTPUT_PATH, f"ChuanDian_faults_3D_X{CONST.Z_SCALE}.vtk")
    writer2 = vtk.vtkDataSetWriter()
    writer2.SetFileName(combined_vtk_file)
    writer2.SetInputData(append_filter.GetOutput())
    writer2.Write()
    
    # Visualize the combined VTK file
    visualize_vtk_file(combined_vtk_file)



def write_vtk(filename, vertices, triangles, depth_attribute):
    """Write the mesh data to a VTK file."""
    
    ugrid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()

    depth = vtk.vtkFloatArray()
    depth.SetNumberOfComponents(1)
    depth.SetName('fault depth(km)')
    
    # Add vertices and depth attribute
    for v, dep in zip(vertices, depth_attribute):
        points.InsertNextPoint(v)
        depth.InsertNextValue(round(-1.0 * dep * CONST.angle_to_kilometers, 3))
    # Add triangle cells
    for t in triangles:
        ugrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, t)
    
    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(depth)
    
    # Write to VTK file
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(filename)
    writer.Update()


if __name__ == "__main__":
    main()
