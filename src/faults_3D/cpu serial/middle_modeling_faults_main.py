import pymeshlab
import os
import vtk
import numpy as np
import trimesh
import random
import pandas as pd
from scipy.spatial import cKDTree
import sys

sys.path.append(r"src")
import CONST

sys.path.append(r'src\faults_3D\lib')
from vis_func import (
    write_obj,
    vtk_append_files,
    visualize_vtk_file
)



input_OBJ_filePath = r'resultData\faults_3D\Initial_modeling_faults\faults_obj'

output_path = r'resultData\faults_3D\remesh_faults'
result_OBJ_filePath = os.path.join(output_path, r'remesh_faults_obj')
result_VTK_filePath = os.path.join(output_path, r'remesh_faults_vtk')
result_CSV_filePath = os.path.join(output_path, r'remesh_faults_points_csv')
os.makedirs(result_OBJ_filePath, exist_ok=True)
os.makedirs(result_VTK_filePath, exist_ok=True)
os.makedirs(result_CSV_filePath, exist_ok=True)

unit_area = 1000.0 ** 2.0

Ball_small_faults = [
    'Anninghe_Zemuhe_Xiaojiang_Faults', 'Bailongjiang_fault_N_1', 'Bailongjiang_fault_N_2', 'Bailongjiang_fault_S',
    'Daliangshan_1_1', 'Daliangshan_2_1', 'Daliangshan_3_1', 'Daliangshan_4_1',
    'Dari_Fault', 'E Kunlun_fault_Maqin_Maqu_N', 'Emei_fault', 'Ganzi-Yushu_Fault-1', 'Ganzi-Yushu_Fault-2',
    'Ganzi-Yushu_Fault-3', 'Ganzi-Yushu_Fault-4', 'Range_Frontal_Thrust_fault', 'Yingxiu_fault',
    'Guanggaishan-Dieshan_fault_S', 'Hanan_fault_N', 'Hanan_fault_S_1', 'Hanyuan_Ganluo_fault_2',
    'Honghe Fault-Range_Front_Fault', 'Huya_fault_N', 'Tazang_fault_5', 'weixi-qiaohou fault',
    'Yuanmou_fault_1_1', 'Huya_fault_S', 'Jiangyou_Guangyuan_fault', 'jinshajiangzhu fault',
    'Jinshajiang_fault_2', 'lancangjiang fault', 'Langmusi_fault', 'Longquanshan_fault_1',
    'Longquanshan_fault_2', 'Longquanshan_fault_3', 'Longquanshan_fault_4',
    'Longriba_fault-Longriqu', 'Longriba_fault-Maoergai', 'Longriba_fault_S_2',
    'Longriba_fault_S_3', 'Maduo-Gande_Fault-1', 'Maduo-Gande_Fault-2-1', 'Maduo-Gande_Fault-2-2',
    'Maowen_Wenchuan_fault', 'Minjiang_fault', 'nujiang fault', 'Pengguan_fault',
    'Pengxian_blind_fault', 'Puduhe_fault_1_0', 'Puduhe_fault_2_0', 'Puduhe_fault_3_0',
    'Puduhe_fault_4_0', 'Qingchuan_fault', 'Sansuchang_fault_1', 'Shiping-Jianshui_fault_1',
    'Songgang-Fubianhe_fault', 'Tanglang-Yimeng_fault_1_0', 'Tanglang-Yimeng_fault_2_0',
    'Tazang_fault_1', 'Tazang_fault_2', 'Tazang_fault_4', 'Xiongpo_fault_2',
    'Wenxian_fault_1', 'Wenxian_fault_2', 'Wudaoliang-Changshagongma_Fault', 'Xianshuihe_Fault-1',
    'Xianshuihe_Fault-2', 'Xianshuihe_Fault-3', 'Xianshuihe_fault_1', 'Xianshuihe_fault_2',
    'zigashi-deqing fault'
]

# List of faults with holes
holes_faults = ['weixi-qiaohou fault', 'lancangjiang fault', 'Hanyuan_Ganluo_fault_1', 'ganzi-litang faults2']



def main():
    # Create a vtkAppendFilter object to combine multiple datasets
    append_filter = vtk.vtkAppendFilter()
    
    for fault in os.listdir(input_OBJ_filePath):
        print(f"Remeshing fault:  {fault} ...")
        
        fault_name = fault.split('.')[0]
        fault_file = os.path.join(input_OBJ_filePath, fault)

        # Calculate the surface area of the triangular mesh to determine the sampling point density
        area_mesh = trimesh.load_mesh(fault_file)
        num_samples = int(round((area_mesh.area) / unit_area, 0))
        
        # Create MeshSet and load mesh
        ms = pymeshlab.MeshSet()               
        ms.load_new_mesh(fault_file)
        
        # Mesh repair 1
        ms.meshing_remove_null_faces()               # Remove faces with zero area
        ms.meshing_remove_duplicate_faces()          # Remove duplicate faces
        ms.meshing_remove_duplicate_vertices()       # Remove duplicate vertices   
        ms.meshing_remove_unreferenced_vertices()    # Remove vertices not referenced by any face
        # ms.meshing_remove_folded_faces()
        if len(ms.current_mesh().vertex_matrix()) >= 1000:
            ms.meshing_close_holes(maxholesize=10)
        
        # Resampling
        ms.generate_sampling_poisson_disk(samplenum=num_samples * 5)
        
        xyz0 = ms.current_mesh().vertex_matrix()
        mean_step = svole_mean_step(xyz0)
        xyz1 = remove_close_points_average(xyz0, mean_step / 2)
        xyz = remove_close_points(xyz1, mean_step / 4)
        
        # Create MeshSet and load mesh
        ms2 = pymeshlab.MeshSet()
        mesh2 = pymeshlab.Mesh(vertex_matrix=xyz)
        ms2.add_mesh(mesh2)

        # Determine the ball radius based on fault name
        if fault_name in ["Lushan_EQ_fault_2", "Yingxiu_Beichuan_fault", 'Honghe Fault-Mid_Valley_Fault', 'Guanggaishan-Dieshan_fault_N_1']:
            radius1 = 8 + 5
        elif fault_name in ['41_1_fault', 'ganzi-litang faults', 'bianba-luolong fault', 'Hanyuan_Ganluo_fault_1']:
            radius1 = 160 + 5
        elif fault_name in Ball_small_faults:
            radius1 = 80 + 5
        else:
            radius1 = 300 + 5
        
        for i in range(1, radius1 + 1):
            if i <= 5:
                i = 0
            else:
                i -= 5
            ballradius_percentage1 = pymeshlab.PercentageValue(0.25 * i)
            ms2.generate_surface_reconstruction_ball_pivoting(ballradius=ballradius_percentage1)
        
        # Surface smoothing 1
        # ms2.apply_coord_laplacian_smoothing()                     # May cause surface shrinkage
        # ms2.apply_coord_hc_laplacian_smoothing()
        ms2.apply_coord_laplacian_smoothing_scale_dependent()
        try:
            ms2.apply_coord_two_steps_smoothing()
            ms2.apply_coord_taubin_smoothing()
        except:
            ms2.apply_coord_taubin_smoothing()
        
        # Mesh repair 2
        if len(ms2.current_mesh().vertex_matrix()) >= 3000:
            ms2.meshing_close_holes()                                           # maxholesize = 30
        elif len(ms2.current_mesh().vertex_matrix()) >= 1000:
            ms2.meshing_close_holes(maxholesize=10)
            
        ms2.generate_sampling_poisson_disk(samplenum=num_samples)
        
        if fault_name in ['Guanggaishan-Dieshan_fault_N_1']:
            ms2.generate_surface_reconstruction_ball_pivoting()
        else:
            ms2.generate_surface_reconstruction_ball_pivoting()
            ms2.generate_surface_reconstruction_ball_pivoting(ballradius=pymeshlab.PercentageValue(2))
            ms2.generate_surface_reconstruction_ball_pivoting(ballradius=pymeshlab.PercentageValue(4))
                    
        # Mesh repair 2
        ms2.meshing_remove_null_faces()               # Remove faces with zero area
        ms2.meshing_remove_duplicate_faces()          # Remove duplicate faces
        ms2.meshing_remove_duplicate_vertices()       # Remove duplicate vertices   
        ms2.meshing_remove_unreferenced_vertices()    # Remove vertices not referenced by any face
        
        vertices = ms2.current_mesh().vertex_matrix()
        faces = ms2.current_mesh().face_matrix()
        
        # Write OBJ file
        out_obj_file = os.path.join(result_OBJ_filePath, fault)
        write_obj(out_obj_file, vertices, faces)
        
        # Write VTK file
        out_vtk_file = os.path.join(result_VTK_filePath, fault.replace('.obj', '.vtk'))
        write_vtk(out_vtk_file, vertices, faces, vertices[:, -1])
        vtk_append_files(out_vtk_file, append_filter)
        
        # Write CSV file
        out_csv_file = os.path.join(result_CSV_filePath, fault.replace('.obj', '.csv'))
        vertices_Data = pd.DataFrame(vertices, columns=['x', 'y', 'z'])
        vertices_Data.to_csv(out_csv_file, index=False)


    # Update and write the combined dataset to a new VTK file
    append_filter.Update()
    
    combined_vtk_file = os.path.normpath(  os.path.join(output_path, f"middle_faults_3D.vtk")  )
    writer2 = vtk.vtkDataSetWriter()
    writer2.SetFileName(combined_vtk_file)
    writer2.SetInputData(append_filter.GetOutput())
    writer2.Write()
    
    # Visualize the combined VTK file
    visualize_vtk_file(combined_vtk_file)



def write_vtk(filename, vertices, triangles, depth_Attribute):
    """Write the mesh data to a VTK file."""
    
    ugrid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()

    depth = vtk.vtkFloatArray()
    depth.SetNumberOfComponents(1)
    depth.SetName('fault depth(km)')
    
    for v,dep in zip(vertices, depth_Attribute):
        points.InsertNextPoint(v)
        depth.InsertNextValue(  dep / 1000.0  )
    for t in triangles:
        ugrid.InsertNextCell(  vtk.VTK_TRIANGLE, 3, t  )
    
    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(depth)
    
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(filename)
    writer.Update()


def remove_close_points_average(points, min_distance):
    """Remove points closer than the minimum distance and replace them with the average."""
    tree = cKDTree(points)
    to_remove = set()
    to_add = []

    for i, p in enumerate(points):
        if i in to_remove:
            continue
        indices = tree.query_ball_point(p, r=min_distance)
        indices = [idx for idx in indices if idx != i]
        if indices:
            cluster_points = points[indices]
            cluster_points = np.vstack((cluster_points, p))
            mean_point = np.mean(cluster_points, axis=0)
            to_add.append(mean_point)
            to_remove.update(indices)
            to_remove.add(i)

    mask = np.ones(len(points), dtype=bool)
    mask[list(to_remove)] = False
    remaining_points = points[mask]
    
    if to_add:
        return np.vstack((remaining_points, to_add))
    else:
        return remaining_points


def remove_close_points(points, min_distance):
    """Remove points closer than the minimum distance."""
    tree = cKDTree(points)
    to_remove = set()

    for i, p in enumerate(points):
        if i in to_remove:
            continue
        indices = tree.query_ball_point(p, r=min_distance)
        indices = [idx for idx in indices if idx != i]
        if indices:
            chosen = np.random.choice(indices)
            to_remove.add(chosen)

    mask = np.ones(len(points), dtype=bool)
    mask[list(to_remove)] = False
    return points[mask]


def svole_mean_step(pointsXYZ):
    """Calculate the average step based on the distance to the nearest neighbor."""
    random_p = random.randint(0, len(pointsXYZ) - 1)
    p_kdtree = cKDTree(pointsXYZ)
    indices = p_kdtree.query_ball_point(pointsXYZ[random_p], r=(unit_area ** 0.5) * 2.5)
    
    xyz_P = pointsXYZ[indices]
    xyz_kdtree = cKDTree(xyz_P)
    step = []
    for xyz in xyz_P:
        _, p0_indice = xyz_kdtree.query(xyz)
        step.append(np.linalg.norm(xyz - xyz_P[p0_indice]))
    
    return sum(step) / len(step)


if __name__ == "__main__":
    main()
