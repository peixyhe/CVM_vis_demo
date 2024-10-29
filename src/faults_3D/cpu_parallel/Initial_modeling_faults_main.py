import os
import sys
import vtk
from concurrent.futures import ProcessPoolExecutor, as_completed

# Add custom module paths
sys.path.append(r'src\faults_3D\lib')
from pca import PCA_mesh
from meshLab_Ball_Pivoting import Ball_Pivoting_mesh
from AABB_Delaunay_mesh import AABB_Delaunay_mesh
from diagonal_line import diagonal_line_mesh
from transData import preprocess_fault_data
from vis_func import (
    write_obj,
    vtk_append_files,
    visualize_vtk_file
)



# Define input and output paths
INPUT_RAW_DATA = os.path.join('rawData', 'CSES_3DFM_V2.0.txt')

OUTPUT_ROOT_PATH = os.path.join('resultData', 'faults_3D', 'Initial_modeling_faults')
os.makedirs(OUTPUT_ROOT_PATH, exist_ok=True)

OUTPUT_TXT_PATH = os.path.join(OUTPUT_ROOT_PATH, 'faults_noDenseData_txt')
os.makedirs(OUTPUT_TXT_PATH, exist_ok=True)

OUTPUT_VTK_PATH = os.path.join(OUTPUT_ROOT_PATH, 'faults_vtk')
os.makedirs(OUTPUT_VTK_PATH, exist_ok=True)

OUTPUT_OBJ_PATH = os.path.join(OUTPUT_ROOT_PATH, 'faults_obj')
os.makedirs(OUTPUT_OBJ_PATH, exist_ok=True)

# Define fault lists for different meshing methods
ABANDON_FAULTS = ['20_fault', 'dege-xiangcheng fault']

PCA_FAULTS = [
    '10_fault', '12_fault', '13_fault', '14_fault', '15_fault', '19_fault', '21_fault', '23_fault',
    '31_fault', '37_3_fault', '41_2_fault', '41_3_fault', '57_fault', '60_fault', '6_fault',
    '9_fault', 'Awancang_Fault-1', 'Awancang_Fault-2', 'Bailongjiang_fault_N_3',
    'Bailongjiang_fault_N_4', 'Baishuxi_fault', 'baqing-leiniaoqiao fault', 'Dachuan_Shuangshi_fault',
    'Daguanbei_fault', 'Dayi_blind_fault', 'Dianlanba_fault', 'Guanggaishan-Dieshan_fault_N_2',
    'Hanan_fault_S_2', 'hanmuba fault', 'heihe fault', 'honghe fault', 'Huayingshan_fault_1',
    'Huayingshan_fault_2', 'Huya_fault_M', 'jiali fault', 'jiali fault2', 'Leibo-Mahu_fault_2',
    'Leibo-Mahu_fault_3', 'Lianfeng_fault_1', 'Lianfeng_fault_2', 'lijiang-xiaojinhe fault',
    'Lijiang-Xiaojinhe_fault-E_fault', 'Lijiang-Xiaojinhe_fault-W_fault', 'Longriba_fault_N_E',
    'Longriba_fault_S_1', 'Longriba_fault_S_4', 'Lushan_EQ_fault_1', 'Lushan_EQ_fault_3',
    'malipo-zengjiaying fault', 'Meigu-Yiping_fault', 'niqu-yuke fault', 'nujiangxizhi fault',
    'puer fault', 'puwen fault', 'Qiaoxi-Weihou-Weishan_fault-3', 'Qiongxi_fault', 'Qujiang_fault_1',
    'Sansuchang_fault_2', 'shidian fault', 'sudian fault', 'Tazang_fault_3', 'Tazang_fault_6',
    'Tianquan_fault', 'tuojiao-tengchong fault', 'Wanwantan_fault', 'Xiaxi_fault', 'Xiezhiba_fault',
    'xinlong fault', 'Xiongpo_fault_1', 'Xiongpo_fault_3', 'Xueshanliangzi_fault', 'Yanjingou_fault',
    'Yazha-Zhuyuanba_fault', 'yinzhidaohu fault2', 'yongkang fault', 'Yunongxi_fault',
    'yupengcun fault', 'zengke-shuoqu fault', 'zengke-shuoqu fault2', 'Zhongcun_fault',
    'Zhongdu_fault', 'mogang fault', 'ninglang fault', 'longling-ruili fault', 'yinzhidaohu fault',
    '54_fault', '55_fault'
]

BALL_SMALL_FAULTS = [
    'Anninghe_Zemuhe_Xiaojiang_Faults', 'Bailongjiang_fault_N_1', 'Bailongjiang_fault_N_2',
    'Bailongjiang_fault_S', 'Daliangshan_1_1', 'Daliangshan_2_1', 'Daliangshan_3_1',
    'Daliangshan_4_1', 'Dari_Fault', 'E Kunlun_fault_Maqin_Maqu_N', 'Emei_fault',
    'Ganzi-Yushu_Fault-1', 'Ganzi-Yushu_Fault-2', 'Ganzi-Yushu_Fault-3', 'Ganzi-Yushu_Fault-4',
    'Guanggaishan-Dieshan_fault_N_1', 'Guanggaishan-Dieshan_fault_S', 'Hanan_fault_N',
    'Hanan_fault_S_1', 'Hanyuan_Ganluo_fault_1', 'Hanyuan_Ganluo_fault_2',
    'Honghe Fault-Mid_Valley_Fault', 'Honghe Fault-Range_Front_Fault', 'Huya_fault_N',
    'Huya_fault_S', 'Jiangyou_Guangyuan_fault', 'jinshajiangzhu fault', 'Jinshajiang_fault_2',
    'lancangjiang fault', 'Langmusi_fault', 'Longquanshan_fault_1', 'Longquanshan_fault_2',
    'Longquanshan_fault_3', 'Longquanshan_fault_4', 'Longriba_fault-Longriqu',
    'Longriba_fault-Maoergai', 'Longriba_fault_S_2', 'Longriba_fault_S_3', 'Lushan_EQ_fault_2',
    'Maduo-Gande_Fault-1', 'Maduo-Gande_Fault-2-1', 'Maduo-Gande_Fault-2-2',
    'Maowen_Wenchuan_fault', 'Minjiang_fault', 'nujiang fault', 'Pengguan_fault',
    'Pengxian_blind_fault', 'Puduhe_fault_1_0', 'Puduhe_fault_2_0', 'Puduhe_fault_3_0',
    'Puduhe_fault_4_0', 'Qingchuan_fault', 'Range_Frontal_Thrust_fault', 'Sansuchang_fault_1',
    'Shiping-Jianshui_fault_1', 'Songgang-Fubianhe_fault', 'Tanglang-Yimeng_fault_1_0',
    'Tanglang-Yimeng_fault_2_0', 'Tazang_fault_1', 'Tazang_fault_2', 'Tazang_fault_4',
    'Tazang_fault_5', 'weixi-qiaohou fault', 'Wenxian_fault_1', 'Wenxian_fault_2',
    'Wudaoliang-Changshagongma_Fault', 'Xianshuihe_Fault-1', 'Xianshuihe_Fault-2',
    'Xianshuihe_Fault-3', 'Xianshuihe_fault_1', 'Xianshuihe_fault_2', 'Xiongpo_fault_2',
    'Yingxiu_Beichuan_fault', 'Yingxiu_fault', 'Yuanmou_fault_1_1', 'zigashi-deqing fault'
]

BALL_BIGGER_FAULTS = [
    '37_1_fault', '41_1_fault', '41_4_fault', '63_fault', 'bianba-luolong fault',
    'E Kunlun_fault_Maqin_Maqu_S', 'ganzi-litang faults', 'ganzi-litang faults2',
    'hanmuba-lancang fault', 'Maisu_fault', 'xisashi fault'
]

DIAGONAL_LINE_FAULTS = [
    'Yinjing_fault', 'Sanhekou_Yanfeng_fault', 'Jinyang_fault', 'Leibo-Mahu_fault_1',
    'Xintian_fault', 'Malao_fault'
]

AABB_DELAUNAY_MESH_FAULTS = [
    'Huangyi_fault', 'Jinping_fault', 'Sanhe_fault', 'wanting-anding fault', 'Xihe-Meigu_fault',
    'Yanjing_Wulong_fault', 'jinshajiang fault', '22_fault', 'Nantinghe_Fault_fault'
]



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
        depth.InsertNextValue(dep / 1000.0)  # Convert depth to km

    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(depth)

    # Create triangle cells
    cells = vtk.vtkCellArray()
    for t in triangles:
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, int(t[0]))
        triangle.GetPointIds().SetId(1, int(t[1]))
        triangle.GetPointIds().SetId(2, int(t[2]))
        cells.InsertNextCell(triangle)
    ugrid.SetCells(vtk.VTK_TRIANGLE, cells)

    # Write to VTK file
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(ugrid)
    writer.Write()


def process_fault_file(fault_file, OUTPUT_TXT_PATH, OUTPUT_VTK_PATH, OUTPUT_OBJ_PATH):
    """Process a single fault file: meshing and writing VTK/OBJ files."""
    try:
        fault_path = os.path.join(OUTPUT_TXT_PATH, fault_file)
        fault_name, _ = os.path.splitext(fault_file)
        print(f"Processing fault:  {fault_name} ...")

        # Determine the meshing method based on fault name
        if fault_name in ["Lushan_EQ_fault_2", "Yingxiu_Beichuan_fault", 'Honghe Fault-Mid_Valley_Fault']:
            vertex, face = Ball_Pivoting_mesh(fault_path, 410)
        elif fault_name in ['41_1_fault', 'ganzi-litang faults', 'bianba-luolong fault']:
            vertex, face = Ball_Pivoting_mesh(fault_path, 8010)
        elif fault_name in BALL_SMALL_FAULTS:
            vertex, face = Ball_Pivoting_mesh(fault_path, 4010)
        elif fault_name in BALL_BIGGER_FAULTS:
            vertex, face = Ball_Pivoting_mesh(fault_path, 15010)
        elif fault_name in PCA_FAULTS:
            vertex, face = PCA_mesh(fault_path)
        elif fault_name in DIAGONAL_LINE_FAULTS:
            vertex, face = diagonal_line_mesh(fault_path)
        elif fault_name in AABB_DELAUNAY_MESH_FAULTS:
            vertex, face = AABB_Delaunay_mesh(fault_path)
        else:
            print(f"Skipping unknown fault:  {fault_name}")
            return None  # Skip unknown faults

        # Write VTK file
        vtk_file_name = os.path.join(OUTPUT_VTK_PATH, f"{fault_name}.vtk")
        write_vtk(vtk_file_name, vertex, face, vertex[:, -1])

        # Write OBJ file
        obj_file_name = os.path.join(OUTPUT_OBJ_PATH, f"{fault_name}.obj")
        write_obj(obj_file_name, vertex, face)

        return vtk_file_name  # Return the path for aggregation

    except Exception as e:
        print(f"Error processing  {fault_file}:  {e}")
        return None


def main():
    print('Converting data...')
    preprocess_fault_data(INPUT_RAW_DATA, OUTPUT_ROOT_PATH)

    # List all fault files to process
    fault_files = [f for f in os.listdir(OUTPUT_TXT_PATH) if os.path.isfile(os.path.join(OUTPUT_TXT_PATH, f))]

    vtk_files = []

    # Determine the number of parallel workers (half the CPU cores to balance performance and memory usage)
    cpu_cores = os.cpu_count() or 1
    max_workers = max(1, cpu_cores // 2)
    
    print(f"Number of parallel processes: {max_workers}")
    
    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_fault = {
            executor.submit(process_fault_file, fault_file, OUTPUT_TXT_PATH, OUTPUT_VTK_PATH, OUTPUT_OBJ_PATH): fault_file
            for fault_file in fault_files
        }

        # Collect results as they complete
        for future in as_completed(future_to_fault):
            fault_file = future_to_fault[future]
            try:
                vtk_file = future.result()
                if vtk_file:
                    vtk_files.append(vtk_file)
            except Exception as exc:
                print(f"{fault_file} generated an exception: {exc}")

    # Create a vtkAppendFilter object to combine multiple datasets
    append_filter = vtk.vtkAppendFilter()

    # Append all generated VTK files
    for vtk_file in vtk_files:
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtk_file)
        reader.Update()
        append_filter.AddInputData(reader.GetOutput())

    # Update and write the combined dataset to a new VTK file
    append_filter.Update()
    combined_vtk_file = os.path.join(OUTPUT_ROOT_PATH, "initial_faults_3D.vtk")
    writer = vtk.vtkDataSetWriter()
    writer.SetFileName(combined_vtk_file)
    writer.SetInputData(append_filter.GetOutput())
    writer.Write()

    # Visualize the combined VTK file
    visualize_vtk_file(combined_vtk_file)


if __name__ == "__main__":
    main()
