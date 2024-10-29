"""
Summary:
    Generates a 3D scatter point mesh based on the Ball Pivoting algorithm.
Returns:
    tuple:
        - numpy.ndarray: 3D scatter points
        - list of lists: Triangle mesh connectivity information
"""

import os
import numpy as np
import pymeshlab



# Specify the number of layers for discrete points for each fault
faults_layer_num = {
    'Yinjing_fault': 5,
    'Sanhekou_Yanfeng_fault': 2,
    'Jinyang_fault': 2,
    'Leibo-Mahu_fault_1': 2,
    'Xintian_fault': 2,
    'Malao_fault': 2
}



def angle_edgeRate_normal_of_triangle(p1, p2, p3):
    """
    Calculates the minimum angle, maximum angle, edge rate, and normal vector of a triangle.

    Parameters:
        p1, p2, p3 (numpy.ndarray): The three vertices of the triangle.

    Returns:
        tuple:
            - float: Minimum angle in degrees.
            - float: Maximum angle in degrees.
            - float: Edge rate (minimum edge length divided by maximum edge length).
            - numpy.ndarray: Normal vector of the triangle.
    """
    # Calculate the vectors representing the edges of the triangle
    p12 = p2 - p1
    p13 = p3 - p1
    p23 = p3 - p2

    # Calculate the lengths of the edges
    lengths = [np.linalg.norm(p12), np.linalg.norm(p13), np.linalg.norm(p23)]

    def angle_between(v1, v2):
        """
        Calculates the angle between two vectors in degrees.

        Parameters:
            v1, v2 (numpy.ndarray): Input vectors.

        Returns:
            float: Angle between v1 and v2 in degrees.
        """
        unit_v1 = v1 / np.linalg.norm(v1)
        unit_v2 = v2 / np.linalg.norm(v2)
        dot_product = np.dot(unit_v1, unit_v2)
        angle = np.arccos(np.clip(dot_product, -1.0, 1.0))
        return np.degrees(angle)

    # Calculate the three angles of the triangle
    angles = [
        angle_between(p12, p13),
        angle_between(-p12, p23),
        angle_between(-p13, -p23)
    ]

    # Calculate the normal vector of the triangle
    normal_vector = np.cross(p12, p13)
    norm = np.linalg.norm(normal_vector)
    if norm != 0:
        normal_vector /= norm
    else:
        normal_vector = np.zeros(3)

    # Return the minimum angle, maximum angle, edge rate, and normal vector
    return min(angles), max(angles), min(lengths) / max(lengths), normal_vector


def mesh_adjustment(cells, p_xyz, control_factor):
    """
    Controls the quality of the mesh based on angle, edge rate, and normal vector criteria.

    Parameters:
        cells (list of lists): Triangle mesh connectivity information.
        p_xyz (numpy.ndarray): Original 3D points.
        control_factor (float): Threshold factor for normal vector alignment.

    Returns:
        list of lists: Filtered triangle mesh based on quality criteria.
    """
    angle_edgesRate = []
    sur_normal = []
    for simplex in cells:
        min_angle, max_angle, edge_rate, normal = angle_edgeRate_normal_of_triangle(
            p_xyz[simplex[0]], p_xyz[simplex[1]], p_xyz[simplex[2]]
        )
        angle_edgesRate.append([min_angle, max_angle, edge_rate])
        sur_normal.append(normal)

    sur_normal = np.array(sur_normal)
    mean_sur_normal = np.mean(sur_normal, axis=0)

    tri_mesh = []
    for i in range(len(cells)):
        if (abs(np.dot(sur_normal[i], mean_sur_normal)) >= control_factor and
            angle_edgesRate[i][0] >= 2.5 and
            angle_edgesRate[i][1] <= 150.0 and
            angle_edgesRate[i][2] >= 0.02):
            tri_mesh.append(cells[i])

    return tri_mesh


def Ball_Pivoting_mesh(fault_file, recons_num):
    """
    Generates a Delaunay triangulated mesh using the Ball Pivoting algorithm.

    Parameters:
        fault_file (str): Path to the file containing fault data.
        recons_num (int): Number of reconstruction iterations.

    Returns:
        tuple:
            - numpy.ndarray: 3D scatter points.
            - list of lists: Triangle mesh connectivity information.
    """
    fault_name = os.path.splitext(os.path.basename(fault_file))[0]

    ms = pymeshlab.MeshSet()  # Create a MeshSet
    ms.load_new_mesh(fault_file, rowtoskip=1, strformat='X Y Z', separator=',')

    # Define faults that require multiple Ball Pivoting iterations
    specific_faults = [
        "Jiangyou_Guangyuan_fault",
        "Guanggaishan-Dieshan_fault_S",
        "Daliangshan_1_1",
        "Dachuan_Shuangshi_fault",
        "Xiezhiba_fault",
        "Sanhekou_Yanfeng_fault",
        "Honghe Fault-Range_Front_Fault",
        "Yazha-Zhuyuanba_fault"
    ]

    if fault_name in specific_faults:
        for i in range(1, 11):
            ms.generate_surface_reconstruction_ball_pivoting()
    else:
        for i in range(1, 1 + recons_num):
            if i <= 10:
                iteration = 0
            else:
                iteration = i - 10

            ballradius_percentage = pymeshlab.PercentageValue(0.005 * iteration)
            ms.generate_surface_reconstruction_ball_pivoting(ballradius=ballradius_percentage)

    vertices = ms.current_mesh().vertex_matrix()
    faces = ms.current_mesh().face_matrix()

    # Apply mesh quality control based on fault name
    if fault_name in ["Anninghe_Zemuhe_Xiaojiang_Faults"]:
        cells = mesh_adjustment(faces, vertices, 0.2)
    elif fault_name in ['Daliangshan_1_1']:
        cells = mesh_adjustment(faces, vertices, 0.5)
    elif fault_name not in ["Jiangyou_Guangyuan_fault"]:
        cells = mesh_adjustment(faces, vertices, 0.1)
    else:
        cells = faces

    return vertices, cells
