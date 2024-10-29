"""
Summary:
    Projects 3D scatter points onto an appropriate face of an Axis-Aligned Bounding Box (AABB),
    then performs Delaunay 2D triangulation.
Returns:
    tuple: 
        - numpy.ndarray: 3D scatter points
        - list of tuples: Triangle mesh connectivity information
"""

import os
import pandas as pd
import numpy as np
from scipy.spatial import Delaunay



def AABB_Delaunay_TriMesh(xy):
    """
    Performs Delaunay triangulation on 2D points.

    Parameters:
        xy (numpy.ndarray): 2D points for triangulation.

    Returns:
        numpy.ndarray: Indices of the points forming the simplices (triangles).
    """
    tri = Delaunay(xy)
    return tri.simplices


def grid_quality_control(mesh, xyz, control_factor):
    """
    Controls the quality of the mesh based on angle and edge rate criteria.

    Parameters:
        mesh (list of tuples): Triangle mesh connectivity information.
        xyz (numpy.ndarray): Original 3D points.
        control_factor (float): Threshold factor for normal vector alignment.

    Returns:
        list of tuples: Filtered triangle mesh based on quality criteria.
    """
    angle_edgesRate = []
    sur_normal = []
    for simplex in mesh:
        min_angle, max_angle, edge_rate, normal, max_edge = angle_edgeRate_normal_of_triangle(
            xyz[simplex[0]], xyz[simplex[1]], xyz[simplex[2]]
        )
        angle_edgesRate.append([min_angle, max_angle, edge_rate])
        sur_normal.append(normal)

    sur_normal = np.array(sur_normal)
    mean_sur_normal = np.mean(sur_normal, axis=0)

    tri_mesh = []
    for i in range(len(mesh)):
        if (abs(np.dot(sur_normal[i], mean_sur_normal)) >= control_factor) and (angle_edgesRate[i][0] >= 1.5):
            tri_mesh.append(mesh[i])

    return tri_mesh


def angle_edgeRate_normal_of_triangle(p1, p2, p3):
    """
    Calculates the minimum angle, maximum angle, edge rate, normal vector, and maximum edge length of a triangle.

    Parameters:
        p1, p2, p3 (numpy.ndarray): The three vertices of the triangle.

    Returns:
        tuple:
            - float: Minimum angle in degrees.
            - float: Maximum angle in degrees.
            - float: Edge rate (minimum edge length divided by maximum edge length).
            - numpy.ndarray: Normal vector of the triangle.
            - float: Maximum edge length.
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

    # Return the minimum angle, maximum angle, edge rate, normal vector, and maximum edge length
    return min(angles), max(angles), min(lengths) / max(lengths), normal_vector, max(lengths)


def AABB_Delaunay_mesh(fault_file):
    """
    Generates a Delaunay triangulated mesh from a fault file.

    Parameters:
        fault_file (str): Path to the CSV file containing fault data with 'x', 'y', 'z' columns.

    Returns:
        tuple:
            - numpy.ndarray: 3D scatter points.
            - list of tuples: Triangle mesh connectivity information.
    """
    fault_name = os.path.splitext(os.path.basename(fault_file))[0]
    df = pd.read_csv(fault_file)

    xyz = df[['x', 'y', 'z']].values
    # Calculate the center of the AABB
    AABB_box_center = (np.min(xyz, axis=0) + np.max(xyz, axis=0)) / 2.0

    # Translate points to center the AABB
    translated_xyz = xyz - AABB_box_center
    sur_normal1 = translated_xyz[0] / np.linalg.norm(translated_xyz[0])
    sur_normal2 = translated_xyz[-1] / np.linalg.norm(translated_xyz[-1])
    sur_normal = (sur_normal1 + sur_normal2) / 2.0

    # Determine the axis with the maximum absolute normal component
    abs_sur_normal = np.abs(sur_normal).tolist()
    axis_id = abs_sur_normal.index(max(abs_sur_normal))

    # Project points onto the dominant plane for 2D triangulation
    if axis_id == 0:
        projected_xy = xyz[:, [1, 2]]
    elif axis_id == 1:
        projected_xy = xyz[:, [0, 2]]
    else:
        projected_xy = xyz[:, [0, 1]]

    faces = AABB_Delaunay_TriMesh(projected_xy)

    # Apply grid quality control based on fault name
    if fault_name in ['Yanjing_Wulong_fault', '22_fault']:
        tri_mesh = grid_quality_control(faces, xyz, 0.55)
    elif fault_name in ['wanting-anding fault', 'Nantinghe_Fault_fault', 'jinshajiang fault']:
        tri_mesh = grid_quality_control(faces, xyz, 0.35)
    elif fault_name in ['Xihe-Meigu_fault']:
        tri_mesh = grid_quality_control(faces, xyz, 0.4)
    elif fault_name in ['Malao_fault']:
        tri_mesh = grid_quality_control(faces, xyz, 0.6)
    else:
        tri_mesh = grid_quality_control(faces, xyz, 0.0)

    return xyz, tri_mesh
