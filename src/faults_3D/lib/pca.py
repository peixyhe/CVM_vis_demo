"""
Summary:
    Generates a 3D scatter point mesh based on Principal Component Analysis (PCA):
        1) Projects the original points onto the plane defined by the two principal components with the largest eigenvalues.
        2) Performs Delaunay triangulation on the projected 2D points.
        3) Maps the triangulated mesh back to the original 3D points.
Returns:
    tuple:
        - numpy.ndarray: 3D scatter points
        - list of lists: Triangle mesh connectivity information
"""

import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.spatial import Delaunay



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


def tri_mesh(fault_name, xy, p_xyz):
    """
    Performs 2D Delaunay triangulation and controls the quality of the triangular mesh.

    Parameters:
        fault_name (str): Name of the fault to determine specific quality control criteria.
        xy (numpy.ndarray): 2D projected points for triangulation.
        p_xyz (numpy.ndarray): Original 3D points.

    Returns:
        list of lists: Filtered triangle mesh connectivity information based on quality criteria.
    """
    # Perform Delaunay triangulation on the 2D points
    tri = Delaunay(xy)
    Delaunay_mesh = tri.simplices

    angle_edgesRate = []
    sur_normal = []
    for simplex in Delaunay_mesh:
        min_angle, max_angle, edge_rate, normal, max_edge = angle_edgeRate_normal_of_triangle(
            p_xyz[simplex[0]], p_xyz[simplex[1]], p_xyz[simplex[2]]
        )
        angle_edgesRate.append([min_angle, max_angle, edge_rate])
        sur_normal.append(normal)

    sur_normal = np.array(sur_normal)
    mean_sur_normal = np.mean(sur_normal, axis=0)

    tri_mesh = []
    # Define specific quality control criteria based on fault name
    if fault_name in ['Qiaoxi-Weihou-Weishan_fault-3']:
        tri_mesh = Delaunay_mesh.tolist()
    elif fault_name == 'longling-ruili fault':
        for i in range(len(Delaunay_mesh)):
            if (abs(np.dot(sur_normal[i], mean_sur_normal)) >= 0.85 and
                angle_edgesRate[i][2] >= 0.015 and
                angle_edgesRate[i][0] >= 1.5 and
                angle_edgesRate[i][1] <= 150.0 and
                max_angle <= 1500.0):
                tri_mesh.append(Delaunay_mesh[i].tolist())
    elif fault_name in ['21_fault', '15_fault', '13_fault', '12_fault', '9_fault']:
        for i in range(len(Delaunay_mesh)):
            if (abs(np.dot(sur_normal[i], mean_sur_normal)) >= 0.85 and
                angle_edgesRate[i][2] >= 0.015 and
                angle_edgesRate[i][0] >= 1.5 and
                angle_edgesRate[i][1] <= 150.0):
                tri_mesh.append(Delaunay_mesh[i].tolist())
    elif fault_name in ["Lijiang-Xiaojinhe_fault-E_fault", 'Leibo-Mahu_fault_3', 'Nantinghe_Fault_fault']:
        for i in range(len(Delaunay_mesh)):
            if (abs(np.dot(sur_normal[i], mean_sur_normal)) >= 0.5 and
                angle_edgesRate[i][2] >= 0.015 and
                angle_edgesRate[i][0] >= 1.5 and
                angle_edgesRate[i][1] <= 150.0):
                tri_mesh.append(Delaunay_mesh[i].tolist())
    elif fault_name == 'mogang fault':
        for i in range(len(Delaunay_mesh)):
            if (abs(np.dot(sur_normal[i], mean_sur_normal)) >= 0.5 and
                angle_edgesRate[i][0] >= 0.5):
                tri_mesh.append(Delaunay_mesh[i].tolist())
    elif fault_name in ['54_fault', '55_fault']:
        for i in range(len(Delaunay_mesh)):
            if (abs(np.dot(sur_normal[i], mean_sur_normal)) >= 0.15 and
                angle_edgesRate[i][0] >= 2.0):
                tri_mesh.append(Delaunay_mesh[i].tolist())
    else:
        for i in range(len(Delaunay_mesh)):
            if (abs(np.dot(sur_normal[i], mean_sur_normal)) >= 0.6 and
                angle_edgesRate[i][2] >= 0.015 and
                angle_edgesRate[i][0] >= 1.5 and
                angle_edgesRate[i][1] <= 150.0):
                tri_mesh.append(Delaunay_mesh[i].tolist())

    return tri_mesh


def lowDim(xyz):
    """
    Performs PCA dimensionality reduction to map 3D points to 2D.

    Parameters:
        xyz (numpy.ndarray): Original 3D points.

    Returns:
        numpy.ndarray: 2D projected points after PCA.
    """
    pca = PCA(n_components=2)
    pca.fit(xyz)
    return pca.transform(xyz)  # Reduced 2D data


def PCA_mesh(fault_file):
    """
    Generates a Delaunay triangulated mesh based on PCA projection of 3D scatter points.

    Parameters:
        fault_file (str): Path to the CSV file containing fault data with 'x', 'y', 'z' columns.

    Returns:
        tuple:
            - numpy.ndarray: Original 3D scatter points.
            - list of lists: Triangle mesh connectivity information.
    """
    fault_name = os.path.splitext(os.path.basename(fault_file))[0]

    # Read point cloud data
    data = pd.read_csv(fault_file)
    xyz = data[['x', 'y', 'z']].values

    # Perform PCA to reduce dimensionality to 2D
    projected_xy = lowDim(xyz)

    # Generate the triangular mesh based on the projected 2D points and original 3D points
    tri_cell = tri_mesh(fault_name, projected_xy, xyz)

    return xyz, tri_cell
