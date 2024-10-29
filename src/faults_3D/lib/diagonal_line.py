"""
Summary:
    Generates a 3D scatter point mesh based on an Axis-Aligned Bounding Box (AABB):
        1) Translates the original points to a coordinate system centered at the AABB.
        2) Separates the scatter points into layers based on the number of layers specified (faults_layer_num).
        3) Projects the scatter points onto the diagonal of the AABB and sorts them.
        4) Performs Delaunay triangulation and maps the mesh back to the original points.
Returns:
    tuple:
        - numpy.ndarray: 3D scatter points
        - list of lists: Triangle mesh connectivity information
"""

import os
import pandas as pd
import numpy as np
from scipy.spatial import Delaunay



# Specify the number of layers for discrete points for each fault
faults_layer_num = {
    'Yinjing_fault': 5,
    'Sanhekou_Yanfeng_fault': 2,
    'Jinyang_fault': 2,
    'Leibo-Mahu_fault_1': 2,
    'Xintian_fault': 2,
    'Malao_fault': 2
}



def find_nearest(array, target):
    """
    Finds the nearest number in a list to the target number.

    Parameters:
        array (list or numpy.ndarray): List of numbers to search.
        target (float): The target number.

    Returns:
        float: The number in the array closest to the target.
    """
    min_distance = float('inf')
    nearest_value = None

    # Iterate through each value in the array
    for value in array:
        # Calculate the distance between the current value and the target
        distance = abs(value - target)
        # Update the nearest value if a smaller distance is found
        if distance < min_distance:
            min_distance = distance
            nearest_value = value
    return nearest_value


def find_element_index(matrix, element):
    """
    Finds the row and column index of an element in a 2D array.

    Parameters:
        matrix (list of lists or 2D numpy.ndarray): The 2D array to search.
        element (list or numpy.ndarray): The element to find.

    Returns:
        tuple:
            - int: Row index of the element.
            - int: Column index of the element.
            If the element is not found, returns (None, None).
    """
    # Iterate through each row and column in the 2D array
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            # Check if the current value matches the target element
            if value == element:
                return i, j
    # Return (None, None) if the element is not found
    return None, None


def diagonal_line_mesh(fault_file):
    """
    Generates a Delaunay triangulated mesh based on the diagonal projection of 3D scatter points within an AABB.

    Parameters:
        fault_file (str): Path to the CSV file containing fault data with 'x', 'y', 'z' columns.

    Returns:
        tuple:
            - numpy.ndarray: Original 3D scatter points.
            - list of lists: Triangle mesh connectivity information.
    """
    fault_name = os.path.splitext(os.path.basename(fault_file))[0]
    df = pd.read_csv(fault_file)

    xyz = df[['x', 'y', 'z']].values
    # Calculate the center of the AABB
    AABB_box_center = (np.min(xyz, axis=0) + np.max(xyz, axis=0)) / 2.0

    # Translate points to center the AABB
    translated_xyz = xyz - AABB_box_center
    # Generate layer boundaries based on the number of layers for the fault
    layer_z = np.linspace(
        np.min(translated_xyz[:, -1]),
        np.max(translated_xyz[:, -1]),
        faults_layer_num.get(fault_name, 1)
    )

    # Assign each point to the nearest layer
    for i in range(len(translated_xyz)):
        z0 = find_nearest(layer_z, translated_xyz[i][-1])
        translated_xyz[i][-1] = z0

    layer_points_id = []
    delaunay_xy = []
    delaunay_xy_for_index = []
    delaunay_y = np.linspace(0, 10, faults_layer_num.get(fault_name, 1))
    
    for z, y0 in zip(layer_z, delaunay_y):
        p_id = []
        delaunay_xy0 = []
        # Select points that belong to the current layer
        layer_xyz = translated_xyz[translated_xyz[:, -1] == z]

        # Define the diagonal line on the current layer's AABB projection
        layer_box_plane = (np.min(layer_xyz, axis=0), np.max(layer_xyz, axis=0))
        p0 = np.array([layer_box_plane[1][0], layer_box_plane[0][1]])
        p1 = np.array([layer_box_plane[0][0], layer_box_plane[1][1]])

        # Calculate the direction and length of the diagonal
        u = p1 - p0
        len_u = np.linalg.norm(u)
        
        for point in layer_xyz:
            p_plane = point[0:2]
            v = p_plane - p0
            # Project the point onto the diagonal line
            projection_length = np.dot(v, u) / len_u
            translated_xyz[np.where((translated_xyz == point).all(axis=1))[0][0], -1] = projection_length

        # Sort points based on their projection length
        sorted_layer_xyz = layer_xyz[np.argsort(layer_xyz[:, -1])]

        for sorted_point in sorted_layer_xyz:
            sorted_point[-1] = z
            # Find the index of the sorted point in the original translated_xyz array
            p_id.append(np.where((translated_xyz == sorted_point).all(axis=1))[0][0])

        layer_points_id.append(p_id)

        # Generate Delaunay points for the current layer
        delaunay_x = np.linspace(0, 20, len(p_id))
        for x0 in delaunay_x:
            delaunay_xy0.append([x0, y0])
            delaunay_xy.append([x0, y0])
        delaunay_xy_for_index.append(delaunay_xy0)

    # Perform Delaunay triangulation on the 2D projected points
    cells = Delaunay(np.array(delaunay_xy)).simplices

    tri_mesh = []
    for cell in cells:
        delaunay_p0 = delaunay_xy[cell[0]]
        delaunay_p1 = delaunay_xy[cell[1]]
        delaunay_p2 = delaunay_xy[cell[2]]

        # Find the indices of the triangle vertices in the original 3D points
        p0_row_index, p0_column_index = find_element_index(delaunay_xy_for_index, delaunay_p0)
        p1_row_index, p1_column_index = find_element_index(delaunay_xy_for_index, delaunay_p1)
        p2_row_index, p2_column_index = find_element_index(delaunay_xy_for_index, delaunay_p2)

        if None in [p0_row_index, p0_column_index, p1_row_index, p1_column_index, p2_row_index, p2_column_index]:
            continue  # Skip if any index is not found

        p0_index = layer_points_id[p0_row_index][p0_column_index]
        p1_index = layer_points_id[p1_row_index][p1_column_index]
        p2_index = layer_points_id[p2_row_index][p2_column_index]

        tri_mesh.append([p0_index, p1_index, p2_index])

    return xyz, tri_mesh
