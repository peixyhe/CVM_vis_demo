import pandas as pd
import numpy as np
import os
from scipy.spatial import cKDTree

import sys
sys.path.append(r"src")
import CONST



def inverse_distance_weighte(p0, xyz):
    """
    Calculate the new point position using inverse distance weighting.
    
    Args:
        p0 (np.ndarray): Original point coordinates.
        xyz (np.ndarray): Neighboring points coordinates.
    
    Returns:
        np.ndarray: New point coordinates after weighting.
    """
    vector_diff = []
    weighte = []
    for xyz0 in xyz:
        if not np.array_equal(xyz0, p0):
            diff = xyz0 - p0
            vector_diff.append(diff)
            weighte.append(1.0 / np.linalg.norm(diff))
    
    weighte = np.array(weighte) / sum(weighte)
    for weighte0, diff0 in zip(weighte, vector_diff):
        p0 += weighte0 * diff0
    
    return np.round(p0, 6)


def remove_dense_points(points, threshold):
    """
    Remove overly dense points.
    
    Args:
        points (np.ndarray): Array of point coordinates.
        threshold (float): Distance threshold to consider points as dense.
    
    Returns:
        np.ndarray: Array of points after removing dense points.
    """
    p_kdtree = cKDTree(points)
    keep_p = []
    
    for p in points:
        indices = p_kdtree.query_ball_point(p, r=threshold)
        
        if len(indices) == 1:
            keep_p.append(p)
        elif len(indices) == 0:
            continue
        else:
            # Apply inverse distance weighting to compute new point
            keep_p.append(inverse_distance_weighte(p, points[indices]))
            
            # Remove the points that have been processed
            points = np.delete(points, indices, axis=0)
            p_kdtree = cKDTree(points)
    
    return np.array(keep_p)


def preprocess_fault_data(input_file, output_root_path):
    """
    Preprocess fault data by converting raw data to CSV, splitting into separate files per fault,
    and removing dense points.
    
    Args:
        input_file (str): Path to the raw input data file.
        output_root_path (str): Directory to save the processed data.
    
    Returns:
        None
    """
    faults_name = set()
    id_num = set()

    # Read and convert raw data to CSV
    with open(input_file, 'r') as r:
        lines = r.readlines()
        with open(os.path.normpath(os.path.join(output_root_path, "CSES_3DFM_V2.csv")), 'w') as w:
            w.write("x,y,z,fault_name,id\n")
            for line in lines:
                lineData = line.split()
                
                z = lineData[-1]
                y = lineData[-2]
                x = lineData[-3]
                name_id = lineData[-4]
                name = ' '.join(lineData[0:-4])
                
                data = [x, y, z, name, name_id]
                faults_name.add(name)
                id_num.add(name_id)
                
                w.write(','.join(data) + '\n')

    # Define output directories
    faults_rawData_txt_path = os.path.join(output_root_path, 'faults_rawData_txt')
    os.makedirs(faults_rawData_txt_path, exist_ok=True)
    
    faults_noDenseDataFile_path = os.path.join(output_root_path, "faults_noDenseData_txt")
    os.makedirs(faults_noDenseDataFile_path, exist_ok=True)

    # Read the consolidated CSV file
    df = pd.read_csv(os.path.normpath(os.path.join(output_root_path, "CSES_3DFM_V2.csv")))
    df_data = df[['x', 'y', 'z', 'fault_name', 'id']].values

    # Split data into separate files per fault
    for fault in faults_name:
        with open(os.path.join(faults_rawData_txt_path, fault + '.txt'), 'w') as w2:
            w2.write("x,y,z,fault_name,id\n")
            fault_data = df_data[df_data[:, -2] == fault]
            for data_w in fault_data:
                w2.write(','.join(map(str, data_w)) + '\n')

    # Remove dense points and save the cleaned data
    for file in os.listdir(faults_rawData_txt_path):
        # Read point cloud data
        data = pd.read_csv(os.path.join(faults_rawData_txt_path, file))
        xyz0 = data[['x', 'y', 'z']].values
        
        xyz = remove_dense_points(xyz0, 1250.0)
        unDenseData = pd.DataFrame(xyz, columns=['x', 'y', 'z'])
        unDenseData.to_csv(os.path.join(faults_noDenseDataFile_path, file), index=False)
