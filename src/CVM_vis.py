import os
import sys
import numpy as np
import pandas as pd
import vtk
from scipy.interpolate import RegularGridInterpolator
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

sys.path.append(r'src')
import CONST    # User-defined constants module



def main():
    """
    Main function to orchestrate the velocity model processing.
    """
    # Coordinate steps
    xyz_step = np.array(  [0.1, 0.1, 1.0]  )

    # Color Transfer Functions' Parameters
    prec0 = 0.254
    prec1 = 0.238
    y0 = -10.0
    y3 = 10.0

    # Input and output paths
    input_rawData_file = r'rawData\SWChinaCVMv2.0.txt'
    output_path = r'resultData\CVM'
    os.makedirs(output_path, exist_ok=True)

    # Step 1: Read raw data and create interpolators of Vp, Vs
    X, Y, Z, vp_interpolator, vs_interpolator = read_and_interpolate_data(input_rawData_file)

    # Step 2: Compute gradients
    grads_df = compute_gradients(X, Y, Z, vp_interpolator, vs_interpolator, xyz_step)

    # Step 3: Compute dVp, dVs, and dVp_Vs and write to CSV
    compute_dvp_dvs(y0, y3, prec0, prec1, grads_df, output_path, xyz_step)

    # Step 4: Create VTK model and visualize the CVM
    create_vtk_model_and_visualize(output_path, xyz_step)



def read_and_interpolate_data(data_file):
    """
    Reads the raw data from a text file and creates interpolators for Vp and Vs.

    Returns:
        X, Y, Z (numpy.ndarray): Sorted unique coordinate values.
        vp_interpolator, vs_interpolator: Interpolator objects for Vp and Vs.
    """
    print("Step 1: Reading raw data and creating interpolators ...")

    # Use pandas to read the data efficiently
    df = pd.read_csv(data_file, skiprows=1, delim_whitespace=True, names=['lon', 'lat', 'dep', 'Vp', 'Vs'])
    data = df.values
    print(f"        - Read {len(data)} data points;")

    # Generate sorted arrays of unique X, Y, Z values
    X = np.unique(data[:, 0])
    Y = np.unique(data[:, 1])
    Z = np.unique(data[:, 2])
    print(f"        - Grid dimensions of the Raw Data: X({len(X)}), Y({len(Y)}), Z({len(Z)});")

    # Initialize Vp and Vs arrays
    Vp = np.full(  (len(X), len(Y), len(Z)), np.nan  )
    Vs = np.full(  (len(X), len(Y), len(Z)), np.nan  )

    # Map indices for X, Y, Z
    x_idx = {x: i for i, x in enumerate(X)}
    y_idy = {y: i for i, y in enumerate(Y)}
    z_idz = {z: i for i, z in enumerate(Z)}

    # Fill Vp and Vs arrays
    for row in data:
        x, y, z, vp_val, vs_val = row
        ix = x_idx[x]
        iy = y_idy[y]
        iz = z_idz[z]
        Vp[ix, iy, iz] = vp_val
        Vs[ix, iy, iz] = vs_val

    # Create RegularGridInterpolators
    vp_interpolator = RegularGridInterpolator(  (X, Y, Z), Vp, bounds_error=False, fill_value=np.nan  )
    vs_interpolator = RegularGridInterpolator(  (X, Y, Z), Vs, bounds_error=False, fill_value=np.nan  )

    return X, Y, Z, vp_interpolator, vs_interpolator


def compute_gradients(X, Y, Z, vp_interpolator, vs_interpolator, xyz_step):
    """
    Computes the gradients of Vp and Vs, and returns the results as a DataFrame.

    Parameters:
        X, Y, Z (numpy.ndarray): Coordinate values.
        vp_interpolator, vs_interpolator: Interpolator functions for Vp and Vs.
        xyz_step (numpy.ndarray): Coordinate steps.

    Returns:
        pandas.DataFrame: DataFrame containing the computed results.
    """
    print("Step 2: Computing gradients ...")

    # Define the coordinate ranges for results
    result_X = np.round(np.arange(X[0] + xyz_step[0], X[-1] - xyz_step[0]*0.5, xyz_step[0]), 5)
    result_Y = np.round(np.arange(Y[0] + xyz_step[1], Y[-1] - xyz_step[1]*0.5, xyz_step[1]), 5)
    result_Z = np.round(np.arange(Z[0] + xyz_step[2], Z[-1] - xyz_step[2]*0.5, xyz_step[2]), 5)

    # Generate all combinations of coordinates
    grid_points = np.array(np.meshgrid(result_X, result_Y, result_Z, indexing='ij')).reshape(3, -1).T
    print(f"        - Grid dimensions of the Interpolator Data: Interpolator_X({len(result_X)}), Interpolator_Y({len(result_Y)}), Interpolator_Z({len(result_Z)});")

    # Prepare arrays to store results
    lon_list = []
    lat_list = []
    dep_list = []
    vp_values = []
    vs_values = []
    vp_grad_norms = []
    vs_grad_norms = []
    vp_divide_vs_norms = []

    # Process in batches to optimize memory usage
    batch_size = 10000  # Adjust based on memory capacity
    total_points = grid_points.shape[0]

    for start_idx in tqdm(range(0, total_points, batch_size), desc="        - Processing gradients"):
        end_idx = min(start_idx + batch_size, total_points)
        batch_points = grid_points[start_idx:end_idx]

        # Compute central differences
        # Create shifted points for gradient computation
        points_x_plus = batch_points + [xyz_step[0], 0, 0]
        points_x_minus = batch_points - [xyz_step[0], 0, 0]
        points_y_plus = batch_points + [0, xyz_step[1], 0]
        points_y_minus = batch_points - [0, xyz_step[1], 0]

        # Interpolate Vp and Vs values at required points
        vp_center = vp_interpolator(batch_points)
        vs_center = vs_interpolator(batch_points)

        vp_x_plus = vp_interpolator(points_x_plus)
        vp_x_minus = vp_interpolator(points_x_minus)
        vp_y_plus = vp_interpolator(points_y_plus)
        vp_y_minus = vp_interpolator(points_y_minus)

        vs_x_plus = vs_interpolator(points_x_plus)
        vs_x_minus = vs_interpolator(points_x_minus)
        vs_y_plus = vs_interpolator(points_y_plus)
        vs_y_minus = vs_interpolator(points_y_minus)

        # Check for NaN values
        invalid_mask = np.isnan(vp_center) | np.isnan(vs_center) | \
                       np.isnan(vp_x_plus) | np.isnan(vp_x_minus) | \
                       np.isnan(vp_y_plus) | np.isnan(vp_y_minus) | \
                       np.isnan(vs_x_plus) | np.isnan(vs_x_minus) | \
                       np.isnan(vs_y_plus) | np.isnan(vs_y_minus)

        # Filter valid data
        valid_indices = np.where(~invalid_mask)[0]

        if valid_indices.size > 0:
            vp_center = vp_center[valid_indices]
            vs_center = vs_center[valid_indices]
            vp_x_plus = vp_x_plus[valid_indices]
            vp_x_minus = vp_x_minus[valid_indices]
            vp_y_plus = vp_y_plus[valid_indices]
            vp_y_minus = vp_y_minus[valid_indices]
            vs_x_plus = vs_x_plus[valid_indices]
            vs_x_minus = vs_x_minus[valid_indices]
            vs_y_plus = vs_y_plus[valid_indices]
            vs_y_minus = vs_y_minus[valid_indices]
            valid_points = batch_points[valid_indices]

            # Compute gradients
            vp_grad_x = (vp_x_plus - vp_x_minus) / (2 * xyz_step[0])
            vp_grad_y = (vp_y_plus - vp_y_minus) / (2 * xyz_step[1])
            vp_grad_norm = np.sqrt(vp_grad_x ** 2 + vp_grad_y ** 2)

            vs_grad_x = (vs_x_plus - vs_x_minus) / (2 * xyz_step[0])
            vs_grad_y = (vs_y_plus - vs_y_minus) / (2 * xyz_step[1])
            vs_grad_norm = np.sqrt(vs_grad_x ** 2 + vs_grad_y ** 2)

            # Compute Vp/Vs and its gradient
            with np.errstate(divide='ignore', invalid='ignore'):
                vp_vs_ratio = vp_center / vs_center
                vp_vs_x_plus = vp_x_plus / vs_x_plus
                vp_vs_x_minus = vp_x_minus / vs_x_minus
                vp_vs_y_plus = vp_y_plus / vs_y_plus
                vp_vs_y_minus = vp_y_minus / vs_y_minus

                vp_vs_grad_x = (vp_vs_x_plus - vp_vs_x_minus) / (2 * xyz_step[0])
                vp_vs_grad_y = (vp_vs_y_plus - vp_vs_y_minus) / (2 * xyz_step[1])
                vp_vs_grad_norm = np.sqrt(vp_vs_grad_x ** 2 + vp_vs_grad_y ** 2)

            # Filter out any remaining NaN or Inf values after calculations
            final_invalid_mask = np.isnan(vp_grad_norm) | np.isnan(vs_grad_norm) | \
                                 np.isnan(vp_vs_grad_norm) | np.isinf(vp_grad_norm) | \
                                 np.isinf(vs_grad_norm) | np.isinf(vp_vs_grad_norm)

            if np.any(final_invalid_mask):
                valid_indices = np.where(~final_invalid_mask)[0]
                valid_points = valid_points[valid_indices]
                vp_center = vp_center[valid_indices]
                vs_center = vs_center[valid_indices]
                vp_grad_norm = vp_grad_norm[valid_indices]
                vs_grad_norm = vs_grad_norm[valid_indices]
                vp_vs_grad_norm = vp_vs_grad_norm[valid_indices]

            # Append results
            lon_list.extend(valid_points[:, 0])
            lat_list.extend(valid_points[:, 1])
            dep_list.extend(valid_points[:, 2])
            vp_values.extend(vp_center)
            vs_values.extend(vs_center)
            vp_grad_norms.extend(vp_grad_norm)
            vs_grad_norms.extend(vs_grad_norm)
            vp_divide_vs_norms.extend(vp_vs_grad_norm)

    # Prepare DataFrame
    data = {
        'lon': lon_list,
        'lat': lat_list,
        'dep': dep_list,
        'Vp': vp_values,
        'Vs': vs_values,
        'vp_grad_norm': vp_grad_norms,
        'vs_grad_norm': vs_grad_norms,
        'vp_divide_vs_norm': vp_divide_vs_norms
    }

    results_df = pd.DataFrame(data)

    # Round numerical columns to 12 decimal places
    numeric_cols = ['lon', 'lat', 'dep', 'Vp', 'Vs', 'vp_grad_norm', 'vs_grad_norm', 'vp_divide_vs_norm']
    results_df[numeric_cols] = results_df[numeric_cols].round(12)
    
    return results_df


def compute_dvp_dvs(y0, y3, prec0, prec1, df, output_path, xyz_step):
    """
    Processes the DataFrame to compute dVp, dVs, and dVp_Vs, and writes the results to a new CSV file.
    """
    print("Step 3: Computing dVp, dVs, and dVp_Vs ...")

    y1 = y0 + prec0 * (y3 - y0)
    y2 = y3 - prec1 * (y3 - y0)

    def trans_functions(x, d):
        """
        Vectorized transformation function.
        """
        res = np.empty_like(x)
        mask1 = x < d[1]
        mask2 = x > d[2]
        mask3 = ~(mask1 | mask2)
        res[mask1] = ((y1 - y0) / (d[1] - d[0])) * (x[mask1] - d[0]) + y0
        res[mask2] = ((y3 - y2) / (d[3] - d[2])) * (x[mask2] - d[2]) + y2
        res[mask3] = ((y2 - y1) / (d[2] - d[1])) * (x[mask3] - d[1]) + y1
        return res

    def compute_d(x_min, x_max):
        D0 = (x_max - x_min) * prec0
        D1 = (x_max - x_min) * prec1

        d0 = x_min
        d1 = x_min + D0
        d2 = x_max - D1
        d3 = x_max
        return [d0, d1, d2, d3]

    # Compute Vp/Vs ratio
    df['Vp/Vs'] = df['Vp'] / df['Vs']

    # Group by 'dep' and process in parallel
    dep_groups = list(df.groupby('dep'))

    def process_dep_group(dep_group):
        dep_value, data0 = dep_group

        Vp_vals = data0['Vp'].values
        Vs_vals = data0['Vs'].values
        vp_vs_ratio = data0['Vp/Vs'].values

        vp_d = compute_d(Vp_vals.min(), Vp_vals.max())
        vs_d = compute_d(Vs_vals.min(), Vs_vals.max())
        vp_vs_d = compute_d(vp_vs_ratio.min(), vp_vs_ratio.max())

        # Apply transformation functions
        data0 = data0.copy()
        data0['dVp'] = trans_functions(Vp_vals, vp_d)
        data0['dVs'] = trans_functions(Vs_vals, vs_d)
        data0['dVp_Vs'] = trans_functions(vp_vs_ratio, vp_vs_d)

        return data0

    # Use ThreadPoolExecutor to parallelize processing
    with ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(process_dep_group, dep_groups), total=len(dep_groups), desc="        - Processing transfer function"))

    # Concatenate results
    final_df = pd.concat(results, ignore_index=True)

    # Round numerical columns to 12 decimal places
    numeric_cols = ['Vp', 'Vs', 'Vp/Vs', 'dVp', 'dVs', 'dVp_Vs']
    final_df[numeric_cols] = final_df[numeric_cols].round(12)

    output_csv_file = os.path.join(output_path, f'compute_grads_dVs_dVp_step{xyz_step[-1]}.csv')
    final_df.to_csv(output_csv_file, index=False)
    print(f"        - Computed dVp, dVs, dVp_Vs and saved to {output_csv_file}.")


def create_vtk_model_and_visualize(output_path, xyz_step):
    """
    Creates a VTK model from the computed data and visualizes it.
    """
    print("Step 4: Creating VTK model and visualizing ...")
    input_file = os.path.join(output_path, f'compute_grads_dVs_dVp_step{xyz_step[-1]}.csv')
    df = pd.read_csv(input_file)

    # Initialize vtk variables
    points = vtk.vtkPoints()
    vtk_arrays = initialize_vtk_arrays()

    # Insert data into VTK arrays
    for idx, row in df.iterrows():
        insert_vtk_point_and_data(points, vtk_arrays, row)

    # Create the unstructured grid with tetrahedral cells
    ugrid = create_unstructured_grid(points, vtk_arrays, df)

    # Write the unstructured grid to a VTK file
    output_vtk_file = os.path.join(output_path, f'compute_grads_dVs_dVp_step{xyz_step[-1]}_X{CONST.Z_SCALE}.vtu')
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(output_vtk_file)
    writer.Write()

    # Visualization (Optional)
    visualize_vtk_model(ugrid)
    
    print(f"        - Visualize the VTK model, and saved to {output_vtk_file}.")


def initialize_vtk_arrays():
    """
    Initializes and returns a dictionary of vtkFloatArray for storing data.
    """
    array_names = [
        ('Vp', 'Vp'),
        ('Vs', 'Vs'),
        ('vp_grad_norm', 'vp_grad_norm'),
        ('vs_grad_norm', 'vs_grad_norm'),
        ('Vp/Vs', 'Vp/Vs'),
        ('vp_divide_vs_norm', 'vp_divide_vs_norm'),
        ('dVp', 'dVp'),
        ('dVs', 'dVs'),
        ('dVp_Vs', 'dVp_Vs')
    ]
    vtk_arrays = {}
    for var_name, array_name in array_names:
        vtk_array = vtk.vtkFloatArray()
        vtk_array.SetNumberOfComponents(1)
        vtk_array.SetName(var_name)
        vtk_arrays[array_name] = vtk_array
    return vtk_arrays


def insert_vtk_point_and_data(points, vtk_arrays, row):
    """
    Inserts a point and associated data into the VTK structures.

    Args:
        points (vtk.vtkPoints): VTK points object.
        vtk_arrays (dict): Dictionary of VTK arrays.
        row (pandas.Series): Data row containing point coordinates and associated values.
    """
    lon, lat, dep = row['lon'], row['lat'], row['dep']
    Vp_val, Vs_val = row['Vp'], row['Vs']
    vp_grad_norm_val, vs_grad_norm_val = row['vp_grad_norm'], row['vs_grad_norm']
    vp_divide_vs_val = row['Vp/Vs']
    vp_divide_vs_norm_val = row['vp_divide_vs_norm']
    dVp_val, dVs_val, dVp_Vs_val = row['dVp'], row['dVs'], row['dVp_Vs']

    points.InsertNextPoint(lon, lat, -1.0 * dep / CONST.angle_to_kilometers)

    vtk_arrays['Vp'].InsertNextValue(Vp_val)
    vtk_arrays['Vs'].InsertNextValue(Vs_val)
    vtk_arrays['vp_grad_norm'].InsertNextValue(vp_grad_norm_val)
    vtk_arrays['vs_grad_norm'].InsertNextValue(vs_grad_norm_val)
    vtk_arrays['Vp/Vs'].InsertNextValue(vp_divide_vs_val)
    vtk_arrays['vp_divide_vs_norm'].InsertNextValue(vp_divide_vs_norm_val)
    vtk_arrays['dVp'].InsertNextValue(dVp_val)
    vtk_arrays['dVs'].InsertNextValue(dVs_val)
    vtk_arrays['dVp_Vs'].InsertNextValue(dVp_Vs_val)


def create_unstructured_grid(points, vtk_arrays, df):
    """
    Creates a vtkUnstructuredGrid with tetrahedral cells.

    Args:
        points (vtk.vtkPoints): VTK points object.
        vtk_arrays (dict): Dictionary of VTK arrays.
        df (pandas.DataFrame): DataFrame used to construct the grid.

    Returns:
        vtk.vtkUnstructuredGrid: The constructed unstructured grid.
    """
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    # Add point data arrays
    for vtk_array in vtk_arrays.values():
        ugrid.GetPointData().AddArray(vtk_array)

    # Prepare for constructing tetrahedral cells
    xyz = df[['lon', 'lat', 'dep']].values
    z_values = np.unique(xyz[:, 2])

    # Get the points at the top layer
    xy = xyz[xyz[:, 2] == z_values[0]][:, 0:2]
    triangle = CONST.triMesh_oneLayer(xy)
    xy_len = len(xy)
    total_layers = len(z_values)

    # Construct tetrahedral cells between layers
    for k in range(total_layers - 1):
        for tri in triangle:
            id_up = (tri + k * xy_len).astype(int)
            id_down = id_up + xy_len

            face_id0 = [id_up[0], id_up[1], id_up[2], id_down[0]]
            face_id1 = [id_up[1], id_up[2], id_down[0], id_down[1]]
            face_id2 = [id_up[2], id_down[0], id_down[1], id_down[2]]

            cell0 = vtk.vtkTetra()
            cell1 = vtk.vtkTetra()
            cell2 = vtk.vtkTetra()

            for i in range(4):
                cell0.GetPointIds().SetId(i, face_id0[i])
                cell1.GetPointIds().SetId(i, face_id1[i])
                cell2.GetPointIds().SetId(i, face_id2[i])

            ugrid.InsertNextCell(cell0.GetCellType(), cell0.GetPointIds())
            ugrid.InsertNextCell(cell1.GetCellType(), cell1.GetPointIds())
            ugrid.InsertNextCell(cell2.GetCellType(), cell2.GetPointIds())

    return ugrid


def visualize_vtk_model(ugrid):
    """
    Visualizes the VTK unstructured grid.

    Args:
        ugrid (vtk.vtkUnstructuredGrid): The unstructured grid to visualize.
    """
    colors = vtk.vtkNamedColors()
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(ugrid)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(colors.GetColor3d('Silver'))

    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.SetWindowName('Velocity Model Visualization')
    render_window.AddRenderer(renderer)
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('Salmon'))
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(30)
    renderer.GetActiveCamera().Elevation(30)
    render_window.Render()
    render_window_interactor.Start()


if __name__ == "__main__":
    main()
