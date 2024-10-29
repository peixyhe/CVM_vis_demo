"""
This script reads elevation data from a NetCDF file, extracts a subset within specified latitude and longitude ranges,
and saves the data to both CSV and VTI files. The elevation data is shifted along the z-axis by a specified amount.

Dependencies:
- xarray
- numpy
- vtk
- CONST module containing the following constants:
    - latitude_down
    - latitude_up
    - longitude_left
    - longitude_right
"""

import os
import xarray as xr
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

import CONST  # Assume CONST module contains necessary constants



def main():
    z_Up = 1.5  # Amount to shift along the z-axis

    # Read the dataset from the NetCDF file
    data = xr.open_dataset(r'rawData\ETOPO_2022_v1_30s_N90W180_surface.nc')

    # Select data within specified latitude and longitude ranges using xarray
    subset = data.sel(
        lat=slice(CONST.latitude_down, CONST.latitude_up),
        lon=slice(CONST.longitude_left, CONST.longitude_right)
    )

    # Extract elevation, latitude, and longitude data
    z = subset.z.values  # Elevation data
    lat = subset.lat.values
    lon = subset.lon.values

    # Ensure data types are double precision to preserve accuracy
    z = z.astype(np.float64)
    lat = lat.astype(np.float64)
    lon = lon.astype(np.float64)

    # Create output directory if it doesn't exist
    output_dem_path = r'resultData\DEM_30s'
    if not os.path.exists(output_dem_path):
        os.makedirs(output_dem_path)

    # Define output CSV file path
    output_csv_file = os.path.join(
        output_dem_path,
        f"dem_lat{CONST.latitude_down}-{CONST.latitude_up}_"
        f"lon{CONST.longitude_left}-{CONST.longitude_right}.csv"
    )

    # Generate meshgrid of longitude and latitude values
    lon_grid, lat_grid = np.meshgrid(lon, lat)

    # Flatten the data arrays to 1D arrays in row-major order
    xyz = np.column_stack(
        (lon_grid.ravel(order='C'), lat_grid.ravel(order='C'), z.ravel(order='C'))
    )

    # Write data to CSV file with high precision
    with open(output_csv_file, 'w') as w:
        w.write('lon,lat,elevation\n')
        np.savetxt(w, xyz, delimiter=',', fmt='%.12g,%.12g,%.12g')

    # Define output VTI file path
    output_vti_file = os.path.join(
        output_dem_path,
        f"dem_lat{CONST.latitude_down}-{CONST.latitude_up}_"
        f"lon{CONST.longitude_left}-{CONST.longitude_right}_zUP{z_Up}.vti"
    )

    # Create VTI file with elevation data
    create_vti_file(lon, lat, z, z_Up, output_vti_file)


def create_vti_file(lon, lat, z, z_Up, output_vti_filename):
    """
    Write elevation data to a VTI format file.

    Parameters:
    - lon: Longitude array
    - lat: Latitude array
    - z: Elevation data array
    - z_Up: Amount to shift along the z-axis
    - output_vti_filename: Output VTI file name
    """
    nx = len(lon)
    ny = len(lat)
    nz = 1  # Number of points in z-direction (since data is 2D)

    # Create vtkImageData object to store image data
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(nx, ny, nz)

    # Set the origin of the data (minimum longitude, minimum latitude, z_Up)
    origin = (lon.min(), lat.min(), z_Up)
    imageData.SetOrigin(origin)

    # Calculate spacing between points
    if nx > 1:
        dx = (lon.max() - lon.min()) / (nx - 1)
    else:
        dx = 1.0  # Default spacing if only one point in x

    if ny > 1:
        dy = (lat.max() - lat.min()) / (ny - 1)
    else:
        dy = 1.0  # Default spacing if only one point in y

    dz = 1.0  # Spacing in z-direction
    imageData.SetSpacing(dx, dy, dz)

    # Flatten elevation data to a 1D array in row-major order
    elevation_flat = z.ravel(order='C')

    # Convert numpy array to VTK array, using double precision to preserve accuracy
    elevation_vtk = numpy_to_vtk(elevation_flat, deep=True, array_type=vtk.VTK_FLOAT)
    elevation_vtk.SetName("Elevation")

    # Add elevation data to vtkImageData's point data
    imageData.GetPointData().SetScalars(elevation_vtk)

    # Write image data to VTI file
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(output_vti_filename)
    writer.SetInputData(imageData)
    writer.Write()


if __name__ == "__main__":
    main()
