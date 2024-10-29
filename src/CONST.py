"""
Constants and Common Functions Definition Script

"""

import math
import numpy as np
from scipy.spatial import KDTree



# Earth's radius in kilometers and meters
EARTH_RADIUS_KM = 6_371.393     # in kilometers
EARTH_RADIUS_M = 6_371_393.0    # in meters

# Geographic region boundaries (in degrees)
# Latitude ranges from -90 to 90 degrees
latitude_up = 33.4      # Northern boundary latitude
latitude_down = 21.8    # Southern boundary latitude

# Longitude ranges from -180 to 180 degrees
longitude_left = 97.8   # Western boundary longitude
longitude_right = 107.0 # Eastern boundary longitude

# Other constants
Z_SCALE = 10            # Scaling factor for Z-axis
MIN_MAGNITUDE = 1       # Minimum magnitude for seismic events



# Calculate the average latitude of the region
average_latitude = (latitude_up + latitude_down) / 2.0

# Calculate the conversion factor from degrees to kilometers
# This formula computes the distance corresponding to one degree difference in coordinates,
# adjusted by the scaling factor and Earth's curvature at the given latitude.
# It accounts for the fact that the length of a degree of longitude varies with latitude.

# Conversion factor from degrees to kilometers
angle_to_kilometers = (
    (2 * math.pi * EARTH_RADIUS_KM) / 360.0 / Z_SCALE
) * (
    math.sqrt(1 + math.cos(math.radians(average_latitude)) ** 2) / math.sqrt(2)
)

# Conversion factor from degrees to meters
angle_to_meters = (
    (2 * math.pi * EARTH_RADIUS_M) / 360.0 / Z_SCALE
) * (
    math.sqrt(1 + math.cos(math.radians(average_latitude)) ** 2) / math.sqrt(2)
)

# print(f"One degree is equal to {angle_to_kilometers} kilometers")
# print(f"Z-axis data scale: {-1.0 / angle_to_kilometers}\n")



def Get_mag2radius(mag0):
    # Convert magnitude to radius using a power function
    return round(math.pow(2.0, (mag0 - 4.56) / 1.96), 6)

def geodetic_to_ecef_m(lon_deg, lat_deg, depth_m):
    # WGS84 ellipsoid parameters
    a = 6378137.0  # Semi-major axis, in meters
    f = 1 / 298.257223563  # Flattening factor
    e2 = 2 * f - f ** 2  # Square of the first eccentricity

    # Convert angles from degrees to radians
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)

    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)

    # Calculate the radius of curvature in the prime vertical N
    N = a / math.sqrt(1 - e2 * sin_lat ** 2)

    # Calculate ECEF coordinates (x_m, y_m, z_m)
    x_m = (N + depth_m) * cos_lat * cos_lon
    y_m = (N + depth_m) * cos_lat * sin_lon
    z_m = (N * (1 - e2) + depth_m) * sin_lat

    return x_m, y_m, z_m

def geodetic_to_ecef_km(lon_deg, lat_deg, depth_m):
    # WGS84 ellipsoid parameters
    a = 6378137.0  # Semi-major axis, in meters
    f = 1 / 298.257223563  # Flattening factor
    e2 = 2 * f - f ** 2  # Square of the first eccentricity

    # Convert angles from degrees to radians
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)

    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)

    # Calculate the radius of curvature in the prime vertical N
    N = a / math.sqrt(1 - e2 * sin_lat ** 2)

    # Calculate ECEF coordinates (x_m, y_m, z_m)
    x_m = (N + depth_m) * cos_lat * cos_lon
    y_m = (N + depth_m) * cos_lat * sin_lon
    z_m = (N * (1 - e2) + depth_m) * sin_lat

    return x_m / 1000.0, y_m / 1000.0, z_m / 1000.0

def triMesh_oneLayer(xy):
    """
    Generates a triangular mesh for a single layer based on input coordinates.

    Parameters:
    xy (numpy.ndarray): A 2D array of shape (n_points, 2) containing x and y coordinates.

    Returns:
    numpy.ndarray: An array of triangles represented by indices of the points in 'xy'.
    """
    # Initialize an empty list to store the triangles
    tri_mesh = []

    # Round coordinates to 5 decimal places to mitigate floating point errors
    xy = np.round(xy, 5)

    # Build a KDTree for efficient nearest-neighbor searches
    xy_kdtree = KDTree(xy)

    # Extract sorted unique x and y coordinates
    x_coords = sorted(set(xy[:, 0]))
    y_coords = sorted(set(xy[:, 1]))

    # Calculate the grid spacing in x and y directions
    dx = round(x_coords[1] - x_coords[0], 5)
    dy = round(y_coords[1] - y_coords[0], 5)

    # Iterate over each point in 'xy'
    for i in range(len(xy)):
        xy0 = xy[i]

        # Define neighboring points to form triangles
        xy1 = np.round(xy0 + np.array([dx, 0.0]), 5)
        xy2 = np.round(xy0 + np.array([dx, dy]), 5)
        xy3 = np.round(xy0 + np.array([0.0, dy]), 5)

        # Find the indices of the closest points in 'xy' using KDTree
        _, id1 = xy_kdtree.query(xy1)
        _, id2 = xy_kdtree.query(xy2)
        _, id3 = xy_kdtree.query(xy3)

        closest_point1 = xy[id1]
        closest_point2 = xy[id2]
        closest_point3 = xy[id3]

        # Check if the closest points match the expected grid points
        if (np.array_equal(closest_point1, xy1) and
            np.array_equal(closest_point2, xy2) and
            np.array_equal(closest_point3, xy3)):

            # Retrieve the indices of the neighboring points
            id0 = i
            id1_array = np.where((xy == xy1).all(axis=1))[0]
            id2_array = np.where((xy == xy2).all(axis=1))[0]
            id3_array = np.where((xy == xy3).all(axis=1))[0]

            # Squeeze the arrays to get scalar indices
            id1 = np.squeeze(id1_array)
            id2 = np.squeeze(id2_array)
            id3 = np.squeeze(id3_array)

            # Add two triangles to the mesh
            tri_mesh.append([id0, id1, id2])
            tri_mesh.append([id0, id2, id3])

    return np.array(tri_mesh)
