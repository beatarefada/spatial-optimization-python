"""
Spatial Optimization & Coordinate Conversion
------------------------------------------------
Project: Residential Location Choice Model
Author: Betania Alves
Description: 
    This script performs a spatial optimization analysis to find the optimal 
    location for a resident based on proximity to specific amenities.
    
    It includes:
    1. Transformation of global coordinates (Lat/Lon) to local Euclidean (km).
    2. Unconstrained optimization of a utility function.
    3. Constrained optimization (Lagrange Multipliers) along a street network.
    4. Inverse transformation of results back to global coordinates.
"""

from math import pi, cos
from sympy import symbols, diff, solve

# ------------------------------------------------------------------------------
# 1. COORDINATE CONVERSION UTILITIES
# ------------------------------------------------------------------------------

def deg_to_kms(point_long, point_lat, origin_long, origin_lat):
    """
    Converts global coordinates of a point (point_lat, point_long) measured in degrees 
    to local coordinates based at point (origin_long, origin_lat) measured in kilometers.
    """
    earth_radius = 6371.0
    
    # Convert longitudes and latitudes from degrees to radians
    long = point_long * pi / 180
    reference_long = origin_long * pi / 180
    lat = point_lat * pi / 180
    reference_lat = origin_lat * pi / 180
    
    # Calculate the difference in latitude and longitude in radians
    long_difference_rad = long - reference_long
    lat_difference_rad = lat - reference_lat
    
    # Calculate the signed distance using the Earth's meridian length formula
    # (Equirectangular projection approximation for small distances)
    diff_long = earth_radius * long_difference_rad * cos(reference_lat)
    diff_lat = earth_radius * lat_difference_rad
    
    return diff_long, diff_lat

def kms_to_deg(local_x, local_y, origin_long, origin_lat):
    """
    Converts local coordinates of a point (local_x, local_y) measured in kms 
    from the origin into global coordinates measured in degrees.
    """
    earth_radius = 6371.0
    
    # Convert origin latitude to radians for the longitude calculation
    reference_lat_rad = origin_lat * pi / 180
    
    # Calculate difference in radians (Inverse operations of the function above)
    long_difference_rad = local_x / (earth_radius * cos(reference_lat_rad))
    lat_difference_rad = local_y / earth_radius
    
    # Convert differences back to degrees
    diff_long_deg = long_difference_rad * 180 / pi
    diff_lat_deg = lat_difference_rad * 180 / pi
    
    # Add differences to origin to get final coordinates
    point_long = origin_long + diff_long_deg
    point_lat = origin_lat + diff_lat_deg
    
    return point_long, point_lat

# ------------------------------------------------------------------------------
# 2. DATA SETUP
# ------------------------------------------------------------------------------

# Define the Origin Point (Reference for local system)
ORIGIN_LONG = -58.37788955179407
ORIGIN_LAT = -34.595228892628455

# Coordinates of Amenities (Global Degrees)
# di: Disco, ob: Obelisk, ga: Galerias Pacifico (Examples)
di_global = (-58.38467378144673, -34.596156182566006)
ob_global = (-58.38094967105265, -34.6034559421601)
ga_global = (-58.37401050128944, -34.598868856938026)

# Coordinates of Street Constraints (Paraguay Street endpoints)
ho_global = (-58.38358725902865, -34.598181576896955)
res_global = (-58.38026132000657, -34.597792990501425)

# Convert all points to Local Coordinates (km)
x1, y1 = deg_to_kms(di_global[0], di_global[1], ORIGIN_LONG, ORIGIN_LAT)
x2, y2 = deg_to_kms(ob_global[0], ob_global[1], ORIGIN_LONG, ORIGIN_LAT)
x3, y3 = deg_to_kms(ga_global[0], ga_global[1], ORIGIN_LONG, ORIGIN_LAT)

ho_x, ho_y = deg_to_kms(ho_global[0], ho_global[1], ORIGIN_LONG, ORIGIN_LAT)
res_x, res_y = deg_to_kms(res_global[0], res_global[1], ORIGIN_LONG, ORIGIN_LAT)

print(f"Local Coords - Amenity 1: ({x1:.4f}, {y1:.4f})")
print(f"Local Coords - Amenity 2: ({x2:.4f}, {y2:.4f})")
print(f"Local Coords - Amenity 3: ({x3:.4f}, {y3:.4f})")

# ------------------------------------------------------------------------------
# 3. OPTIMIZATION MODEL (SymPy)
# ------------------------------------------------------------------------------

# Define symbolic variables
x, y, L = symbols('x y L')

# Define Utility Function U(x,y)
# Minimizing weighted squared Euclidean distance to amenities
# Weights: 1 for Amenity 1, 2 for Amenity 2, 3 for Amenity 3
U = 1*((x - x1)**2 + (y - y1)**2) + \
    2*((x - x2)**2 + (y - y2)**2) + \
    3*((x - x3)**2 + (y - y3)**2)

# --- PART A: Unconstrained Optimization ---
print("\nSolving Unconstrained Optimization...")
eq_1 = diff(U, x)
eq_2 = diff(U, y)

# Solve system of equations (Gradient = 0)
unconstrained_sol = solve([eq_1, eq_2], (x, y))
print(f"Optimal Point (Unconstrained): {unconstrained_sol}")

# --- PART B: Constrained Optimization (Street Network) ---
print("\nSolving Constrained Optimization...")

# Define the constraint line (y = mx + b) based on Paraguay Street points
m_slope = (res_y - ho_y) / (res_x - ho_x)
b_intercept = ho_y - m_slope * ho_x

# Constraint function g(x,y) = 0 => y - mx - b = 0
g = y - m_slope*x - b_intercept

# Lagrange Multiplier Method: grad(U) = L * grad(g)
grad_g1 = diff(g, x)
grad_g2 = diff(g, y)

# System of equations for Lagrange
# eq_1 is dU/dx, eq_2 is dU/dy
lagrange_eq1 = eq_1 - L * grad_g1
lagrange_eq2 = eq_2 - L * grad_g2
constraint_eq = g

constrained_sol = solve([lagrange_eq1, lagrange_eq2, constraint_eq], (x, y, L))
print(f"Optimal Point (Constrained): {constrained_sol}")

# ------------------------------------------------------------------------------
# 4. RESULTS CONVERSION
# ------------------------------------------------------------------------------

print("\n--- Final Results (Global Coordinates) ---")

# Extract numeric values from SymPy solution (converting from dictionary/tuple)
# Unconstrained
unc_x_local = float(unconstrained_sol[x])
unc_y_local = float(unconstrained_sol[y])

# Constrained (SymPy often returns a list of tuples, we take the first valid one)
if isinstance(constrained_sol, list):
    cons_x_local = float(constrained_sol[0][0])
    cons_y_local = float(constrained_sol[0][1])
elif isinstance(constrained_sol, dict):
    cons_x_local = float(constrained_sol[x])
    cons_y_local = float(constrained_sol[y])

# Convert back to Degrees
unc_long, unc_lat = kms_to_deg(unc_x_local, unc_y_local, ORIGIN_LONG, ORIGIN_LAT)
cons_long, cons_lat = kms_to_deg(cons_x_local, cons_y_local, ORIGIN_LONG, ORIGIN_LAT)

print(f"Unconstrained Optimal Location: {unc_long}, {unc_lat}")
print(f"Constrained Optimal Location:   {cons_long}, {cons_lat}")
