# Spatial Optimization & Coordinate Transformation (Python)

## Overview
This repository contains a Python-based analytical engine for solving spatial optimization problems. It was developed to determine optimal residential locations by minimizing weighted Euclidean distances to key amenities, both in an unconstrained 2D plane and subject to linear network constraints (e.g., specific street segments).

The project demonstrates the intersection of **geospatial analysis** and **mathematical optimization**, implementing custom projection algorithms to bridge the gap between global geodetic coordinates (Latitude/Longitude) and local Cartesian systems (Kilometers).

## Key Features
* **Geodetic Projections:** Implements the Equirectangular approximation to convert between spherical coordinates (degrees) and a local tangent plane (km) without relying on external GIS libraries.
* **Symbolic Optimization:** Uses `SymPy` to define and solve multivariate utility functions symbolically.
* **Constrained Optimization:** Solves for local maxima/minima along linear constraints using **Lagrange Multipliers**.
* **Inverse Transformation:** Automatically maps local optimized solutions back to real-world global coordinates.

## Mathematical Framework

### 1. Coordinate Transformation
To perform Euclidean optimization, global coordinates are projected onto a local flat surface using the standard approximation for small distances:
$$x = R \cdot \Delta\lambda \cdot \cos(\phi_{ref})$$
$$y = R \cdot \Delta\phi$$
Where $R$ is the Earth's radius (6371 km).

### 2. The Optimization Problem
The script minimizes a weighted cost function representing the "disutility" of distance from specific amenities (e.g., Disco, Obelisk, Galleries):

$$\min_{x,y} U(x,y) = \sum_{i=1}^{n} w_i \left[ (x - x_i)^2 + (y - y_i)^2 \right]$$

### 3. The Constraint
The constrained solution restricts the optimal point to lie along a vector defined by two points (representing a street segment):
$$g(x,y) = y - mx - b = 0$$
The system is solved using the Lagrangian:
$$\mathcal{L}(x, y, \lambda) = U(x,y) - \lambda \cdot g(x,y)$$

## Requirements
* Python 3.x
* `sympy`

## Usage
1. Install dependencies:
   ```bash
   pip install sympy
